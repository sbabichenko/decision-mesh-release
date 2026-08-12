#include <array>
#include <string>
#include "vertex.h"
#include "edge.h"
#include "face.h"
#include "hierarchical_design.h"
#include "mesh.h"
#include <cmath>
#include <map>
#include <limits>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <numeric>
#include "units.h"

// Propagate a committed height change through geometry-only descendants before
// recomputing any local regressions.  This is the minimal hierarchical-state
// implementation: a constrained midpoint has no independent coefficient and
// therefore must satisfy h_v = mu_lin(v) exactly as either parent moves.
//
// The full hierarchical-basis solver will eventually make this propagation
// implicit in the design matrix.  For now, doing it explicitly keeps the
// existing nodal quadratic solver intact while preserving the same surface.
static inline void refresh_after_height_change(Vertex* source) {
    static thread_local std::vector<Vertex*> locked_queue;
    static thread_local std::vector<Vertex*> to_update;
    static thread_local std::set<Vertex*> seen_locked;

    locked_queue.clear();
    to_update.clear();
    seen_locked.clear();

    to_update.insert(to_update.end(), source->affected_vertices.begin(),
                     source->affected_vertices.end());
    to_update.insert(to_update.end(), source->prior_children.begin(),
                     source->prior_children.end());

    for (Vertex* child : source->prior_children) {
        if (child && child->active && child->parent_edge &&
            !child->free_coefficient) {
            if (seen_locked.insert(child).second) locked_queue.push_back(child);
        }
    }

    // Birth depth is a topological order for midpoint dependencies.  Sorting
    // makes a constrained child's target use already-synchronized parents and
    // also keeps the operation deterministic.
    std::sort(locked_queue.begin(), locked_queue.end(),
              [](const Vertex* a, const Vertex* b) {
                  if (a->depth != b->depth) return a->depth < b->depth;
                  return a->_id < b->_id;
              });

    for (std::size_t qi = 0; qi < locked_queue.size(); ++qi) {
        Vertex* vertex = locked_queue[qi];
        if (!vertex || !vertex->active || vertex->free_coefficient ||
            !vertex->parent_edge) {
            continue;
        }

        vertex->height = vertex->mu_lin();
        vertex->new_height = vertex->height;
        vertex->delta_pooled = 0.0;
        vertex->lambda_v = 0.0;
        vertex->sigma_pooled = 0.0;
        vertex->loss_reduction = 0.0;
        vertex->mesh->heap_set(vertex, 0.0);

        to_update.insert(to_update.end(), vertex->affected_vertices.begin(),
                         vertex->affected_vertices.end());
        to_update.insert(to_update.end(), vertex->prior_children.begin(),
                         vertex->prior_children.end());

        for (Vertex* child : vertex->prior_children) {
            if (child && child->active && child->parent_edge &&
                !child->free_coefficient) {
                if (seen_locked.insert(child).second) locked_queue.push_back(child);
            }
        }
    }

    std::sort(to_update.begin(), to_update.end());
    to_update.erase(std::unique(to_update.begin(), to_update.end()),
                    to_update.end());
    for (Vertex* vertex : to_update) {
        if (vertex && vertex != source) vertex->update_info();
    }
}


extern CascadeCounters g_counters;

int Vertex::_seq = 0;

double Vertex::minimum_effective_support() {
    static const double value = [] {
        const char* raw = std::getenv("DMESH_MIN_CURRENT_KEFF");
        if (!raw) return 0.5;
        const double parsed = std::atof(raw);
        return (std::isfinite(parsed) && parsed >= 0.0) ? parsed : 0.5;
    }();
    return value;
}

Vertex::Vertex(DecisionMesh* mesh, double x, double y, bool active,
               Edge* parent_edge, bool real_vertex)
    : _id(Vertex::_seq++), x(x), y(y), active(active),
      real_vertex(real_vertex), mesh(mesh), parent_edge(parent_edge)
{
    // Initial/root vertices are the unthresholded coarse model and are free.
    // Inactive midpoint candidates are not free until explicitly admitted.
    free_coefficient = active;
    if (mesh) mesh->coord_map[DecisionMesh::coord_key(x, y)] = this;
    if (real_vertex) {
        mesh->vertices.insert(this);
    }
    if (parent_edge) {
        height = (parent_edge->vertex0->height + parent_edge->vertex1->height) / 2.0;
        depth = std::max(parent_edge->vertex0->depth, parent_edge->vertex1->depth) + 1;
        // Register with parent edge endpoints for prior propagation
        if (real_vertex) {
            parent_edge->vertex0->prior_children.insert(this);
            parent_edge->vertex1->prior_children.insert(this);
        }
    } else {
        height = 0.0;
        depth = 0;
    }
}

std::string Vertex::sid() const {
    char buf[16];
    snprintf(buf, sizeof(buf), "V%03d", _id);
    return buf;
}

void Vertex::add_edge(Edge* e) {
    edges.insert(e);
    Vertex* other = e->other_vertex(this);
    if (other) neighbors.insert(other);
}

void Vertex::remove_edge(Edge* e) {
    edges.erase(e);
    Vertex* other = e->other_vertex(this);
    if (other) neighbors.erase(other);
}

void Vertex::activate(bool make_free_coefficient) {
    height_frozen = height;
    if (active) {
        if (make_free_coefficient && !retired_coefficient) {
            free_coefficient = true;
            dormant_coefficient = false;
        }
        return;
    }
    g_counters.cascade_activations++;
    active = true;
    free_coefficient = make_free_coefficient;
    dormant_coefficient = false;
    retired_coefficient = false;
    if (!make_free_coefficient) conformity_origin = true;
    // A conformity-only midpoint must reproduce the pre-split P1 surface.
    // Never inherit the candidate's potentially enormous local proposal.
    height = free_coefficient ? new_height : mu_lin();

    // Set EB posterior fields before split (new midpoints need these during split)
    if (mesh->use_eb && parent_edge != nullptr) {
        delta_pooled = height - mu_lin();
        double inv_sigma = (sigma_sq < 1e300) ? (1.0 / sigma_sq) : 0.0;
        double xTx_reg = inv_sigma + lambda_v;
        sigma_pooled = (xTx_reg > 0) ? (1.0 / xTx_reg) : 1e300;
    }

    parent_edge->split();

    refresh_after_height_change(this);
    loss_reduction = 0;
    mesh->heap_set(this, 0);
}

void Vertex::update_info() {
    g_counters.update_info_calls++;
    auto result = loc_regress();
    affected_vertices = result.affected;

    double xTx = result.xTx;
    double xTr = result.xTr;
    double rTr = result.rTr;

    sigma_sq = (xTx > 0) ? (1.0 / xTx) : 1e300;
    current_n_points = result.n_points;
    current_fit_points = result.fit_points;
    current_xTx = xTx;
    current_fit_xTx = result.fit_xTx;
    current_keff = (result.xT2x > 0.0) ? (xTx * xTx / result.xT2x) : 0.0;

    // During adaptive refinement, exact loss of nodal support causes dormancy,
    // not retirement. Preserve the last accepted height so that an empty-region
    // reset cannot alter later interpolation, gate scores, or topology. Because
    // topology can still change, a dormant coefficient wakes automatically when
    // positive training curvature reappears.
    if (active && free_coefficient && parent_edge != nullptr) {
        if (!(current_fit_xTx > 0.0)) {
            dormant_coefficient = true;
            static const bool CONTAIN_H2 =
                std::getenv("DMESH_CONTAIN_HEIGHT") != nullptr;
            new_height = CONTAIN_H2 ? height
                                    : mu_prior() + contained_surplus;
            loss_reduction = 0.0;
            if (!disqualified) mesh->heap_set(this, 0.0);
            return;
        }
        dormant_coefficient = false;
    }

    // Geometric completion and finally retired vertices are not statistical degrees of
    // freedom. Their zero-surplus constraint is dynamic because parent heights
    // may move.
    if (active && parent_edge != nullptr && !free_coefficient) {
        const double target = mu_lin();
        new_height = target;
        delta_pooled = 0.0;
        lambda_v = 0.0;
        sigma_pooled = 0.0;
        const double diff = height - target;
        loss_reduction = diff * diff;
        if (!disqualified) mesh->heap_set(this, -loss_reduction);
        return;
    }

    // Admission-time support is not enough in a nodal basis: deeper
    // refinement can collapse a formerly broad star to only a handful of
    // training observations.  Until the full hierarchical design is in
    // place, freeze such a coefficient at its last accepted value rather than
    // letting a near-unidentified coordinate produce an extreme update.
    if (active && free_coefficient &&
        (current_fit_points < kMinimumFitPoints ||
         current_keff < minimum_effective_support())) {
        static const bool CONTAIN_HEIGHT =
            std::getenv("DMESH_CONTAIN_HEIGHT") != nullptr;  // legacy A/B
        if (CONTAIN_HEIGHT) {
            new_height = height;
        } else if (parent_edge == nullptr) {
            // Fix 1b: starved ROOT. No parents to track and no usable data;
            // free-at-initialization roots otherwise leak their init value
            // into any predictions their faces touch. Continue the local
            // field instead.
            new_height = neighbor_mean_height();
        } else {
            // Hold surplus: height tracks live parents plus the surplus at
            // this vertex's last committed fit.
            new_height = mu_prior() + contained_surplus;
        }
        loss_reduction = 0.0;
        if (!disqualified) mesh->heap_set(this, 0.0);
        return;
    }
    contained_state = false;  // support returned; normal fitting resumes

    if (mesh->use_eb) {
        auto [lv, mu_v] = compute_eb_params(xTx);
        lambda_v = lv;

        if (lv > 0 && xTx > 0) {
            double xTx_reg = xTx + lv;
            double xTr_reg = xTr + lv * mu_v;
            double beta_reg = xTr_reg / xTx_reg;

            double beta_orig = height;
            double orig_loss = beta_orig * beta_orig * xTx
                             - 2.0 * beta_orig * xTr + rTr
                             + lv * (beta_orig - mu_v) * (beta_orig - mu_v);
            double sse_post = rTr + lv * mu_v * mu_v - xTr_reg * xTr_reg / xTx_reg;

            new_height = beta_reg;
            loss_reduction = std::max(0.0, orig_loss - sse_post);
            sigma_pooled = 1.0 / xTx_reg;
            delta_pooled = beta_reg - mu_lin();
            // Optional legacy per-sweep displacement cap. The former
            // saturation-triggered exact fallback was removed after B1: all
            // admitted coordinates are re-profiled and committed on the exact
            // hierarchical objective, and the fallback was regression-inert.
            static const double COMMIT_CAP = []{
                const char* e = std::getenv("DMESH_COMMIT_CAP");
                return (e ? std::atof(e) : 0.0) * kHeightScale;
            }();
            if (COMMIT_CAP > 0.0) {
                new_height = height +
                    std::clamp(new_height - height, -COMMIT_CAP, COMMIT_CAP);
            }
            // DMESH_EXACT_SCORE: when the quadratic's own proposal exceeds the
            // trust radius (~0.5 logits) or the star is boundary-clipped, refit
            // the commit and heap gain from the exact tempered likelihood.
            static const bool EXACT = std::getenv("DMESH_EXACT_SCORE") != nullptr;
            if (EXACT) {
                const double raw_disp = std::fabs(xTr / xTx - mu_v);
                const double eps = 1e-9;
                const bool onb = x < eps || x > 1 - eps || y < eps || y > 1 - eps;
                if (raw_disp > 0.5 * kHeightScale || onb) {
                    auto es = exact_score(mu_v, lv);
                    if (es.ok) {
                        // Heap gain and gate displacement: always honest.
                        loss_reduction = std::max(0.0, es.W);
                        delta_pooled = es.b_star_units - mu_lin();
                        // Commit override only under DMESH_EXACT_COMMIT: the
                        // three-attempt record shows exact commits fight the
                        // working-scale sweep semantics.
                        static const bool EXC = std::getenv("DMESH_EXACT_COMMIT") != nullptr;
                        if (EXC) new_height = es.b_tax_units;
                    }
                }
            }
        } else {
            new_height = result.beta_opt;
            loss_reduction = result.loss_reduction;
            if (parent_edge != nullptr && xTx > 0) {
                delta_pooled = result.beta_opt - mu_lin();
                sigma_pooled = sigma_sq;
            }
        }
    } else {
        new_height = result.beta_opt;
        loss_reduction = result.loss_reduction;
    }

    if (!disqualified) {
        mesh->heap_set(this, -loss_reduction);
    }
}

namespace {
// DMESH_HEIGHT_LIMIT_LOGIT: reinstated safety clamp at a configurable
// threshold (legacy semantics: RESET to parent midpoint, not a clip).
// Default OFF (0). At 12 logits it cannot bind legitimate structure (the
// deepest law-marginal dot is -7.95, heights ~ -9.3) but catches the
// split-dependent runaway eruptions the lifecycle mechanisms miss.
double height_limit_units() {
    static const double L = []{
        const char* e = std::getenv("DMESH_HEIGHT_LIMIT_LOGIT");
        return e ? std::atof(e) * 100.0 : 0.0;
    }();
    return L;
}
}  // namespace

void Vertex::clamp_height() {
    const double L = height_limit_units();
    if (!std::isfinite(height)) {
        height = parent_edge ? mu_lin() : 0.0;
        new_height = height;
        return;
    }
    if (L <= 0.0 || std::fabs(height) <= L) return;
    // DMESH_HEIGHT_LIMIT_MODE=reset restores the legacy trapdoor; default is
    // CLIP (saturate at the wall). The reset injects a discontinuous kick on
    // every violation, which feeds the oscillation it is meant to contain;
    // a clip is an absorbing boundary.
    static const bool RESET = []{
        const char* e = std::getenv("DMESH_HEIGHT_LIMIT_MODE");
        return e && std::string(e) == "reset";
    }();
    if (RESET) height = parent_edge ? mu_lin() : 0.0;
    else height = (height > 0.0) ? L : -L;
    new_height = height;
}

void Vertex::update_height() {
    height = new_height;
    clamp_height();
    // Fix 1 (v2): maintain the last-committed surplus continuously. Holds
    // (containment, dormancy) then track mu_prior(live parents) + this
    // surplus, so parent motion after a vertex's last real fit propagates
    // instead of embalming a stale absolute height.
    if (parent_edge != nullptr) contained_surplus = height - mu_prior();
    if (parent_edge != nullptr) {
        delta_pooled = height - mu_lin();
    }
    refresh_after_height_change(this);
    loss_reduction = 0;
    mesh->heap_set(this, 0);
}

Vertex::FacesResult Vertex::get_faces(bool want_neighbors) const {
    FacesResult result;
    result.faces.reserve(8);
    for (Edge* e : edges) {
        for (int s = 0; s < 2; ++s) {
            Face* f = e->faces[s];
            if (f) {
                bool has_self = false;
                for (int i = 0; i < 3; ++i) {
                    if (f->vertices[i] == this) { has_self = true; break; }
                }
                if (has_self) result.faces.push_back(f);
            }
        }
    }
    std::sort(result.faces.begin(), result.faces.end());
    result.faces.erase(std::unique(result.faces.begin(), result.faces.end()),
                       result.faces.end());
    if (want_neighbors) {
        result.neighbors.reserve(result.faces.size() * 2);
        for (Face* f : result.faces) {
            for (int i = 0; i < 3; ++i)
                if (f->vertices[i] != this) result.neighbors.push_back(f->vertices[i]);
        }
        std::sort(result.neighbors.begin(), result.neighbors.end());
        result.neighbors.erase(std::unique(result.neighbors.begin(), result.neighbors.end()),
                               result.neighbors.end());
    }
    return result;
}

Vertex::FacesResult Vertex::get_sim_faces(bool want_neighbors) const {
    FacesResult result;
    if (!parent_edge) return result;

    result.faces.reserve(4);
    for (int s = 0; s < 2; ++s) {
        Face* face = parent_edge->faces[s];
        if (!face) continue;
        int idx = -1;
        for (int i = 0; i < 3; ++i) {
            if (face->edges[i] == parent_edge) { idx = i; break; }
        }
        if (idx < 0) continue;
        if (face->sub_divisions[idx].plus_face)
            result.faces.push_back(face->sub_divisions[idx].plus_face);
        if (face->sub_divisions[idx].minus_face)
            result.faces.push_back(face->sub_divisions[idx].minus_face);
    }
    std::sort(result.faces.begin(), result.faces.end());
    result.faces.erase(std::unique(result.faces.begin(), result.faces.end()),
                       result.faces.end());
    if (want_neighbors) {
        result.neighbors.reserve(result.faces.size() * 2);
        for (Face* f : result.faces) {
            for (int i = 0; i < 3; ++i)
                if (f->vertices[i] != this) result.neighbors.push_back(f->vertices[i]);
        }
        std::sort(result.neighbors.begin(), result.neighbors.end());
        result.neighbors.erase(std::unique(result.neighbors.begin(), result.neighbors.end()),
                               result.neighbors.end());
    }
    return result;
}

std::set<Vertex*> Vertex::neighbors_after_split() const {
    std::set<Vertex*> nbrs;
    if (!parent_edge) return nbrs;
    nbrs.insert(const_cast<Vertex*>(this));
    nbrs.insert(parent_edge->vertex0);
    nbrs.insert(parent_edge->vertex1);
    if (parent_edge->opposing_vertices[0])
        nbrs.insert(parent_edge->opposing_vertices[0]);
    if (parent_edge->opposing_vertices[1])
        nbrs.insert(parent_edge->opposing_vertices[1]);
    return nbrs;
}

Vertex::LocRegressResult Vertex::loc_regress(bool want_affected) {
    LocRegressResult res;
    res.beta_opt = height;
    res.loss_reduction = 0.0;
    res.n_points = 0;
    res.fit_points = 0;
    res.fit_xTx = 0.0;
    res.xTx = 0.0;
    res.xTr = 0.0;
    res.rTr = 0.0;

    FacesResult fr;
    if (active) {
        fr = get_faces(want_affected);
    } else {
        fr = get_sim_faces(want_affected);
    }
    if (want_affected) res.affected = std::move(fr.neighbors);

    if (fr.faces.empty()) {
        return res;
    }

    double xTx = 0.0, xTr = 0.0, rTr = 0.0;
    int total_points = 0;

    for (Face* f : fr.faces) {
        int si = -1;
        for (int i = 0; i < 3; ++i) {
            if (f->vertices[i] == this) { si = i; break; }
        }
        if (si < 0 || f->n_covered == 0) continue;

        total_points += f->n_covered;
        res.fit_points += f->n_fit_covered;

        double h[3];
        for (int i = 0; i < 3; ++i)
            h[i] = f->vertices[i]->height;

        xTx += f->S_ww[si][si];
        res.fit_xTx += f->S_fit_xx[si];
        res.xT2x += f->S_w2[si];
        res.xT2xs += f->S_w2s[si];
        res.xT4k += f->S_w4k[si];

        double xTr_face = f->S_wy[si];
        for (int j = 0; j < 3; ++j) {
            if (j == si) continue;
            xTr_face -= h[j] * f->S_ww[si][j];
        }
        xTr += xTr_face;

        double rTr_face = f->S_yy;
        for (int j = 0; j < 3; ++j) {
            if (j == si) continue;
            rTr_face -= 2.0 * h[j] * f->S_wy[j];
            for (int k = 0; k < 3; ++k) {
                if (k == si) continue;
                rTr_face += h[j] * h[k] * f->S_ww[j][k];
            }
        }
        rTr += rTr_face;
    }

    if (total_points == 0) return res;

    double beta_orig = height;
    double orig_loss = rTr - 2.0 * beta_orig * xTr + beta_orig * beta_orig * xTx;

    double beta_opt, sse_post;
    if (xTx > 0.0) {
        beta_opt = xTr / xTx;
        sse_post = rTr - (xTr * xTr) / xTx;
    } else {
        beta_opt = beta_orig;
        sse_post = rTr;
    }

    res.beta_opt = beta_opt;
    res.loss_reduction = orig_loss - sse_post;
    res.n_points = total_points;
    res.xTx = xTx;
    res.xTr = xTr;
    res.rTr = rTr;
    return res;
}

// --- Empirical Bayes helpers ---

double Vertex::mu_lin() const {
    if (parent_edge == nullptr) return 0.0;
    return (parent_edge->vertex0->height + parent_edge->vertex1->height) / 2.0;
}

double Vertex::mu_prior() const {
    // Shrink target. Default: parent midpoint (zero-second-difference prior,
    // which taxes curvature). DMESH_QUAD_PRIOR=1: locally quadratic target
    // along the parent edge (Dubuc-Deslauriers 4-point where both collinear
    // extensions exist; one-sided quadratic at boundaries), so the prior
    // carries convexity and surplus prices deviation-from-quadratic.
    if (parent_edge == nullptr) return 0.0;
    static const bool QUAD = std::getenv("DMESH_QUAD_PRIOR") != nullptr;
    // DMESH_GHOST_BC=1|2: implicit boundary condition on the shrink target.
    // For vertices ON the boundary, extrapolate along the INWARD NORMAL
    // (ghost-cell style), not along the parent edge (which runs along the
    // boundary and is itself dragged). 1 = first-derivative continuation
    // mu = 2 v(h) - v(2h); 2 = second-derivative mu = 3 v(h) - 3 v(2h) + v(3h).
    static const char* gbc_env = std::getenv("DMESH_GHOST_BC");
    static const int GBC = gbc_env ? std::atoi(gbc_env) : 0;
    if (GBC == 1 || GBC == 2) {
        const double eps = 1e-9;
        double nx = 0.0, ny = 0.0;
        if (x < eps) nx = 1.0; else if (x > 1.0 - eps) nx = -1.0;
        if (y < eps) ny = 1.0; else if (y > 1.0 - eps) ny = -1.0;
        // one boundary coordinate only (corners fall through to default)
        if ((nx != 0.0) != (ny != 0.0)) {
            for (int j = 10; j >= 1; --j) {
                const double h = std::ldexp(1.0, -j);
                const Vertex* v1 = mesh->vertex_at(x + h * nx, y + h * ny);
                const Vertex* v2 = mesh->vertex_at(x + 2 * h * nx, y + 2 * h * ny);
                if (!v1 || !v2) continue;
                if (GBC == 1) return 2.0 * v1->height - v2->height;
                const Vertex* v3 = mesh->vertex_at(x + 3 * h * nx, y + 3 * h * ny);
                if (!v3) continue;
                return 3.0 * v1->height - 3.0 * v2->height + v3->height;
            }
        }
    }
    if (!QUAD) {
        const Vertex* p0 = parent_edge->vertex0;
        const Vertex* p1 = parent_edge->vertex1;
        return 0.5 * (p0->height + p1->height);
    }
    return mu_quad();
}

double Vertex::mu_quad() const {
    // Locally quadratic target along the parent edge.  The symmetric stencil
    // is the midpoint value implied by the average left/right second
    // difference (equivalently the interpolatory four-point rule); a
    // one-sided quadratic is used when only one extension exists.
    double value = 0.0;
    for (const auto& [vertex, weight] : prior_stencil())
        value += weight * vertex->height;
    return value;
}

std::vector<std::pair<const Vertex*, double>> Vertex::prior_stencil() const {
    std::vector<std::pair<const Vertex*, double>> stencil;
    if (parent_edge == nullptr) return stencil;
    const Vertex* p0 = parent_edge->vertex0;
    const Vertex* p1 = parent_edge->vertex1;
    static const bool QUAD = std::getenv("DMESH_QUAD_PRIOR") != nullptr;
    if (!QUAD) {
        stencil.push_back({p0, 0.5});
        stencil.push_back({p1, 0.5});
        return stencil;
    }
    const Vertex* a = mesh->vertex_at(2 * p0->x - p1->x, 2 * p0->y - p1->y);
    const Vertex* b = mesh->vertex_at(2 * p1->x - p0->x, 2 * p1->y - p0->y);
    if (a && b) {
        stencil.push_back({a, -1.0 / 16.0});
        stencil.push_back({p0, 9.0 / 16.0});
        stencil.push_back({p1, 9.0 / 16.0});
        stencil.push_back({b, -1.0 / 16.0});
    } else if (a) {
        stencil.push_back({a, -1.0 / 8.0});
        stencil.push_back({p0, 6.0 / 8.0});
        stencil.push_back({p1, 3.0 / 8.0});
    } else if (b) {
        stencil.push_back({p0, 3.0 / 8.0});
        stencil.push_back({p1, 6.0 / 8.0});
        stencil.push_back({b, -1.0 / 8.0});
    } else {
        stencil.push_back({p0, 0.5});
        stencil.push_back({p1, 0.5});
    }
    return stencil;
}



double Vertex::neighbor_mean_height() {
    // Mean height of adjacent active vertices (face-weighted). Used to hold
    // starved roots: a coefficient with no usable data and no parents must
    // continue the local field, not its initialization value.
    FacesResult fr = get_faces(false);
    double sum = 0.0; int cnt = 0;
    for (Face* f : fr.faces) {
        for (Vertex* v : f->vertices) {
            if (v != this && v->active) { sum += v->height; ++cnt; }
        }
    }
    return cnt > 0 ? sum / cnt : height;
}

Vertex::ExactScore Vertex::exact_score(double mu_units, double lam_prec_units,
                                       bool use_live_neighbors) {
    // Profile the exact (tempered) binomial log-likelihood in this vertex's
    // own coefficient, holding neighbors fixed. Units: heights are
    // centi-log-odds; the Newton runs in logits. Temper c_i reproduces the
    // law weights' per-pool information discount exactly where the quadratic
    // path applies it, so mega-hostage protection carries over. Held-out
    // pools enter with weight ~1e-12 and are excluded by the same route.
    ExactScore es;
    if (mesh->obs_store == nullptr || mesh->refresh_epoch == 0) return es;
    const auto& obs = *mesh->obs_store;
    FacesResult fr = active ? get_faces(false) : get_sim_faces(false);
    if (fr.faces.empty()) return es;

    std::vector<double> phi, eta0, cc, kk, nn;
    for (Face* f : fr.faces) {
        int si = -1;
        for (int i = 0; i < 3; ++i) if (f->vertices[i] == this) { si = i; break; }
        if (si < 0 || f->n_covered == 0) continue;
        const Vertex* v0 = f->vertices[0];
        const Vertex* v1 = f->vertices[1];
        const Vertex* v2 = f->vertices[2];
        const double ax = v1->x - v0->x, ay = v1->y - v0->y;
        const double bx = v2->x - v0->x, by = v2->y - v0->y;
        const double d00 = ax*ax + ay*ay, d11 = bx*bx + by*by;
        if (d00 < 1e-14 || d11 < 1e-14) continue;
        const double h[3] = {
            use_live_neighbors ? v0->height : v0->height_frozen,
            use_live_neighbors ? v1->height : v1->height_frozen,
            use_live_neighbors ? v2->height : v2->height_frozen };
        for (int idx : f->coords_indices) {
            const Observation& o = obs[idx];
            const double px = o.x - v0->x, py = o.y - v0->y;
            const double u = (px*ax + py*ay) / d00;
            const double vv = (px*bx + py*by) / d11;
            const double ph3[3] = { 1.0 - u - vv, u, vv };
            const double phi_i = ph3[si];
            if (phi_i <= 1e-9) continue;
            // eta at the frozen epoch, rebuilt from frozen HEIGHTS through the
            // CURRENT geometry (stored per-pool eta goes stale when faces split
            // after a refresh; heights frozen at last touch cannot mismatch).
            const double eta_epoch =
                o.offset + (h[0]*ph3[0] + h[1]*ph3[1] + h[2]*ph3[2]) / kHeightScale;
            double eta_other = eta_epoch - h[si] * ph3[si] / kHeightScale;
            const double eta_cur = eta_epoch;
            const double p_cur = 1.0 / (1.0 + std::exp(-eta_cur));
            const double w_bin = o.exposure * p_cur * (1.0 - p_cur);
            if (w_bin <= 0) continue;
            const double c = (idx < (int)mesh->likelihood_temper.size())
                ? mesh->likelihood_temper[idx]
                : std::clamp(mesh->weights[idx] * kVarianceScale / w_bin,
                             0.0, 1.0);
            if (c < 1e-8) continue;   // held-out or fully discounted
            phi.push_back(phi_i); eta0.push_back(eta_other);
            cc.push_back(c); kk.push_back((double)o.events); nn.push_back((double)o.exposure);
        }
    }
    if (phi.size() < 10) return es;

    auto softplus = [](double x) {
        return x > 30 ? x : (x < -30 ? std::exp(x) : std::log1p(std::exp(x)));
    };
    auto ell = [&](double b) {
        double L = 0;
        for (size_t i = 0; i < phi.size(); ++i) {
            const double e = eta0[i] + b * phi[i];
            L += cc[i] * (kk[i] * e - nn[i] * softplus(e));
        }
        return L;
    };
    auto newton = [&](double b, double lam_logit, double mu_logit) {
        for (int it = 0; it < 60; ++it) {
            double g = -lam_logit * (b - mu_logit), H = -lam_logit;
            for (size_t i = 0; i < phi.size(); ++i) {
                const double e = eta0[i] + b * phi[i];
                const double p = 1.0 / (1.0 + std::exp(-e));
                g += cc[i] * phi[i] * (kk[i] - nn[i] * p);
                H -= cc[i] * phi[i] * phi[i] * nn[i] * p * (1.0 - p);
            }
            if (H >= -1e-12) break;
            double step = std::clamp(g / (-H), -1.5, 1.5);
            b += step;
            if (std::fabs(step) < 1e-7) break;
        }
        return b;
    };
    const double mu_logit = mu_units / kHeightScale;
    const double b_start = height / kHeightScale;
    const double b_star = newton(b_start, 0.0, 0.0);
    const double lam_logit = std::max(lam_prec_units, 0.0) * kVarianceScale;
    const double b_tax = (lam_logit > 0.0) ? newton(b_star, lam_logit, mu_logit) : b_star;
    es.ok = true;
    es.b_star_units = b_star * kHeightScale;
    es.b_tax_units = b_tax * kHeightScale;
    es.W = 2.0 * (ell(b_star) - ell(mu_logit));
    return es;
}

Vertex::ExactCandidateGain Vertex::exact_candidate_gain(
    double null_units,
    double prior_mean_units,
    double lam_prec_units,
    bool validate_grid) const {
    ExactCandidateGain empty;
    empty.null_units = null_units;
    empty.prior_mean_units = prior_mean_units;
    empty.lambda_logit = std::max(lam_prec_units, 0.0) * kVarianceScale;
    if (mesh == nullptr || mesh->obs_store == nullptr) return empty;

    std::vector<std::pair<int, double>> column;
    double current_units = null_units;
    if (!active) {
        // Exact hypothetical midpoint column on each CURRENT parent face.
        // Face::coords_indices can retain historical entries after later
        // splits, so point_face/point_bary are the canonical membership and
        // barycentric coordinates here as in build_active_design.
        if (parent_edge == nullptr) return empty;
        for (Face* face : parent_edge->faces) {
            if (face == nullptr || face->coords_indices.empty()) continue;
            int ia = -1, ib = -1;
            for (int j = 0; j < 3; ++j) {
                if (face->vertices[j] == parent_edge->vertex0) ia = j;
                if (face->vertices[j] == parent_edge->vertex1) ib = j;
            }
            if (ia < 0 || ib < 0) continue;
            for (int idx : face->coords_indices) {
                if (idx < 0 || idx >= mesh->n_points ||
                    mesh->point_face[idx] != face) continue;
                const auto& bc = mesh->point_bary[idx];
                const double value = 2.0 * std::min(bc[ia], bc[ib]);
                if (std::fabs(value) > 1e-12) column.push_back({idx, value});
            }
        }
    } else {
        hierarchical_design::BuildOptions options;
        options.training_only = false;
        column = hierarchical_design::build_active_column(*mesh, this, options);
        current_units = height;
    }

    // Shared profiling path for hypothetical candidates, active promotions,
    // commit validation, and the frozen-topology fixed-point audit.
    return exact_gain_on_column(column, current_units, null_units,
                                prior_mean_units, lam_prec_units,
                                validate_grid);
}


namespace {

struct HistogramLawConfig {
    double total_variance[4] = {0.35, 0.35, 0.35, 0.35};
    double mixture_weight[4] = {0.0, 0.0, 0.0, 0.0};
    double variance_ratio[4] = {1.0, 1.0, 1.0, 1.0};
    int explicit_rows = 8;
    bool use_point_variance = true;
};

const HistogramLawConfig& histogram_law_config() {
    static const HistogramLawConfig config = [] {
        HistogramLawConfig out;
        if (const char* gev = std::getenv("DMESH_GATE_EXTRA_VAR")) {
            const double v = std::atof(gev);
            if (std::isfinite(v) && v > 0.0)
                for (double& x : out.total_variance) x = v;
        }
        if (const char* tv = std::getenv("DMESH_HIST_TOTAL_VAR")) {
            double a, b, c, d;
            if (std::sscanf(tv, "%lf,%lf,%lf,%lf", &a, &b, &c, &d) == 4) {
                const double vals[4] = {a, b, c, d};
                for (int i = 0; i < 4; ++i)
                    if (std::isfinite(vals[i]) && vals[i] > 0.0)
                        out.total_variance[i] = vals[i];
            }
        }
        // Per stratum: mixture probability and wide/core variance ratio.
        // Example: "0.27,5.85;0.16,5.92;0.07,7.79;0.10,91.6".
        if (const char* mix = std::getenv("DMESH_HIST_MIX_SHAPE")) {
            double p0, r0, p1, r1, p2, r2, p3, r3;
            if (std::sscanf(mix,
                    "%lf,%lf;%lf,%lf;%lf,%lf;%lf,%lf",
                    &p0, &r0, &p1, &r1, &p2, &r2, &p3, &r3) == 8) {
                const double ps[4] = {p0, p1, p2, p3};
                const double rs[4] = {r0, r1, r2, r3};
                for (int i = 0; i < 4; ++i) {
                    out.mixture_weight[i] = std::clamp(ps[i], 0.0, 0.5);
                    out.variance_ratio[i] = std::max(rs[i], 1.0);
                }
            }
        }
        if (const char* top = std::getenv("DMESH_HIST_TOP")) {
            const int parsed = std::atoi(top);
            if (parsed >= 0) out.explicit_rows = std::min(parsed, 14);
        }
        if (const char* ignore = std::getenv("DMESH_HIST_IGNORE_POINT_VAR")) {
            out.use_point_variance = std::atoi(ignore) == 0;
        }
        return out;
    }();
    return config;
}

inline double stable_logistic(double x) {
    if (x >= 0.0) return 1.0 / (1.0 + std::exp(-x));
    const double ex = std::exp(x);
    return ex / (1.0 + ex);
}

inline double normal_cdf_hist(double z) {
    return 0.5 * std::erfc(-z / std::sqrt(2.0));
}

inline double upper_probit_hist(double p) {
    p = std::clamp(p, 1e-300, 0.5);
    const double t = std::sqrt(-2.0 * std::log(p));
    return t - (2.515517 + 0.802853*t + 0.010328*t*t)
             / (1.0 + 1.432788*t + 0.189269*t*t + 0.001308*t*t*t);
}

struct HistogramComponentMoment {
    double mean = 0.0;
    double variance = 0.0;
};

struct MarginalComponentDerivatives {
    bool ok = false;
    double log_marginal = 0.0;
    double score = 0.0;
    double second = 0.0;
};

inline double log1pexp_hist(double x) {
    if (x > 0.0) return x + std::log1p(std::exp(-x));
    return std::log1p(std::exp(x));
}

MarginalComponentDerivatives marginal_component_derivatives(
        double eta, double events, double exposure, double effect_variance) {
    MarginalComponentDerivatives out;
    const double s2 = std::max(effect_variance, 1e-10);
    // Posterior mode under the Gaussian component.  The score equation is
    // k - n p(eta+u) - u/s2 = 0 and is globally monotone.
    double u = 0.0;
    const double empirical_p = std::clamp((events + 0.5) / (exposure + 1.0),
                                          1e-10, 1.0 - 1e-10);
    u = std::clamp(std::log(empirical_p / (1.0 - empirical_p)) - eta,
                   -12.0, 12.0);
    for (int it = 0; it < 40; ++it) {
        const double p = stable_logistic(eta + u);
        const double curvature = exposure * p * (1.0 - p) + 1.0 / s2;
        const double gradient = events - exposure * p - u / s2;
        const double step = std::clamp(gradient / std::max(curvature, 1e-12),
                                       -2.0, 2.0);
        u += step;
        u = std::clamp(u, -20.0, 20.0);
        if (std::fabs(step) < 1e-10) break;
    }
    const double p_mode = stable_logistic(eta + u);
    const double curvature = exposure * p_mode * (1.0 - p_mode) + 1.0 / s2;
    const double scale = 1.0 / std::sqrt(std::max(curvature, 1e-12));

    // Nine-node mode-centered Gauss-Hermite integration.  These are the
    // probabilists' nodes with weights normalized to E_{N(0,1)}.
    static constexpr double nodes[9] = {
        -4.5127458633997826, -3.2054290028564698,
        -2.0768479786778302, -1.0232556637891326,
         0.0,
         1.0232556637891326,  2.0768479786778302,
         3.2054290028564698,  4.5127458633997826
    };
    static constexpr double weights[9] = {
        2.2345844007746607e-05, 0.0027891413212317692,
        0.049916406765217823,  0.24409750289493953,
        0.40634920634920635,
        0.24409750289493953,  0.049916406765217823,
        0.0027891413212317692, 2.2345844007746607e-05
    };
    double log_terms[9];
    double d1[9], d2[9];
    double max_log = -std::numeric_limits<double>::infinity();
    for (int j = 0; j < 9; ++j) {
        const double uj = u + scale * nodes[j];
        const double theta = eta + uj;
        const double pj = stable_logistic(theta);
        const double log_kernel = events * theta
            - exposure * log1pexp_hist(theta)
            - 0.5 * uj * uj / s2;
        // Change of variables around the posterior mode, with the x^2/2
        // correction because the quadrature weights integrate against phi(x).
        log_terms[j] = std::log(weights[j]) + 0.5 * nodes[j] * nodes[j]
            + log_kernel + std::log(scale) - 0.5 * std::log(s2);
        max_log = std::max(max_log, log_terms[j]);
        d1[j] = events - exposure * pj;
        d2[j] = -exposure * pj * (1.0 - pj);
    }
    double norm = 0.0, e1 = 0.0, e2 = 0.0, e11 = 0.0;
    for (int j = 0; j < 9; ++j) {
        const double w = std::exp(std::clamp(log_terms[j] - max_log,
                                             -745.0, 0.0));
        norm += w;
        e1 += w * d1[j];
        e2 += w * d2[j];
        e11 += w * d1[j] * d1[j];
    }
    if (!(norm > 0.0) || !std::isfinite(max_log)) return out;
    e1 /= norm;
    e2 /= norm;
    e11 /= norm;
    out.ok = true;
    out.log_marginal = max_log + std::log(norm);
    out.score = e1;
    out.second = e2 + e11 - e1 * e1;
    return out;
}

HistogramComponentMoment score_component_moment(double eta_null,
                                                 double p_null,
                                                 double exposure,
                                                 double score_loading,
                                                 double effect_variance) {
    // Five-node Gauss-Hermite rule for a standard normal expectation.
    static constexpr double nodes[5] = {
        -2.8569700138728056, -1.355626179974266,
         0.0,
         1.355626179974266,  2.8569700138728056
    };
    static constexpr double weights[5] = {
        0.01125741132772069, 0.2220759220056126,
        0.5333333333333333,
        0.2220759220056126, 0.01125741132772069
    };
    const double sd = std::sqrt(std::max(effect_variance, 0.0));
    double eq = 0.0, eq2 = 0.0;
    if (sd <= 1e-14) {
        eq = p_null;
        eq2 = p_null * p_null;
    } else {
        for (int j = 0; j < 5; ++j) {
            const double q = stable_logistic(eta_null + sd * nodes[j]);
            eq += weights[j] * q;
            eq2 += weights[j] * q * q;
        }
    }
    const double mean_k = exposure * eq;
    const double var_k = std::max(
        exposure * (eq - eq2)
        + exposure * exposure * (eq2 - eq * eq), 0.0);
    HistogramComponentMoment out;
    out.mean = score_loading * (mean_k - exposure * p_null);
    out.variance = score_loading * score_loading * var_k;
    return out;
}

}  // namespace

Vertex::HistogramNullResult Vertex::histogram_null_score(
    double null_units,
    double prior_mean_units,
    double lam_prec_units,
    const ExactCandidateGain& exact) const {
    (void)prior_mean_units;
    (void)lam_prec_units;
    HistogramNullResult out;
    if (!exact.ok || mesh == nullptr || mesh->obs_store == nullptr) return out;

    std::vector<std::pair<int, double>> column;
    double current_units = null_units;
    if (!active) {
        if (parent_edge == nullptr) return out;
        for (Face* face : parent_edge->faces) {
            if (face == nullptr || face->coords_indices.empty()) continue;
            int ia = -1, ib = -1;
            for (int j = 0; j < 3; ++j) {
                if (face->vertices[j] == parent_edge->vertex0) ia = j;
                if (face->vertices[j] == parent_edge->vertex1) ib = j;
            }
            if (ia < 0 || ib < 0) continue;
            for (int idx : face->coords_indices) {
                if (idx < 0 || idx >= mesh->n_points ||
                    mesh->point_face[idx] != face) continue;
                const auto& bc = mesh->point_bary[idx];
                const double value = 2.0 * std::min(bc[ia], bc[ib]);
                if (std::fabs(value) > 1e-12) column.push_back({idx, value});
            }
        }
    } else {
        hierarchical_design::BuildOptions options;
        options.training_only = false;
        column = hierarchical_design::build_active_column(*mesh, this, options);
        current_units = height;
    }
    if (column.empty()) return out;
    std::sort(column.begin(), column.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });
    std::vector<std::pair<int, double>> merged;
    merged.reserve(column.size());
    for (const auto& entry : column) {
        if (!merged.empty() && merged.back().first == entry.first)
            merged.back().second += entry.second;
        else
            merged.push_back(entry);
    }

    struct RowMoment {
        double pi = 0.0;
        HistogramComponentMoment core;
        HistogramComponentMoment wide;
        double mean = 0.0;
        double variance = 0.0;
        double importance = 0.0;
    };
    std::vector<RowMoment> moments;
    moments.reserve(merged.size());
    const auto& observations = *mesh->obs_store;
    const auto& law = histogram_law_config();
    const double active_height_adjustment = null_units - current_units;

    static const bool MARGINAL_SCORE = [] {
        const char* e = std::getenv("DMESH_HIST_MARGINAL");
        return e != nullptr && std::atoi(e) != 0;
    }();
    if (MARGINAL_SCORE) {
        double score = 0.0;
        double information = 0.0;
        double information_sq = 0.0;
        int used = 0;
        for (const auto& [idx, value] : merged) {
            if (idx < 0 || idx >= mesh->n_points || std::fabs(value) <= 1e-12)
                continue;
            const double train_temper = (idx < (int)mesh->likelihood_temper.size())
                ? mesh->likelihood_temper[idx] : 1.0;
            // The marginal law already prices pool heterogeneity.  Tempering
            // it again would double-discount the same nuisance variance; use
            // the frozen temper only as the train/held-out indicator.
            if (train_temper < 1e-8) continue;
            const Observation& observation = observations[idx];
            const double prediction_null = mesh->predict_index(idx)
                + value * active_height_adjustment;
            const double eta_null = observation.offset
                + prediction_null / kHeightScale;
            const double n = observation.exposure;
            const int stratum = n < 128.0 ? 0 : n < 512.0 ? 1
                                : n < 4096.0 ? 2 : 3;
            double total_s2 = law.total_variance[stratum];
            if (law.use_point_variance && observation.pool_variance > 0.0
                && std::isfinite(observation.pool_variance))
                total_s2 = observation.pool_variance;
            total_s2 = std::max(total_s2, 1e-10);
            const double pi = law.mixture_weight[stratum];
            const double ratio = law.variance_ratio[stratum];
            const double core_s2 = total_s2
                / std::max((1.0 - pi) + pi * ratio, 1e-12);
            const double wide_s2 = ratio * core_s2;
            const auto core = marginal_component_derivatives(
                eta_null, observation.events, observation.exposure, core_s2);
            const auto wide = marginal_component_derivatives(
                eta_null, observation.events, observation.exposure, wide_s2);
            if (!core.ok || !wide.ok) continue;
            const double la = std::log(std::max(1.0 - pi, 1e-300))
                            + core.log_marginal;
            const double lb = std::log(std::max(pi, 1e-300))
                            + wide.log_marginal;
            const double lm = std::max(la, lb);
            const double wa = std::exp(la - lm);
            const double wb = std::exp(lb - lm);
            const double norm = wa + wb;
            const double qa = wa / norm, qb = wb / norm;
            const double d1 = qa * core.score + qb * wide.score;
            const double d2 = qa * (core.second + core.score * core.score)
                            + qb * (wide.second + wide.score * wide.score)
                            - d1 * d1;
            const double row_score = value * d1;
            const double row_info = std::max(-value * value * d2, 0.0);
            score += row_score;
            information += row_info;
            information_sq += row_info * row_info;
            ++used;
        }
        if (!(information > 1e-12) || used == 0) return out;
        const double z = score / std::sqrt(information);
        const double p_value = std::clamp(
            std::erfc(std::fabs(z) / std::sqrt(2.0)), 1e-300, 1.0);
        out.ok = true;
        out.p_value = p_value;
        out.gaussian_equivalent = z;
        out.observed_data_score = score;
        out.null_mean = 0.0;
        out.null_sd = std::sqrt(information);
        out.tail_keff = information_sq > 0.0
            ? information * information / information_sq : 0.0;
        out.rows = used;
        out.explicit_mixture_rows = 0;
        return out;
    }
    double observed_data_score = 0.0;
    double variance_sq_sum = 0.0;
    double total_mean = 0.0;
    double total_variance = 0.0;

    for (const auto& [idx, value] : merged) {
        if (idx < 0 || idx >= mesh->n_points || std::fabs(value) <= 1e-12)
            continue;
        const Observation& observation = observations[idx];
        const double prediction_null = mesh->predict_index(idx)
            + value * active_height_adjustment;
        const double eta_null = observation.offset
            + prediction_null / kHeightScale;
        const double p_null = stable_logistic(eta_null);
        const double information = observation.exposure * p_null * (1.0 - p_null);
        if (!(information > 0.0)) continue;
        const double temper = (idx < (int)mesh->likelihood_temper.size())
            ? mesh->likelihood_temper[idx]
            : std::clamp(mesh->weights[idx] * kVarianceScale / information,
                         0.0, 1.0);
        if (temper < 1e-8) continue;
        const double loading = temper * value;
        observed_data_score += loading
            * (static_cast<double>(observation.events)
               - observation.exposure * p_null);

        const double n = observation.exposure;
        const int stratum = n < 128.0 ? 0 : n < 512.0 ? 1
                            : n < 4096.0 ? 2 : 3;
        double total_s2 = law.total_variance[stratum];
        if (law.use_point_variance && observation.pool_variance > 0.0
            && std::isfinite(observation.pool_variance)) {
            total_s2 = observation.pool_variance;
        }
        total_s2 = std::max(total_s2, 1e-10);
        const double pi = law.mixture_weight[stratum];
        const double ratio = law.variance_ratio[stratum];
        const double core_s2 = total_s2
            / std::max((1.0 - pi) + pi * ratio, 1e-12);
        const double wide_s2 = ratio * core_s2;

        RowMoment row;
        row.pi = pi;
        row.core = score_component_moment(
            eta_null, p_null, observation.exposure, loading, core_s2);
        row.wide = score_component_moment(
            eta_null, p_null, observation.exposure, loading, wide_s2);
        row.mean = (1.0 - pi) * row.core.mean + pi * row.wide.mean;
        const double second = (1.0 - pi)
            * (row.core.variance + row.core.mean * row.core.mean)
            + pi * (row.wide.variance + row.wide.mean * row.wide.mean);
        row.variance = std::max(second - row.mean * row.mean, 1e-18);
        const double mean_gap = row.wide.mean - row.core.mean;
        const double sd_gap = std::sqrt(std::max(row.wide.variance, 0.0))
                            - std::sqrt(std::max(row.core.variance, 0.0));
        row.importance = pi * (1.0 - pi)
            * (mean_gap * mean_gap + sd_gap * sd_gap);
        moments.push_back(row);
        total_mean += row.mean;
        total_variance += row.variance;
        variance_sq_sum += row.variance * row.variance;
    }
    if (moments.empty() || !(total_variance > 0.0)) return out;

    std::vector<int> order(moments.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        if (moments[a].importance != moments[b].importance)
            return moments[a].importance > moments[b].importance;
        return moments[a].variance > moments[b].variance;
    });
    std::vector<int> selected_rows;
    selected_rows.reserve(std::min<int>(law.explicit_rows, order.size()));
    for (int index : order) {
        if ((int)selected_rows.size() >= law.explicit_rows) break;
        const RowMoment& row = moments[index];
        if (row.pi <= 1e-8 || row.importance <= 1e-18) continue;
        selected_rows.push_back(index);
    }
    const int explicit_rows = static_cast<int>(selected_rows.size());

    const std::size_t states = std::size_t{1} << explicit_rows;
    double cdf = 0.0;
    double weight_sum = 0.0;
    for (std::size_t mask = 0; mask < states; ++mask) {
        double state_weight = 1.0;
        double state_mean = total_mean;
        double state_variance = total_variance;
        for (int j = 0; j < explicit_rows; ++j) {
            const RowMoment& row = moments[selected_rows[j]];
            const bool wide = ((mask >> j) & 1U) != 0;
            const double component_weight = wide ? row.pi : (1.0 - row.pi);
            const auto& component = wide ? row.wide : row.core;
            state_weight *= component_weight;
            state_mean += component.mean - row.mean;
            state_variance += component.variance - row.variance;
        }
        if (!(state_weight > 0.0)) continue;
        state_variance = std::max(state_variance, 1e-18);
        const double z = (observed_data_score - state_mean)
            / std::sqrt(state_variance);
        cdf += state_weight * normal_cdf_hist(z);
        weight_sum += state_weight;
    }
    if (!(weight_sum > 0.0)) return out;
    cdf = std::clamp(cdf / weight_sum, 0.0, 1.0);
    const double p_value = std::clamp(2.0 * std::min(cdf, 1.0 - cdf),
                                      1e-300, 1.0);
    const double centered = observed_data_score - total_mean;
    const double z_equiv = centered == 0.0 ? 0.0
        : (centered < 0.0 ? -1.0 : 1.0) * upper_probit_hist(0.5 * p_value);

    out.ok = true;
    out.p_value = p_value;
    out.gaussian_equivalent = z_equiv;
    out.observed_data_score = observed_data_score;
    out.null_mean = total_mean;
    out.null_sd = std::sqrt(total_variance);
    out.tail_keff = variance_sq_sum > 0.0
        ? total_variance * total_variance / variance_sq_sum : 0.0;
    out.rows = static_cast<int>(moments.size());
    out.explicit_mixture_rows = explicit_rows;
    return out;
}

Vertex::ExactCandidateGain Vertex::exact_gain_on_column(
    const std::vector<std::pair<int, double>>& input_column,
    double current_units,
    double null_units,
    double prior_mean_units,
    double lam_prec_units,
    bool validate_grid) const {
    ExactCandidateGain out;
    out.null_units = null_units;
    out.prior_mean_units = prior_mean_units;
    out.lambda_logit = std::max(lam_prec_units, 0.0) * kVarianceScale;
    if (mesh == nullptr || mesh->obs_store == nullptr || input_column.empty())
        return out;

    // Defensive merge: an observation on a shared face boundary can occur in
    // more than one local acceleration list.  The canonical active design is
    // unique by observation index, so exact profiling must be unique too.
    std::vector<std::pair<int, double>> column = input_column;
    std::sort(column.begin(), column.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });
    std::vector<std::pair<int, double>> merged;
    merged.reserve(column.size());
    for (const auto& entry : column) {
        if (!merged.empty() && merged.back().first == entry.first)
            merged.back().second += entry.second;
        else
            merged.push_back(entry);
    }

    const auto& obs = *mesh->obs_store;
    struct Row {
        int index = -1;
        double column = 0.0;
        double eta_null = 0.0;
        double temper = 0.0;
        double events = 0.0;
        double exposure = 0.0;
    };
    std::vector<Row> rows;
    rows.reserve(merged.size());
    const double active_height_adjustment = null_units - current_units;
    for (const auto& [idx, value] : merged) {
        if (idx < 0 || idx >= mesh->n_points || std::fabs(value) <= 1e-12)
            continue;
        const Observation& o = obs[idx];
        ++out.structural_points;
        if (mesh->weights[idx] > 1e-9) ++out.training_structural_points;
        const double prediction_null = mesh->predict_index(idx)
            + value * active_height_adjustment;
        const double eta_null = o.offset + prediction_null / kHeightScale;
        const double p = eta_null >= 0.0
            ? 1.0 / (1.0 + std::exp(-eta_null))
            : std::exp(eta_null) / (1.0 + std::exp(eta_null));
        const double binomial_information = o.exposure * p * (1.0 - p);
        if (!(binomial_information > 0.0)) continue;
        ++out.positive_information_points;
        const double temper = (idx < (int)mesh->likelihood_temper.size())
            ? mesh->likelihood_temper[idx]
            : std::clamp(mesh->weights[idx] * kVarianceScale
                             / binomial_information,
                         0.0, 1.0);
        if (temper < 1e-8) continue;
        rows.push_back({idx, value, eta_null, temper,
                        static_cast<double>(o.events),
                        static_cast<double>(o.exposure)});
    }

    // Optimization and certification are separate. Any nonempty likelihood
    // column has a finite MAP under positive EB precision; replication is
    // priced later by the gate calibration.
    if (rows.empty()) return out;
    out.fit_points = static_cast<int>(rows.size());

    auto softplus = [](double x) {
        if (x > 30.0) return x;
        if (x < -30.0) return std::exp(x);
        return std::log1p(std::exp(x));
    };
    auto logistic = [](double x) {
        if (x >= 0.0) return 1.0 / (1.0 + std::exp(-x));
        const double ex = std::exp(x);
        return ex / (1.0 + ex);
    };
    const double null_logit = null_units / kHeightScale;
    const double prior_logit = prior_mean_units / kHeightScale;
    const double lambda = out.lambda_logit;

    auto data_objective = [&](double b) {
        double value = 0.0;
        const double shift = b - null_logit;
        for (const Row& row : rows) {
            const double eta = row.eta_null + row.column * shift;
            value += row.temper
                * (row.events * eta - row.exposure * softplus(eta));
        }
        return value;
    };
    auto penalty_objective = [&](double b) {
        return -0.5 * lambda * (b - prior_logit) * (b - prior_logit);
    };
    auto objective = [&](double b) {
        return data_objective(b) + penalty_objective(b);
    };
    auto derivatives = [&](double b) {
        double gradient = -lambda * (b - prior_logit);
        double hessian = -lambda;
        const double shift = b - null_logit;
        for (const Row& row : rows) {
            const double eta = row.eta_null + row.column * shift;
            const double p = logistic(eta);
            gradient += row.temper * row.column
                * (row.events - row.exposure * p);
            hessian -= row.temper * row.column * row.column
                * row.exposure * p * (1.0 - p);
        }
        return std::pair<double, double>{gradient, hessian};
    };

    const auto [g_null, h_null] = derivatives(null_logit);
    out.score_at_null = g_null;
    out.information_at_null = -h_null;
    out.data_objective_null = data_objective(null_logit);
    out.objective_null = out.data_objective_null + penalty_objective(null_logit);
    out.min_eta_null = std::numeric_limits<double>::infinity();
    out.max_eta_null = -std::numeric_limits<double>::infinity();
    out.min_p_null = 1.0;
    out.max_p_null = 0.0;
    for (const Row& row : rows) {
        out.min_eta_null = std::min(out.min_eta_null, row.eta_null);
        out.max_eta_null = std::max(out.max_eta_null, row.eta_null);
        const double p0 = logistic(row.eta_null);
        out.min_p_null = std::min(out.min_p_null, p0);
        out.max_p_null = std::max(out.max_p_null, p0);
    }

    // Strict concavity gives a monotone score. Bracket its root before Newton
    // so saturated stars cannot generate a global quadratic jump.
    double lo = null_logit, hi = null_logit;
    double g_lo = g_null, g_hi = g_null;
    double width = 0.5;
    if (g_null > 0.0) {
        for (int it = 0; it < 20 && g_hi > 0.0; ++it) {
            hi = null_logit + width;
            g_hi = derivatives(hi).first;
            width *= 2.0;
        }
    } else if (g_null < 0.0) {
        for (int it = 0; it < 20 && g_lo < 0.0; ++it) {
            lo = null_logit - width;
            g_lo = derivatives(lo).first;
            width *= 2.0;
        }
    }
    if (g_lo < 0.0 || g_hi > 0.0) return out;

    double b = std::clamp(prior_logit, lo, hi);
    if (!(b > lo && b < hi)) b = 0.5 * (lo + hi);
    for (int it = 0; it < 100; ++it) {
        const auto [g, h] = derivatives(b);
        if (g > 0.0) lo = b; else hi = b;
        double proposal = 0.5 * (lo + hi);
        if (h < -1e-14) {
            const double newton = b - g / h;
            if (newton > lo && newton < hi && std::isfinite(newton))
                proposal = newton;
        }
        if (std::fabs(proposal - b) < 1e-10 || hi - lo < 1e-10) {
            b = proposal;
            break;
        }
        b = proposal;
    }

    out.data_objective_optimum = data_objective(b);
    out.objective_optimum = out.data_objective_optimum + penalty_objective(b);
    out.raw_gain_nats = out.objective_optimum - out.objective_null;
    if (out.raw_gain_nats < -1e-8 || !std::isfinite(out.raw_gain_nats))
        return out;
    if (out.raw_gain_nats < 0.0) {
        b = null_logit;
        out.data_objective_optimum = out.data_objective_null;
        out.objective_optimum = out.objective_null;
        out.raw_gain_nats = 0.0;
    }
    out.optimum_units = b * kHeightScale;
    out.gain_nats = out.raw_gain_nats;
    out.data_gain_nats = out.data_objective_optimum - out.data_objective_null;
    out.penalty_gain_nats = penalty_objective(b) - penalty_objective(null_logit);
    const double decomposed_gain = out.data_gain_nats + out.penalty_gain_nats;
    if (std::fabs(decomposed_gain - out.gain_nats)
        > 1e-9 * (1.0 + std::fabs(out.gain_nats)))
        return ExactCandidateGain{};
    const double displacement = b - null_logit;
    const double direction = (displacement > 0.0) - (displacement < 0.0);
    out.signed_root_gain = direction
        * std::sqrt(std::max(0.0, 2.0 * out.gain_nats));
    out.data_information_optimum =
        std::max(0.0, -derivatives(b).second - lambda);

    if (validate_grid) {
        const double left = std::min(null_logit, b) - 2.0;
        const double right = std::max(null_logit, b) + 2.0;
        double best_grid = -std::numeric_limits<double>::infinity();
        constexpr int kGrid = 2000;
        for (int j = 0; j <= kGrid; ++j) {
            const double q = left + (right - left) * j / kGrid;
            best_grid = std::max(best_grid, objective(q));
        }
        out.max_grid_excess = best_grid - out.objective_optimum;
        if (out.max_grid_excess > 2e-6) return ExactCandidateGain{};
    }

    out.ok = true;
    return out;
}

std::pair<double, double> Vertex::compute_eb_params(double xTx) const {
    if (parent_edge == nullptr || xTx <= 0) return {0.0, 0.0};

    double tau = tau_override_sq > 0.0
                     ? std::max(tau_override_sq, mesh->tau_floor_sq)
                     : mesh->tau_for_depth(depth);
    // DMESH_EDGE_TAU_SCALE=<f>: vertices with boundary-clipped stars carry
    // one-sided (roughly halved) information, and their ancestral shrink
    // targets are themselves boundary-dragged; scaling tau up by f weakens
    // the ancestral grip exactly there, letting boundary evidence speak.
    static const char* ets = std::getenv("DMESH_EDGE_TAU_SCALE");
    static const double edge_scale = ets ? std::atof(ets) : 1.0;
    if (edge_scale != 1.0) {
        const double eps = 1e-9;
        const bool on_boundary = x < eps || x > 1.0 - eps
                              || y < eps || y > 1.0 - eps;
        if (on_boundary && tau > 0.0) tau *= edge_scale;
    }
    const double precision = (tau <= 0.0) ? 1e12 : 1.0 / tau;
    return {precision, mu_prior()};
}
