#include "mesh.h"
#include <cstdlib>
#include "tree.h"
#include "vertex.h"
#include "edge.h"
#include "face.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <sstream>
#include <numeric>
#include <limits>

bool EdgeCoordLess::operator()(const Edge* a, const Edge* b) const {
    if (a == b) return false;
    if (a->vertex0->x != b->vertex0->x) return a->vertex0->x < b->vertex0->x;
    if (a->vertex0->y != b->vertex0->y) return a->vertex0->y < b->vertex0->y;
    if (a->vertex1->x != b->vertex1->x) return a->vertex1->x < b->vertex1->x;
    if (a->vertex1->y != b->vertex1->y) return a->vertex1->y < b->vertex1->y;
    return a->_id < b->_id;
}


CascadeCounters g_counters;

DecisionMesh::DecisionMesh(const std::vector<double>& x_data,
                           const std::vector<double>& y_data,
                           const std::vector<double>& z_data,
                           bool use_eb_flag)
    : use_eb(use_eb_flag), rng(42)
{
    if (const char* s = std::getenv("DMESH_TAU0_LOGIT_SD")) {
        const double sd = std::atof(s);
        if (std::isfinite(sd) && sd > 0) tau_default_sq = 1.0e4 * sd * sd;
    }
    if (const char* s = std::getenv("DMESH_TAU_FLOOR_LOGIT_SD")) {
        const double sd = std::atof(s);
        if (std::isfinite(sd) && sd > 0) tau_floor_sq = 1.0e4 * sd * sd;
    }
    n_points = (int)x_data.size();
    X.resize(n_points * 2);
    values = z_data;
    weights.assign(values.size(), 1.0);
    likelihood_temper.assign(values.size(), 1.0);
    point_face.assign(n_points, nullptr);
    point_bary.assign(n_points, {0.0, 0.0, 0.0});
    point_basis.assign(n_points, {0.0, 0.0, 0.0});
    for (int i = 0; i < n_points; ++i) {
        X[i * 2 + 0] = x_data[i];
        X[i * 2 + 1] = y_data[i];
    }

    root = make_node();
    leaves.insert(root);

    create_outer_vertices();
    create_outer_edges();

    split_edge = make_edge(outer_vertices["bottom_left"],
                           outer_vertices["top_right"], true);

    // Compute mask: X @ split_edge.normal >= split_edge.intercept
    std::vector<bool> mask_pos(n_points), mask_neg(n_points);
    for (int i = 0; i < n_points; ++i) {
        double val = split_edge->normal[0] * get_x(i, 0) +
                     split_edge->normal[1] * get_x(i, 1);
        mask_pos[i] = (val >= split_edge->intercept);
        mask_neg[i] = !mask_pos[i];
    }

    top_face = make_face(split_edge, outer_edges["left"], outer_edges["top"],
                         mask_pos, true, "+");
    bottom_face = make_face(split_edge, outer_edges["bottom"], outer_edges["right"],
                            mask_neg, true, "-");
    root->split(split_edge, top_face, bottom_face);

    for (Vertex* v : std::set<Vertex*>(vertices)) {
        v->update_info();
    }

    // Bootstrap empirical Bayes: compute tau_sq from initial estimates,
    // then re-update non-corner vertices with regularization
    if (use_eb) {
        recompute_tau_sq();
        for (Vertex* v : std::set<Vertex*>(vertices)) {
            if (v->parent_edge != nullptr) {
                v->update_info();
            }
        }
    }
}

DecisionMesh::~DecisionMesh() {
    for (auto* p : all_nodes) delete p;
    for (auto* p : all_vertices) delete p;
    for (auto* p : all_edges) delete p;
    for (auto* p : all_faces) delete p;
}

TreeNode* DecisionMesh::make_node(TreeNode* parent, Face* face) {
    auto* n = new TreeNode(this, parent, face);
    all_nodes.push_back(n);
    return n;
}

Vertex* DecisionMesh::make_vertex(double x, double y, bool active,
                                   Edge* parent_edge, bool real_vertex) {
    auto* v = new Vertex(this, x, y, active, parent_edge, real_vertex);
    all_vertices.push_back(v);
    return v;
}

Edge* DecisionMesh::make_edge(Vertex* v0, Vertex* v1, bool active) {
    g_counters.edges_created++;
    auto* e = new Edge(this, v0, v1, active);
    all_edges.push_back(e);
    return e;
}

Face* DecisionMesh::make_face(Edge* e0, Edge* e1, Edge* e2,
                               const std::vector<bool>& mask, bool active,
                               const std::string& path, bool skip_coords) {
    g_counters.faces_created++;
    auto* f = new Face(this, e0, e1, e2, mask, active, path, skip_coords);
    all_faces.push_back(f);
    return f;
}

// --- Heap operations using std::map<double, std::set<Vertex*>> ---

void DecisionMesh::heap_set(Vertex* v, double neg_loss) {
    if (!maintain_heap) return;
    // O(log n) removal of old entry using cached key
    if (v->in_heap) {
        auto it = loss_heap.find(v->heap_key);
        if (it != loss_heap.end()) {
            it->second.erase(v);
            if (it->second.empty()) {
                loss_heap.erase(it);
            }
        }
        v->in_heap = false;
    }
    // O(log n) insertion at new key
    loss_heap[neg_loss].insert(v);
    v->heap_key = neg_loss;
    v->in_heap = true;
}

void DecisionMesh::heap_pop(Vertex* v) {
    if (!v->in_heap) return;
    auto it = loss_heap.find(v->heap_key);
    if (it != loss_heap.end()) {
        it->second.erase(v);
        if (it->second.empty()) {
            loss_heap.erase(it);
        }
    }
    v->in_heap = false;
}

std::pair<Vertex*, double> DecisionMesh::heap_peek() const {
    if (loss_heap.empty()) return {nullptr, 0.0};
    auto it = loss_heap.begin();  // smallest key = most negative = best
    Vertex* v = *it->second.begin();
    return {v, it->first};
}

Face* DecisionMesh::random_face() {
    std::vector<Face*> faces(active_faces.begin(), active_faces.end());
    if (faces.empty()) return nullptr;

    std::vector<double> weights(faces.size());
    for (size_t i = 0; i < faces.size(); ++i) {
        double a = faces[i]->area();
        int count = faces[i]->n_covered;
        weights[i] = (a > 0 && std::isfinite(a)) ? a * count : 0.0;
    }

    double total = 0;
    for (double w : weights) total += w;

    if (total <= 0) {
        std::uniform_int_distribution<int> dist(0, (int)faces.size() - 1);
        return faces[dist(rng)];
    }

    std::uniform_real_distribution<double> dist(0.0, total);
    double r = dist(rng);
    double cum = 0;
    for (size_t i = 0; i < faces.size(); ++i) {
        cum += weights[i];
        if (r <= cum) return faces[i];
    }
    return faces.back();
}

void DecisionMesh::update_best_vertex(double) {
    if (use_eb) {
        maybe_recompute_tau_sq();
    }

    auto [best, _] = heap_peek();
    if (!best) return;

    if (best->active) {
        best->update_height();
    } else {
        best->activate(true);
    }
}

TimingRecord DecisionMesh::update_best_vertex_timed(double /*random_prob*/, int iteration) {
    using clock = std::chrono::high_resolution_clock;
    TimingRecord rec{};
    rec.iteration = iteration;

    if (use_eb) {
        maybe_recompute_tau_sq();
    }

    // Phase 1: Find best vertex (heap peek is O(1))
    auto t0 = clock::now();
    auto [best, _] = heap_peek();
    auto t1 = clock::now();
    rec.find_best_us = std::chrono::duration<double, std::micro>(t1 - t0).count();

    if (!best) {
        rec.total_us = rec.find_best_us;
        rec.active_faces = (int)active_faces.size();
        rec.n_vertices = (int)vertices.size();
        return rec;
    }

    // Phase 2: Split/activate — reset cascade counters
    g_counters = CascadeCounters{};
    auto t2 = clock::now();
    if (best->active) {
        best->height = best->new_height;
        if (best->parent_edge != nullptr) {
            best->delta_pooled = best->height - best->mu_lin();
        }
        auto t2b = clock::now();
        std::vector<Vertex*> to_update(best->affected_vertices.begin(),
                                       best->affected_vertices.end());
        if (use_eb) {
            to_update.insert(to_update.end(), best->prior_children.begin(),
                             best->prior_children.end());
            std::sort(to_update.begin(), to_update.end());
            to_update.erase(std::unique(to_update.begin(), to_update.end()),
                            to_update.end());
        }
        for (Vertex* v : to_update) {
            v->update_info();
        }
        best->loss_reduction = 0;
        heap_set(best, 0);
        auto t3 = clock::now();
        rec.split_us = std::chrono::duration<double, std::micro>(t2b - t2).count();
        rec.update_info_us = std::chrono::duration<double, std::micro>(t3 - t2b).count();
    } else {
        // activate: set height, split parent edge (creates new faces), then update_info
        best->active = true;
        best->free_coefficient = true;
        best->height = best->new_height;

        // Set EB posterior fields before split
        if (use_eb && best->parent_edge != nullptr) {
            best->delta_pooled = best->height - best->mu_lin();
            double inv_sigma = (best->sigma_sq < 1e300) ? (1.0 / best->sigma_sq) : 0.0;
            double xTx_reg = inv_sigma + best->lambda_v;
            best->sigma_pooled = (xTx_reg > 0) ? (1.0 / xTx_reg) : 1e300;
        }

        best->parent_edge->split();
        auto t2b = clock::now();
        std::vector<Vertex*> to_update(best->affected_vertices.begin(),
                                       best->affected_vertices.end());
        if (use_eb) {
            to_update.insert(to_update.end(), best->prior_children.begin(),
                             best->prior_children.end());
            std::sort(to_update.begin(), to_update.end());
            to_update.erase(std::unique(to_update.begin(), to_update.end()),
                            to_update.end());
        }
        for (Vertex* v : to_update) {
            v->update_info();
        }
        best->loss_reduction = 0;
        heap_set(best, 0);
        auto t3 = clock::now();
        rec.split_us = std::chrono::duration<double, std::micro>(t2b - t2).count();
        rec.update_info_us = std::chrono::duration<double, std::micro>(t3 - t2b).count();
    }

    auto tend = clock::now();
    rec.total_us = std::chrono::duration<double, std::micro>(tend - t0).count();
    rec.active_faces = (int)active_faces.size();
    rec.n_vertices = (int)vertices.size();
    rec.cascade_activations = g_counters.cascade_activations;
    rec.faces_created = g_counters.faces_created;
    rec.edges_created = g_counters.edges_created;
    rec.add_face_calls = g_counters.add_face_calls;
    rec.update_info_calls = g_counters.update_info_calls;
    return rec;
}

// --- Empirical Bayes ---

void DecisionMesh::recompute_tau_sq() {
    // Estimate slab variance from discoveries only. Coarse pre-refinement
    // and NVB completion vertices are structural, not selected slab draws.
    std::map<int, std::vector<Vertex*>> by_depth;
    std::set<int> target_depths;
    for (Vertex* v : vertices) {
        if (v->parent_edge == nullptr) continue;
        if (v->sigma_sq < 1e300 && v->sigma_sq > 0) target_depths.insert(v->depth);
        if (!v->active || !v->gate_admitted || v->dormant_coefficient ||
            v->retired_coefficient) continue;
        if (v->sigma_sq >= 1e300 || v->sigma_sq <= 0) continue;
        // DMESH_TAU_EXCLUDE_CONTAINED=1: support-starved vertices are frozen
        // at their prior by update_height (containment) yet still vote here
        // with stale admission-time sigma and noise surpluses, deflating tau
        // where their deltas are contained and inflating it where they are
        // wild. Excluded, they stay fully pooled without a hyperparameter vote.
        static const bool kExcludeContained =
            std::getenv("DMESH_TAU_EXCLUDE_CONTAINED") != nullptr;
        if (kExcludeContained &&
            (v->current_fit_points < Vertex::kMinimumFitPoints ||
             v->current_keff < Vertex::minimum_effective_support())) continue;
        by_depth[v->depth].push_back(v);
    }

    // Rebuild the map rather than silently carrying a depth estimate from a
    // different discovery set.  Unrepresented depths borrow below.
    tau_sq.clear();
    tau_sq_vertex_counts.clear();

    // DMESH_INLOOP_EB: use the principled-EB estimator (selection-corrected
    // marginal maximum likelihood + hierarchical borrow toward a global fit)
    // for the IN-LOOP prior, so it governs admission (the candidate exact gain
    // is scored against prior precision 1/tau) and growth, not just the final
    // reshrink pass. Zero-centered marginal N(0, tau + s) on the shrunk
    // surplus delta_pooled, truncated to the admission region |delta| > a_v.
    static const bool INLOOP_EB = std::getenv("DMESH_INLOOP_EB") != nullptr;
    if (INLOOP_EB) {
        auto mle_tau = [&](const std::vector<Vertex*>& vs) {
            if (vs.size() < 3) return 0.0;
            auto logL = [&](double t2) {
                double s = 0.0;
                for (Vertex* v : vs) {
                    double T = t2 + v->sigma_sq;
                    if (T <= 0.0) return -1e300;
                    double q = std::erfc((v->admit_a > 0.0 ? v->admit_a : 0.0) /
                                         std::sqrt(2.0 * T));
                    if (q < 1e-300) q = 1e-300;
                    s += -0.5 * std::log(T) -
                         0.5 * v->delta_pooled * v->delta_pooled / T - std::log(q);
                }
                return s;
            };
            double l0 = logL(0.0), hi = 1.0;
            while (hi < 1e8 && logL(2.0 * hi) > logL(hi)) hi *= 2.0;
            double a = 0.0, b = hi;
            const double gr = 0.6180339887498949;
            double c = b - gr * (b - a), d2 = a + gr * (b - a);
            double fc = logL(c), fd = logL(d2);
            for (int it = 0; it < 70; ++it) {
                if (fc < fd) { a = c; c = d2; fc = fd; d2 = a + gr * (b - a); fd = logL(d2); }
                else         { b = d2; d2 = c; fd = fc; c = b - gr * (b - a); fc = logL(c); }
            }
            double t = 0.5 * (a + b);
            return logL(t) > l0 ? t : 0.0;
        };
        static const double BORROW = [] {
            const char* e = std::getenv("DMESH_EB_BORROW");
            return e ? std::atof(e) : 8.0;
        }();
        std::vector<Vertex*> all;
        for (auto& [d, verts] : by_depth) {
            tau_sq_vertex_counts[d] = (int)verts.size();
            for (Vertex* v : verts) all.push_back(v);
        }
        const double tau_g = mle_tau(all);
        for (auto& [d, verts] : by_depth) {
            if ((int)verts.size() < 3) continue;
            const double td = mle_tau(verts);
            const double m = (double)verts.size();
            const double t = (m * td + BORROW * tau_g) / (m + BORROW);
            tau_sq[d] = std::max(t, tau_floor_sq);
        }
    } else {
    for (auto& [d, verts] : by_depth) {
        int m = (int)verts.size();
        tau_sq_vertex_counts[d] = m;

        if (m < 3) continue;

        double sum_w = 0.0;
        double sum_wd = 0.0;
        for (Vertex* v : verts) {
            double w = 1.0 / v->sigma_sq;
            sum_w += w;
            sum_wd += v->delta_pooled * w;
        }
        if (sum_w <= 0) continue;

        double mu_d = sum_wd / sum_w;

        double chi_sq = 0.0;
        for (Vertex* v : verts) {
            double diff = v->delta_pooled - mu_d;
            chi_sq += diff * diff / v->sigma_sq;
        }
        // Truncation-corrected DL: admitted delta-hats are selected for
        // |delta_hat| > a_v, so E[(delta-mu)^2/s | sel] = kappa_v (tau^2+s_v)/s_v
        // with kappa(alpha) = 1 + alpha phi(alpha)/(1-Phi(alpha)),
        // alpha = a_v/sqrt(tau^2+s_v). Solve chi^2 = sum kappa (tau^2+s)/s - 1
        // by bisection; reduces to plain DL when all a_v = 0. DMESH_NOTRUNC
        // restores the uncorrected moment.
        double tau;
        static const bool NOTRUNC = std::getenv("DMESH_TRUNC") == nullptr;   // correction OPT-IN: see writeup
        bool any_a = false;
        for (Vertex* v : verts) if (v->admit_a > 0) { any_a = true; break; }
        if (NOTRUNC || !any_a) {
            tau = std::max(0.0, (chi_sq - (m - 1)) / sum_w);
        } else {
            auto kap = [](double alpha) {
                if (alpha <= 1e-12) return 1.0;
                double ph = std::exp(-0.5 * alpha * alpha) / 2.5066282746310002;
                double Q = 0.5 * std::erfc(alpha / 1.4142135623730951);
                if (Q < 1e-300) return alpha * alpha + 1.0;
                return 1.0 + alpha * ph / Q;
            };
            auto target = [&](double t2) {
                double acc = 0;
                for (Vertex* v : verts) {
                    double sv = v->sigma_sq;
                    double al = v->admit_a / std::sqrt(t2 + sv);
                    acc += kap(al) * (t2 + sv) / sv;
                }
                return acc - 1.0;
            };
            if (chi_sq <= target(0.0)) tau = 0.0;
            else {
                double lo = 0.0, hi = 1.0;
                while (target(hi) < chi_sq && hi < 1e8) hi *= 2;
                for (int it = 0; it < 60; ++it) {
                    double mid = 0.5 * (lo + hi);
                    if (target(mid) < chi_sq) lo = mid; else hi = mid;
                }
                tau = 0.5 * (lo + hi);
            }
        }
        tau_sq[d] = std::max(tau, tau_floor_sq);
    }
    }

    // For depths with fewer than three discoveries, borrow from the nearest
    // depth with an estimable slab variance.  Before the first discoveries,
    // use a fixed weak prior rather than estimating the slab from structural
    // vertices.  This keeps coarse scales regularized without contaminating
    // the admitted-only hyperparameter estimate.
    std::vector<int> depths_with_tau;
    for (auto& [d, _] : tau_sq) depths_with_tau.push_back(d);
    std::sort(depths_with_tau.begin(), depths_with_tau.end());

    for (int d : target_depths) {
        if (tau_sq.find(d) != tau_sq.end()) continue;
        if (depths_with_tau.empty()) {
            tau_sq[d] = std::max(tau_default_sq, tau_floor_sq);
        } else {
            int nearest = depths_with_tau[0];
            for (int d2 : depths_with_tau) {
                if (std::abs(d2 - d) < std::abs(nearest - d)) nearest = d2;
            }
            tau_sq[d] = tau_sq[nearest];
        }
    }

    steps_since_tau_recompute = 0;
}

void DecisionMesh::maybe_recompute_tau_sq() {
    steps_since_tau_recompute++;

    if (steps_since_tau_recompute >= tau_sq_recompute_interval) {
        recompute_tau_sq();
        return;
    }

    // Check if any depth has grown by 50%
    std::map<int, int> by_depth;
    for (Vertex* v : vertices) {
        if (v->parent_edge != nullptr && v->sigma_sq < 1e300 &&
            v->active && v->gate_admitted && !v->dormant_coefficient &&
            !v->retired_coefficient) {
            by_depth[v->depth]++;
        }
    }

    for (auto& [d, count] : by_depth) {
        auto it = tau_sq_vertex_counts.find(d);
        int old_count = (it != tau_sq_vertex_counts.end()) ? it->second : 0;
        if (old_count > 0 && count >= (int)(old_count * 1.5)) {
            recompute_tau_sq();
            return;
        }
    }
}

double DecisionMesh::tau_for_depth(int depth) const {
    auto it = tau_sq.find(depth);
    if (it != tau_sq.end()) return std::max(it->second, tau_floor_sq);
    if (tau_sq.empty()) return std::max(tau_default_sq, tau_floor_sq);

    auto hi = tau_sq.lower_bound(depth);
    if (hi == tau_sq.begin()) return std::max(hi->second, tau_floor_sq);
    if (hi == tau_sq.end()) return std::max(std::prev(hi)->second, tau_floor_sq);
    auto lo = std::prev(hi);
    return std::max((depth - lo->first <= hi->first - depth) ? lo->second : hi->second,
                    tau_floor_sq);
}

// --- Prediction: walk tree to leaf face, then affine interpolate ---

double DecisionMesh::predict(double px, double py) const {
    // Walk the binary tree from root to leaf: O(depth) = O(log faces)
    const TreeNode* node = root;
    while (node->has_split) {
        double val = node->split_normal[0] * px + node->split_normal[1] * py;
        node = (val >= node->split_intercept) ? node->p : node->m;
    }

    // node->face is the leaf face containing (px, py)
    Face* f = node->face;
    if (!f) return 0.0;

    Vertex* v0 = f->vertices[0];
    Vertex* v1 = f->vertices[1];
    Vertex* v2 = f->vertices[2];
    double h0 = v0->height, h1 = v1->height, h2 = v2->height;

    // Compute barycentric coordinates for the query point
    double ax = v1->x - v0->x, ay = v1->y - v0->y;
    double bx = v2->x - v0->x, by = v2->y - v0->y;
    double det = ax * by - ay * bx;

    if (std::abs(det) < 1e-14) {
        // Degenerate triangle: return average height
        return (h0 + h1 + h2) / 3.0;
    }

    double inv = 1.0 / det;
    double dpx = px - v0->x, dpy = py - v0->y;
    double u = (dpx * by - dpy * bx) * inv;  // weight for v1
    double v = (ax * dpy - ay * dpx) * inv;   // weight for v2
    double w = 1.0 - u - v;                   // weight for v0

    // Clamp barycentric coords to [0,1] to avoid extrapolation on thin triangles
    // when the query point falls slightly outside due to numerical issues
    double cw = std::max(0.0, std::min(1.0, w));
    double cu = std::max(0.0, std::min(1.0, u));
    double cv = std::max(0.0, std::min(1.0, v));
    double total = cw + cu + cv;
    if (total > 0) {
        cw /= total; cu /= total; cv /= total;
    } else {
        cw = cu = cv = 1.0 / 3.0;
    }

    return cw * h0 + cu * h1 + cv * h2;
}


double DecisionMesh::predict_index(int i) const {
    if (i < 0 || i >= n_points) return 0.0;
    Face* f = point_face[i];
    if (!f) return predict(get_x(i, 0), get_x(i, 1));
    const auto& b = point_bary[i];
    return (double)b[0] * f->vertices[0]->height
         + (double)b[1] * f->vertices[1]->height
         + (double)b[2] * f->vertices[2]->height;
}

std::vector<double> DecisionMesh::predict_batch(const std::vector<double>& px,
                                                 const std::vector<double>& py) const {
    int n = (int)px.size();
    std::vector<double> result(n);
    for (int i = 0; i < n; ++i) {
        result[i] = predict(px[i], py[i]);
    }
    return result;
}

void DecisionMesh::create_outer_vertices() {
    xmin = X[0]; xmax = X[0];
    ymin = X[1]; ymax = X[1];
    for (int i = 0; i < n_points; ++i) {
        double xi = get_x(i, 0), yi = get_x(i, 1);
        if (xi < xmin) xmin = xi;
        if (xi > xmax) xmax = xi;
        if (yi < ymin) ymin = yi;
        if (yi > ymax) ymax = yi;
    }
    outer_vertices["bottom_left"]  = make_vertex(xmin, ymin, true);
    outer_vertices["top_left"]     = make_vertex(xmin, ymax, true);
    outer_vertices["bottom_right"] = make_vertex(xmax, ymin, true);
    outer_vertices["top_right"]    = make_vertex(xmax, ymax, true);
}

void DecisionMesh::create_outer_edges() {
    outer_edges["left"]   = make_edge(outer_vertices["top_left"],    outer_vertices["bottom_left"], false);
    outer_edges["bottom"] = make_edge(outer_vertices["bottom_right"],outer_vertices["bottom_left"], false);
    outer_edges["right"]  = make_edge(outer_vertices["top_right"],   outer_vertices["bottom_right"], false);
    outer_edges["top"]    = make_edge(outer_vertices["top_right"],   outer_vertices["top_left"], false);
}

void DecisionMesh::write_svg(const std::string& filename, double w, double h) const {
    // Find height range for coloring
    double hmin = 1e300, hmax = -1e300;
    for (Vertex* v : vertices) {
        if (v->height < hmin) hmin = v->height;
        if (v->height > hmax) hmax = v->height;
    }
    if (hmax - hmin < 1e-12) { hmin -= 1; hmax += 1; }

    double padx = 20, pady = 20;
    double sx = (w - 2 * padx) / (xmax - xmin);
    double sy = (h - 2 * pady) / (ymax - ymin);
    double scale = std::min(sx, sy);

    auto tx = [&](double x) { return padx + (x - xmin) * scale; };
    auto ty = [&](double y) { return h - pady - (y - ymin) * scale; };

    // Simple blue-white-red colormap
    auto color = [&](double val) -> std::string {
        double t = (val - hmin) / (hmax - hmin);
        t = std::max(0.0, std::min(1.0, t));
        int r, g, b;
        if (t < 0.5) {
            double s = t * 2;
            r = (int)(68 + s * (255 - 68));
            g = (int)(1 + s * (255 - 1));
            b = (int)(84 + s * (255 - 84));
        } else {
            double s = (t - 0.5) * 2;
            r = (int)(255 - s * (255 - 253));
            g = (int)(255 - s * (255 - 231));
            b = (int)(255 - s * (255 - 37));
        }
        char buf[16];
        snprintf(buf, sizeof(buf), "#%02x%02x%02x", r, g, b);
        return buf;
    };

    std::ofstream out(filename);
    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << w
        << "\" height=\"" << h << "\">\n";
    out << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    // Draw filled triangles
    for (Face* f : active_faces) {
        double avg_h = (f->vertices[0]->height + f->vertices[1]->height + f->vertices[2]->height) / 3.0;
        out << "<polygon points=\""
            << tx(f->vertices[0]->x) << "," << ty(f->vertices[0]->y) << " "
            << tx(f->vertices[1]->x) << "," << ty(f->vertices[1]->y) << " "
            << tx(f->vertices[2]->x) << "," << ty(f->vertices[2]->y) << "\" "
            << "fill=\"" << color(avg_h) << "\" stroke=\"black\" stroke-width=\"0.5\"/>\n";
    }

    out << "</svg>\n";
    out.close();
}

void DecisionMesh::write_mesh_csv(const std::string& prefix) const {
    // Assign stable IDs to vertices
    std::map<const Vertex*, int> vid;
    int next_id = 0;
    for (const Vertex* v : vertices) {
        vid[v] = next_id++;
    }

    // Write vertices: id, x, y, height, shrinkage, depth
    {
        std::ofstream out(prefix + "_vertices.csv");
        out << "id,x,y,height,shrinkage,depth\n";
        for (const Vertex* v : vertices) {
            // shrinkage = lambda_v / (lambda_v + xTx), where xTx = 1/sigma_sq
            double xTx = (v->sigma_sq > 0 && v->sigma_sq < 1e200) ? 1.0 / v->sigma_sq : 0.0;
            double shrinkage = (v->lambda_v + xTx > 0) ? v->lambda_v / (v->lambda_v + xTx) : 0.0;
            out << vid[v] << "," << v->x << "," << v->y << "," << v->height
                << "," << shrinkage << "," << v->depth << "\n";
        }
    }

    // Write triangles: v0, v1, v2 (vertex ids)
    {
        std::ofstream out(prefix + "_triangles.csv");
        out << "v0,v1,v2\n";
        for (const Face* f : active_faces) {
            out << vid[f->vertices[0]] << ","
                << vid[f->vertices[1]] << ","
                << vid[f->vertices[2]] << "\n";
        }
    }
}
