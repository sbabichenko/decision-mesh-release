#include <unordered_map>
#include "pipeline.h"

#include "diagnostics.h"
#include "edge.h"
#include "face.h"
#include "hierarchical_design.h"
#include "lfdr.h"
#include "mesh.h"
#include "observation.h"
#include "output.h"
#include "units.h"
#include "vertex.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <stdexcept>
#include <tuple>
#include <vector>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>

namespace {

double now_seconds() {
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
}

struct PipelineTiming {
    double height_sweeps_seconds = 0.0;
    double candidate_scoring_seconds = 0.0;
};

// Diagnostic state for tracing clamp-induced trajectory bifurcations.
static int g_diag_stage = -1;
static int g_diag_round = -1;
static int g_diag_sweep = -1;
static const char* g_diag_phase = "init";
static long long g_diag_event = 0;

struct B1GainDiagnostics {
    long long attempted = 0;
    long long scored = 0;
    long long unscoreable = 0;
    long long unscoreable_no_training_column = 0;
    long long unscoreable_tempered_out = 0;
    long long negative = 0;
    long long grid_checked = 0;
    long long grid_failed = 0;
    double max_grid_excess = 0.0;
    double nearest_regime_distance = std::numeric_limits<double>::infinity();
    int nearest_stage = -1;
    int nearest_round = -1;
    int nearest_id = -1;
    double nearest_x = 0.0;
    double nearest_y = 0.0;
    double nearest_gain = 0.0;
    double nearest_null = 0.0;
    double nearest_optimum = 0.0;
    long long live_reprofile_attempted = 0;
    long long live_reprofile_scored = 0;
    long long live_reprofile_rejected = 0;
    long long commit_validated = 0;
    long long commit_failed = 0;
    double max_commit_gain_error = 0.0;
    double max_commit_optimum_error_units = 0.0;
};

static B1GainDiagnostics g_b1_diag;

static std::ofstream& b1_gain_stream() {
    static std::ofstream stream;
    static bool initialized = false;
    if (!initialized) {
        initialized = true;
        // B2 name is preferred; retain the B1 name so closure scripts remain
        // backwards compatible.
        const char* path = std::getenv("DMESH_B2_SCORE_DUMP");
        if (!(path && *path)) path = std::getenv("DMESH_B1_GAIN_DUMP");
        if (path && *path) {
            stream.open(path);
            stream << std::setprecision(17);
            stream << "stage,round,id,x,y,depth,active_conformity,legacy_z,"
                      "signed_root_gain,gate_score,legacy_abs_rank_input,exact_abs_rank_input,"
                      "legacy_center_units,legacy_displacement_units,legacy_sd_units,"
                      "null_units,prior_units,optimum_units,exact_displacement_units,"
                      "gain_nats,data_gain_nats,penalty_gain_nats,raw_gain_nats,"
                      "score_null,information_null,data_information_optimum,lambda_logit,"
                      "structural_points,training_structural_points,positive_information_points,"
                      "fit_points,keff,min_eta_null,max_eta_null,min_p_null,max_p_null,"
                      "phi_disp,gate_var,grid_excess,b3_null_scale,b3_calibrated_score,"
                      "b3_stage0_legacy,hist_p,hist_z,hist_score,hist_null_mean,"
                      "hist_null_sd,hist_tail_keff,hist_rows,hist_explicit_rows\n";
        }
    }
    return stream;
}

static double height_limit_units() {
    const char* e = std::getenv("DMESH_HEIGHT_LIMIT_LOGIT");
    if (!e) return std::numeric_limits<double>::infinity();
    const double v = std::atof(e);
    return (v > 0 && std::isfinite(v)) ? v * kHeightScale
                                      : std::numeric_limits<double>::infinity();
}

static std::ofstream& diag_stream() {
    static std::ofstream f;
    static bool initialized = false;
    if (!initialized) {
        initialized = true;
        const char* path = std::getenv("DMESH_BRANCH_DIAG");
        if (path) {
            f.open(path);
            f << std::setprecision(17);
            f << "event,stage,round,sweep,phase,kind,id,x,y,depth,gate,old_h,proposal_h,post_h,mu_lin,mu_prior,raw_beta,xTx,xTr,n_points,fit_points,lambda,sigma_sq,loss,parent0,parent1,p0h,p1h,clamped\n";
        }
    }
    return f;
}

static bool in_diag_region(const Vertex& v) {
    // WALA roughly 10--50 and WAC 7.7--8.9 in normalized real-data coordinates.
    static const double* box = [] {
        static double b[4] = {0.025, 0.15, 0.80, 0.99};
        if (const char* e = std::getenv("DMESH_DIAG_REGION"))
            std::sscanf(e, "%lf,%lf,%lf,%lf", &b[0], &b[1], &b[2], &b[3]);
        return b;
    }();
    return v.x >= box[0] && v.x <= box[1] && v.y >= box[2] && v.y <= box[3];
}

static void log_vertex_event(const char* kind, Vertex& v, double old_h,
                             double proposal_h, double post_h, bool clamped) {
    auto& f = diag_stream();
    if (!f || (!in_diag_region(v) && !clamped)) return;
    auto lr = v.loc_regress(false);
    const double raw_beta = lr.xTx > 0 ? lr.xTr / lr.xTx : v.height;
    int p0 = -1, p1 = -1; double p0h = 0, p1h = 0;
    if (v.parent_edge) {
        p0 = v.parent_edge->vertex0->_id; p1 = v.parent_edge->vertex1->_id;
        p0h = v.parent_edge->vertex0->height; p1h = v.parent_edge->vertex1->height;
    }
    f << (++g_diag_event) << ',' << g_diag_stage << ',' << g_diag_round << ','
      << g_diag_sweep << ',' << g_diag_phase << ',' << kind << ',' << v._id << ','
      << v.x << ',' << v.y << ',' << v.depth << ',' << (v.gate_admitted?1:0) << ','
      << old_h << ',' << proposal_h << ',' << post_h << ',' << v.mu_lin() << ','
      << v.mu_prior() << ',' << raw_beta << ',' << lr.xTx << ',' << lr.xTr << ','
      << lr.n_points << ',' << lr.fit_points << ',' << v.lambda_v << ',' << v.sigma_sq << ','
      << v.loss_reduction << ',' << p0 << ',' << p1 << ',' << p0h << ',' << p1h
      << ',' << (clamped?1:0) << '\n';
}

bool proposal_is_valid(double proposed_height) {
    const double limit = height_limit_units();
    return std::isfinite(proposed_height) && std::fabs(proposed_height) <= limit;
}

// Validate before mutating any height. This avoids the old failure mode where
// a rejected state had already propagated stale proposals to its neighbors.
bool prepare_vertex_proposal(Vertex& vertex, double& proposal) {
    if (proposal_is_valid(proposal)) return true;
    const double rejected = proposal;
    proposal = vertex.parent_edge ? vertex.mu_lin() : vertex.height;
    vertex.new_height = proposal;
    log_vertex_event("reject_precommit", vertex, vertex.height, rejected,
                     vertex.height, true);
    return false;
}



// ---------------------------------------------------------------------------
// Phase A: the exact stage solver (DMESH_PL_SOLVER=1).
// Coordinate ascent on the exact tempered binomial likelihood plus the EB
// prior (mu_prior, lambda = 1/tau(depth)), on hierarchical columns, with a
// live per-pool eta maintained incrementally. Every accepted step increases
// a bounded-above objective: eruption is impossible by construction, and
// one-sided stars are handled natively by the prior term (MAP is the solver,
// not a special pass). The gate is untouched in this phase: scores and
// admissions still come from the quadratic refresh snapshot.
namespace pl_solver {

struct StaticRow {
    int index;
    double x;
    double events;
    double exposure;
    double fixed_temper;  // NaN means derive from the current working weight.
};

struct CoordinateRows {
    Vertex* vertex;
    const hierarchical_design::SparseColumn* column;
    std::vector<StaticRow> rows;
};

struct VisitRow {
    double eta_other;
    double x;
    double events;
    double exposure;
    double temper;
};

struct CoupledPenalty {
    double residual_logit = 0.0;
    double lambda_logit = 0.0;
};

struct CoupledPenaltyTerm {
    int penalty_index = -1;
    double derivative = 0.0;
};

void run(DecisionMesh& mesh, double tol_units, int max_sweeps) {
    if (mesh.obs_store == nullptr) return;
    const auto& obs = *mesh.obs_store;

    hierarchical_design::BuildOptions design_options;
    // The exact objective still decides which rows carry likelihood weight;
    // build all geometric rows so live predictions remain synchronized for
    // both training and held-out observations.
    design_options.training_only = false;
    auto design = hierarchical_design::build_active_design(mesh, design_options);
    auto& vs = design.active_vertices;
    auto& cols = design.columns;
    std::vector<double> pred = std::move(design.predictions);

    std::vector<Vertex*> free_vs;
    for (Vertex* v : vs) {
        const auto found = cols.find(v);
        if (v->free_coefficient && v->parent_edge != nullptr &&
            found != cols.end()) {
            free_vs.push_back(v);
        }
    }
    std::sort(free_vs.begin(), free_vs.end(), [](const Vertex* a, const Vertex* b) {
        return a->depth != b->depth ? a->depth < b->depth : a->_id < b->_id;
    });

    // The sparse column, counts, and frozen likelihood temper are constant for
    // an entire exact solve.  Materialize them once instead of rebuilding five
    // temporary arrays at every coordinate visit.  eta_other remains visit-
    // specific because earlier Gauss--Seidel coordinates update pred live.
    std::vector<CoordinateRows> coordinates;
    coordinates.reserve(free_vs.size());
    std::size_t maximum_rows = 0;
    for (Vertex* v : free_vs) {
        auto found = cols.find(v);
        if (found == cols.end()) continue;
        CoordinateRows coordinate{v, &found->second, {}};
        coordinate.rows.reserve(found->second.size());
        for (const auto& entry : found->second) {
            const Observation& o = obs[entry.first];
            double fixed_temper = std::numeric_limits<double>::quiet_NaN();
            if (entry.first < static_cast<int>(mesh.likelihood_temper.size())) {
                fixed_temper = mesh.likelihood_temper[entry.first];
                if (fixed_temper < 1e-8) continue;
            }
            coordinate.rows.push_back({entry.first, entry.second,
                                       static_cast<double>(o.events),
                                       static_cast<double>(o.exposure),
                                       fixed_temper});
        }
        maximum_rows = std::max(maximum_rows, coordinate.rows.size());
        coordinates.push_back(std::move(coordinate));
    }

    std::vector<VisitRow> visit_rows;
    visit_rows.reserve(maximum_rows);

    // The declared hierarchical Gaussian prior is a joint penalty on every
    // free surplus. Moving an ancestor changes the prior residuals of free
    // descendants whose parents inherit that ancestor. The old solver priced
    // only the visited coordinate's own penalty, so it was not coordinate
    // ascent on the stated penalized objective. Build the exact sparse penalty
    // derivatives once; the residuals are updated incrementally with each
    // Gauss--Seidel move, just like the likelihood predictions.
    // Linear and quadratic midpoint priors are both fixed linear stencils in
    // the active vertex heights, so both admit the exact sparse coupled
    // penalty construction. Ghost-boundary targets are location-dependent
    // searches and retain their specialized diagonal treatment.
    const bool coupled_stencil_prior =
        std::getenv("DMESH_GHOST_BC") == nullptr;
    std::vector<CoupledPenalty> coupled_penalties;
    std::unordered_map<Vertex*, std::vector<CoupledPenaltyTerm>>
        penalty_terms_by_coordinate;
    if (coupled_stencil_prior) {
        for (Vertex* w : vs) {
            if (!w->free_coefficient || w->retired_coefficient ||
                w->parent_edge == nullptr) continue;
            const auto [lambda_units, prior_units] = w->compute_eb_params(1.0);
            const double lambda_logit = std::max(lambda_units, 0.0) * kVarianceScale;
            if (!(lambda_logit > 0.0)) continue;
            const int penalty_index = static_cast<int>(coupled_penalties.size());
            coupled_penalties.push_back({
                (w->height - prior_units) / kHeightScale,
                lambda_logit
            });

            std::unordered_map<Vertex*, double> derivative;
            derivative[w] += 1.0;
            for (const auto& [stencil_vertex, stencil_weight] :
                 w->prior_stencil()) {
                const auto found = design.parameter_map.find(stencil_vertex);
                if (found == design.parameter_map.end()) continue;
                for (const auto& [coefficient, weight] : found->second)
                    derivative[coefficient] -= stencil_weight * weight;
            }
            for (const auto& [coefficient, value] : derivative) {
                if (std::fabs(value) <= 1e-14) continue;
                penalty_terms_by_coordinate[coefficient].push_back(
                    {penalty_index, value});
            }
        }
    }

    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        double max_move = 0.0;
        for (CoordinateRows& coordinate : coordinates) {
            Vertex* v = coordinate.vertex;
            static const int TRACE_VERTEX_ID = [] {
                const char* e = std::getenv("DMESH_TRACE_VERTEX_ID");
                return e ? std::atoi(e) : -1;
            }();
            auto [lv, mu_v] = v->compute_eb_params(1.0);
            const double lam_logit = std::max(lv, 0.0) * kVarianceScale;
            const double mu_l = mu_v / kHeightScale;
            double b = v->height / kHeightScale;
            const double b_entry = b;

            visit_rows.clear();
            for (const StaticRow& row : coordinate.rows) {
                const Observation& o = obs[row.index];
                const double eta_cur = o.offset + pred[row.index] / kHeightScale;
                const double p_cur = 1.0 / (1.0 + std::exp(-eta_cur));
                const double wb = row.exposure * p_cur * (1.0 - p_cur);
                if (wb <= 0.0) continue;
                const double c = std::isfinite(row.fixed_temper)
                    ? row.fixed_temper
                    : std::clamp(mesh.weights[row.index] * kVarianceScale / wb,
                                 0.0, 1.0);
                if (c < 1e-8) continue;
                visit_rows.push_back({eta_cur - row.x * b, row.x,
                                      row.events, row.exposure, c});
            }
            if (v->_id == TRACE_VERTEX_ID) {
                std::fprintf(stderr,
                    "[pl trace] sweep %d id %d old %.12g prior %.12g lambda %.12g col %zu usable %zu\n",
                    sweep, v->_id, v->height, mu_v, lam_logit,
                    coordinate.column->size(), visit_rows.size());
            }
            const auto penalty_found = penalty_terms_by_coordinate.find(v);
            const bool have_coupled_penalty = coupled_stencil_prior &&
                penalty_found != penalty_terms_by_coordinate.end();
            if (visit_rows.empty() && !have_coupled_penalty && lam_logit <= 0.0)
                continue;
            for (int itn = 0; itn < 40; ++itn) {
                double g = 0.0, H = 0.0;
                if (have_coupled_penalty) {
                    const double shift = b - b_entry;
                    for (const CoupledPenaltyTerm& term : penalty_found->second) {
                        const CoupledPenalty& penalty =
                            coupled_penalties[term.penalty_index];
                        const double residual = penalty.residual_logit +
                            term.derivative * shift;
                        g -= penalty.lambda_logit * term.derivative * residual;
                        H -= penalty.lambda_logit * term.derivative * term.derivative;
                    }
                } else {
                    g = -lam_logit * (b - mu_l);
                    H = -lam_logit;
                }
                for (const VisitRow& row : visit_rows) {
                    const double e2 = row.eta_other + b * row.x;
                    const double p = 1.0 / (1.0 + std::exp(-e2));
                    g += row.temper * row.x * (row.events - row.exposure * p);
                    H -= row.temper * row.x * row.x * row.exposure * p * (1.0 - p);
                }
                if (H >= -1e-12) break;
                const double step = std::clamp(g / (-H), -1.5, 1.5);
                b += step;
                if (std::fabs(step) < 1e-7) break;
            }
            const double old_height = v->height;
            const double proposed_height = b * kHeightScale;
            if (!std::isfinite(proposed_height)) continue;
            v->height = proposed_height;
            v->new_height = v->height;
            v->clamp_height();
            if (v->_id == TRACE_VERTEX_ID) {
                std::fprintf(stderr,
                    "[pl trace] sweep %d id %d raw-proposal %.12g committed %.12g\n",
                    sweep, v->_id, proposed_height, v->height);
            }
            const double committed_delta = v->height - old_height;
            if (committed_delta == 0.0) continue;
            for (const auto& entry : *coordinate.column) {
                pred[entry.first] += entry.second * committed_delta;
            }
            if (have_coupled_penalty) {
                const double delta_logit = committed_delta / kHeightScale;
                for (const CoupledPenaltyTerm& term : penalty_found->second)
                    coupled_penalties[term.penalty_index].residual_logit +=
                        term.derivative * delta_logit;
            }
            max_move = std::max(max_move, std::fabs(committed_delta));
        }
        if (max_move < tol_units) break;
    }

    if (std::getenv("DMESH_COUPLED_FIXED_POINT_AUDIT")) {
        double max_abs_gradient = 0.0;
        double sum_squared_gradient = 0.0;
        double max_abs_interior_gradient = 0.0;
        double sum_squared_interior_gradient = 0.0;
        int audited = 0;
        int interior_audited = 0;
        int boundary_kkt_violations = 0;
        int max_id = -1;
        int max_interior_id = -1;
        for (const CoordinateRows& coordinate : coordinates) {
            Vertex* v = coordinate.vertex;
            double gradient = 0.0;
            const auto penalty_found = penalty_terms_by_coordinate.find(v);
            if (coupled_stencil_prior &&
                penalty_found != penalty_terms_by_coordinate.end()) {
                for (const CoupledPenaltyTerm& term : penalty_found->second) {
                    const CoupledPenalty& penalty =
                        coupled_penalties[term.penalty_index];
                    gradient -= penalty.lambda_logit * term.derivative *
                                penalty.residual_logit;
                }
            } else {
                const auto [lv, mu_v] = v->compute_eb_params(1.0);
                const double lambda = std::max(lv, 0.0) * kVarianceScale;
                gradient -= lambda *
                    ((v->height - mu_v) / kHeightScale);
            }
            for (const StaticRow& row : coordinate.rows) {
                const Observation& o = obs[row.index];
                const double eta = o.offset + pred[row.index] / kHeightScale;
                const double probability = 1.0 / (1.0 + std::exp(-eta));
                const double current_information =
                    row.exposure * probability * (1.0 - probability);
                if (!(current_information > 0.0)) continue;
                const double temper = std::isfinite(row.fixed_temper)
                    ? row.fixed_temper
                    : std::clamp(mesh.weights[row.index] * kVarianceScale /
                                 current_information, 0.0, 1.0);
                if (temper < 1e-8) continue;
                gradient += temper * row.x *
                    (row.events - row.exposure * probability);
            }
            ++audited;
            sum_squared_gradient += gradient * gradient;
            if (std::fabs(gradient) > max_abs_gradient) {
                max_abs_gradient = std::fabs(gradient);
                max_id = v->_id;
            }
            const double limit_units = 12.0 * kHeightScale;
            const bool at_lower = v->height <= -limit_units + 1e-7;
            const bool at_upper = v->height >= limit_units - 1e-7;
            if (at_lower || at_upper) {
                const bool kkt_ok = (at_lower && gradient <= 1e-6) ||
                                    (at_upper && gradient >= -1e-6);
                if (!kkt_ok) ++boundary_kkt_violations;
            } else {
                ++interior_audited;
                sum_squared_interior_gradient += gradient * gradient;
                if (std::fabs(gradient) > max_abs_interior_gradient) {
                    max_abs_interior_gradient = std::fabs(gradient);
                    max_interior_id = v->_id;
                }
            }
        }
        const double rms_gradient = audited > 0
            ? std::sqrt(sum_squared_gradient / audited) : 0.0;
        const double rms_interior_gradient = interior_audited > 0
            ? std::sqrt(sum_squared_interior_gradient / interior_audited) : 0.0;
        std::fprintf(stderr,
            "[coupled fixed-point audit] audited %d max/rms gradient "
            "%.9g/%.9g (max id %d); interior %d max/rms %.9g/%.9g "
            "(max id %d); boundary KKT violations %d\n",
            audited, max_abs_gradient, rms_gradient, max_id,
            interior_audited, max_abs_interior_gradient,
            rms_interior_gradient, max_interior_id,
            boundary_kkt_violations);
    }

    for (Vertex* v : vs) {
        if (v->parent_edge == nullptr && v->free_coefficient &&
            v->current_fit_points < Vertex::kMinimumFitPoints) {
            v->height = v->neighbor_mean_height();
            v->new_height = v->height;
        }
    }
    hierarchical_design::materialize_constraints(mesh);
}

}  // namespace pl_solver

// ---------------------------------------------------------------------------
// Fix 2: hierarchical design columns (DMESH_HIER_FIT=1).
// After the nodal sweeps, re-fit every free coefficient against its
// HIERARCHICAL column: the derivative of the fitted surface with respect to
// that coefficient, propagated through constrained and held descendants
// (h_child = 0.5 h_a + 0.5 h_b + const). A coefficient's identifiability then
// depends on the structure it introduced, not on how much of its territory
// remains unrefined: the shrinking-nodal-star orphaning cannot occur.
// Coordinate descent with lambda = 0 (matching the tree's admitted-coefficient
// EB) on the frozen working responses is monotone in weighted SSE.
// Constrained heights are materialized in birth order afterwards.
namespace hier_fit {

void run(DecisionMesh& mesh, double tol_units, int max_sweeps) {
    if (mesh.obs_store == nullptr) return;
    const auto& obs = *mesh.obs_store;

    std::vector<Vertex*> vs;
    for (Vertex* v : mesh.vertices) if (v->active) vs.push_back(v);
    std::sort(vs.begin(), vs.end(), [](const Vertex* a, const Vertex* b) {
        return a->_id < b->_id;   // birth order: parents precede children
    });

    // Param maps: vertex -> {free coefficient -> weight}.
    std::unordered_map<const Vertex*, std::unordered_map<Vertex*, double>> pmap;
    for (Vertex* v : vs) {
        auto& m = pmap[v];
        const bool starved_root =
            v->parent_edge == nullptr && v->current_fit_points < 10;
        if (starved_root) continue;              // held constant (fix 1b pass)
        if (v->parent_edge == nullptr) { m[v] = 1.0; continue; }
        // Held-free vertices MOVE WITH THEIR PARENTS (mu_prior + frozen
        // surplus), so their columns must propagate like constrained ones;
        // classifying them as free coefficients would let the fit optimize
        // under an assumption the materialization step then violates.
        static const bool HIER_HELD = []{
            const char* e = std::getenv("DMESH_HIER_HELD");
            return e != nullptr && std::atoi(e) != 0;
        }();
        const bool held_free = HIER_HELD && v->free_coefficient &&
            (v->current_fit_points < 10 ||
             v->current_keff < Vertex::minimum_effective_support());
        if (v->free_coefficient && !held_free) { m[v] = 1.0; continue; }
        for (const Vertex* p : { v->parent_edge->vertex0, v->parent_edge->vertex1 }) {
            auto it = pmap.find(p);
            if (it == pmap.end()) continue;
            for (auto& kv : it->second) m[kv.first] += 0.5 * kv.second;
        }
    }

    // Columns: free coefficient -> [(pool index, x weight)].
    std::unordered_map<Vertex*, std::vector<std::pair<int, double>>> cols;
    std::vector<double> pred(mesh.n_points, 0.0);
    for (int i = 0; i < mesh.n_points; ++i) {
        pred[i] = mesh.predict_index(i);
        Face* f = mesh.point_face[i];
        if (f == nullptr) continue;
        const Vertex* c0 = f->vertices[0];
        const Vertex* c1 = f->vertices[1];
        const Vertex* c2 = f->vertices[2];
        const double ax = c1->x - c0->x, ay = c1->y - c0->y;
        const double bx = c2->x - c0->x, by = c2->y - c0->y;
        const double den = ax * by - ay * bx;
        if (std::fabs(den) < 1e-15) continue;
        const Observation& o = obs[i];
        const double px = o.x - c0->x, py = o.y - c0->y;
        const double u = (px * by - py * bx) / den;
        const double w2 = (ax * py - ay * px) / den;
        const double bc[3] = { 1.0 - u - w2, u, w2 };
        const Vertex* corner[3] = { c0, c1, c2 };
        for (int s = 0; s < 3; ++s) {
            if (bc[s] <= 1e-9) continue;
            auto it = pmap.find(corner[s]);
            if (it == pmap.end()) continue;
            for (auto& kv : it->second) {
                if (mesh.weights[i] < 1e-9) continue;   // held-out half
                cols[kv.first].push_back({ i, bc[s] * kv.second });
            }
        }
    }
    // merge duplicate pool entries per column
    for (auto& kv : cols) {
        auto& c = kv.second;
        std::sort(c.begin(), c.end());
        std::vector<std::pair<int, double>> merged;
        for (auto& e : c) {
            if (!merged.empty() && merged.back().first == e.first) merged.back().second += e.second;
            else merged.push_back(e);
        }
        c.swap(merged);
    }

    // DMESH_HIER_MAP=1: weak coordinates (few column pools or collapsed
    // keff) are not held; they are fitted to their POSTERIOR MODE: exact
    // tempered binomial likelihood (one-sided evidence handled honestly)
    // plus a Gaussian prior centered at mu_prior(live parents) with sd
    // DMESH_MAP_SD logits (default 1.0, matching the tree's EB default).
    // The prior floor comes from the hierarchy; the data shifts it.
    static const bool HIER_MAP = []{
        // Default ON (adoption decision): weak coordinates take posterior
        // modes, not holds. Frozen holds were history-dependent values:
        // accidentally good where sweep history was kind, silently wrong
        // where it wasn't, and unauditable either way. DMESH_HIER_MAP=0
        // restores holds for A/B.
        const char* e = std::getenv("DMESH_HIER_MAP");
        return e == nullptr || std::atoi(e) != 0;
    }();
    static const double MAP_SD = []{
        const char* e = std::getenv("DMESH_MAP_SD");
        return e ? std::atof(e) : 1.0;
    }();
    // Coordinate descent, coarse to fine.
    std::vector<Vertex*> free_vs;
    std::unordered_map<Vertex*, bool> weak;
    for (Vertex* v : vs) {
        if (!v->free_coefficient) continue;
        auto it = cols.find(v);
        const std::size_t csz = (it == cols.end()) ? 0 : it->second.size();
        const bool w = csz < 10 ||
                       v->current_keff < Vertex::minimum_effective_support();
        if (w && !HIER_MAP) continue;            // legacy: weak stay held
        if (csz == 0 && !HIER_MAP) continue;
        // Weak ROOTS are excluded from MAP: mu_prior() is 0 for roots, which
        // is not a prior but an artifact. The starved-root neighbor-mean
        // pass (fix 1b) governs them.
        if (w && v->parent_edge == nullptr) continue;
        weak[v] = w;
        free_vs.push_back(v);
    }
    std::sort(free_vs.begin(), free_vs.end(), [](const Vertex* a, const Vertex* b) {
        return a->depth != b->depth ? a->depth < b->depth : a->_id < b->_id;
    });
    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        double max_move = 0.0;
        for (Vertex* v : free_vs) {
            double delta = 0.0;
            if (weak[v]) {
                // Exact 1-D MAP in logits: prior N(mu_prior, MAP_SD^2).
                const auto& obs = *mesh.obs_store;
                static const bool MAP_QUAD = []{
                    const char* e = std::getenv("DMESH_MAP_QUAD");
                    return e != nullptr && std::atoi(e) != 0;
                }();
                const double mu_l =
                    (MAP_QUAD ? v->mu_quad() : v->mu_prior()) / kHeightScale;
                const double lam = 1.0 / (MAP_SD * MAP_SD);
                double b = v->height / kHeightScale;
                // Hoist per-pool constants at ENTRY, matching exact_score:
                // the temper c must be FROZEN for the Newton, otherwise the
                // law's information discount dissolves as the vertex dives
                // (p -> 0 makes w_bin -> 0 and the discount saturate to 1),
                // and the objective is no longer the tempered likelihood.
                std::vector<double> q_eo, q_x, q_k, q_n, q_c;
                for (auto& e : cols[v]) {
                    const Observation& o = obs[e.first];
                    const double eta_cur = o.offset + pred[e.first] / kHeightScale;
                    const double p_cur = 1.0 / (1.0 + std::exp(-eta_cur));
                    const double wbin = o.exposure * p_cur * (1.0 - p_cur);
                    if (wbin <= 0.0) continue;
                    const double c_frozen =
                        (e.first < (int)mesh.likelihood_temper.size())
                        ? mesh.likelihood_temper[e.first]
                        : std::clamp(mesh.weights[e.first] * kVarianceScale / wbin,
                                     0.0, 1.0);
                    if (c_frozen < 1e-8) continue;
                    q_eo.push_back(eta_cur - e.second * b);
                    q_x.push_back(e.second);
                    q_k.push_back((double)o.events);
                    q_n.push_back((double)o.exposure);
                    q_c.push_back(c_frozen);
                }
                for (int it2 = 0; it2 < 40; ++it2) {
                    double gsum = -lam * (b - mu_l), H = -lam;
                    for (std::size_t q = 0; q < q_x.size(); ++q) {
                        const double e2 = q_eo[q] + b * q_x[q];
                        const double p = 1.0 / (1.0 + std::exp(-e2));
                        gsum += q_c[q] * q_x[q] * (q_k[q] - q_n[q] * p);
                        H -= q_c[q] * q_x[q] * q_x[q] * q_n[q] * p * (1.0 - p);
                    }
                    if (H >= -1e-12) break;
                    const double step = std::clamp(gsum / (-H), -1.5, 1.5);
                    b += step;
                    if (std::fabs(step) < 1e-7) break;
                }
                delta = b * kHeightScale - v->height;
            } else {
                double num = 0.0, den2 = 0.0;
                for (auto& e : cols[v]) {
                    const double r = mesh.values[e.first] - pred[e.first];
                    num += mesh.weights[e.first] * e.second * r;
                    den2 += mesh.weights[e.first] * e.second * e.second;
                }
                if (den2 <= 1e-12) continue;
                delta = num / den2;
            }
            if (!std::isfinite(delta) || delta == 0.0) continue;
            v->height += delta;
            v->new_height = v->height;
            v->clamp_height();
            for (auto& e : cols[v]) pred[e.first] += e.second * delta;
            max_move = std::max(max_move, std::fabs(delta));
        }
        if (max_move < tol_units) break;
    }

    // Materialize constrained AND held heights in birth order, so holds
    // re-evaluate against post-pass parent positions (a held vertex whose
    // parent just moved must move with it).
    for (Vertex* v : vs) {
        if (v->parent_edge == nullptr) {
            if (v->free_coefficient && v->current_fit_points < 10) {
                v->height = v->neighbor_mean_height();
                v->new_height = v->height;
            }
            continue;
        }
        if (!v->free_coefficient) {
            v->height = v->mu_lin();
            v->new_height = v->height;
        } else if (v->current_fit_points < 10 ||
                   v->current_keff < Vertex::minimum_effective_support()) {
            static const bool HIER_HELD2 = []{
                const char* e = std::getenv("DMESH_HIER_HELD");
                return e != nullptr && std::atoi(e) != 0;
            }();
            static const bool LEGACY_H =
                std::getenv("DMESH_CONTAIN_HEIGHT") != nullptr;
            if (HIER_HELD2 && !LEGACY_H) {
                v->height = v->mu_prior() + v->contained_surplus;
                v->new_height = v->height;
            }
        }
    }
}

}  // namespace hier_fit

void run_height_sweeps(DecisionMesh& mesh,
                       double tolerance,
                       int maximum_sweeps,
                       PipelineTiming& timing) {
    static const bool PL_SOLVER = []{
        const char* e = std::getenv("DMESH_PL_SOLVER");
        return e != nullptr && std::atoi(e) != 0;
    }();
    if (PL_SOLVER) {
        // Honor the caller's convergence contract.  Phase A historically
        // hard-coded 0.5 centi-logit / 30 sweeps, which was adequate for the
        // legacy topology but left newly admitted B3 coordinates materially
        // off their exact fixed point.  Keep at least 30 sweeps during adaptive
        // rounds and allow the fixed-topology closure to request 120.
        pl_solver::run(mesh, tolerance, std::max(maximum_sweeps, 30));
        // Refresh quadratic scores from the new heights so the gate and heap
        // stay fed (phase A contract: gate machinery untouched).
        for (Vertex* vertex : mesh.vertices) {
            if (vertex->active || vertex->parent_edge) vertex->update_info();
        }
        return;
    }

    const double started = now_seconds();
    std::vector<Vertex*> active_vertices;
    active_vertices.reserve(mesh.vertices.size());
    for (Vertex* vertex : mesh.vertices) {
        if (vertex->active) active_vertices.push_back(vertex);
    }
    std::sort(active_vertices.begin(), active_vertices.end(), [](const Vertex* left, const Vertex* right) {
        return left->x != right->x ? left->x < right->x : left->y < right->y;
    });

    for (int sweep = 0; sweep < maximum_sweeps; ++sweep) {
        g_diag_sweep = sweep;
        double largest_gain = 0.0;
        for (Vertex* vertex : active_vertices) {
            // A geometry-only midpoint is a constraint, not a coordinate.  Do
            // not subject its synchronization to the solver tolerance: even a
            // tiny stale gap would create an unintended nonzero surplus.
            if (vertex->parent_edge && !vertex->free_coefficient) {
                const double target = vertex->mu_lin();
                if (vertex->height != target) {
                    vertex->new_height = target;
                    vertex->update_height();
                }
                continue;
            }
            if (vertex->loss_reduction <= tolerance) continue;
            largest_gain = std::max(largest_gain, vertex->loss_reduction);
            const double old_h = vertex->height;
            double proposal_h = vertex->new_height;
            const bool valid = prepare_vertex_proposal(*vertex, proposal_h);
            if (in_diag_region(*vertex))
                log_vertex_event("pre_update", *vertex, old_h, proposal_h, old_h, false);
            if (!valid) {
                vertex->loss_reduction = 0.0;
                mesh.heap_set(vertex, 0.0);
                continue;
            }
            vertex->update_height();
            if (in_diag_region(*vertex))
                log_vertex_event("post_update", *vertex, old_h, proposal_h,
                                 vertex->height, false);
        }
        if (largest_gain <= tolerance) break;
    }
    g_diag_sweep = -1;

    // Fix 1b (starved roots): the sweep loop never visits vertices whose
    // loss_reduction is zero, so a dataless free root keeps its
    // initialization value forever while its faces still shape predictions.
    // Deterministic pass: any active free root with almost no usable
    // training support continues the local field instead.
    static const bool CONTAIN_HEIGHT_R =
        std::getenv("DMESH_CONTAIN_HEIGHT") != nullptr;
    if (!CONTAIN_HEIGHT_R) {
        for (Vertex* vertex : mesh.vertices) {
            if (vertex->active && vertex->free_coefficient &&
                vertex->parent_edge == nullptr &&
                vertex->current_fit_points < 10) {
                vertex->height = vertex->neighbor_mean_height();
                vertex->new_height = vertex->height;
            }
        }
    }
    static const int HIER_FIT = []{
        const char* e = std::getenv("DMESH_HIER_FIT");
        return e ? std::atoi(e) : 0;
    }();
    if (HIER_FIT == 1) hier_fit::run(mesh, 0.5, 30);
    timing.height_sweeps_seconds += now_seconds() - started;
}

struct FinalRetirementStats {
    int retired = 0;
    int batches = 0;
    int woken = 0;
    double max_immediate_train_change = 0.0;
    double max_immediate_heldout_change = 0.0;
    double max_train_change = 0.0;
    double max_heldout_change = 0.0;
    double rms_train_change = 0.0;
    double rms_heldout_change = 0.0;
};

namespace {

struct PredictionChange {
    double max_train = 0.0;
    double max_heldout = 0.0;
    double rms_train = 0.0;
    double rms_heldout = 0.0;
};

std::vector<double> current_predictions(const DecisionMesh& mesh) {
    std::vector<double> result(mesh.n_points);
    for (int i = 0; i < mesh.n_points; ++i) result[i] = mesh.predict_index(i);
    return result;
}

PredictionChange prediction_change(const DecisionMesh& mesh,
                                   const std::vector<double>& before,
                                   bool heldout_split) {
    PredictionChange result;
    double train_sq = 0.0;
    double heldout_sq = 0.0;
    std::size_t train_count = 0;
    std::size_t heldout_count = 0;
    for (int i = 0; i < mesh.n_points; ++i) {
        const double change = mesh.predict_index(i) - before[i];
        if (hierarchical_design::is_training_row(mesh, i, heldout_split)) {
            result.max_train = std::max(result.max_train, std::fabs(change));
            train_sq += change * change;
            ++train_count;
        } else {
            result.max_heldout =
                std::max(result.max_heldout, std::fabs(change));
            heldout_sq += change * change;
            ++heldout_count;
        }
    }
    result.rms_train = train_count ? std::sqrt(train_sq / train_count) : 0.0;
    result.rms_heldout =
        heldout_count ? std::sqrt(heldout_sq / heldout_count) : 0.0;
    return result;
}

void refresh_hierarchical_support_diagnostics(DecisionMesh& mesh,
                                                bool heldout_split) {
    hierarchical_design::BuildOptions options;
    options.heldout_split = heldout_split;
    options.training_only = true;
    const auto design = hierarchical_design::build_active_design(mesh, options);
    constexpr double kStructuralSupportTolerance = 1e-24;
    for (Vertex* vertex : mesh.vertices) {
        if (!vertex || !vertex->active || !vertex->free_coefficient ||
            vertex->parent_edge == nullptr || vertex->retired_coefficient) {
            continue;
        }
        const auto found = design.columns.find(vertex);
        const bool supported = found != design.columns.end() &&
            hierarchical_design::structural_norm_sq(found->second) >
                kStructuralSupportTolerance;
        vertex->dormant_coefficient = !supported;
        if (!supported) {
            vertex->current_n_points = 0;
            vertex->current_fit_points = 0;
            vertex->current_xTx = 0.0;
            vertex->current_fit_xTx = 0.0;
            vertex->current_keff = 0.0;
            vertex->sigma_sq = 1e300;
            vertex->loss_reduction = 0.0;
            mesh.heap_set(vertex, 0.0);
            continue;
        }

        double xTx = 0.0;
        double xT2x = 0.0;
        int points = 0;
        for (const auto& [index, column] : found->second) {
            if (std::fabs(column) <= options.coefficient_tolerance) continue;
            ++points;
            const double weight = std::max(mesh.weights[index], 0.0);
            xTx += weight * column * column;
            xT2x += weight * weight * column * column;
        }
        vertex->current_n_points = points;
        vertex->current_fit_points = points;
        vertex->current_xTx = xTx;
        vertex->current_fit_xTx = xTx;
        vertex->current_keff = xT2x > 0.0 ? xTx * xTx / xT2x : 0.0;
        vertex->sigma_sq = xTx > 0.0 ? 1.0 / xTx : 1e300;
    }
}

}  // namespace

// Fixed-topology cleanup in the SAME hierarchical coordinate system used by
// the exact stage solver. A coefficient is obsolete only if its full
// hierarchical basis function is zero on the training design. Nodal-star
// curvature is deliberately not consulted: a nodally empty ancestor may still
// act through constrained descendants.
FinalRetirementStats retire_zero_support_at_fixed_topology(
    DecisionMesh& mesh,
    bool heldout_split,
    PipelineTiming& timing) {
    FinalRetirementStats stats;
    const std::vector<double> before_cleanup = current_predictions(mesh);

    constexpr double kStructuralSupportTolerance = 1e-24;
    constexpr double kTrainingInvarianceTolerance = 1e-7;  // centi-logits

    while (true) {
        hierarchical_design::BuildOptions options;
        options.heldout_split = heldout_split;
        options.training_only = true;
        const auto design = hierarchical_design::build_active_design(mesh, options);

        std::vector<Vertex*> unsupported;
        int deepest = -1;
        for (Vertex* vertex : mesh.vertices) {
            if (!vertex || !vertex->active || !vertex->free_coefficient ||
                vertex->parent_edge == nullptr || vertex->retired_coefficient) {
                continue;
            }
            const auto found = design.columns.find(vertex);
            const double norm_sq = found == design.columns.end()
                ? 0.0
                : hierarchical_design::structural_norm_sq(found->second);
            if (norm_sq <= kStructuralSupportTolerance) {
                deepest = std::max(deepest, vertex->depth);
                unsupported.push_back(vertex);
            } else if (vertex->dormant_coefficient) {
                // The nodal star may be empty while the hierarchical column is
                // supported through constrained descendants. Wake it without
                // calling update_info(), which still measures nodal support.
                vertex->dormant_coefficient = false;
                ++stats.woken;
            }
        }
        if (unsupported.empty()) break;

        std::vector<Vertex*> batch;
        for (Vertex* vertex : unsupported) {
            if (vertex->depth == deepest) batch.push_back(vertex);
        }
        std::sort(batch.begin(), batch.end(),
                  [](const Vertex* left, const Vertex* right) {
                      return left->_id < right->_id;
                  });

        const std::vector<double> before_batch = current_predictions(mesh);
        for (Vertex* vertex : batch) {
            vertex->dormant_coefficient = false;
            vertex->retired_coefficient = true;
            vertex->free_coefficient = false;
            vertex->lambda_v = 0.0;
            vertex->sigma_pooled = 0.0;
            vertex->delta_pooled = 0.0;
            vertex->loss_reduction = 0.0;
            mesh.heap_set(vertex, 0.0);
        }
        // Retired vertices and every constrained descendant are now exactly at
        // zero surplus, evaluated against already-materialized parents.
        hierarchical_design::materialize_constraints(mesh);

        const PredictionChange immediate =
            prediction_change(mesh, before_batch, heldout_split);
        stats.max_immediate_train_change =
            std::max(stats.max_immediate_train_change, immediate.max_train);
        stats.max_immediate_heldout_change =
            std::max(stats.max_immediate_heldout_change,
                     immediate.max_heldout);
        if (immediate.max_train > kTrainingInvarianceTolerance) {
            std::fprintf(stderr,
                         "[hierarchical retirement] INVARIANT FAILURE at depth %d: "
                         "batch=%zu train max/rms %.12g/%.12g centi-logits\n",
                         deepest, batch.size(), immediate.max_train,
                         immediate.rms_train);
            throw std::runtime_error(
                "retiring a zero hierarchical column changed training predictions");
        }

        stats.retired += static_cast<int>(batch.size());
        ++stats.batches;
        std::fprintf(stderr,
                     "[hierarchical retirement] depth %d retired %zu; "
                     "immediate prediction change train max/rms %.6g/%.6g, "
                     "heldout max/rms %.6g/%.6g centi-logits\n",
                     deepest, batch.size(), immediate.max_train,
                     immediate.rms_train, immediate.max_heldout,
                     immediate.rms_heldout);
        // Rebuild on the next loop. Retiring a deepest descendant can expose
        // training support for an ancestor that was previously blocked by the
        // descendant's independent coordinate.
    }

    if (stats.retired > 0) {
        mesh.recompute_tau_sq();
        run_height_sweeps(mesh, 0.0005, 120, timing);
    }
    refresh_hierarchical_support_diagnostics(mesh, heldout_split);

    // Final fixed-point invariant: every remaining free non-root coefficient
    // has a nonempty training hierarchical column.
    {
        hierarchical_design::BuildOptions options;
        options.heldout_split = heldout_split;
        options.training_only = true;
        const auto design = hierarchical_design::build_active_design(mesh, options);
        for (Vertex* vertex : mesh.vertices) {
            if (!vertex || !vertex->active || !vertex->free_coefficient ||
                vertex->parent_edge == nullptr || vertex->retired_coefficient) {
                continue;
            }
            const auto found = design.columns.find(vertex);
            const double norm_sq = found == design.columns.end()
                ? 0.0
                : hierarchical_design::structural_norm_sq(found->second);
            if (norm_sq <= kStructuralSupportTolerance) {
                throw std::runtime_error(
                    "fixed-topology retirement left an unsupported free coefficient");
            }
        }
    }

    const PredictionChange final_change =
        prediction_change(mesh, before_cleanup, heldout_split);
    stats.max_train_change = final_change.max_train;
    stats.max_heldout_change = final_change.max_heldout;
    stats.rms_train_change = final_change.rms_train;
    stats.rms_heldout_change = final_change.rms_heldout;

    std::fprintf(stderr,
                 "[hierarchical retirement] retired %d coefficients in %d batches, "
                 "woke %d nodally-empty supported coefficients; "
                 "max immediate train/heldout %.6g/%.6g; "
                 "post-refit change train max/rms %.6g/%.6g, "
                 "heldout max/rms %.6g/%.6g centi-logits\n",
                 stats.retired, stats.batches, stats.woken,
                 stats.max_immediate_train_change,
                 stats.max_immediate_heldout_change,
                 stats.max_train_change, stats.rms_train_change,
                 stats.max_heldout_change, stats.rms_heldout_change);
    return stats;
}


struct B1FixedPointAuditStats {
    int audited = 0;
    int unscoreable = 0;
    double max_abs_move_units = 0.0;
    double rms_move_units = 0.0;
    int max_id = -1;
};

B1FixedPointAuditStats run_b1_fixed_point_audit(DecisionMesh& mesh) {
    B1FixedPointAuditStats stats;
    const char* dump_path = std::getenv("DMESH_B1_STALE_DUMP");
    const bool validate = [] {
        const char* e = std::getenv("DMESH_B1_STALE_VALIDATE");
        return e != nullptr && std::atoi(e) != 0;
    }();
    std::ofstream dump;
    if (dump_path && *dump_path) {
        dump.open(dump_path);
        dump << std::setprecision(17)
             << "id,x,y,depth,height_units,prior_units,optimum_units,"
                "move_units,move_logits,z_move,gain_nats,score,information,"
                "lambda_logit,fit_points,structural_points,gate_admitted,"
                "conformity_origin\n";
    }

    hierarchical_design::BuildOptions audit_options;
    audit_options.training_only = false;
    const auto audit_design =
        hierarchical_design::build_active_design(mesh, audit_options);

    double sum_sq = 0.0;
    for (Vertex* vertex : mesh.vertices) {
        if (!vertex || !vertex->active || !vertex->free_coefficient ||
            vertex->retired_coefficient || vertex->parent_edge == nullptr) {
            continue;
        }
        const auto found_column = audit_design.columns.find(vertex);
        if (found_column == audit_design.columns.end()) {
            ++stats.unscoreable;
            continue;
        }
        const auto [lambda, prior] = vertex->compute_eb_params(1.0);
        const auto profile = vertex->exact_gain_on_column(
            found_column->second, vertex->height, vertex->height,
            prior, lambda, false);
        if (!profile.ok) {
            ++stats.unscoreable;
            continue;
        }
        const double move_units = profile.optimum_units - vertex->height;
        const double move_logits = move_units / kHeightScale;
        const double posterior_information =
            std::max(0.0, profile.data_information_optimum +
                              profile.lambda_logit);
        const double z_move = move_logits * std::sqrt(posterior_information);
        ++stats.audited;
        sum_sq += move_units * move_units;
        if (std::fabs(move_units) > stats.max_abs_move_units) {
            stats.max_abs_move_units = std::fabs(move_units);
            stats.max_id = vertex->_id;
        }
        if (dump) {
            dump << vertex->_id << ',' << vertex->x << ',' << vertex->y << ','
                 << vertex->depth << ',' << vertex->height << ',' << prior << ','
                 << profile.optimum_units << ',' << move_units << ','
                 << move_logits << ',' << z_move << ',' << profile.gain_nats << ','
                 << profile.score_at_null << ',' << profile.information_at_null
                 << ',' << profile.lambda_logit << ',' << profile.fit_points << ','
                 << profile.training_structural_points << ','
                 << (vertex->gate_admitted ? 1 : 0) << ','
                 << (vertex->conformity_origin ? 1 : 0) << '\n';
        }
    }
    stats.rms_move_units = stats.audited > 0
        ? std::sqrt(sum_sq / stats.audited) : 0.0;
    std::fprintf(stderr,
                 "[B1 fixed-point audit] audited %d unscoreable %d "
                 "max/rms move %.9g/%.9g centi-logits (max id %d)\n",
                 stats.audited, stats.unscoreable, stats.max_abs_move_units,
                 stats.rms_move_units, stats.max_id);

    // Named regression: the pre-refined boundary coordinate at normalized
    // (0, 0.375), i.e. raw WALA 0 and WAC approximately 4.52. Price the
    // exact hierarchical alternative against its P1 null in the final live
    // state, independently of the legacy Wald score.
    Vertex* regime = nullptr;
    double regime_distance = std::numeric_limits<double>::infinity();
    for (Vertex* vertex : mesh.vertices) {
        if (!vertex || !vertex->active) continue;
        const double distance = std::hypot(vertex->x, vertex->y - 0.375);
        if (distance < regime_distance) {
            regime_distance = distance;
            regime = vertex;
        }
    }
    if (regime && regime_distance < 1e-10) {
        const auto [lambda, prior] = regime->compute_eb_params(1.0);
        const auto profile = regime->exact_candidate_gain(
            regime->mu_lin(), prior, lambda, true);
        std::fprintf(stderr,
                     "[B1 regime-straddle pin] id %d xy=(%.8f,%.8f) "
                     "current %.9f null %.9f optimum %.9f gain %.12g nats "
                     "grid_excess %.3g ok %d\n",
                     regime->_id, regime->x, regime->y,
                     regime->height / kHeightScale,
                     regime->mu_lin() / kHeightScale,
                     profile.optimum_units / kHeightScale,
                     profile.gain_nats, profile.max_grid_excess,
                     profile.ok ? 1 : 0);
        if (validate && (!profile.ok || profile.gain_nats < -1e-10)) {
            throw std::runtime_error(
                "regime-straddle exact profile regression failed");
        }
    } else if (validate) {
        throw std::runtime_error(
            "regime-straddle coordinate (0,0.375) is missing");
    }

    // The final exact coordinate sweep uses a 5e-4 centi-logit stopping
    // threshold. Allow a small multiple for the independent re-profile audit.
    if (validate && (stats.unscoreable != 0 ||
                     stats.max_abs_move_units > 0.01)) {
        throw std::runtime_error(
            "final hierarchical coefficients are not at exact fixed point");
    }
    return stats;
}

void pre_refine_mesh(DecisionMesh& mesh, int target_face_count) {
    while (static_cast<int>(mesh.active_faces.size()) < target_face_count) {
        std::vector<Vertex*> frontier;
        for (Edge* edge : mesh.active_edges) {
            Vertex* midpoint = edge->midpoint;
            if (midpoint && !midpoint->active && !midpoint->disqualified) {
                frontier.push_back(midpoint);
            }
        }
        std::sort(frontier.begin(), frontier.end(), [](const Vertex* left, const Vertex* right) {
            return left->x != right->x ? left->x < right->x : left->y < right->y;
        });
        if (frontier.empty()) break;
        for (Vertex* vertex : frontier) {
            if (vertex->active || vertex->disqualified) continue;
            double proposal = vertex->new_height;
            prepare_vertex_proposal(*vertex, proposal);
            vertex->new_height = proposal;
            vertex->activate(true);
        }
    }

    // The deliberately unthresholded coarse mesh is one statistical model,
    // including vertices that happened to be introduced by NVB completion.
    // Promote every vertex present at the end of pre-refinement; only later
    // conformity vertices remain locked at zero surplus.
    for (Vertex* vertex : mesh.vertices) {
        if (!vertex->active || vertex->free_coefficient) continue;
        vertex->retired_coefficient = false;
        vertex->free_coefficient = true;
        vertex->update_info();
    }
}

static double env_gate_extra_var() {
    static const double v = std::getenv("DMESH_GATE_EXTRA_VAR")
        ? std::atof(std::getenv("DMESH_GATE_EXTRA_VAR")) : 0.0;
    return v;
}

double estimate_dispersion(const DecisionMesh& mesh) {
    int active_vertex_count = 0;
    for (const Vertex* vertex : mesh.vertices) {
        if (vertex->active) ++active_vertex_count;
    }
    // With an explicit per-candidate extra-variance model (DMESH_GATE_EXTRA_VAR,
    // logit^2), divide each Pearson term by its modeled inflation so phi does not
    // double-count the pool-level dispersion the gate already prices.
    const double GEV = env_gate_extra_var();
    const bool onelaw = !mesh.gate_s2.empty();
    double pearson_sum = 0.0;
    int contributing = 0;  // exclude held-out points (weight ~ 0) from the
                           // denominator; counting them halved phi under
                           // DMESH_SPLIT and inflated default-gate z by sqrt(2)
    for (int i = 0; i < mesh.n_points; ++i) {
        if (mesh.weights[i] < 1e-9) continue;
        ++contributing;
        const double residual = to_log_odds(mesh.values[i] - mesh.predict_index(i));
        double term = mesh.weights[i] * kVarianceScale * residual * residual;
        if (onelaw) term /= 1.0 + mesh.weights[i] * kVarianceScale * mesh.gate_s2[i];
        else if (GEV > 0.0) term /= 1.0 + mesh.weights[i] * kVarianceScale * GEV;
        pearson_sum += term;
    }
    return std::max(pearson_sum / std::max(contributing - active_vertex_count, 1), 1.0);
}

void refresh_irls_working_data(DecisionMesh& mesh,
                               const std::vector<Observation>& observations,
                               double extra_variance,
                               bool heldout_split) {
    mesh.eta_frozen.resize(mesh.n_points);
    mesh.likelihood_temper.resize(mesh.n_points);
    mesh.refresh_epoch++;
    for (Vertex* v : mesh.vertices) v->height_frozen = v->height;
    for (int i = 0; i < mesh.n_points; ++i) {
        const double fitted = mesh.predict_index(i);
        const double eta = observations[i].offset + to_log_odds(fitted);
        mesh.eta_frozen[i] = eta;
        const double probability = std::clamp(1.0 / (1.0 + std::exp(-eta)), 1e-4, 1.0 - 1e-4);
        const double information = observations[i].exposure * probability * (1.0 - probability);
        // Safeguarded IRLS: cap the per-pool working step. At saturated
        // probabilities the raw Newton step (k - n p)/(n p (1-p)) reaches
        // hundreds of logits and a single saturated vertex poisons every
        // pool in its faces, propagating a global eruption (the
        // split-dependent runaway). Capping preserves every fixed point
        // (converged steps are ~0) and bounds the divergence dynamics.
        // DMESH_IRLS_STEP_CAP in logits; 0 = legacy uncapped.
        static const double STEP_CAP = []{
            const char* e = std::getenv("DMESH_IRLS_STEP_CAP");
            return e ? std::atof(e) : 4.0;
        }();
        double work_step =
            (observations[i].events - observations[i].exposure * probability) / information;
        if (STEP_CAP > 0.0) work_step = std::clamp(work_step, -STEP_CAP, STEP_CAP);
        mesh.values[i] = fitted + from_log_odds(work_step);
        const double base_weight = information / kVarianceScale;
        static const bool REFRESH_LAW_FIX = std::getenv("DMESH_REFRESH_LAW_FIX") != nullptr;
        const double point_extra = (REFRESH_LAW_FIX && !mesh.gate_s2.empty())
            ? kVarianceScale * mesh.gate_s2[i] : extra_variance;
        mesh.weights[i] = point_extra < 0.0
            ? base_weight
            : 1.0 / (1.0 / base_weight + point_extra);
        static const bool HELDOUT_ZERO = std::getenv("DMESH_HELDOUT_ZERO") != nullptr;
        if (heldout_split && !hierarchical_design::is_training_row(mesh, i, heldout_split))
            mesh.weights[i] = HELDOUT_ZERO ? 0.0 : 1e-12;
        // Freeze the exact-likelihood temper at this IRLS epoch.  Since
        // base_weight = information/kVarianceScale, weights/base_weight is
        // exactly the declared law/overdispersion discount c_i.
        mesh.likelihood_temper[i] = (base_weight > 0.0)
            ? std::clamp(mesh.weights[i] / base_weight, 0.0, 1.0)
            : 0.0;
    }

    for (Face* face : mesh.active_faces) face->refresh_suff_stats_from_cache();
    for (Vertex* vertex : mesh.vertices) {
        if (vertex->active || vertex->parent_edge) vertex->update_info();
    }
}

}  // namespace

using Obs = Observation;

int run_pipeline(const RunConfig& config, Dataset& dataset, const SurfaceSpec& surface) {
    std::vector<Obs>& obs = dataset.observations;
    const int pid = dataset.pool_count;
    const double sigma_u = dataset.truth_sigma_u;
    const int N = static_cast<int>(obs.size());
    std::printf("pool-months %d, pools %d\n", N, pid);

    auto t_start = std::chrono::steady_clock::now();
    auto elapsed = [&]() {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
    };
    PipelineTiming timing;

    std::vector<double> sx(N), sy(N), y0(N), wbin(N);
    for (int i = 0; i < N; ++i) { sx[i] = obs[i].x; sy[i] = obs[i].y; y0[i] = obs[i].working_response; wbin[i] = obs[i].working_weight; }

    // --- learned axis metrics (DMESH_METRIC=1): after the stage-0 pilot fit,
    // estimate per-band directional curvature energy of the fitted surface
    // (data-density weighted), and equidistribute each axis so that refinement
    // budget per unit of transformed coordinate is constant. This incorporates
    // the local scale of each axis into the algorithm instead of fixing it by
    // the global rescale to the unit square. Piecewise-linear monotone maps.
    const bool use_metric = config.fit.learned_metric;
    std::vector<double> sx0, sy0;   // pre-transform coords (empty if metric off)
    const int MB = 24;                       // bands per axis
    std::vector<double> phiX(MB + 1), phiY(MB + 1);
    for (int b = 0; b <= MB; ++b) phiX[b] = phiY[b] = b / (double)MB;
    auto apply_map = [&](const std::vector<double>& phi, double t) {
        double u = std::clamp(t, 0.0, 1.0) * MB;
        int b = std::min((int)u, MB - 1);
        return phi[b] + (u - b) * (phi[b + 1] - phi[b]);
    };
    double sig_cur = 0;
    std::unique_ptr<DecisionMesh> mesh;
    const int n_stages = config.diagnostics.stage_count;
    int nv_final = 0;
    for (int stage = 0; stage < n_stages; ++stage) {
        const double stage_input_sigma = sig_cur;
        g_diag_stage = stage;
        g_diag_round = -1;
        g_diag_phase = "stage_init";
        const bool oracle_mode = config.fit.oracle_mode();
        const bool use_eb = config.fit.use_eb();
        if (oracle_mode) for (int i = 0; i < N; ++i) y0[i] = obs[i].true_surface_scaled;   // noiseless truth
        mesh = std::make_unique<DecisionMesh>(sx, sy, y0, use_eb);
        mesh->obs_store = &obs;   // pool-level data for exact-likelihood scoring
        mesh->training_mask = dataset.training_mask;
        mesh->maintain_heap = false;  // lfdr gate never uses the greedy heap
        if (stage == 0 && use_eb) {
            std::printf("EB: zero-mean surplus prior, gate-admitted tau, "
                        "default-sd=%.3f, floor-sd=%.3f log-odds\n",
                        std::sqrt(mesh->tau_default_sq) / kHeightScale,
                        std::sqrt(mesh->tau_floor_sq) / kHeightScale);
        }
        const bool heldout_split = config.fit.heldout_split;
        // One-law gate: per-point pool-effect variance and kurtosis from the
        // measured mixture (strata <128 / <512 / <4096 / rest). Filled once at
        // fit entry so every face accumulation, from stage 0 on, sees the law.
        if (const char* ol = std::getenv("DMESH_ONELAW")) {
            double s2[4], k4[4];
            if (std::sscanf(ol, "%lf,%lf,%lf,%lf;%lf,%lf,%lf,%lf",
                            &s2[0],&s2[1],&s2[2],&s2[3],&k4[0],&k4[1],&k4[2],&k4[3]) == 8) {
                mesh->gate_s2.assign(mesh->n_points, 0.0);
                mesh->gate_k4.assign(mesh->n_points, 0.0);
                for (int i = 0; i < mesh->n_points; ++i) {
                    double n = obs[i].exposure;
                    int b = n < 128 ? 0 : n < 512 ? 1 : n < 4096 ? 2 : 3;
                    mesh->gate_s2[i] = s2[b];
                    mesh->gate_k4[i] = s2[b] * s2[b] * k4[b];
                }
                std::printf("one-law gate: strata s2 %.3f/%.3f/%.3f/%.3f kurt %.2f/%.2f/%.2f/%.2f\n",
                            s2[0],s2[1],s2[2],s2[3],k4[0],k4[1],k4[2],k4[3]);
                // Optional per-point pool-effect variances (one float per obs
                // line), e.g. the mixture-posterior gamma-weighted variances
                // from the joint fixed point. Overrides the stratum variance
                // per point; the kurtosis table still supplies shape.
                if (const char* pv = std::getenv("DMESH_POINT_VAR")) {
                    std::FILE* f = std::fopen(pv, "r");
                    if (f) {
                        for (int i = 0; i < mesh->n_points; ++i) {
                            double v;
                            if (std::fscanf(f, "%lf", &v) != 1) break;
                            double n = obs[i].exposure;
                            int b = n < 128 ? 0 : n < 512 ? 1 : n < 4096 ? 2 : 3;
                            mesh->gate_s2[i] = v;
                            mesh->gate_k4[i] = v * v * k4[b];
                        }
                        std::fclose(f);
                        std::printf("one-law gate: per-point variances loaded from %s\n", pv);
                    }
                }
            }
        }
        // DMESH_EXTRA_VAR_FLOOR (logit^2): lower bound on the extra variance in
        // the DL working weights, applied from STAGE 0. On single-snapshot data
        // sig_cur floors near zero (no within-pool pairs), weights collapse to
        // binomial, and overdispersed mega pools (w*sigma >> 1) dominate the
        // fit and blow up the candidate null: measured sd(z) 3.0-4.3 at deep
        // depths on the matched-null benchmark vs the sigma0 <= 1.4 box. The
        // floor caps effective weight at 1/sigma and makes candidate variances
        // consistent with the actual residual law.
        static const double EXTRA_FLOOR = std::getenv("DMESH_EXTRA_VAR_FLOOR")
            ? std::atof(std::getenv("DMESH_EXTRA_VAR_FLOOR")) : 0.0;
        double extra;
        if (EXTRA_FLOOR > 0.0)
            extra = std::max((stage == 0) ? 0.0 : sig_cur, EXTRA_FLOOR) * kVarianceScale;
        else
            extra = (stage == 0) ? -1 : sig_cur * kVarianceScale;  // sigma_w demoted to diagnostic (inconsistent; see BUILD)
        mesh->last_extra = extra;
        for (int i = 0; i < N; ++i) {
            static const bool LAW_W = std::getenv("DMESH_LAW_WEIGHTS") != nullptr;
            double extra_i = (LAW_W && !mesh->gate_s2.empty())
                ? kVarianceScale * mesh->gate_s2[i] : extra;
            mesh->weights[i] = extra_i < 0 ? wbin[i] : 1.0 / (1.0 / wbin[i] + extra_i);
            static const bool HELDOUT_ZERO_INIT = std::getenv("DMESH_HELDOUT_ZERO") != nullptr;
            if (heldout_split && !hierarchical_design::is_training_row(*mesh, i, heldout_split))
                mesh->weights[i] = HELDOUT_ZERO_INIT ? 0.0 : 1e-12;  // held out
            // Initial-stage likelihood temper, before the first IRLS refresh.
            // wbin is the binomial working precision in the same units as
            // mesh.weights, so their ratio is c_i.
            mesh->likelihood_temper[i] = (wbin[i] > 0.0)
                ? std::clamp(mesh->weights[i] / wbin[i], 0.0, 1.0)
                : 0.0;
        }
        for (Face* fc : std::vector<Face*>(mesh->active_faces.begin(), mesh->active_faces.end()))
            fc->refresh_suff_stats_from_cache();
        // The constructor initializes geometry before split membership exists.
        // Recompute root support after installing the explicit mask so the
        // pre-refinement pass cannot inherit held-out sufficient statistics.
        for (Vertex* v : mesh->vertices) v->update_info();

        g_diag_phase = "pre_refine";
        pre_refine_mesh(*mesh, 48);
        mesh->recompute_tau_sq();
        g_diag_phase = "pre_round_sweeps";
        run_height_sweeps(*mesh, 0.01, 20, timing);
        {
            int n_lr = 0, n_act = 0; double mh = 0;
            for (Vertex* v : mesh->vertices) {
                if (!v->active) continue;
                n_act++;
                if (v->loss_reduction > 0.01) n_lr++;
                mh += std::fabs(v->height);
            }
            std::printf("  [diag pre-rounds] active %d, with loss_red>tol %d, mean|h| %.2f, pred(0.5,0.5)=%.1f\n",
                        n_act, n_lr, mh / std::max(n_act, 1), mesh->predict(0.5, 0.5));
        }
        // Stage 0 needs only a rough sigma_u^2 (small cap). Later stages cap at a
        // data-density budget: ~1 face per 12 points prevents runaway refinement
        // when sigma_u is small (DL weights ~ binomial) without touching the fit
        // in the frailty regime, where the gate stops well below the cap anyway.
        int vertex_cap = (stage == 0) ? 1200 : std::max(600, N / config.fit.points_per_face);
        // Quasi-likelihood dispersion for candidate-score calibration: without it,
        // overdispersion (white pool noise not in the weights) inflates every z,
        // blows past the empirical-null sigma0<=1.4 cap, and the gate mass-admits.
        double phi_disp = estimate_dispersion(*mesh);
        int stalls = 0;
        for (int rnd = 0; rnd < config.diagnostics.max_gate_rounds; ++rnd) {
            g_diag_round = rnd;
            if ((int)mesh->active_faces.size() > vertex_cap) break;
            const bool oracle_scoring = config.fit.oracle_mode();
            g_diag_phase = "round_sweeps";
            if (!oracle_scoring) { refresh_hierarchical_support_diagnostics(*mesh, heldout_split); mesh->recompute_tau_sq(); run_height_sweeps(*mesh, 0.01, 12, timing); }
            g_diag_phase = "candidate_score";

            // --- Exact-binomial likelihood-ratio candidate score (DMESH_LR_SCORE=1) ---
            // Signed root of 2[l(delta_hat) - l(0)] where l profiles the candidate
            // surplus through the exact binomial pmf with the pool effect
            // integrated out per observation (Gauss-Hermite). The mixing law is a
            // drop-in: Gaussian at the measured size profile sigma^2(n), plus an
            // scale 4x, size-decaying like the measured kurtosis). Binomial
            // coefficients cancel in the ratio; k=n is exact, no Gaussian tails.
            static const bool LR_ON = std::getenv("DMESH_LR_SCORE") != nullptr;
            auto sig2_of_n = [](double nl) {
                double l2 = std::log2(std::max(nl, 32.0));
                const double xs[5] = {5, 8, 10, 13, 21}, ys[5] = {0.47, 0.36, 0.30, 0.38, 0.38};
                if (l2 <= xs[0]) return ys[0];
                for (int i = 1; i < 5; ++i)
                    if (l2 <= xs[i]) return ys[i-1] + (ys[i]-ys[i-1])*(l2-xs[i-1])/(xs[i]-xs[i-1]);
                return ys[4];
            };
            auto lr_score = [&](Vertex* v, double sp_centi2) -> double {
                auto frr = v->active ? v->get_faces(false) : v->get_sim_faces(false);
                std::vector<double> phv, kv, nv, ev, s2v;
                for (Face* f : frr.faces) {
                    int si = -1;
                    for (int i = 0; i < 3; ++i) if (f->vertices[i] == v) { si = i; break; }
                    if (si < 0 || f->coords_indices.empty()) continue;
                    Vertex *A = f->vertices[0], *B = f->vertices[1], *C = f->vertices[2];
                    double det = (B->x - A->x)*(C->y - A->y) - (C->x - A->x)*(B->y - A->y);
                    if (std::fabs(det) < 1e-14) continue;
                    // H0 baseline from THIS face's plane with the candidate's height
                    // at its parent interpolation (exactly the pre-split surface):
                    // predict_index is cache-based and goes stale after mid-round
                    // subdivisions, which leaked just-fitted signal into the score.
                    double hv[3];
                    for (int j = 0; j < 3; ++j)
                        hv[j] = (j == si) ? v->mu_prior() : f->vertices[j]->height;
                    for (int i : f->coords_indices) {
                        double l0 = ((B->x - sx[i])*(C->y - sy[i]) - (C->x - sx[i])*(B->y - sy[i]))/det;
                        double l1 = ((C->x - sx[i])*(A->y - sy[i]) - (A->x - sx[i])*(C->y - sy[i]))/det;
                        double lb2 = 1.0 - l0 - l1;
                        double lam = (si == 0) ? l0 : (si == 1) ? l1 : lb2;
                        if (lam <= 1e-6) continue;
                        double pred = l0*hv[0] + l1*hv[1] + lb2*hv[2];
                        phv.push_back(lam);
                        kv.push_back((double)obs[i].events);
                        nv.push_back(obs[i].exposure);
                        ev.push_back(obs[i].offset + pred / kHeightScale);
                        double nl = obs[i].exposure;
                        s2v.push_back(sig2_of_n(nl)
                            + config.fit.parent_variance_scale * sp_centi2 / kVarianceScale);
                    }
                }
                if (phv.size() < 6) return std::nan("");
                // Retired: heavy-tail LR statistic (superseded; 407-543
                // certified-false on its own benchmark). lr_score now
                // serves only the permutation gate.
                return std::nan("");
            };
            double _tsc = now_seconds();
            std::vector<Vertex*> cands;
            std::vector<double> zs;
            std::vector<double> zsd;
            std::vector<double> sps;
            std::vector<Vertex::ExactCandidateGain> exact_profiles;
            std::vector<double> histogram_pvalues;

            // Candidate family = new inactive edge midpoints plus conformity-only
            // vertices that already exist geometrically but have zero surplus.
            // The latter may be promoted only after passing the same evidence gate.
            std::vector<Vertex*> frontier_vertices;
            frontier_vertices.reserve(mesh->active_edges.size() + mesh->vertices.size() / 8);
            for (Edge* edge : mesh->active_edges) {
                Vertex* midpoint = edge->midpoint;
                if (midpoint && !midpoint->active && !midpoint->disqualified)
                    frontier_vertices.push_back(midpoint);
            }
            // Geometry-only vertices remain members of future candidate
            // families.  Promotion is part of the model semantics rather than
            // an experimental switch: geometry may be required today and earn
            // an independent coefficient only after later evidence appears.
            if (!config.fit.oracle_mode()) {
                for (Vertex* vertex : mesh->vertices) {
                    if (vertex->active && vertex->parent_edge &&
                        !vertex->free_coefficient && !vertex->retired_coefficient &&
                        !vertex->disqualified)
                        frontier_vertices.push_back(vertex);
                }
            }
            std::sort(frontier_vertices.begin(), frontier_vertices.end(),
                      [](const Vertex* left, const Vertex* right) {
                          if (left->x != right->x) return left->x < right->x;
                          if (left->y != right->y) return left->y < right->y;
                          return left->_id < right->_id;
                      });
            frontier_vertices.erase(std::unique(frontier_vertices.begin(), frontier_vertices.end()),
                                    frontier_vertices.end());

            for (Vertex* v : frontier_vertices) {
                if (!v || v->retired_coefficient || v->disqualified ||
                    (v->active && v->free_coefficient)) continue;
                Edge* parent = v->parent_edge;
                if (!parent) continue;
                auto lr = v->loc_regress(false);
                const double candidate_keff = lr.xT2x > 0.0
                    ? lr.xTx * lr.xTx / lr.xT2x : 0.0;
                if (lr.xTx <= 0 || lr.fit_points < Vertex::kMinimumFitPoints ||
                    candidate_keff < Vertex::minimum_effective_support()) continue;
                double sp = 0;
                for (Vertex* p : {parent->vertex0, parent->vertex1})
                    if (p->sigma_pooled < 1e290) sp += p->sigma_pooled;
                const double GEV = env_gate_extra_var();
                // Test-side null: Var(beta_hat) = 1/xTx + sigma_pool^2 * xT2x/xTx^2.
                // Fit weights stay as configured; only the ADMISSION null prices
                // pool-level variance, so admitted structure means cross-pool
                // replicated beyond pool noise (the identified question on a
                // single snapshot). A star dominated by one pool has
                // xT2x/xTx^2 ~ 1/w_dominant and its null variance blows up:
                // the continuous form of a min-effective-pools gate.
                static const bool ONELAW = std::getenv("DMESH_ONELAW") != nullptr;
                double gate_var = 1.0 / lr.xTx
                    + (ONELAW ? kVarianceScale * lr.xT2xs / (lr.xTx * lr.xTx)
                              : GEV * kVarianceScale * lr.xT2x / (lr.xTx * lr.xTx));
                double sd_cal = std::sqrt(
                    (gate_var + config.fit.parent_variance_scale * sp) * phi_disp);
                // B1 computes the exact penalized alternative against the true
                // P1-preserving null for every candidate.  The legacy gate
                // statistic is deliberately left unchanged until B2, so its
                // historical center remains separate from the exact null.
                const double exact_null = v->mu_lin();
                auto [candidate_lambda, candidate_prior] =
                    v->compute_eb_params(lr.xTx);
                static const bool B1_VALIDATE = [] {
                    const char* e = std::getenv("DMESH_B1_VALIDATE");
                    return e != nullptr && std::atoi(e) != 0;
                }();
                const bool near_regime =
                    v->x <= 0.015 && std::fabs(v->y - 0.3750871809178406) <= 0.03;
                const bool grid_check = B1_VALIDATE
                    && (near_regime || g_b1_diag.grid_checked < 24);
                ++g_b1_diag.attempted;
                const auto exact = v->exact_candidate_gain(
                    exact_null, candidate_prior, candidate_lambda, grid_check);
                if (grid_check) {
                    ++g_b1_diag.grid_checked;
                    g_b1_diag.max_grid_excess =
                        std::max(g_b1_diag.max_grid_excess, exact.max_grid_excess);
                }
                if (!exact.ok) {
                    ++g_b1_diag.unscoreable;
                    if (exact.training_structural_points == 0)
                        ++g_b1_diag.unscoreable_no_training_column;
                    else
                        ++g_b1_diag.unscoreable_tempered_out;
                    if (grid_check) ++g_b1_diag.grid_failed;
                    continue;
                }
                ++g_b1_diag.scored;
                if (exact.raw_gain_nats < -1e-8) {
                    ++g_b1_diag.negative;
                    continue;
                }
                // A zero-gain alternative contains no evidence beyond the
                // constrained model.  Do not allow a stale Wald statistic to
                // admit it while B2 calibration is still pending.
                if (!(exact.gain_nats > 1e-12)) continue;

                const double regime_distance = std::hypot(
                    v->x, v->y - 0.3750871809178406);
                if (regime_distance < g_b1_diag.nearest_regime_distance) {
                    g_b1_diag.nearest_regime_distance = regime_distance;
                    g_b1_diag.nearest_stage = stage;
                    g_b1_diag.nearest_round = rnd;
                    g_b1_diag.nearest_id = v->_id;
                    g_b1_diag.nearest_x = v->x;
                    g_b1_diag.nearest_y = v->y;
                    g_b1_diag.nearest_gain = exact.gain_nats;
                    g_b1_diag.nearest_null = exact.null_units;
                    g_b1_diag.nearest_optimum = exact.optimum_units;
                }

                const double legacy_center = (v->active && !v->free_coefficient)
                    ? v->mu_lin() : v->mu_prior();
                double cand_disp = lr.xTr / lr.xTx - legacy_center;
                // DMESH_EXACT_SCORE: replace the Wald displacement with the
                // exact-likelihood profile MLE displacement when the quadratic
                // proposal leaves its trust region or the star is boundary-
                // clipped. sd_cal and the kurtosis-matched null are unchanged:
                // only the biased bounded displacement estimate is repaired.
                static const bool EXACT_G = std::getenv("DMESH_EXACT_SCORE") != nullptr;
                if (EXACT_G && !LR_ON) {
                    const double eps = 1e-9;
                    const bool onb = v->x < eps || v->x > 1 - eps || v->y < eps || v->y > 1 - eps;
                    if (std::fabs(cand_disp) > 0.5 * kHeightScale || onb) {
                        auto es = v->exact_score(legacy_center, 0.0);
                        if (es.ok) cand_disp = es.b_star_units - legacy_center;
                    }
                }
                double legacy_zc = LR_ON ? lr_score(v, sp)
                                         : cand_disp / sd_cal;
                if (ONELAW && !LR_ON && std::isfinite(legacy_zc) && lr.xT4k > 0) {
                    // Monotone probability match through the declared law's own
                    // tail: the null is modeled as a two-component Gaussian scale
                    // mixture matched to the candidate's variance and fourth
                    // cumulant (from the measured per-stratum mixture). Strictly
                    // monotone in z by construction; no series expansion.
                    double x4 = lr.xTx * lr.xTx * lr.xTx * lr.xTx;
                    double vtot = gate_var + config.fit.parent_variance_scale * sp;
                    double lam4 = kVarianceScale * kVarianceScale * lr.xT4k / x4
                                  / (vtot * vtot);
                    if (lam4 > 1e-3) {
                        double w = 1.8 / (lam4 + 3.0);            // < 3/(lam4+3): va > 0
                        double dd = std::sqrt(lam4 / (3.0 * w * (1.0 - w)));
                        double va = 1.0 - w * dd, vb = 1.0 + (1.0 - w) * dd;
                        double az = std::fabs(legacy_zc);
                        double S = (1.0 - w) * 0.5 * std::erfc(az / std::sqrt(2.0 * va))
                                 + w * 0.5 * std::erfc(az / std::sqrt(2.0 * vb));
                        // upper-tail probit (Acklam rational approximation)
                        auto probit_upper = [](double p) {
                            p = std::clamp(p, 1e-300, 0.5);
                            double t = std::sqrt(-2.0 * std::log(p));
                            return t - (2.515517 + 0.802853*t + 0.010328*t*t)
                                       / (1.0 + 1.432788*t + 0.189269*t*t + 0.001308*t*t*t);
                        };
                        double zt = probit_upper(S);
                        legacy_zc = (legacy_zc < 0 ? -zt : zt);
                    }
                }
                // B3 calibration.  The raw signed-root gain is exact for the
                // tempered penalized objective, but its null variance is not one
                // when the objective is tempered, the EB penalty is active, or
                // pool-level noise exceeds the working likelihood.  In the local
                // quadratic limit,
                //
                //   r_exact ~= delta_hat * sqrt(I_exact),
                //
                // while sd_cal is the pool-noise/replication-aware null standard
                // error of delta_hat in centi-logit units.  Their product gives
                // the candidate-specific null scale of r_exact.  Dividing by it
                // is the sandwich-calibrated exact score.  This retains the exact
                // nonlinear direction and gain while reusing only the declared
                // null variance model from the legacy gate.
                const double b3_null_scale =
                    (sd_cal / kHeightScale)
                    * std::sqrt(std::max(exact.information_at_null, 0.0));
                const double b3_calibrated_score =
                    (b3_null_scale > 1e-12 && std::isfinite(b3_null_scale))
                    ? exact.signed_root_gain / b3_null_scale
                    : std::numeric_limits<double>::quiet_NaN();

                static const bool B3_USE_SCORE = [] {
                    const char* e = std::getenv("DMESH_B3_USE_SCORE");
                    return e != nullptr && std::atoi(e) != 0;
                }();
                static const bool B3_STAGE0_EXACT = [] {
                    const char* e = std::getenv("DMESH_B3_STAGE0_EXACT");
                    return e != nullptr && std::atoi(e) != 0;
                }();
                // Stage 0 is a pilot used to identify the pool-variance scale.
                // Before that scale exists, using an exact-score pivot would be
                // circular.  Keep the established pilot calibration there by
                // default; stages 1+ use the B3 exact score.  The stage-0 exact
                // path remains available only as a diagnostic.
                const bool b3_stage0_legacy = B3_USE_SCORE && stage == 0
                    && !B3_STAGE0_EXACT;

                static const bool HIST_NULL = [] {
                    const char* e = std::getenv("DMESH_HIST_NULL");
                    return e != nullptr && std::atoi(e) != 0;
                }();
                Vertex::HistogramNullResult histogram_null;
                if (HIST_NULL) {
                    histogram_null = v->histogram_null_score(
                        exact_null, candidate_prior, candidate_lambda, exact);
                    if (!histogram_null.ok || !std::isfinite(histogram_null.p_value)
                        || !std::isfinite(histogram_null.gaussian_equivalent)) {
                        continue;
                    }
                }

                // B2's raw-score switch remains as a backwards-compatible
                // diagnostic.  B3 takes precedence when both are set.
                static const bool B2_USE_SCORE = [] {
                    const char* e = std::getenv("DMESH_B2_USE_SCORE");
                    return e != nullptr && std::atoi(e) != 0;
                }();
                const double zc = HIST_NULL
                    ? histogram_null.gaussian_equivalent
                    : (B3_USE_SCORE
                        ? (b3_stage0_legacy ? legacy_zc : b3_calibrated_score)
                        : (B2_USE_SCORE ? exact.signed_root_gain : legacy_zc));
                if (std::isfinite(zc) && std::isfinite(legacy_zc)
                    && std::isfinite(exact.signed_root_gain)) {
                    cands.push_back(v);
                    zs.push_back(zc);
                    zsd.push_back(sd_cal);
                    sps.push_back(sp);
                    exact_profiles.push_back(exact);
                    histogram_pvalues.push_back(HIST_NULL ? histogram_null.p_value : 1.0);
                    auto& gain_dump = b1_gain_stream();
                    if (gain_dump) {
                        gain_dump << stage << ',' << rnd << ',' << v->_id << ','
                                  << v->x << ',' << v->y << ',' << v->depth << ','
                                  << ((v->active && !v->free_coefficient) ? 1 : 0) << ','
                                  << legacy_zc << ',' << exact.signed_root_gain << ','
                                  << zc << ',' << std::fabs(legacy_zc) << ','
                                  << std::fabs(exact.signed_root_gain) << ','
                                  << legacy_center << ',' << cand_disp << ',' << sd_cal << ','
                                  << exact.null_units << ',' << exact.prior_mean_units << ','
                                  << exact.optimum_units << ','
                                  << (exact.optimum_units - exact.null_units) << ','
                                  << exact.gain_nats << ',' << exact.data_gain_nats << ','
                                  << exact.penalty_gain_nats << ',' << exact.raw_gain_nats << ','
                                  << exact.score_at_null << ',' << exact.information_at_null << ','
                                  << exact.data_information_optimum << ',' << exact.lambda_logit << ','
                                  << exact.structural_points << ','
                                  << exact.training_structural_points << ','
                                  << exact.positive_information_points << ','
                                  << exact.fit_points << ',' << candidate_keff << ','
                                  << exact.min_eta_null << ',' << exact.max_eta_null << ','
                                  << exact.min_p_null << ',' << exact.max_p_null << ','
                                  << phi_disp << ',' << gate_var << ','
                                  << exact.max_grid_excess << ','
                                  << b3_null_scale << ',' << b3_calibrated_score << ','
                                  << (b3_stage0_legacy ? 1 : 0) << ','
                                  << (HIST_NULL ? histogram_null.p_value : 1.0) << ','
                                  << (HIST_NULL ? histogram_null.gaussian_equivalent : 0.0) << ','
                                  << (HIST_NULL ? histogram_null.observed_data_score : 0.0) << ','
                                  << (HIST_NULL ? histogram_null.null_mean : 0.0) << ','
                                  << (HIST_NULL ? histogram_null.null_sd : 0.0) << ','
                                  << (HIST_NULL ? histogram_null.tail_keff : 0.0) << ','
                                  << (HIST_NULL ? histogram_null.rows : 0) << ','
                                  << (HIST_NULL ? histogram_null.explicit_mixture_rows : 0)
                                  << '\n';
                    }
                    if (config.diagnostics.dump_score_diagnostics)
                        {
                            double keff = (lr.xT2x > 0 ? lr.xTx*lr.xTx/lr.xT2x : 0.0);
                            std::fprintf(stderr, "ZD %d %d %d %.5f %.5f %.5f %.2f\n", stage, rnd, v->depth, zc, v->x, v->y, keff);
                        }
                }
            }
            // Retired: Monte Carlo null calibration (superseded by the
            // composite/one-law analytic null and the permutation gate).
            timing.candidate_scoring_seconds += now_seconds() - _tsc;
            if (cands.empty()) break;
            const int oracle_vertex_budget = config.fit.oracle_vertex_budget;
            if (oracle_vertex_budget > 0) {
                // greedy best-N-term surrogate: admit by exact projection gain
                {
                    int na = 0;
                    for (Vertex* v : mesh->vertices) if (v->active) na++;
                    if (na >= oracle_vertex_budget) break;
                }
                std::vector<int> og(cands.size());
                for (size_t i = 0; i < og.size(); ++i) og[i] = (int)i;
                // scores: xTr^2/xTx already reflected in |z| ordering? recompute gain:
                std::vector<double> gain(cands.size());
                for (size_t i = 0; i < cands.size(); ++i) {
                    auto lr = cands[i]->loc_regress(false);
                    gain[i] = lr.xTx > 0 ? lr.xTr * lr.xTr / lr.xTx : 0.0;
                }
                std::sort(og.begin(), og.end(), [&](int a, int b) { return gain[a] > gain[b]; });
                int nact = 0;
                for (Vertex* v : mesh->vertices) if (v->active) nact++;
                int kk2 = std::max(1, (int)(0.15 * cands.size()));
                kk2 = std::min(kk2, std::max(1, oracle_vertex_budget - nact));
                for (int i = 0; i < kk2; ++i) {
                    Vertex* v = cands[og[i]];
                    if (v->active || v->disqualified) continue;
                    double proposal = v->new_height;
                    prepare_vertex_proposal(*v, proposal);
                    v->new_height = proposal;
                    v->activate(true);
                }
                // nodal interpolation oracle: heights are exact truth at vertices;
                // no fitting, so no under-determined junk between data points
                for (Vertex* v : mesh->vertices)
                    if (v->active) v->height = from_log_odds(surface.evaluate(-4 + 8 * v->x, -4 + 8 * v->y));
                nact = 0;
                for (Vertex* v : mesh->vertices) if (v->active) nact++;
                if (nact >= oracle_vertex_budget) break;
                continue;
            }
            std::vector<double> lf;
            LfdrResult lfdr_result;
            std::vector<int> direct_bh_order;
            int direct_bh_k = -1;
            bool lindsey_fallback = false;
            static const bool HIST_NULL_GATE = [] {
                const char* e = std::getenv("DMESH_HIST_NULL");
                return e != nullptr && std::atoi(e) != 0;
            }();
            static const bool CANDIDATE_BH = [] {
                const char* e = std::getenv("DMESH_CANDIDATE_BH");
                return e != nullptr && std::atoi(e) != 0;
            }();
            auto run_direct_bh = [&]() {
                std::vector<int> op(zs.size());
                std::vector<double> pv(zs.size());
                for (size_t i = 0; i < zs.size(); ++i) {
                    pv[i] = HIST_NULL_GATE ? histogram_pvalues[i]
                                           : 2 * normal_survival(std::fabs(zs[i]));
                    op[i] = static_cast<int>(i);
                }
                std::sort(op.begin(), op.end(), [&](int a, int b) { return pv[a] < pv[b]; });
                int kbh = 0;
                for (size_t i = 0; i < op.size(); ++i)
                    if (pv[op[i]] <= config.fit.fdr_level * (i + 1) / op.size())
                        kbh = static_cast<int>(i) + 1;
                direct_bh_order = std::move(op);
                direct_bh_k = kbh;
                lf.assign(zs.size(), 1.0);
                for (int i = 0; i < kbh; ++i) lf[direct_bh_order[i]] = 0.0;
            };

            if (!HIST_NULL_GATE && !CANDIDATE_BH && zs.size() >= 50) {
                // B3 default: after the stage-0 legacy pilot, let central
                // matching estimate the residual family-level scale outside
                // the obsolete [1,1.4] box.  Candidate-specific sandwich
                // scaling handles heterogeneity; this second scale captures
                // adaptive fit-state contraction/expansion.  Keep an explicit
                // fixed-box ablation for regression comparisons.
                const bool b3_score_enabled = [] {
                    const char* e = std::getenv("DMESH_B3_USE_SCORE");
                    return e != nullptr && std::atoi(e) != 0;
                }();
                const bool fixed_empirical_null = [] {
                    const char* e = std::getenv("DMESH_B3_FIXED_EMPIRICAL_NULL");
                    return e != nullptr && std::atoi(e) != 0;
                }();
                const bool free_empirical_null =
                    b3_score_enabled && stage > 0 && !fixed_empirical_null;
                lfdr_result = lindsey_lfdr(zs, free_empirical_null);
                if (lfdr_result.used_lindsey) {
                    lf = lfdr_result.lfdr;
                } else {
                    // The local-density fit is optional calibration machinery,
                    // not authorization to fail open.  If it is nonconverged,
                    // non-normalizable, or lacks a concave central peak, use
                    // the declared candidate scores with direct BH instead.
                    lindsey_fallback = true;
                    run_direct_bh();
                }
            } else {
                // This branch already performed a frequentist BH decision.
                // Do not encode pass/fail as pseudo-lfdr values and then run
                // the mean-lfdr prefix a second time: that would append about
                // q/(1-q) arbitrary failed hypotheses after the BH discoveries.
                run_direct_bh();
            }
            std::vector<int> order;
            int kk = 0;
            static const bool ADMIT_ALL = std::getenv("DMESH_ADMIT_ALL") != nullptr;
            if (direct_bh_k >= 0) {
                order = direct_bh_order;
                kk = ADMIT_ALL ? static_cast<int>(order.size()) : direct_bh_k;
            } else {
                order.resize(lf.size());
                for (size_t i = 0; i < lf.size(); ++i) order[i] = (int)i;
                std::sort(order.begin(), order.end(), [&](int a, int b) { return lf[a] < lf[b]; });
                if (ADMIT_ALL) {
                    // clean-measurement mode: refine everything scoreable; no selection
                    kk = (int)order.size();
                } else {
                    double cum = 0.0;
                    for (size_t i = 0; i < order.size(); ++i) {
                        cum += lf[order[i]];
                        if (cum / (i + 1) <= config.fit.fdr_level)
                            kk = (int)i + 1;
                        else
                            break;
                    }
                }
            }
            {
                double zmean = 0.0, zsq = 0.0;
                for (double z : zs) zmean += z;
                zmean /= std::max<std::size_t>(zs.size(), 1);
                for (double z : zs) zsq += (z - zmean) * (z - zmean);
                const double zsd = zs.size() > 1
                    ? std::sqrt(zsq / (zs.size() - 1)) : 0.0;
                std::fprintf(stderr,
                    "[B3 calibration] stage %d round %d M %zu score mean/sd %.6g/%.6g "
                    "emp-null mean/sd/pi0 %.6g/%.6g/%.6g method %s admitted-prefix %d "
                    "lindsey-status %s iter %d mass %.9g curvature %.9g\n",
                    stage, rnd, zs.size(), zmean, zsd,
                    lfdr_result.empirical_mean, lfdr_result.empirical_sd,
                    lfdr_result.pi0,
                    HIST_NULL_GATE ? "histogram-BH"
                        : (CANDIDATE_BH ? "candidate-BH"
                            : (lfdr_result.used_lindsey ? "lindsey"
                                : (lindsey_fallback ? "BH-fallback" : "BH"))), kk,
                    lindsey_status_name(lfdr_result.status),
                    lfdr_result.iterations, lfdr_result.fitted_mass,
                    lfdr_result.central_curvature);
            }
            if (kk == 0) {
                stalls++;
                if (stalls >= 3) break;
                if (stalls == 1) { g_diag_phase = "stall_sweeps"; run_height_sweeps(*mesh, 0.003, 25, timing); continue; }
                // Outer IRLS refresh at the current fit.
                g_diag_phase = "stall_refresh";
                refresh_irls_working_data(*mesh, obs, extra, heldout_split);
                run_height_sweeps(*mesh, 0.003, 25, timing);
                phi_disp = estimate_dispersion(*mesh);
                continue;
            }
            stalls = 0;
            double c_round = 1e300;
            for (int i = 0; i < kk; ++i) c_round = std::min(c_round, std::fabs(zs[order[i]]));
            for (int i = 0; i < kk; ++i) {
                const int candidate_index = order[i];
                Vertex* v = cands[candidate_index];
                const auto& scored_exact = exact_profiles[candidate_index];
                if (v->retired_coefficient || v->disqualified ||
                    (v->active && v->free_coefficient)) continue;
                if (!scored_exact.ok || scored_exact.raw_gain_nats < -1e-8 ||
                    !(scored_exact.gain_nats > 1e-12)) {
                    // B1 invariant: no candidate without a positive exact
                    // penalized gain may cross the admission boundary.
                    continue;
                }

                const bool was_active = v->active;
                g_diag_phase = was_active ? "gate_promote" : "gate_activate";
                const double old_h_diag = v->height;

                // For a previously inactive midpoint, first realize the NVB
                // closure at exactly zero surplus. The statistical decision is
                // then profiled on the *actual* active hierarchical column.
                // This separates irreversible geometry from statistical
                // freedom and avoids trusting an approximate simulated cascade.
                if (!was_active) {
                    v->new_height = v->mu_lin();
                    v->activate(false);
                }

                // Candidates are scored as a batch, but commits are sequential.
                // Earlier admissions can move parents and predictors. Re-profile
                // at the actual commit state, after any conformity closure, so a
                // stale positive batch score can never authorize a negative live
                // coefficient.
                ++g_b1_diag.live_reprofile_attempted;
                const double live_null = v->mu_lin();
                const auto [live_lambda, live_prior] = v->compute_eb_params(1.0);
                const auto exact = v->exact_candidate_gain(
                    live_null, live_prior, live_lambda, false);
                if (!exact.ok || exact.raw_gain_nats < -1e-8 ||
                    !(exact.gain_nats > 1e-12)) {
                    ++g_b1_diag.live_reprofile_rejected;
                    v->gate_admitted = false;
                    continue;
                }
                ++g_b1_diag.live_reprofile_scored;

                const double scored_h_diag = exact.optimum_units;
                if (in_diag_region(*v))
                    log_vertex_event(was_active ? "pre_promote" : "pre_activate",
                                     *v, old_h_diag, scored_h_diag, v->height, false);

                bool accepted = false;
                static const double kMinKeff = [] {
                    const char* e = std::getenv("DMESH_MIN_KEFF");
                    return (e && *e) ? std::atof(e) : 0.0;
                }();
                if (kMinKeff > 0.0 && v->current_keff < kMinKeff) {
                    v->gate_admitted = false;
                    continue;
                }
                double proposal = exact.optimum_units;
                const double data_precision_centi =
                    exact.data_information_optimum / kVarianceScale;
                v->lambda_v = exact.lambda_logit / kVarianceScale;
                v->sigma_sq = data_precision_centi > 0.0
                    ? 1.0 / data_precision_centi : 1e300;
                v->sigma_pooled = data_precision_centi + v->lambda_v > 0.0
                    ? 1.0 / (data_precision_centi + v->lambda_v) : 1e300;
                v->loss_reduction = exact.gain_nats;

                // At this point every candidate is active geometry. Promotion
                // changes only statistical freedom and commits the exact live
                // optimum on the realized column.
                v->retired_coefficient = false;
                accepted = prepare_vertex_proposal(*v, proposal);
                if (accepted) {
                    v->free_coefficient = true;
                    v->dormant_coefficient = false;
                    v->new_height = proposal;
                    v->update_height();
                } else {
                    v->free_coefficient = false;
                    v->new_height = v->height;
                    v->loss_reduction = 0.0;
                    mesh->heap_set(v, 0.0);
                }
                v->gate_admitted = accepted;

                // Optional expensive regression invariant: after activation
                // and its conformity cascade, rebuild the profile on the
                // actual active hierarchical column. It must reproduce the
                // live pre-commit optimum and gain.
                static const bool B1_COMMIT_VALIDATE = [] {
                    const char* e = std::getenv("DMESH_B1_COMMIT_VALIDATE");
                    return e != nullptr && std::atoi(e) != 0;
                }();
                static const bool B1_COMMIT_STRICT = [] {
                    const char* e = std::getenv("DMESH_B1_COMMIT_VALIDATE");
                    return e != nullptr && std::atoi(e) >= 2;
                }();
                if (accepted && B1_COMMIT_VALIDATE) {
                    const auto realized = v->exact_candidate_gain(
                        live_null, live_prior, live_lambda, false);
                    ++g_b1_diag.commit_validated;
                    const double gain_error = realized.ok
                        ? std::fabs(realized.gain_nats - exact.gain_nats)
                        : std::numeric_limits<double>::infinity();
                    const double optimum_error = realized.ok
                        ? std::fabs(realized.optimum_units - v->height)
                        : std::numeric_limits<double>::infinity();
                    g_b1_diag.max_commit_gain_error =
                        std::max(g_b1_diag.max_commit_gain_error, gain_error);
                    g_b1_diag.max_commit_optimum_error_units =
                        std::max(g_b1_diag.max_commit_optimum_error_units,
                                 optimum_error);
                    const double gain_tolerance =
                        2e-7 * (1.0 + std::fabs(exact.gain_nats));
                    const double optimum_tolerance_units = 2e-5;
                    if (!realized.ok || gain_error > gain_tolerance ||
                        optimum_error > optimum_tolerance_units) {
                        ++g_b1_diag.commit_failed;
                        std::fprintf(stderr,
                                     "[B1 commit invariant] FAILED id %d "
                                     "live gain %.12g realized %.12g "
                                     "height %.12g realized optimum %.12g "
                                     "rows %d/%d structural %d/%d was_active %d\n",
                                     v->_id, exact.gain_nats,
                                     realized.gain_nats, v->height,
                                     realized.optimum_units,
                                     exact.fit_points, realized.fit_points,
                                     exact.training_structural_points,
                                     realized.training_structural_points,
                                     was_active ? 1 : 0);
                        if (B1_COMMIT_STRICT) {
                            throw std::runtime_error(
                                "exact candidate commit changed its profiled objective");
                        }
                    }
                }
                if (in_diag_region(*v))
                    log_vertex_event(accepted
                                         ? (was_active ? "post_promote" : "post_activate")
                                         : (was_active ? "post_promote_rejected" : "post_activate_rejected"),
                                     *v, old_h_diag, proposal, v->height, !accepted);
                if (accepted) {
                    v->admit_a = c_round * zsd[candidate_index];
                    v->admit_sd = zsd[candidate_index];
                }
            }
            if (rnd % 3 == 2) {   // periodic re-expansion: keep the target honest
                g_diag_phase = "periodic_refresh";
                refresh_irls_working_data(*mesh, obs, extra, heldout_split);
                run_height_sweeps(*mesh, 0.005, 20, timing);
                phi_disp = estimate_dispersion(*mesh);
            }
        }
        // DMESH_TAX_DEST=1: tax the destination, not the journey. Before the
        // taxed finalization, iterate {relinearize, untaxed sweeps} so raw
        // proposals converge to where the evidence actually points (the
        // bounded working-response step otherwise limits descent to ~1 logit
        // per relinearization at k=0-saturated edges, and per-sweep EB tax
        // turns that into Zeno's descent). The EB tax is then applied to the
        // converged proposal by the existing taxed sweeps.
        static const bool TAX_DEST = std::getenv("DMESH_TAX_DEST") != nullptr;
        if (TAX_DEST && !config.fit.oracle_mode()) {
            const bool eb_saved = mesh->use_eb;
            mesh->use_eb = false;
            for (int cyc = 0; cyc < 6; ++cyc) {
                refresh_irls_working_data(*mesh, obs, extra, heldout_split);
                run_height_sweeps(*mesh, 0.003, 15, timing);
            }
            mesh->use_eb = eb_saved;
            refresh_irls_working_data(*mesh, obs, extra, heldout_split);
        }
        g_diag_phase = "final_sweeps";
        if (!config.fit.oracle_mode()) run_height_sweeps(*mesh, 0.005, 15, timing);
        // Unsupported free coefficients have already been frozen at their
        // last accepted value by Vertex::update_info().  Do not overwrite
        // them here: a direct reset would again violate descendant constraints
        // and would make the final model depend on cleanup order.
        nv_final = 0;
        for (Vertex* v : mesh->vertices) if (v->active) nv_final++;

        // pairwise sigma_u^2 at this stage
        std::vector<double> Sr(pid, 0), Sr2(pid, 0), Sw(pid, 0), Sw2(pid, 0);
        for (int i = 0; i < N; ++i) {
            double g = mesh->predict_index(i);
            double ph = 1 / (1 + std::exp(-(obs[i].offset + g / kHeightScale)));
            double inf = obs[i].exposure * ph * (1 - ph);
            double r = (obs[i].events - obs[i].exposure * ph) / std::sqrt(inf);
            Sr[obs[i].pool_id] += r; Sr2[obs[i].pool_id] += r * r;
            Sw[obs[i].pool_id] += std::sqrt(inf); Sw2[obs[i].pool_id] += inf;
        }
        double nu = 0, de = 0;
        for (int p = 0; p < pid; ++p) { nu += Sr[p] * Sr[p] - Sr2[p]; de += Sw[p] * Sw[p] - Sw2[p]; }
        // Lag-matched attribution. Within-pool pairs sit at month-lags along a
        // trajectory; shared spatial miss contributes to cross-pool pairs at the
        // SAME displacement equally, frailty only within. Bucket obs by
        // (WAC band, WALA month); per band b and lag l:
        //   within: sum_P sum_w Rp[b,P,w] * Rp[b,P,w+l]
        //   cross : sum_w  R[b,w] * R[b,w+l]  -  (within part)
        // idio = weighted avg over lags of (within/dw - cross/dc).
        double sig_shared = 0, sig_idio = 0, sig_white = 0;
        double sr2_d = 0, sw2_d = 0; long n_d = 0;
        {
            const int NB2 = config.diagnostics.lag_match_bands;
            const int MAXLAG = config.diagnostics.max_lag;
            // maps keyed by (band, month) and (band, pool, month)
            std::map<std::pair<int,int>, std::pair<double,double>> cell;      // R, S
            std::map<std::tuple<int,int,int>, std::pair<double,double>> pcell;
            for (int i = 0; i < N; ++i) {
                double g = mesh->predict_index(i);
                double ph = 1 / (1 + std::exp(-(obs[i].offset + g / kHeightScale)));
                double inf = obs[i].exposure * ph * (1 - ph);
                double r = (obs[i].events - obs[i].exposure * ph) / std::sqrt(inf);
                double sq = std::sqrt(inf);
                double gx = sx0.empty() ? sx[i] : sx0[i];
                double gy = sy0.empty() ? sy[i] : sy0[i];
                int b = std::clamp((int)(gy * NB2), 0, NB2 - 1);
                int w = (int)std::lround(gx * 60.0);
                auto& c = cell[{b, w}]; c.first += r; c.second += sq;
                auto& pc = pcell[{b, obs[i].pool_id, w}]; pc.first += r; pc.second += sq;
                sr2_d += r * r; sw2_d += inf; n_d += 1;
            }
            // Zero-lag contrast: diagonal (r^2 - 1, the within-pool lag-0 'pair')
            // minus co-located cross-pool moment. Co-located pairs share the
            // spatial miss AND the local smoothing error of g-hat, so the
            // subtraction cancels both to first order: no trace, no edf, no
            // selection accounting. contrast(0) - contrast(l>=1) = sigma_w^2.
            double c0_within = 0, c0_cross = 0;
            {
                double nw0 = 0, dw0 = 0, nc0 = 0, dc0 = 0;
                for (int i = 0; i < N; ++i) {
                    double g = mesh->predict_index(i);
                    double ph = 1 / (1 + std::exp(-(obs[i].offset + g / kHeightScale)));
                    double inf = obs[i].exposure * ph * (1 - ph);
                    double r = (obs[i].events - obs[i].exposure * ph) / std::sqrt(inf);
                    nw0 += r * r - 1.0; dw0 += inf;
                }
                for (auto& kv : cell) {
                    nc0 += kv.second.first * kv.second.first;
                    dc0 += kv.second.second * kv.second.second;
                }
                // remove same-pool (diagonal) part from the cell products
                for (auto& kv : pcell) {
                    nc0 -= kv.second.first * kv.second.first;
                    dc0 -= kv.second.second * kv.second.second;
                }
                c0_within = dw0 > 1e-3 ? nw0 / dw0 : 0.0;
                c0_cross  = dc0 > 1e-3 ? nc0 / dc0 : 0.0;
            }
            double num_i = 0, den_i = 0;
            std::vector<double> lag_val, lag_idio, lag_w;
            for (int l = 1; l <= MAXLAG; ++l) {
                double nw = 0, dw = 0, nc = 0, dc = 0;
                for (auto& kv : pcell) {
                    auto [b, P, w] = kv.first;
                    auto it = pcell.find({b, P, w + l});
                    if (it == pcell.end()) continue;
                    nw += kv.second.first * it->second.first;
                    dw += kv.second.second * it->second.second;
                }
                for (auto& kv : cell) {
                    auto [b, w] = kv.first;
                    auto it = cell.find({b, w + l});
                    if (it == cell.end()) continue;
                    nc += kv.second.first * it->second.first;
                    dc += kv.second.second * it->second.second;
                }
                nc -= nw; dc -= dw;                     // cross-pool only
                if (dw < 1e-3 || dc < 1e-3) continue;
                num_i += nw - dw * (nc / dc);           // within minus matched shared
                den_i += dw;
                sig_shared += dw * (nc / dc);           // shared profile, dw-weighted
                lag_val.push_back(l); lag_idio.push_back(nw / dw - nc / dc); lag_w.push_back(dw);
            }
            sig_idio = den_i > 1e-3 ? num_i / den_i : 0.0;
            sig_shared = den_i > 1e-3 ? sig_shared / den_i : 0.0;
            // White-plus-permanent model: C(l) = sigma_a^2 flat for all l >= 1,
            // so the lag-matched flat level IS sigma_a^2 (no attenuation, and the
            // Sec. 9.6 leverage correction applies to it directly). The white
            // component is the diagonal contrast: E[r^2]-1 per unit info reads
            // (miss + sigma_a^2 + sigma_w^2); subtracting the pairwise level and
            // the shared (spatial-miss) term leaves sigma_w^2, which absorbs any
            // residual sub-resolution miss (conservative, per the attribution
            // semantics of the writeup).
            double diag_excess = sw2_d > 1e-9 ? (sr2_d - (n_d - nv_final)) / sw2_d : 0.0;
            sig_white = std::max(diag_excess - std::max(sig_idio, 0.0)
                                 - std::max(sig_shared, 0.0), 0.0);
            double sig_w0 = std::max((c0_within - c0_cross) - sig_idio, 0.0);
            std::fprintf(stderr, "[pool state] sigma_a^2 %.4f (persistent) + sigma_w^2 %.4f (white) | "
                         "diag excess %.4f, shared %.4f\n",
                         sig_idio, sig_white, diag_excess, sig_shared);
            std::fprintf(stderr, "[sigma_w zero-lag] contrast0 %.4f (within %.4f, cross %.4f) - lag>=1 %.4f "
                         "=> sigma_w^2 %.4f\n",
                         c0_within - c0_cross, c0_within, c0_cross, sig_idio, sig_w0);
        }

        if (config.diagnostics.oracle_shared) {
            double nu_o = 0, de_o = 0;
            std::map<int, double> pr, pw;
            for (Face* fc : mesh->active_faces) {
                if (fc->coords_indices.size() < 2) continue;
                pr.clear(); pw.clear();
                double fr = 0, fw = 0;
                for (int i : fc->coords_indices) {
                    double ph = 1 / (1 + std::exp(-(obs[i].offset + obs[i].true_surface_scaled / kHeightScale + obs[i].latent_effect)));
                    double inf = obs[i].exposure * ph * (1 - ph);
                    double r = (obs[i].events - obs[i].exposure * ph) / std::sqrt(inf);
                    double sq = std::sqrt(inf);
                    fr += r; fw += sq;
                    pr[obs[i].pool_id] += r; pw[obs[i].pool_id] += sq;
                }
                double sr = 0, sw = 0;
                for (auto& kv : pr) sr += kv.second * kv.second;
                for (auto& kv : pw) sw += kv.second * kv.second;
                nu_o += fr * fr - sr; de_o += fw * fw - sw;
            }
            std::fprintf(stderr, "[oracle shared] %.4f (nu %.3g de %.3g)\n",
                         de_o > 1e-3 ? nu_o / de_o : 0.0, nu_o, de_o);
        }
        // one-obs-per-pool data has no within-pool pairs: de == 0 and sigma_u
        // is not identifiable; hold at the floor rather than propagating NaN.
        sig_cur = (de > 1e-3 && std::isfinite(nu / de)) ? std::max(nu / de, 1e-4) : 1e-4;
        std::printf("stage %d: vertices %d, sigma_u^2 %.4f (shared %.4f, idio %.4f)  [%.1f s]\n",
                    stage + 1, nv_final, sig_cur, sig_shared, sig_idio, elapsed());
        if (config.fit.oracle_mode()) break;   // oracle is single-stage: later
                                                  // stages' DL weights corrupt the projection
        // Once a post-pilot stage maps sigma_u^2 to exactly the same value,
        // the next stage is a deterministic repeat: it receives identical
        // weights, starts from the same coarse mesh, and follows the same
        // deterministic gate/solver path.  Stop rather than rebuilding it.
        if (stage >= 1 && stage + 1 < n_stages &&
            sig_cur == stage_input_sigma) {
            std::printf("stage convergence: sigma_u^2 unchanged; skipping %d repeated stage(s)\n",
                        n_stages - stage - 1);
            break;
        }
        if (use_metric && stage == 0) {
            // directional curvature energy on a grid, weighted by data density
            const int GG = 96;
            std::vector<double> F(GG * GG), dens(GG * GG, 0.0);
            for (int a = 0; a < GG; ++a) for (int b = 0; b < GG; ++b)
                F[a * GG + b] = mesh->predict((a + 0.5) / GG, (b + 0.5) / GG);
            for (int i = 0; i < N; ++i) {
                int a = std::clamp((int)(sx[i] * GG), 0, GG - 1);
                int b = std::clamp((int)(sy[i] * GG), 0, GG - 1);
                dens[a * GG + b] += mesh->weights[i];
            }
            std::vector<double> ex(MB, 0.0), ey(MB, 0.0);
            for (int a = 1; a < GG - 1; ++a) for (int b = 1; b < GG - 1; ++b) {
                double w = dens[a * GG + b];
                if (w <= 0) continue;
                double fxx = F[(a+1)*GG+b] - 2*F[a*GG+b] + F[(a-1)*GG+b];
                double fyy = F[a*GG+b+1] - 2*F[a*GG+b] + F[a*GG+b-1];
                ex[std::min(a * MB / GG, MB - 1)] += w * fxx * fxx;
                ey[std::min(b * MB / GG, MB - 1)] += w * fyy * fyy;
            }
            // equidistribute sqrt(sqrt(energy)) (P1 error ~ h^2 |f''| -> h ~ |f''|^{-1/2},
            // so knot density ~ |f''|^{1/2} ~ energy^{1/4}); floor keeps empty bands mapped
            auto equi = [&](std::vector<double>& e, std::vector<double>& phi) {
                double mx = 1e-12;
                for (double v : e) mx = std::max(mx, v);
                double tot = 0;
                std::vector<double> d(MB);
                for (int b = 0; b < MB; ++b) { d[b] = std::pow(e[b] / mx, 0.25) + 0.08; tot += d[b]; }
                phi[0] = 0;
                for (int b = 0; b < MB; ++b) phi[b + 1] = phi[b] + d[b] / tot;
            };
            equi(ex, phiX); equi(ey, phiY);
            sx0 = sx; sy0 = sy;   // originals for dumps/figures
            for (int i = 0; i < N; ++i) { sx[i] = apply_map(phiX, sx[i]); sy[i] = apply_map(phiY, sy[i]); }
            std::fprintf(stderr, "[metric] axis maps learned; y-band widths min/max %.3f/%.3f\n",
                         *std::min_element(phiY.begin()+1, phiY.end()) - phiY[0],
                         phiY[MB] - phiY[MB-1]);
        }
    }
    if (config.diagnostics.stop_after_adaptation) {
        // Diagnostic-only exit used for null-mechanism decomposition.  The
        // selected topology and gate-admission ledger are already determined;
        // skip fixed-topology retirement, exact final refits, leverage, and PM
        // diagnostics so adaptive amplification can be measured cheaply at a
        // predeclared round/stage boundary.
        refresh_hierarchical_support_diagnostics(*mesh, config.fit.heldout_split);
        write_figure_dumps(config.output, *mesh, obs, sx, sy, sx0, sy0,
                           phiX, phiY, surface, pid, sig_cur);
        std::fprintf(stderr,
                     "[diagnostic stop] after adaptation: stages %d max-rounds %d\n",
                     n_stages, config.diagnostics.max_gate_rounds);
        std::printf("TOTAL %.1f s (diagnostic stop after adaptation)\n", elapsed());
        return 0;
    }
    // --- Final post-selection re-shrinkage at fixed topology (DMESH_FINAL_RESHRINK=1).
    // Inside the refinement loop tau is load-bearing damping: lowering the floor
    // mid-loop triggers a shrink -> residual-signal -> re-admission spiral
    // (measured: floor sd 0.01 grows the mesh 2.2x and worsens held-out deviance
    // 20.7 -> 27.8). After the gate closes that feedback channel is gone, so the
    // honest prior is safe. Recover RAW surpluses for admitted vertices from
    // local regressions at the final linearization (undoing shrinkage without a
    // global refit), re-estimate tau per depth with the truncation correction,
    // zero-mean chi^2, floored at the configured tau floor (a tau = 0 moment
    // read means "shrink hard", not "weld": the historical no-floor variant
    // produced lambda_v = 1e8 wholesale welds at whole depths and is kept
    // only under DMESH_RESHRINK_NO_FLOOR=1), decay-borrow sparse depths at 4^{-dd}
    // (the C^2 surplus scaling), then re-sweep heights once at fixed topology.
    // Raw-surplus variance omits the 0.25*sp parent term that admit_a includes,
    // which overstates alpha slightly and biases tau DOWN: conservative.
    if (std::getenv("DMESH_FINAL_RESHRINK")) {
        struct RawSurplus { double draw, svar, athr, keff; Vertex* v; };
        std::map<int, std::vector<RawSurplus>> by_depth;
        std::map<int, std::vector<double>> before_prof;
        for (Vertex* v : mesh->vertices) {
            if (!v->active || !v->gate_admitted || v->dormant_coefficient ||
                v->retired_coefficient || !v->parent_edge) continue;
            before_prof[v->depth].push_back(std::fabs(v->height - v->mu_prior()));
            auto lr = v->loc_regress(false);
            if (!(lr.xTx > 0) || !(v->admit_sd > 0)) continue;
            // variance = admission-time calibrated sd^2: alpha(tau=0) is then the
            // actual selection cutoff, making the truncation term self-consistent
            by_depth[v->depth].push_back(
                {lr.xTr / lr.xTx - v->mu_prior(), v->admit_sd * v->admit_sd,
                 v->admit_a, v->current_keff, v});
        }
        auto kap = [](double alpha) {
            if (alpha <= 1e-12) return 1.0;
            double ph = std::exp(-0.5 * alpha * alpha) / 2.5066282746310002;
            double Q = 0.5 * std::erfc(alpha / 1.4142135623730951);
            if (Q < 1e-300) return alpha * alpha + 1.0;
            return 1.0 + alpha * ph / Q;
        };
        // Truncation-corrected moment solve; w = nullptr for the unweighted
        // per-depth read, or per-surplus responsibilities for the two-groups
        // M-step. Returns 0 when the (weighted) chi^2 sits below its
        // selection-inflated null expectation.
        auto solve_tau = [&](const std::vector<RawSurplus>& rs,
                             const std::vector<double>* w) {
            double chi = 0.0;
            for (size_t i = 0; i < rs.size(); ++i)
                chi += (w ? (*w)[i] : 1.0) * rs[i].draw * rs[i].draw / rs[i].svar;
            auto target = [&](double t2) {
                double acc = 0.0;
                for (size_t i = 0; i < rs.size(); ++i) {
                    const auto& r = rs[i];
                    double al = r.athr > 0 ? r.athr / std::sqrt(t2 + r.svar) : 0.0;
                    acc += (w ? (*w)[i] : 1.0) * kap(al) * (t2 + r.svar) / r.svar;
                }
                return acc;
            };
            if (chi <= target(0.0)) return 0.0;
            double lo = 0.0, hi = 1.0;
            while (target(hi) < chi && hi < 1e8) hi *= 2;
            for (int it = 0; it < 60; ++it) {
                double mid = 0.5 * (lo + hi);
                if (target(mid) < chi) lo = mid; else hi = mid;
            }
            return 0.5 * (lo + hi);
        };
        std::map<int, double> tau_new;
        for (auto& [d, rs] : by_depth) {
            if ((int)rs.size() < 3) continue;
            tau_new[d] = solve_tau(rs, nullptr);
        }
        if (!tau_new.empty()) {
            for (auto& [d, rs] : by_depth) {
                if (tau_new.count(d)) continue;
                int nearest = tau_new.begin()->first;
                for (auto& [d0, _] : tau_new)
                    if (std::abs(d0 - d) < std::abs(nearest - d)) nearest = d0;
                tau_new[d] = tau_new[nearest] * std::pow(4.0, -(double)(d - nearest));
            }
            // Improvement-only direction: never shrink LESS than the loop did.
            // The loop's (pooled, floored) tau under-reads at signal depths and
            // that under-read is load-bearing regularization there; the honest
            // moment's contribution is where it reads BELOW the loop, i.e. where
            // raw-surplus evidence says a depth is thinner than the floor let it be.
            for (auto& [d, t] : tau_new) t = std::min(t, mesh->tau_for_depth(d));
            for (auto& [d, rs] : by_depth) {
                double old_tau = mesh->tau_for_depth(d);
                auto& bp = before_prof[d];
                std::nth_element(bp.begin(), bp.begin() + bp.size() / 2, bp.end());
                std::fprintf(stderr,
                             "[reshrink] depth %2d: m=%4zu tau %.1f -> %.1f  median|surplus| %.2f\n",
                             d, rs.size(), old_tau, tau_new[d], bp[bp.size() / 2]);
            }
            // --- Ladder variants (DMESH_RESHRINK_LADDER=cells|twogroups). ---
            // Depth is the wrong exchangeability unit when one tier spans
            // starved corridor vertices and well-supported workhorses: a
            // single pooled moment read then shrinks the whole tier to one
            // tau (writeup sec:aug6-ladder, EB audit defect 1). Both variants
            // keep the truncation-corrected moment machinery and the
            // improvement-only clamp against the loop tau, but assign
            // per-vertex prior variances (Vertex::tau_override_sq) consumed
            // by compute_eb_params in the fixed-topology re-solve below.
            // Vertices outside the raw-surplus set keep the pooled depth tau.
            const char* lme = std::getenv("DMESH_RESHRINK_LADDER");
            const std::string ladder_mode = lme ? lme : "depth";
            if (ladder_mode == "cells") {
                // Evidence-geometry cells: split each depth at a K_eff
                // threshold and give each side its own moment read, so a
                // noise-dominated starved cell cannot weld the workhorse
                // cell (and vice versa).
                static const double KEFF_SPLIT = [] {
                    const char* e = std::getenv("DMESH_TAU_CELL_KEFF");
                    return e ? std::atof(e) : 1.0;
                }();
                for (auto& [d, rs] : by_depth) {
                    std::vector<RawSurplus> lo_cell, hi_cell;
                    for (auto& r : rs)
                        (r.keff < KEFF_SPLIT ? lo_cell : hi_cell).push_back(r);
                    const double pooled = std::max(tau_new[d], mesh->tau_floor_sq);
                    double cell_tau[2];
                    const std::vector<RawSurplus>* cells[2] = {&lo_cell, &hi_cell};
                    for (int c = 0; c < 2; ++c) {
                        if ((int)cells[c]->size() < 3) {
                            cell_tau[c] = pooled;
                            continue;
                        }
                        double t = std::min(solve_tau(*cells[c], nullptr),
                                            mesh->tau_for_depth(d));
                        cell_tau[c] = std::max(t, mesh->tau_floor_sq);
                    }
                    for (int c = 0; c < 2; ++c)
                        for (auto& r : *cells[c])
                            r.v->tau_override_sq = cell_tau[c];
                    std::fprintf(stderr,
                                 "[laddercells] depth %2d: keff<%.2f m=%zu tau %.1f | "
                                 "keff>=%.2f m=%zu tau %.1f\n",
                                 d, KEFF_SPLIT, lo_cell.size(), cell_tau[0],
                                 KEFF_SPLIT, hi_cell.size(), cell_tau[1]);
                }
            } else if (ladder_mode == "twogroups") {
                // Spike-and-slab: per depth, a two-component EM on the raw
                // surpluses -- spike at the configured floor, slab variance
                // from the responsibility-weighted truncation-corrected
                // moment, selection entering each component's likelihood
                // through the admission normalizer P(|draw| clears athr).
                // Each coefficient gets the posterior-mixture prior variance
                // gamma*tau1 + (1-gamma)*tau0: slab-probable coefficients
                // keep their evidence, spike-probable ones shrink to the
                // floor, and no single pooled read can weld a mixed tier.
                const double t0 = std::max(mesh->tau_floor_sq, 1e-8);
                auto logtrunc = [](double x, double v, double a) {
                    double q = 0.5 * std::erfc((a > 0.0 ? a : 0.0) /
                                               std::sqrt(v) / 1.4142135623730951);
                    if (q < 1e-300) q = 1e-300;
                    return -0.5 * std::log(v) - 0.5 * x * x / v - std::log(2.0 * q);
                };
                // DMESH_RESHRINK_2G_FREE=1 lifts the improvement-only cap on
                // the slab variance: the responsibility-weighted moment may
                // then exceed the loop tau, letting slab-probable
                // coefficients shrink LESS than the loop did. The pooled
                // clamp's rationale (the loop under-read is load-bearing at
                // signal depths) was measured for the all-or-nothing pooled
                // estimator; with per-coefficient responsibilities the spike
                // still absorbs the noise mass.
                static const bool TG_FREE =
                    std::getenv("DMESH_RESHRINK_2G_FREE") != nullptr;
                for (auto& [d, rs] : by_depth) {
                    if ((int)rs.size() < 4) continue;  // depth tau stands
                    const double old_tau = mesh->tau_for_depth(d);
                    const double tau1_cap = TG_FREE ? 1e8 : old_tau;
                    double mean_pos = 0.0;
                    for (auto& r : rs)
                        mean_pos += std::max(r.draw * r.draw - r.svar, 0.0);
                    mean_pos /= rs.size();
                    double tau1 = std::clamp(std::max(tau_new[d], mean_pos),
                                             4.0 * t0, tau1_cap);
                    if (tau1 <= 2.0 * t0) continue;  // no separable slab
                    double pi = 0.5;
                    std::vector<double> gam(rs.size(), 0.5);
                    for (int it = 0; it < 30; ++it) {
                        for (size_t i = 0; i < rs.size(); ++i) {
                            const auto& r = rs[i];
                            double ls = std::log(pi) +
                                        logtrunc(r.draw, tau1 + r.svar, r.athr);
                            double l0 = std::log(1.0 - pi) +
                                        logtrunc(r.draw, t0 + r.svar, r.athr);
                            double mx = std::max(ls, l0);
                            gam[i] = std::exp(ls - mx) /
                                     (std::exp(ls - mx) + std::exp(l0 - mx));
                        }
                        double s = 0.0;
                        for (double g : gam) s += g;
                        pi = std::clamp(s / rs.size(), 0.01, 0.99);
                        tau1 = std::clamp(solve_tau(rs, &gam), 4.0 * t0, tau1_cap);
                    }
                    int nslab = 0;
                    for (size_t i = 0; i < rs.size(); ++i) {
                        double tv = gam[i] * tau1 + (1.0 - gam[i]) * t0;
                        rs[i].v->tau_override_sq = std::max(tv, mesh->tau_floor_sq);
                        if (gam[i] > 0.5) ++nslab;
                    }
                    std::fprintf(stderr,
                                 "[ladder2g] depth %2d: m=%4zu pi=%.2f tau1=%.1f "
                                 "slab(g>.5)=%d\n",
                                 d, rs.size(), pi, tau1, nslab);
                }
            }
            mesh->tau_sq = tau_new;
            // The configured floor (DMESH_TAU_FLOOR_LOGIT_SD) stays in force:
            // a depth whose truncation-corrected moment reads tau = 0 shrinks
            // to the floor, not to a 1e-8 hard weld (lambda_v ~ 1e8, T = 0),
            // which silently freezes every admitted coefficient at that depth
            // regardless of its individual information. DMESH_RESHRINK_NO_FLOOR=1
            // restores the historical no-floor behavior for reproducibility.
            if (std::getenv("DMESH_RESHRINK_NO_FLOOR"))
                mesh->tau_floor_sq = 1e-8;
            run_height_sweeps(*mesh, 0.0005, 120, timing);
            std::map<int, std::vector<double>> after_prof;
            for (Vertex* v : mesh->vertices)
                if (v->active && v->gate_admitted && !v->dormant_coefficient &&
                    !v->retired_coefficient && v->parent_edge)
                    after_prof[v->depth].push_back(std::fabs(v->height - v->mu_prior()));
            for (auto& [d, ap] : after_prof) {
                std::vector<double> a2 = ap;
                std::nth_element(a2.begin(), a2.begin() + a2.size() / 2, a2.end());
                std::fprintf(stderr, "[reshrink] depth %2d: median|surplus| after %.2f\n",
                             d, a2[a2.size() / 2]);
            }
        } else {
            std::fprintf(stderr, "[reshrink] no estimable depth: pass skipped\n");
        }
    }

    // Finalize exact obsolescence only after all adaptive topology decisions
    // and optional fixed-topology re-shrinkage are complete. No admission or
    // refinement pass follows this cleanup.
    if (!std::getenv("DMESH_DISABLE_FINAL_RETIREMENT")) {
        retire_zero_support_at_fixed_topology(*mesh, config.fit.heldout_split, timing);
    }

    {
        // DMESH_HIER_FIT=2 historically requested a final hierarchical WLS
        // compatibility pass. Under the Phase-A exact solver that pass is
        // redundant and can move the fit away from the exact likelihood fixed
        // point. Preserve the flag as a request for a final frozen-topology
        // hierarchical refit, but execute it in the same exact objective and
        // columns as the stage solver.
        const char* e = std::getenv("DMESH_HIER_FIT");
        if (e && std::atoi(e) == 2) {
            refresh_irls_working_data(*mesh, obs, mesh->last_extra,
                                      config.fit.heldout_split);
            // Two-currency experiment (DMESH_LAW_WEIGHTS_FINAL): topology and
            // gate were built on binomial-information weights; re-temper the
            // frozen-topology final height solve by the pool law (GLS), so the
            // law discounts fitting information exactly once, here.
            if (std::getenv("DMESH_LAW_WEIGHTS_FINAL") && !mesh->gate_s2.empty()) {
                for (int i = 0; i < mesh->n_points; ++i) {
                    const double eta = mesh->eta_frozen[i];
                    const double pr = std::clamp(1.0 / (1.0 + std::exp(-eta)), 1e-4, 1.0 - 1e-4);
                    const double bw = obs[i].exposure * pr * (1.0 - pr) / kVarianceScale;
                    const double law_extra = kVarianceScale * mesh->gate_s2[i];
                    double w = (bw > 0.0) ? 1.0 / (1.0 / bw + law_extra) : 0.0;
                    if (config.fit.heldout_split &&
                        !hierarchical_design::is_training_row(*mesh, i, config.fit.heldout_split))
                        w = 1e-12;
                    mesh->weights[i] = w;
                    mesh->likelihood_temper[i] = (bw > 0.0)
                        ? std::clamp(w / bw, 0.0, 1.0) : 0.0;
                }
                for (Face* fc : mesh->active_faces) fc->refresh_suff_stats_from_cache();
                for (Vertex* v : mesh->vertices)
                    if (v->active || v->parent_edge) v->update_info();
                std::printf("law-weights-final: frozen-topology GLS temper applied\n");
            }
            const char* p = std::getenv("DMESH_PL_SOLVER");
            const bool exact_line = p != nullptr && std::atoi(p) != 0;
            if (exact_line) {
                int final_cycles = 30;
                if (const char* fc = std::getenv("DMESH_B1_FINAL_CYCLES"))
                    final_cycles = std::max(1, std::atoi(fc));
                for (int cycle = 0; cycle < final_cycles; ++cycle)
                    pl_solver::run(*mesh, 0.0005, 120);
            }
            else hier_fit::run(*mesh, 0.5, 60);
        }
    }
    // Restore lifecycle flags and information in the fitted hierarchical
    // coordinate system before dumps or downstream estimators consume them.
    refresh_hierarchical_support_diagnostics(*mesh, config.fit.heldout_split);
    run_b1_fixed_point_audit(*mesh);
    write_edf_dump(config.output, *mesh, N);

    // The fitted topology, face memberships, and final weights are all that the
    // diagnostics need. Release duplicate training coordinates and responses
    // before allocating the leverage/PM workspaces.
    std::vector<double>().swap(mesh->X);
    std::vector<double>().swap(mesh->values);

    // final estimators
    double sumr2 = 0, suminfo = 0, energy_truth_num = 0, lvl_num = 0;
    for (int i = 0; i < N; ++i) {
        double g = mesh->predict_index(i);
        double ph = 1 / (1 + std::exp(-(obs[i].offset + g / kHeightScale)));
        double inf = obs[i].exposure * ph * (1 - ph);
        double r = (obs[i].events - obs[i].exposure * ph) / std::sqrt(inf);
        sumr2 += r * r; suminfo += inf;
        double D = (obs[i].true_surface_scaled - g) / kHeightScale + obs[i].latent_effect;
        energy_truth_num += inf * D * D;
        lvl_num += std::sqrt(inf) * r;
    }
    double energy = (sumr2 - (N - nv_final)) / suminfo;
    std::printf("energy %.4f (truth %.4f) | sigma_u^2 %.4f (truth %.4f) | level %+.5f\n",
                energy, energy_truth_num / suminfo, sig_cur, sigma_u * sigma_u, lvl_num / suminfo);
    double t_before_corr = elapsed();
    double dsr = 1, cnr = 0, sig_corr = sig_cur;
    if (sig_cur > 5e-4) {
        // Only correct when there is frailty to correct. At the floor the DL
        // weights are ~binomial, the mesh refines fully into surface signal, and
        // V (hence the O(V^2) probe solves) is large for no benefit.
        sig_corr = leverage_correction(mesh.get(), obs, sx, sy, pid, sig_cur, dsr, cnr, config.diagnostics.leverage_probes);
    }
    std::printf("leverage-corrected sigma_u^2 %.4f (naive %.4f, truth %.4f) | "
                "D_sig/D_naive %.3f, C_noise/T %+.4f  [corr %.1f s]\n",
                sig_corr, sig_cur, sigma_u * sigma_u, dsr, cnr, elapsed() - t_before_corr);
    {
        double t_pm = elapsed();
        double pm_tot = pm_sigma_total(mesh.get(), obs, sx, sy);
        // Zero-lag miss moment: co-located (same-face) cross-pool pairs share the
        // local miss m(z) and nothing else (no a_P: different pools; no eta: i.i.d.),
        // so Z = kappa * M^2 with kappa the within-face miss coherence (15/16 for
        // the isotropic P1 interpolation bubble; regime-measurable at zero signal).
        double Z = 0;
        // Per-size-bin pool variance sigma2(n): within each face the local surface
        // bias is COMMON across pools; the face-level cross-pool moment measures it
        // (m2f), and subtracting w_i*m2f from each bin's diagonal excess isolates
        // the pool-level variance bias-free. Bins: log2(loans), clamped 5..13.
        double bin_num[14] = {0}, bin_den[14] = {0}; long bin_pools[14] = {0};
        {
            double nu = 0, de = 0;
            std::map<int, double> pr, pw;
            for (Face* fc : mesh->active_faces) {
                if (fc->coords_indices.size() < 2) continue;
                pr.clear(); pw.clear();
                double fr = 0, fw = 0;
                for (int i : fc->coords_indices) {
                    double g = mesh->predict_index(i);
                    double ph = std::clamp(1 / (1 + std::exp(-(obs[i].offset + g / kHeightScale))), 1e-9, 1 - 1e-9);
                    double inf = obs[i].exposure * ph * (1 - ph);
                    double r = (obs[i].events - obs[i].exposure * ph) / std::sqrt(inf);
                    double sq = std::sqrt(inf);
                    fr += r; fw += sq;
                    pr[obs[i].pool_id] += r; pw[obs[i].pool_id] += sq;
                }
                double sr = 0, sw = 0;
                for (auto& kv : pr) sr += kv.second * kv.second;
                for (auto& kv : pw) sw += kv.second * kv.second;
                nu += fr * fr - sr;
                // local shared (surface-bias) level for this face, per unit info
                double de_f = fw * fw - sw;
                double m2f = de_f > 1e-3 ? std::max((fr * fr - sr) / de_f, 0.0) : 0.0;
                for (int i : fc->coords_indices) {
                    double g = mesh->predict_index(i);
                    double ph = std::clamp(1 / (1 + std::exp(-(obs[i].offset + g / kHeightScale))), 1e-9, 1 - 1e-9);
                    double inf = obs[i].exposure * ph * (1 - ph);
                    double r = (obs[i].events - obs[i].exposure * ph) / std::sqrt(inf);
                    int b = std::clamp((int)std::floor(std::log2((double)std::max(obs[i].exposure, 1))), 5, 13);
                    bin_num[b] += r * r - 1.0 - inf * m2f;
                    bin_den[b] += inf; bin_pools[b]++;
                }
                de += fw * fw - sw;
            }
            Z = de > 1e-3 ? nu / de : 0.0;
        }
        for (int b = 5; b <= 13; ++b)
            if (bin_pools[b] > 200)
                std::fprintf(stderr, "[poolvar] n in [2^%d,2^%d): pools %ld  sigma2(n) %.4f\n",
                             b, b + 1, bin_pools[b], bin_num[b] / std::max(bin_den[b], 1e-9));
        std::printf("PM: total %.4f | Z(zero-lag miss, kappa-scaled) %.4f | "
                    "sigma_w^2 = total - a_corr - Z = %.4f  [pm %.1f s]\n",
                    pm_tot, Z, std::max(pm_tot - sig_corr - std::max(Z, 0.0), 0.0),
                    elapsed() - t_pm);
    }
    if (config.fit.heldout_split) {
        double dev = 0, x2 = 0, wsum = 0; long m = 0;
        for (int i = 0; i < N; ++i) {
            if (hierarchical_design::is_training_row(*mesh, i, true)) continue;
            double g = mesh->predict_index(i);
            double ph = std::clamp(1 / (1 + std::exp(-(obs[i].offset + g / kHeightScale))), 1e-6, 1 - 1e-6);
            double k = obs[i].events, n = obs[i].exposure;
            if (k > 0) dev += 2 * k * std::log(k / (n * ph));
            if (n - k > 0) dev += 2 * (n - k) * std::log((n - k) / (n * (1 - ph)));
            double inf = n * ph * (1 - ph);
            x2 += (k - n * ph) * (k - n * ph) / inf; wsum += inf; m++;
        }
        std::printf("HELDOUT: mean deviance/pool %.3f | X2/info %.4f | %ld pools\n",
                    dev / m, x2 / wsum, m);
    }
    if (!config.simulation.data_file) {
        // Empirical FDP: an admitted vertex is FALSE if the true hierarchical
        // surplus of f at v (relative to its parent-edge endpoints) is below
        // eps. Only lfdr-gate admissions count; pre-refinement does not.
        auto dtrue = [&](Vertex* v) {
            auto F = [&](double gx, double gy) {
                return from_log_odds(surface.evaluate(-4 + 8 * gx, -4 + 8 * gy));
            };
            Vertex* a = v->parent_edge->vertex0; Vertex* b = v->parent_edge->vertex1;
            return F(v->x, v->y) - 0.5 * (F(a->x, a->y) + F(b->x, b->y));
        };
        int G = 0; int F1 = 0, F05 = 0, F2 = 0;
        for (Vertex* v : mesh->vertices) {
            if (!v->active || !v->gate_admitted || !v->parent_edge) continue;
            G++;
            double d = std::fabs(dtrue(v));
            if (d < 1.0) F1++;
            if (d < 0.5) F05++;
            if (d < 2.0) F2++;
        }
        std::printf("FDP: admitted %d, false %d/%d/%d (eps 0.5/1.0/2.0), FDP %.3f/%.3f/%.3f\n",
                    G, F05, F1, F2, G ? (double)F05 / G : 0.0, G ? (double)F1 / G : 0.0,
                    G ? (double)F2 / G : 0.0);
    }
    if (config.diagnostics.run_height_test) {
        int shown = 0;
        for (Vertex* v : mesh->vertices) {
            if (!v->active || !v->parent_edge || v->depth < 6) continue;
            std::fprintf(stderr, "HT depth %d height %.2f predict %.2f\n",
                         v->depth, v->height, mesh->predict(v->x, v->y));
            if (++shown >= 6) break;
        }
    }
    std::fprintf(stderr,
                 "[B1 exact candidate gains] attempted %lld scored %lld "
                 "unscoreable %lld (no-training-column %lld, tempered-out %lld) "
                 "negative %lld grid_checked %lld grid_failed %lld "
                 "max_grid_excess %.3g | live-reprofile %lld/%lld rejected %lld "
                 "commit-check %lld failed %lld max gain/height error %.3g/%.3g\n",
                 g_b1_diag.attempted, g_b1_diag.scored,
                 g_b1_diag.unscoreable,
                 g_b1_diag.unscoreable_no_training_column,
                 g_b1_diag.unscoreable_tempered_out, g_b1_diag.negative,
                 g_b1_diag.grid_checked, g_b1_diag.grid_failed,
                 g_b1_diag.max_grid_excess,
                 g_b1_diag.live_reprofile_scored,
                 g_b1_diag.live_reprofile_attempted,
                 g_b1_diag.live_reprofile_rejected,
                 g_b1_diag.commit_validated,
                 g_b1_diag.commit_failed,
                 g_b1_diag.max_commit_gain_error,
                 g_b1_diag.max_commit_optimum_error_units);
    if (g_b1_diag.nearest_id >= 0) {
        std::fprintf(stderr,
                     "[B1 regime-straddle nearest] stage %d round %d id %d "
                     "xy=(%.8f,%.8f) distance %.6g null %.6f optimum %.6f "
                     "gain %.9g nats\n",
                     g_b1_diag.nearest_stage, g_b1_diag.nearest_round,
                     g_b1_diag.nearest_id, g_b1_diag.nearest_x,
                     g_b1_diag.nearest_y, g_b1_diag.nearest_regime_distance,
                     g_b1_diag.nearest_null / kHeightScale,
                     g_b1_diag.nearest_optimum / kHeightScale,
                     g_b1_diag.nearest_gain);
    }
    std::printf("TOTAL %.1f s\n", elapsed());

    write_figure_dumps(config.output, *mesh, obs, sx, sy, sx0, sy0,
                       phiX, phiY, surface, pid, sig_corr);
    return 0;
}
