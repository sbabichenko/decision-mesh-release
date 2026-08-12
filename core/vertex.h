#pragma once
#include <cmath>
#include <set>
#include <string>
#include <vector>
#include <utility>

struct DecisionMesh;
struct Edge;
struct Face;

struct Vertex {
    static constexpr int kMinimumFitPoints = 10;
    static double minimum_effective_support();
    bool gate_admitted = false;   // set only at lfdr-gate admissions (FDP accounting)
    // Geometry and statistical freedom are distinct. NVB conformity may need
    // an active midpoint even when no coefficient was admitted there. Such a
    // midpoint stays at zero hierarchical surplus (height == mu_lin) until it
    // is explicitly promoted by an unthresholded coarse activation or the
    // lfdr gate.
    bool free_coefficient = false;
    // A previously admitted coefficient whose current nodal column has no
    // training support during adaptive refinement. Dormant coefficients keep
    // their last accepted height, are skipped by coordinate sweeps, and wake
    // automatically if later topology gives them support again.
    bool dormant_coefficient = false;
    // Final fixed-topology state. Retired coefficients are constrained to the
    // regularization-only optimum (zero surplus) only after refinement ends.
    bool retired_coefficient = false;
    // Fix 1: containment holds SURPLUS, not nodal height. While contained,
    // the vertex's height tracks mu_prior(live parents) + frozen surplus,
    // so ancestor motion propagates instead of embalming a stale absolute
    // value mid-slope (the (5.6,5.30) distortion class).
    bool contained_state = false;
    double contained_surplus = 0.0;
    // Fix 1b: data-anchored hold for starved ROOTS (no parents to track).
    double neighbor_mean_height();
    // Historical origin is kept separate from current coefficient status.  A
    // midpoint introduced only to complete a conforming split may later pass
    // the gate and become free without losing the fact that it began as a
    // geometric completion vertex.
    bool conformity_origin = false;
    double admit_a = 0.0;         // selection threshold in delta-hat units (tau^2 truncation correction)
    double admit_sd = 0.0;        // calibrated candidate sd at admission (same scale as admit_a)
    // Self-child mode (DMESH_SELF_CHILD): frozen-increment bookkeeping.
    // sc_ref is the frozen reference height (the delta's zero point);
    // sc_thawed marks a coefficient admitted in the CURRENT gate round, the
    // only ones the stage solver may move while the mesh is frozen.
    bool sc_thawed = false;
    double sc_ref = 0.0;
    // DMESH_SC_SMOOTH: gated-measurement record. sc_obs is the height
    // committed at this coefficient's last admission (or base fit) and
    // sc_obs_var its data-only variance; the between-round Gaussian smoother
    // relaxes rendered heights under the hierarchical prior plus these
    // measurements only - raw data never moves values outside the gate.
    bool sc_has_obs = false;
    double sc_obs = 0.0;
    double sc_obs_var = 1e300;
    static int _seq;

    int _id;
    double x, y;
    double height = 0.0;
    double new_height = 0.0;
    double height_frozen = 0.0;  // height at last IRLS refresh (scoring epoch)
    double loss_reduction = 0.0;
    double heap_key = 0.0;       // current key in loss_heap (for O(log n) removal)
    bool in_heap = false;        // whether this vertex is in the heap
    bool active = false;
    bool disqualified = false;
    bool real_vertex = true;

    DecisionMesh* mesh = nullptr;
    Edge* parent_edge = nullptr;

    std::set<Edge*> edges;
    std::set<Vertex*> neighbors;
    std::vector<Vertex*> affected_vertices;

    // Empirical Bayes partial pooling fields
    int depth = 0;
    double sigma_sq = 1e300;       // data-only variance (1/xTx)
    double delta_pooled = 0.0;     // height - mu_lin
    double sigma_pooled = 1e300;   // posterior variance
    double lambda_v = 0.0;        // prior precision (1/sigma_prior)
    double prior_mean = 0.0;      // EB prior mean
    std::set<Vertex*> prior_children;  // vertices on this vertex's edges (for prior propagation)

    // Current-star diagnostics. Unlike admission-time support, these are
    // recomputed after every topology change and therefore expose coefficients
    // whose support has collapsed under deeper refinement.
    int current_n_points = 0;
    int current_fit_points = 0;
    double current_xTx = 0.0;
    double current_fit_xTx = 0.0;  // training-only coordinate curvature
    double current_keff = 0.0;

    Vertex() : _id(-1), x(0), y(0) {}
    Vertex(DecisionMesh* mesh, double x, double y, bool active = false,
           Edge* parent_edge = nullptr, bool real_vertex = true);

    std::string sid() const;

    // Vector arithmetic (returns temp vertex, not registered)
    double norm() const { return std::hypot(x, y); }
    double dot(const Vertex& o) const { return x * o.x + y * o.y; }
    double dot2(double a, double b) const { return x * a + y * b; }

    void add_edge(Edge* e);
    void remove_edge(Edge* e);

    void activate(bool make_free_coefficient = false);
    void update_info();
    void update_height();

    // Face collection
    struct FacesResult {
        std::vector<Vertex*> neighbors;
        std::vector<Face*> faces;
    };
    FacesResult get_faces(bool want_neighbors = true) const;
    FacesResult get_sim_faces(bool want_neighbors = true) const;

    std::set<Vertex*> neighbors_after_split() const;

    // Local regression
    struct LocRegressResult {
        double beta_opt;
        double loss_reduction;
        std::vector<Vertex*> affected;
        int n_points;
        int fit_points = 0;
        double fit_xTx = 0.0;  // training-only sum w*phi^2 (heldout excluded)
        double xTx, xTr, rTr;  // sufficient statistics for EB
        double xT2x = 0.0;     // sum w^2 phi^2 over the star (gate variance under extra dispersion)
        double xT2xs = 0.0;    // sum w^2 phi^2 s2_i (stratified one-law variance)
        double xT4k = 0.0;     // sum w^4 phi^4 s2_i^2 kurt_i (one-law CF term)
    };
    LocRegressResult loc_regress(bool want_affected = true);
    // Exact-likelihood 1-D profile over this vertex's coefficient
    // (DMESH_EXACT_SCORE): signed-root LR displacement + taxed argmax.
    struct ExactScore { bool ok = false; double b_star_units = 0.0;
                        double W = 0.0; double b_tax_units = 0.0; };
    ExactScore exact_score(double mu_units, double lam_prec_units,
                            bool use_live_neighbors = false);

    // Phase B1 candidate profile.  The null is the exact P1-preserving
    // candidate value (normally mu_lin), while the penalty may be centered at
    // a different EB target (mu_prior).  The alternative is optimized on the
    // exact tempered binomial objective along the candidate's hierarchical
    // column.  Because the null is feasible, gain_nats must be nonnegative up
    // to numerical tolerance.
    struct ExactCandidateGain {
        bool ok = false;
        double null_units = 0.0;
        double prior_mean_units = 0.0;
        double optimum_units = 0.0;
        double gain_nats = 0.0;
        // B2 decomposition of the exact penalized gain.  data_gain_nats
        // is the tempered-binomial log-likelihood change evaluated at the
        // penalized optimum; penalty_gain_nats is the change in the Gaussian
        // EB log-prior. Their sum is gain_nats.
        double data_gain_nats = 0.0;
        double penalty_gain_nats = 0.0;
        // Direction-preserving exact-gain coordinate proposed for the gate.
        // This is not assumed to be N(0,1); B2 measures its null law.
        double signed_root_gain = 0.0;
        double raw_gain_nats = 0.0;
        double objective_null = 0.0;
        double objective_optimum = 0.0;
        double data_objective_null = 0.0;
        double data_objective_optimum = 0.0;
        double score_at_null = 0.0;
        double information_at_null = 0.0;
        double data_information_optimum = 0.0;
        double lambda_logit = 0.0;
        double min_eta_null = 0.0;
        double max_eta_null = 0.0;
        double min_p_null = 0.0;
        double max_p_null = 0.0;
        double max_grid_excess = 0.0;
        int structural_points = 0;
        int training_structural_points = 0;
        int positive_information_points = 0;
        int fit_points = 0;
    };
    ExactCandidateGain exact_candidate_gain(double null_units,
                                             double prior_mean_units,
                                             double lam_prec_units,
                                             bool validate_grid = false) const;

    // Deterministic candidatewise null based on the measured pool-effect law.
    // The exact score column is used; each pool's binomial response is
    // integrated over a size-stratified two-component approximation to the
    // supplied effect histogram.  The few largest contributors are retained
    // as an explicit Gaussian mixture and the remainder is moment-collapsed.
    struct HistogramNullResult {
        bool ok = false;
        double p_value = 1.0;
        double gaussian_equivalent = 0.0;
        double observed_data_score = 0.0;
        double null_mean = 0.0;
        double null_sd = 0.0;
        double tail_keff = 0.0;
        int rows = 0;
        int explicit_mixture_rows = 0;
    };
    HistogramNullResult histogram_null_score(
        double null_units, double prior_mean_units,
        double lam_prec_units, const ExactCandidateGain& exact) const;
    // Profile a caller-supplied canonical hierarchical column.  This is used
    // by the fixed-point audit after building the full active design once, so
    // the audit and the exact stage solver evaluate literally the same sparse
    // column without an O(vertices x observations) targeted rebuild.
    ExactCandidateGain exact_gain_on_column(
        const std::vector<std::pair<int, double>>& column,
        double current_units,
        double null_units,
        double prior_mean_units,
        double lam_prec_units,
        bool validate_grid = false) const;

    // Empirical Bayes helpers
    double mu_lin() const;
    double mu_prior() const;
    double mu_quad() const;
    std::vector<std::pair<const Vertex*, double>> prior_stencil() const;
    void clamp_height();
    std::pair<double, double> compute_eb_params(double xTx) const;

    int degree() const { return (int)edges.size(); }
};
