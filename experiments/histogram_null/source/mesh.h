#pragma once
#include "observation.h"
#include <vector>
#include <array>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>

struct TreeNode;
struct Vertex;
struct Edge;
struct Face;

struct EdgeCoordLess {
    bool operator()(const Edge* a, const Edge* b) const;
};

struct TimingRecord {
    int iteration;
    double find_best_us;   // microseconds spent finding best vertex
    double split_us;       // microseconds spent splitting/activating
    double update_info_us; // microseconds spent in update_info propagation
    double total_us;       // total for this refinement
    int active_faces;
    int n_vertices;
    // Diagnostic counters for understanding cascade behavior
    int cascade_activations;  // how many recursive vertex activations (completions)
    int faces_created;        // faces created during this refinement
    int edges_created;        // edges created during this refinement
    int add_face_calls;       // Edge::add_face calls (each does coord computation)
    int update_info_calls;    // Vertex::update_info calls
};

// Thread-local counters for cascade instrumentation
struct CascadeCounters {
    int cascade_activations = 0;
    int faces_created = 0;
    int edges_created = 0;
    int add_face_calls = 0;
    int update_info_calls = 0;
};

// Global counters (single-threaded, so this is fine)
extern CascadeCounters g_counters;

struct DecisionMesh {
    double max_aspect_ratio = 5.0;
    bool maintain_heap = true;  // heap only needed by the greedy updater, not the lfdr gate
    bool use_eb = false;  // empirical Bayes regularization
    const std::vector<Observation>* obs_store = nullptr;
    double last_extra = 0.0;
    // Frozen linearization epoch: per-pool eta and per-vertex heights captured
    // at the last IRLS refresh. Every scorer (quadratic via face stats, exact
    // via pool loops) answers the same question: "given the LAST REFRESH
    // state, where should this coefficient go" -- keeping sweep dynamics
    // identical across scoring strategies.
    std::vector<double> eta_frozen;
    int refresh_epoch = 0;

    // Data: n rows, 2 columns for X, 1 column for values
    int n_points = 0;
    std::vector<double> X;       // n_points * 2, row-major
    std::vector<double> values;  // n_points
    std::vector<double> weights; // per-point precision (default 1)
    // Fixed likelihood temper for the current fitting epoch.  The exact
    // objective is sum_i c_i ell_i; c_i is set when the stage/IRLS weights are
    // formed and must not be recomputed from the live probability during
    // coordinate ascent.  Recomputing it made the alleged exact objective
    // state-dependent and produced a two-cycle on thin columns.
    std::vector<double> likelihood_temper;
    std::vector<double> gate_s2;  // per-point pool-effect variance (one-law gate); empty = off
    std::vector<double> gate_k4;  // per-point sigma^4 * excess kurtosis (one-law CF term)

    // Training-point geometry cache. Each active face owns a disjoint subset
    // of observations, and updates these entries when it is created/refreshed.
    // This turns repeated in-sample prediction into three height reads and a
    // dot product rather than a tree walk plus barycentric recomputation.
    std::vector<Face*> point_face;
    std::vector<std::array<double, 3>> point_bary;
    std::vector<std::array<double, 3>> point_basis;

    // Incremented only when topology changes. Vertex star caches are valid
    // throughout height sweeps and candidate rescoring at a fixed topology.
    std::size_t topology_version = 0;

    // EB specification: independent, zero-mean hierarchical surpluses
    // delta_v = h_v - (h_p0+h_p1)/2 ~ N(0, tau_depth^2).
    // tau_depth is estimated only from lfdr-gate admissions.
    double tau_default_sq = 1.0e4;  // weak default: prior sd = 1 log-odds (heights are x100)
    double tau_floor_sq = 100.0;    // numerical floor: prior sd = 0.10 log-odds

    double xmin, xmax, ymin, ymax;

    // Object ownership (heap-allocated, cleaned up in destructor)
    std::vector<TreeNode*> all_nodes;
    std::vector<Vertex*> all_vertices;
    // coordinate -> vertex lookup (quantized unit coords) for prior stencils
    std::unordered_map<long long, Vertex*> coord_map;
    static long long coord_key(double x, double y) {
        return (static_cast<long long>(std::llround(x * 1048576.0)) << 21)
             ^ static_cast<long long>(std::llround(y * 1048576.0));
    }
    Vertex* vertex_at(double x, double y) const {
        if (x < -1e-9 || x > 1 + 1e-9 || y < -1e-9 || y > 1 + 1e-9) return nullptr;
        auto it = coord_map.find(coord_key(x, y));
        return it == coord_map.end() ? nullptr : it->second;
    }
    std::vector<Edge*> all_edges;
    std::vector<Face*> all_faces;

    TreeNode* root = nullptr;
    std::set<TreeNode*> leaves;
    std::set<Face*> active_faces;
    std::set<Vertex*> vertices;      // "real" vertices
    std::set<Edge*, EdgeCoordLess> active_edges;

    // Loss heap: map from negative loss_reduction -> vertex
    std::map<double, std::set<Vertex*>> loss_heap;

    // Outer geometry
    std::map<std::string, Vertex*> outer_vertices;
    std::map<std::string, Edge*> outer_edges;
    Edge* split_edge = nullptr;
    Face* top_face = nullptr;
    Face* bottom_face = nullptr;

    std::mt19937 rng;

    // Empirical Bayes state
    std::map<int, double> tau_sq;       // depth -> tau_sq
    std::map<int, int> tau_sq_vertex_counts;  // depth -> count at last recomputation
    int tau_sq_recompute_interval = 20;
    int steps_since_tau_recompute = 0;

    DecisionMesh(const std::vector<double>& x_data,
                 const std::vector<double>& y_data,
                 const std::vector<double>& z_data,
                 bool use_eb = false);
    ~DecisionMesh();

    // Factory methods (own the memory)
    TreeNode* make_node(TreeNode* parent = nullptr, Face* face = nullptr);
    Vertex* make_vertex(double x, double y, bool active = false,
                        Edge* parent_edge = nullptr, bool real_vertex = true);
    Edge* make_edge(Vertex* v0, Vertex* v1, bool active);
    Face* make_face(Edge* e0, Edge* e1, Edge* e2,
                    const std::vector<bool>& mask, bool active = true,
                    const std::string& path = "", bool skip_coords = false);

    // Heap operations
    void heap_set(Vertex* v, double neg_loss);
    void heap_pop(Vertex* v);
    std::pair<Vertex*, double> heap_peek() const;

    // Core
    void create_outer_vertices();
    void create_outer_edges();

    Face* random_face();
    void update_best_vertex(double random_prob = 0.0);
    TimingRecord update_best_vertex_timed(double random_prob, int iteration);

    // Empirical Bayes
    void recompute_tau_sq();
    void maybe_recompute_tau_sq();
    double tau_for_depth(int depth) const;

    double get_x(int i, int col) const { return X[i * 2 + col]; }

    // Prediction: walk tree to leaf in O(log n), then affine interpolate
    double predict(double px, double py) const;
    double predict_index(int i) const;
    // Batch prediction
    std::vector<double> predict_batch(const std::vector<double>& px,
                                      const std::vector<double>& py) const;

    // Export mesh to SVG (simple visualization)
    void write_svg(const std::string& filename, double width = 800, double height = 800) const;

    // Export mesh data for Python plotting (vertices + triangles)
    void write_mesh_csv(const std::string& prefix) const;
};
