#pragma once

#include <unordered_map>
#include <utility>
#include <vector>

struct DecisionMesh;
struct Vertex;

namespace hierarchical_design {

using SparseColumn = std::vector<std::pair<int, double>>;
using ParameterMap =
    std::unordered_map<const Vertex*, std::unordered_map<Vertex*, double>>;

struct BuildOptions {
    bool heldout_split = false;
    bool training_only = true;
    bool hold_starved_roots = true;
    double coefficient_tolerance = 1e-14;
};

struct ActiveDesign {
    std::vector<Vertex*> active_vertices;
    ParameterMap parameter_map;
    std::unordered_map<Vertex*, SparseColumn> columns;
    std::vector<double> predictions;
};

bool is_training_row(int index, bool heldout_split);

// Build the hierarchical columns used by the exact stage solver. A free
// coefficient starts a new column; a constrained midpoint propagates one half
// of each parent's derivative; another free descendant stops propagation.
ActiveDesign build_active_design(DecisionMesh& mesh,
                                 const BuildOptions& options = {});

// Build one active hierarchical column with the same birth-order recurrence
// as build_active_design. This is the canonical targeted column used by
// exact candidate promotion and fixed-point audits.
SparseColumn build_active_column(DecisionMesh& mesh,
                                 const Vertex* coefficient,
                                 const BuildOptions& options = {});

// Sum of squared raw hierarchical-column values on the selected design rows.
// This is geometric support, not Fisher information and not gate replication.
double structural_norm_sq(const SparseColumn& column);

// Materialize all constrained (geometry-only or retired) midpoint heights in
// birth order so every one has exactly zero hierarchical surplus.
void materialize_constraints(DecisionMesh& mesh);

}  // namespace hierarchical_design
