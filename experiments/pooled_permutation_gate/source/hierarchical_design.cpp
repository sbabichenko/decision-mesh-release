#include "hierarchical_design.h"

#include "edge.h"
#include "face.h"
#include "mesh.h"
#include "observation.h"
#include "vertex.h"

#include <algorithm>
#include <cmath>
#include <set>

namespace hierarchical_design {
namespace {

void merge_column(SparseColumn& column, double tolerance) {
    std::sort(column.begin(), column.end(),
              [](const auto& left, const auto& right) {
                  return left.first < right.first;
              });
    SparseColumn merged;
    merged.reserve(column.size());
    for (const auto& entry : column) {
        if (!merged.empty() && merged.back().first == entry.first) {
            merged.back().second += entry.second;
        } else {
            merged.push_back(entry);
        }
    }
    merged.erase(
        std::remove_if(merged.begin(), merged.end(),
                       [tolerance](const auto& entry) {
                           return std::fabs(entry.second) <= tolerance;
                       }),
        merged.end());
    column.swap(merged);
}

}  // namespace

bool is_training_row(const DecisionMesh& mesh, int index, bool heldout_split) {
    if (!heldout_split) return true;
    if (index < 0 || index >= mesh.n_points) return false;
    if (mesh.training_mask.size() == static_cast<std::size_t>(mesh.n_points))
        return mesh.training_mask[index] != 0;
    // Defensive compatibility fallback. New real-data runs always install an
    // explicit mask before fitting.
    return index % 2 == 0;
}

ActiveDesign build_active_design(DecisionMesh& mesh,
                                 const BuildOptions& options) {
    ActiveDesign design;
    design.active_vertices.reserve(mesh.vertices.size());
    for (Vertex* vertex : mesh.vertices) {
        if (vertex && vertex->active) design.active_vertices.push_back(vertex);
    }
    std::sort(design.active_vertices.begin(), design.active_vertices.end(),
              [](const Vertex* left, const Vertex* right) {
                  return left->_id < right->_id;
              });

    // Birth order is a topological order for midpoint-parent dependencies.
    for (Vertex* vertex : design.active_vertices) {
        auto& map = design.parameter_map[vertex];
        const bool starved_root =
            options.hold_starved_roots && vertex->parent_edge == nullptr &&
            vertex->current_fit_points < Vertex::kMinimumFitPoints;
        if (starved_root) continue;
        if (vertex->parent_edge == nullptr || vertex->free_coefficient) {
            map[vertex] = 1.0;
            continue;
        }
        for (const Vertex* parent : {vertex->parent_edge->vertex0,
                                     vertex->parent_edge->vertex1}) {
            const auto parent_it = design.parameter_map.find(parent);
            if (parent_it == design.parameter_map.end()) continue;
            for (const auto& [coefficient, weight] : parent_it->second) {
                map[coefficient] += 0.5 * weight;
            }
        }
    }

    design.predictions.resize(mesh.n_points, 0.0);
    const bool have_observations = mesh.obs_store != nullptr;
    for (int index = 0; index < mesh.n_points; ++index) {
        design.predictions[index] = mesh.predict_index(index);
        if (options.training_only &&
            !is_training_row(mesh, index, options.heldout_split)) {
            continue;
        }
        if (have_observations && (*mesh.obs_store)[index].exposure <= 0.0) {
            continue;
        }
        Face* face = mesh.point_face[index];
        if (face == nullptr) continue;
        const auto& barycentric = mesh.point_bary[index];
        for (int corner = 0; corner < 3; ++corner) {
            const double basis = barycentric[corner];
            if (std::fabs(basis) <= options.coefficient_tolerance) continue;
            const auto map_it =
                design.parameter_map.find(face->vertices[corner]);
            if (map_it == design.parameter_map.end()) continue;
            for (const auto& [coefficient, derivative] : map_it->second) {
                const double value = basis * derivative;
                if (std::fabs(value) <= options.coefficient_tolerance) continue;
                design.columns[coefficient].push_back({index, value});
            }
        }
    }

    for (auto& [coefficient, column] : design.columns) {
        (void)coefficient;
        merge_column(column, options.coefficient_tolerance);
    }
    return design;
}


SparseColumn build_active_column(DecisionMesh& mesh,
                                 const Vertex* coefficient,
                                 const BuildOptions& options) {
    SparseColumn column;
    if (coefficient == nullptr || !coefficient->active) return column;

    std::vector<Vertex*> active;
    active.reserve(mesh.vertices.size());
    for (Vertex* vertex : mesh.vertices) {
        if (vertex && vertex->active) active.push_back(vertex);
    }
    std::sort(active.begin(), active.end(),
              [](const Vertex* left, const Vertex* right) {
                  return left->_id < right->_id;
              });

    std::unordered_map<const Vertex*, double> derivative;
    derivative.reserve(active.size());
    for (Vertex* vertex : active) {
        if (vertex == coefficient) {
            derivative[vertex] = 1.0;
            continue;
        }
        const bool starved_root =
            options.hold_starved_roots && vertex->parent_edge == nullptr &&
            vertex->current_fit_points < Vertex::kMinimumFitPoints;
        if (starved_root || vertex->parent_edge == nullptr ||
            vertex->free_coefficient) {
            derivative[vertex] = 0.0;
            continue;
        }
        double value = 0.0;
        for (const Vertex* parent : {vertex->parent_edge->vertex0,
                                     vertex->parent_edge->vertex1}) {
            const auto found = derivative.find(parent);
            if (found != derivative.end()) value += 0.5 * found->second;
        }
        derivative[vertex] = value;
    }

    std::set<Face*> support_faces;
    for (Vertex* vertex : active) {
        const auto found = derivative.find(vertex);
        if (found == derivative.end() ||
            std::fabs(found->second) <= options.coefficient_tolerance) {
            continue;
        }
        const auto faces = vertex->get_faces(false);
        support_faces.insert(faces.faces.begin(), faces.faces.end());
    }

    const bool have_observations = mesh.obs_store != nullptr;
    for (Face* face : support_faces) {
        if (face == nullptr) continue;
        double corner_derivative[3] = {0.0, 0.0, 0.0};
        for (int corner = 0; corner < 3; ++corner) {
            const auto found = derivative.find(face->vertices[corner]);
            if (found != derivative.end()) corner_derivative[corner] = found->second;
        }
        for (int index : face->coords_indices) {
            if (index < 0 || index >= mesh.n_points) continue;
            // Face::coords_indices is an acceleration cache and may retain
            // historical memberships after later NVB splits.  The canonical
            // live design is mesh.point_face/point_bary; filtering here makes
            // the targeted column exactly the same column used by
            // build_active_design and the stage solver.  Without this check a
            // very thin promoted coefficient could acquire rows from an
            // ancestor face that no longer owns them, so the fixed-point audit
            // profiled a different objective than the fitter optimized.
            if (mesh.point_face[index] != face) continue;
            if (options.training_only &&
                !is_training_row(mesh, index, options.heldout_split)) {
                continue;
            }
            if (have_observations &&
                (*mesh.obs_store)[index].exposure <= 0.0) {
                continue;
            }
            const auto& barycentric = mesh.point_bary[index];
            const double value =
                barycentric[0] * corner_derivative[0] +
                barycentric[1] * corner_derivative[1] +
                barycentric[2] * corner_derivative[2];
            if (std::fabs(value) > options.coefficient_tolerance) {
                column.push_back({index, value});
            }
        }
    }
    merge_column(column, options.coefficient_tolerance);
    return column;
}

double structural_norm_sq(const SparseColumn& column) {
    double norm_sq = 0.0;
    for (const auto& [index, value] : column) {
        (void)index;
        norm_sq += value * value;
    }
    return norm_sq;
}

void materialize_constraints(DecisionMesh& mesh) {
    std::vector<Vertex*> active;
    active.reserve(mesh.vertices.size());
    for (Vertex* vertex : mesh.vertices) {
        if (vertex && vertex->active) active.push_back(vertex);
    }
    std::sort(active.begin(), active.end(),
              [](const Vertex* left, const Vertex* right) {
                  return left->_id < right->_id;
              });
    for (Vertex* vertex : active) {
        if (vertex->parent_edge == nullptr || vertex->free_coefficient) continue;
        vertex->height = vertex->mu_lin();
        vertex->new_height = vertex->height;
        vertex->delta_pooled = 0.0;
        vertex->lambda_v = 0.0;
        vertex->sigma_pooled = 0.0;
        vertex->loss_reduction = 0.0;
        mesh.heap_set(vertex, 0.0);
    }
}

}  // namespace hierarchical_design
