#include "output.h"

#include "face.h"
#include "edge.h"
#include "mesh.h"
#include "units.h"
#include "vertex.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::ofstream open_output(const std::string& filename) {
    std::ofstream stream(filename);
    if (!stream) throw std::runtime_error("unable to open output file: " + filename);
    stream << std::setprecision(10);
    return stream;
}

double apply_axis_map(const std::vector<double>& map, double value) {
    const int band_count = static_cast<int>(map.size()) - 1;
    if (band_count <= 0) return value;
    const double position = std::clamp(value, 0.0, 1.0) * band_count;
    const int band = std::min(static_cast<int>(position), band_count - 1);
    return map[band] + (position - band) * (map[band + 1] - map[band]);
}

}  // namespace

void write_edf_dump(const OutputConfig& output,
                    const DecisionMesh& mesh,
                    int observation_count) {
    if (!output.edf_dump_file) return;
    auto file = open_output(*output.edf_dump_file);
    for (int i = 0; i < observation_count; ++i) {
        const double fitted = mesh.predict_index(i);
        file << to_log_odds(mesh.values[i]) << ','
             << to_log_odds(fitted) << ','
             << mesh.weights[i] * kVarianceScale << '\n';
    }
}

void write_figure_dumps(const OutputConfig& output,
                        const DecisionMesh& mesh,
                        const std::vector<Observation>& observations,
                        const std::vector<double>& x,
                        const std::vector<double>& y,
                        const std::vector<double>& original_x,
                        const std::vector<double>& original_y,
                        const std::vector<double>& x_map,
                        const std::vector<double>& y_map,
                        const SurfaceSpec& surface,
                        int pool_count,
                        double corrected_pool_variance) {
    if (!output.prefix) return;
    const std::string& base = *output.prefix;

    {
        auto file = open_output(base + "_surface.csv");
        file << "gx,gy,wala,wac,g_hat,f_true\n";
        constexpr int grid_size = 60;
        for (int a = 0; a < grid_size; ++a) {
            for (int b = 0; b < grid_size; ++b) {
                const double gx = (a + 0.5) / grid_size;
                const double gy = (b + 0.5) / grid_size;
                const double fitted = mesh.predict(apply_axis_map(x_map, gx), apply_axis_map(y_map, gy));
                const double xx = -4.0 + 8.0 * gx;
                const double yy = -4.0 + 8.0 * gy;
                file << gx << ',' << gy << ',' << gx * 60.0 << ',' << 2.6 + gy * 5.0 << ','
                     << fitted << ',' << from_log_odds(surface.evaluate(xx, yy)) << '\n';
            }
        }
    }

    {
        auto file = open_output(base + "_obs.csv");
        file << "wala,wac,pid,n,k,p_hat,resid,info,g_hat,u,is_train\n";
        for (std::size_t i = 0; i < observations.size(); ++i) {
            const double fitted = mesh.predict_index(static_cast<int>(i));
            const double probability = 1.0 / (1.0 + std::exp(-(observations[i].offset + to_log_odds(fitted))));
            const double information = observations[i].exposure * probability * (1.0 - probability);
            const double residual = (observations[i].events - observations[i].exposure * probability) /
                                    std::sqrt(information);
            const double raw_x = original_x.empty() ? x[i] : original_x[i];
            const double raw_y = original_y.empty() ? y[i] : original_y[i];
            file << raw_x * 60.0 << ',' << 2.6 + raw_y * 5.0 << ',' << observations[i].pool_id << ','
                 << observations[i].exposure << ',' << observations[i].events << ',' << probability << ','
                 << residual << ',' << information << ',' << fitted << ',' << observations[i].latent_effect << ','
                 << ((mesh.training_mask.size() == observations.size()) ? int(mesh.training_mask[i]) : 1) << '\n';
        }
    }

    {
        auto file = open_output(base + "_hier_vertices.csv");
        file << "id,x,y,height,new_height,depth,parent0,parent1,active,gate_admitted,free_coefficient,dormant_coefficient,retired_coefficient,conformity_origin,"
                "sigma_sq,delta_pooled,sigma_pooled,lambda_v,prior_mean,admit_a,"
                "current_n_points,current_fit_points,current_xTx,current_fit_xTx,current_keff\n";
        std::vector<Vertex*> verts(mesh.vertices.begin(), mesh.vertices.end());
        std::sort(verts.begin(), verts.end(), [](const Vertex* a, const Vertex* b) { return a->_id < b->_id; });
        for (const Vertex* v : verts) {
            int p0 = -1, p1 = -1;
            if (v->parent_edge != nullptr) {
                p0 = v->parent_edge->vertex0->_id;
                p1 = v->parent_edge->vertex1->_id;
            }
            file << v->_id << ',' << v->x << ',' << v->y << ',' << v->height << ','
                 << v->new_height << ',' << v->depth << ',' << p0 << ',' << p1 << ','
                 << (v->active ? 1 : 0) << ',' << (v->gate_admitted ? 1 : 0) << ','
                 << (v->free_coefficient ? 1 : 0) << ','
                 << (v->dormant_coefficient ? 1 : 0) << ','
                 << (v->retired_coefficient ? 1 : 0) << ','
                 << (v->conformity_origin ? 1 : 0) << ','
                 << v->sigma_sq << ',' << v->delta_pooled << ',' << v->sigma_pooled << ','
                 << v->lambda_v << ',' << v->prior_mean << ',' << v->admit_a << ','
                 << v->current_n_points << ',' << v->current_fit_points << ','
                 << v->current_xTx << ',' << v->current_fit_xTx << ','
                 << v->current_keff << '\n';
        }
    }

    {
        auto file = open_output(base + "_mesh.csv");
        file << "x0,y0,h0,x1,y1,h1,x2,y2,h2,n_pts\n";
        for (const Face* face : mesh.active_faces) {
            const Vertex* v0 = face->vertices[0];
            const Vertex* v1 = face->vertices[1];
            const Vertex* v2 = face->vertices[2];
            file << v0->x << ',' << v0->y << ',' << v0->height << ','
                 << v1->x << ',' << v1->y << ',' << v1->height << ','
                 << v2->x << ',' << v2->y << ',' << v2->height << ','
                 << face->coords_indices.size() << '\n';
        }
    }

    {
        std::vector<double> weighted_error(pool_count, 0.0);
        std::vector<double> information(pool_count, 0.0);
        std::vector<double> truth(pool_count, 0.0);
        const double scaled_variance = std::max(corrected_pool_variance, 1e-4) * kVarianceScale;
        for (std::size_t i = 0; i < observations.size(); ++i) {
            const double fitted = mesh.predict_index(static_cast<int>(i));
            const double probability = 1.0 / (1.0 + std::exp(-(observations[i].offset + to_log_odds(fitted))));
            const double obs_information = observations[i].exposure * probability * (1.0 - probability);
            const double working_residual =
                (observations[i].events - observations[i].exposure * probability) /
                obs_information * kHeightScale;
            const int pool = observations[i].pool_id;
            weighted_error[pool] += (obs_information / kVarianceScale) * working_residual;
            information[pool] += obs_information / kVarianceScale;
            truth[pool] = from_log_odds(observations[i].latent_effect);
        }
        auto file = open_output(base + "_pools.csv");
        file << "pid,u_hat,u_true,W\n";
        for (int pool = 0; pool < pool_count; ++pool) {
            if (information[pool] <= 0.0) continue;
            const double estimate = weighted_error[pool] /
                                    (information[pool] + 1.0 / scaled_variance);
            file << pool << ',' << estimate << ',' << truth[pool] << ',' << information[pool] << '\n';
        }
    }
}
