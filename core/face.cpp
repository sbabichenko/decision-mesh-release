#include "face.h"
#include "mesh.h"
#include "edge.h"
#include "vertex.h"
#include "tree.h"
#include <cmath>
#include <cstring>
#include <algorithm>

int Face::_seq = 0;

Face::Face(DecisionMesh* mesh, Edge* e0, Edge* e1, Edge* e2,
           const std::vector<bool>& mask, bool active_flag,
           const std::string& path, bool skip_coords)
    : _id(Face::_seq++), path(path), mesh(mesh), mask(mask)
{
    edges[0] = e0;
    edges[1] = e1;
    edges[2] = e2;

    // Resolve vertices: vertex0 is opposite edge0
    Vertex* v1 = e0->vertex0;
    Vertex* v2 = e0->vertex1;
    Vertex* v0 = e2->other_vertex(v1);
    if (!v0) {
        std::swap(v1, v2);
        v0 = e2->other_vertex(v1);
    }
    vertices[0] = v0;
    vertices[1] = v1;
    vertices[2] = v2;

    if (!skip_coords) {
        update_coords();
    }

    if (active_flag) {
        activate();
    }
}

std::string Face::sid() const {
    char buf[16];
    snprintf(buf, sizeof(buf), "F%03d", _id);
    return buf;
}

void Face::activate() {
    if (active) return;
    active = true;
    mesh->active_faces.insert(this);
    refresh_point_cache();
    for (int i = 0; i < 3; ++i) {
        edges[i]->activate();
        edges[i]->add_face(this);
    }
}

void Face::deactivate() {
    if (!active) return;
    active = false;
    mesh->active_faces.erase(this);
    for (int i = 0; i < 3; ++i) {
        edges[i]->remove_face(this);
    }
}

void Face::split(Edge* edge) {
    if (edge != edges[0]) {
        // Completion: must split along refinement edge first
        edges[0]->midpoint->activate(false);
        return;
    }

    auto& sd = sub_divisions[0];
    if (!sd.e) {
        return;  // Edge not activated yet
    }

    deactivate();
    // chord separates the children, not the bisected edge
    node->split(sd.e, sd.plus_face, sd.minus_face);
    sd.plus_face->activate();
    sd.minus_face->activate();

    // This face will never become active again.  Its children now own the
    // point partition, so retaining the parent's indices makes memory grow
    // like N times mesh depth.  Release the capacity, not just the size.
    release_point_storage();
}

double Face::area() const {
    double det = (vertices[0]->x - vertices[2]->x) * (vertices[1]->y - vertices[0]->y) -
                 (vertices[0]->x - vertices[1]->x) * (vertices[2]->y - vertices[0]->y);
    return 0.5 * std::abs(det);
}

double Face::aspect_ratio() const {
    double max_len = 0, min_len = 1e300;
    for (int i = 0; i < 3; ++i) {
        double l = edges[i]->length;
        if (l > max_len) max_len = l;
        if (l < min_len) min_len = l;
    }
    if (min_len <= 0) return 1e300;
    return max_len / min_len;
}

// Accumulate sufficient statistics for one point with barycentric weights w
// and response value y. Inlined to avoid function call overhead in hot loop.
static inline void accum_suff_stats(double S_ww[3][3], double S_wy[3],
                                     double& S_yy, double S_w2[3], double om, double w0, double u, double v,
                                     double y, double S_w2s[3] = nullptr, double S_w4k[3] = nullptr,
                                     double s2i = 0.0, double k4i = 0.0) {
    S_w2[0] += om * om * w0 * w0;
    S_w2[1] += om * om * u * u;
    S_w2[2] += om * om * v * v;
    if (S_w2s) {
        double o2 = om * om, o4 = o2 * o2;
        S_w2s[0] += o2 * w0 * w0 * s2i;
        S_w2s[1] += o2 * u * u * s2i;
        S_w2s[2] += o2 * v * v * s2i;
        S_w4k[0] += o4 * w0 * w0 * w0 * w0 * k4i;
        S_w4k[1] += o4 * u * u * u * u * k4i;
        S_w4k[2] += o4 * v * v * v * v * k4i;
    }
    // S_ww: only upper triangle, symmetric fill at end
    S_ww[0][0] += om * w0 * w0;
    S_ww[0][1] += om * w0 * u;
    S_ww[0][2] += om * w0 * v;
    S_ww[1][1] += om * u * u;
    S_ww[1][2] += om * u * v;
    S_ww[2][2] += om * v * v;
    S_wy[0] += om * w0 * y;
    S_wy[1] += om * u * y;
    S_wy[2] += om * v * y;
    S_yy += om * y * y;
}

static inline void fill_symmetric(double S_ww[3][3]) {
    S_ww[1][0] = S_ww[0][1];
    S_ww[2][0] = S_ww[0][2];
    S_ww[2][1] = S_ww[1][2];
}

void Face::update_coords(double eps) {
    coords_indices.clear();

    // Zero sufficient statistics
    std::memset(S_ww, 0, sizeof(S_ww));
    std::memset(S_wy, 0, sizeof(S_wy));
    std::memset(S_w2, 0, sizeof(S_w2));
    std::memset(S_w2s, 0, sizeof(S_w2s));
    std::memset(S_w4k, 0, sizeof(S_w4k));
    std::memset(S_fit_xx, 0, sizeof(S_fit_xx));
    S_yy = 0.0;
    n_covered = 0;
    n_fit_covered = 0;

    Vertex* v0 = vertices[0];
    Vertex* v1 = vertices[1];
    Vertex* v2 = vertices[2];

    // Right-triangle property: v0 is the right-angle vertex (opposite hypotenuse).
    // Legs a = v1-v0, b = v2-v0 are perpendicular, so a·b = 0.
    // Barycentric coords simplify: u = (p·a)/|a|², v = (p·b)/|b|²
    double ax = v1->x - v0->x, ay = v1->y - v0->y;
    double bx = v2->x - v0->x, by = v2->y - v0->y;

    double d00 = ax * ax + ay * ay;
    double d11 = bx * bx + by * by;

    bool degenerate = (d00 < eps || d11 < eps);

    if (!degenerate) {
        double inv_d00 = 1.0 / d00;
        double inv_d11 = 1.0 / d11;
        double v0x = v0->x, v0y = v0->y;

        for (int i = 0; i < mesh->n_points; ++i) {
            if (!mask[i]) continue;

            double px = mesh->get_x(i, 0) - v0x;
            double py = mesh->get_x(i, 1) - v0y;
            double u = (px * ax + py * ay) * inv_d00;
            double v = (px * bx + py * by) * inv_d11;
            double w0 = 1.0 - u - v;

            coords_indices.push_back(i);
            accum_suff_stats(S_ww, S_wy, S_yy, S_w2, mesh->weights[i], w0, u, v, mesh->values[i],
                             S_w2s, S_w4k,
                             mesh->gate_s2.empty() ? 0.0 : mesh->gate_s2[i],
                             mesh->gate_k4.empty() ? 0.0 : mesh->gate_k4[i]);
            n_covered++;
            if (mesh->weights[i] > 1e-9) {
                n_fit_covered++;
                S_fit_xx[0] += mesh->weights[i] * w0 * w0;
                S_fit_xx[1] += mesh->weights[i] * u * u;
                S_fit_xx[2] += mesh->weights[i] * v * v;
            }
        }
    } else {
        for (int i = 0; i < mesh->n_points; ++i) {
            if (!mask[i]) continue;
            coords_indices.push_back(i);
        }
    }

    fill_symmetric(S_ww);
}

void Face::update_coords_from_indices(const std::vector<int>& indices, double eps) {
    // Guard against self-aliasing: callers pass coords_indices into this method,
    // and clearing it first would empty the input before the loop reads it.
    const std::vector<int>* src = &indices;
    std::vector<int> snapshot;
    if (&indices == &coords_indices) { snapshot = indices; src = &snapshot; }
    const std::vector<int>& idx = *src;

    coords_indices.clear();

    std::memset(S_ww, 0, sizeof(S_ww));
    std::memset(S_wy, 0, sizeof(S_wy));
    std::memset(S_w2, 0, sizeof(S_w2));
    std::memset(S_w2s, 0, sizeof(S_w2s));
    std::memset(S_w4k, 0, sizeof(S_w4k));
    std::memset(S_fit_xx, 0, sizeof(S_fit_xx));
    S_yy = 0.0;
    n_covered = 0;
    n_fit_covered = 0;

    Vertex* v0 = vertices[0];
    Vertex* v1 = vertices[1];
    Vertex* v2 = vertices[2];

    double ax = v1->x - v0->x, ay = v1->y - v0->y;
    double bx = v2->x - v0->x, by = v2->y - v0->y;

    double d00 = ax * ax + ay * ay;
    double d11 = bx * bx + by * by;

    bool degenerate = (d00 < eps || d11 < eps);

    coords_indices.reserve(idx.size());

    if (!degenerate) {
        double inv_d00 = 1.0 / d00;
        double inv_d11 = 1.0 / d11;
        double v0x = v0->x, v0y = v0->y;

        for (int i : idx) {
            double px = mesh->get_x(i, 0) - v0x;
            double py = mesh->get_x(i, 1) - v0y;
            double u = (px * ax + py * ay) * inv_d00;
            double v = (px * bx + py * by) * inv_d11;
            double w0 = 1.0 - u - v;

            coords_indices.push_back(i);
            accum_suff_stats(S_ww, S_wy, S_yy, S_w2, mesh->weights[i], w0, u, v, mesh->values[i],
                             S_w2s, S_w4k,
                             mesh->gate_s2.empty() ? 0.0 : mesh->gate_s2[i],
                             mesh->gate_k4.empty() ? 0.0 : mesh->gate_k4[i]);
            n_covered++;
            if (mesh->weights[i] > 1e-9) {
                n_fit_covered++;
                S_fit_xx[0] += mesh->weights[i] * w0 * w0;
                S_fit_xx[1] += mesh->weights[i] * u * u;
                S_fit_xx[2] += mesh->weights[i] * v * v;
            }
        }
    } else {
        for (int i : idx) {
            coords_indices.push_back(i);
        }
    }

    fill_symmetric(S_ww);
}



void Face::refresh_point_cache(double eps) {
    Vertex* v0 = vertices[0];
    Vertex* v1 = vertices[1];
    Vertex* v2 = vertices[2];
    const double ax = v1->x - v0->x, ay = v1->y - v0->y;
    const double bx = v2->x - v0->x, by = v2->y - v0->y;
    const double d00 = ax * ax + ay * ay;
    const double d11 = bx * bx + by * by;
    const double det = ax * by - ay * bx;
    if (std::abs(det) < eps) {
        for (int i : coords_indices) {
            mesh->point_face[i] = this;
            mesh->point_bary[i] = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
            mesh->point_basis[i] = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        }
        return;
    }
    const double inv = 1.0 / det;
    for (int i : coords_indices) {
        const double dpx = mesh->get_x(i, 0) - v0->x;
        const double dpy = mesh->get_x(i, 1) - v0->y;
        const double su = (dpx * ax + dpy * ay) / d00;
        const double sv = (dpx * bx + dpy * by) / d11;
        mesh->point_basis[i] = {1.0 - su - sv, su, sv};
        double u = (dpx * by - dpy * bx) * inv;
        double v = (ax * dpy - ay * dpx) * inv;
        double w = 1.0 - u - v;
        // Match predict() exactly at numerical boundaries.
        double cw = std::max(0.0, std::min(1.0, w));
        double cu = std::max(0.0, std::min(1.0, u));
        double cv = std::max(0.0, std::min(1.0, v));
        const double total = cw + cu + cv;
        if (total > 0.0) { cw /= total; cu /= total; cv /= total; }
        else { cw = cu = cv = 1.0 / 3.0; }
        mesh->point_face[i] = this;
        mesh->point_bary[i] = {cw, cu, cv};
    }
}


void Face::refresh_suff_stats_from_cache() {
    std::memset(S_ww, 0, sizeof(S_ww));
    std::memset(S_wy, 0, sizeof(S_wy));
    std::memset(S_w2, 0, sizeof(S_w2));
    std::memset(S_w2s, 0, sizeof(S_w2s));
    std::memset(S_w4k, 0, sizeof(S_w4k));
    std::memset(S_fit_xx, 0, sizeof(S_fit_xx));
    S_yy = 0.0;
    n_covered = 0;
    n_fit_covered = 0;
    for (int i : coords_indices) {
        const auto& b = mesh->point_basis[i];
        const double w0 = b[0], u = b[1], v = b[2];
        accum_suff_stats(S_ww, S_wy, S_yy, S_w2, mesh->weights[i], w0, u, v, mesh->values[i],
                             S_w2s, S_w4k,
                             mesh->gate_s2.empty() ? 0.0 : mesh->gate_s2[i],
                             mesh->gate_k4.empty() ? 0.0 : mesh->gate_k4[i]);
        ++n_covered;
        if (mesh->weights[i] > 1e-9) {
            ++n_fit_covered;
            S_fit_xx[0] += mesh->weights[i] * w0 * w0;
            S_fit_xx[1] += mesh->weights[i] * u * u;
            S_fit_xx[2] += mesh->weights[i] * v * v;
        }
    }
    fill_symmetric(S_ww);
}

void Face::release_point_storage() {
    std::vector<int>().swap(coords_indices);
    std::vector<bool>().swap(mask);
    std::string().swap(path);
}

void Face::add_sub_division(Edge* edge, SubDivision sd) {
    for (int i = 0; i < 3; ++i) {
        if (edges[i] == edge) {
            sub_divisions[i] = sd;
            return;
        }
    }
}
