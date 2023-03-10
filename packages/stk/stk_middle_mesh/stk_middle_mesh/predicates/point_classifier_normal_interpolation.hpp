#ifndef POINT_CLASSIFIER_NONLINEAR_H
#define POINT_CLASSIFIER_NONLINEAR_H

#include <stk_middle_mesh/matrix.hpp>
#include <stk_middle_mesh/newton.hpp>
#include "point_classifier_impl.hpp" //TODO: only need this for the struct definitions
#include "point_classifier_normal_interpolation_reverse_solver.hpp"
#include "point_classifier_opts.hpp"

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

using stk::middle_mesh::impl::PointClassifierNormalInterpolationTolerances;

// use the nearly-orthogonal projection from
// "Efficient and Robust Algorithm for Overlaying Nonmatching Surface Meshes" by
// Jiao and Heath
class PointClassifierNormalInterpolation
{
  public:
    explicit PointClassifierNormalInterpolation(const PointClassifierNormalInterpolationTolerances& tolerances =
                                                    PointClassifierNormalInterpolationTolerances(1e-13))
      : m_xiClassifier(tolerances.pointClassifierTol)
      , m_mat3x3(3, 3)
      , m_x0(3)
      , m_newton(3, tolerances.newtonTol, tolerances.newtonItermax)
      , m_tolerances(tolerances)
    {}

    // project pt with normal vector normal onto element el
    PointRecordForTriangle classify(mesh::MeshEntityPtr el, utils::Point pt, utils::Point normal);

    // projection point pt onto element el.  normals contains the normals at the vertices
    // of el
    PointRecordForTriangle classify_reverse(mesh::MeshEntityPtr el, utils::Point pt,
                                            std::array<utils::Point, 4>& normals, bool logicalResultOnly = false);

  private:
    std::array<double, 3> solve_reverse3x3(const std::array<mesh::MeshEntityPtr, 3>& verts, const utils::Point& pt,
                                           const std::array<utils::Point, 4>& normals);

    void check_solution_reverse(const std::array<mesh::MeshEntityPtr, 3>& verts, const utils::Point& pt,
                                const std::array<utils::Point, 4>& normals, const std::array<double, 3>& solution);

    ElementXiClassifier m_xiClassifier;
    utils::impl::Matrix<double> m_mat3x3;
    std::vector<double> m_x0;
    utils::impl::Newton m_newton;
    PointClassifierNormalInterpolationTolerances m_tolerances;
};

} // namespace impl
} // namespace predicates

} // namespace middle_mesh
} // namespace stk
#endif
