#ifndef KRINO_KRINO_GEOMETRY_AKRI_SIMPLEXCURVATURE_HPP_
#define KRINO_KRINO_GEOMETRY_AKRI_SIMPLEXCURVATURE_HPP_
#include <vector>
#include <stk_math/StkVector.hpp>

namespace krino {

double laplacian_of_field_at_node_of_patch_of_triangles(
    const unsigned node,
    const std::vector<std::array<unsigned,3>> & tris,
    const std::vector<stk::math::Vector3d> & coords,
    const std::vector<double> & u);

double laplacian_of_field_at_node_of_patch_of_tetrahedra(
    const unsigned node,
    const std::vector<std::array<unsigned,4>> & tets,
    const std::vector<stk::math::Vector3d> & coords,
    const std::vector<double> & u);

std::pair<double,double> tri_area_and_divergence_of_node_normals(
    const stk::math::Vector3d & coords0,
    const stk::math::Vector3d & coords1,
    const stk::math::Vector3d & coords2,
    const stk::math::Vector3d & nodeNormals0,
    const stk::math::Vector3d & nodeNormals1,
    const stk::math::Vector3d & nodeNormals2);

std::pair<double,double> tet_volume_and_divergence_of_node_normals(
    const stk::math::Vector3d & coords0,
    const stk::math::Vector3d & coords1,
    const stk::math::Vector3d & coords2,
    const stk::math::Vector3d & coords3,
    const stk::math::Vector3d & nodeNormals0,
    const stk::math::Vector3d & nodeNormals1,
    const stk::math::Vector3d & nodeNormals2,
    const stk::math::Vector3d & nodeNormals3);

template<size_t N>
double divergence_of_node_normals_on_patch_of_elements(
    const std::vector<std::array<unsigned,N>> & elemConn,
    const std::vector<stk::math::Vector3d> & coords,
    const std::vector<stk::math::Vector3d> & nodeNormals);

double finite_difference_divergence_of_normal_2d(const double delta, const double f0,
    const double fxm, const double fxp,
    const double fym, const double fyp,
    const double fxmym, const double fxmyp, const double fxpym, const double fxpyp);

double finite_difference_divergence_of_normal_3d(const double delta, const double f0,
    const double fxm, const double fxp,
    const double fym, const double fyp,
    const double fzm, const double fzp,
    const double fxmym, const double fxmyp, const double fxpym, const double fxpyp,
    const double fxmzm, const double fxmzp, const double fxpzm, const double fxpzp,
    const double fymzm, const double fymzp, const double fypzm, const double fypzp);

}

#endif /* KRINO_KRINO_GEOMETRY_AKRI_SIMPLEXCURVATURE_HPP_ */
