/*
 * Akri_SimplexCurvature.cpp
 *
 *  Created on: Jan 23, 2026
 *      Author: drnoble
 */
#include <Akri_SimplexCurvature.hpp>
#include <vector>
#include <stk_math/StkVector.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <tuple>

namespace krino {

template<size_t NPE>
unsigned find_index_of_node_in_element(const std::array<unsigned,NPE> & origTet, const unsigned node)
{
  for (unsigned i=0; i<NPE; ++i)
    if (origTet[i] == node)
      return i;
  STK_ThrowRequireMsg(false, "Node not found in element.");
  return 0;
}

static std::array<unsigned,3> permute_tri_so_that_node_is_first(const std::array<unsigned,3> & origTri, const unsigned node)
{
  static const std::array<std::array<unsigned,3>,3> nodePermute{{ {{0,1,2}}, {{1,2,0}}, {{2,0,1}} }};
  const std::array<unsigned,3> & perm = nodePermute[find_index_of_node_in_element(origTri, node)];
  return {{origTri[perm[0]], origTri[perm[1]], origTri[perm[2]]}};
}

static std::array<unsigned,4> permute_tet_so_that_node_is_first(const std::array<unsigned,4> & origTet, const unsigned node)
{
  static const std::array<std::array<unsigned,4>,4> nodePermute{{ {{0,1,2,3}}, {{1,2,0,3}}, {{2,0,1,3}}, {{3,1,0,2}} }};
  const std::array<unsigned,4> & perm = nodePermute[find_index_of_node_in_element(origTet, node)];
  return {{origTet[perm[0]], origTet[perm[1]], origTet[perm[2]], origTet[perm[3]]}};
}

double laplacian_of_field_at_node_of_patch_of_tetrahedra(
    const unsigned node,
    const std::vector<std::array<unsigned,4>> & tets,
    const std::vector<stk::math::Vector3d> & coords,
    const std::vector<double> & u)
{
  double sumVolW = 0.0;
  double sumContrib = 0.0;

  for(auto & origTet : tets)
  {
    const auto & tet = permute_tet_so_that_node_is_first(origTet, node);

    const double u10 = u[tet[1]] - u[tet[0]];
    const double u20 = u[tet[2]] - u[tet[0]];
    const double u30 = u[tet[3]] - u[tet[0]];
    const stk::math::Vector3d x10 = coords[tet[1]] - coords[tet[0]];
    const stk::math::Vector3d x20 = coords[tet[2]] - coords[tet[0]];
    const stk::math::Vector3d x30 = coords[tet[3]] - coords[tet[0]];

    const stk::math::Vector3d x10_x_x20 = Cross(x10,x20);
    const stk::math::Vector3d x20_x_x30 = Cross(x20,x30);
    const stk::math::Vector3d x30_x_x10 = Cross(x30,x10);

    const double detJ = Dot(x30,x10_x_x20);
    const double vol = detJ/6.;
    const stk::math::Vector3d gradU = (1./detJ)*(u10*x20_x_x30 + u20*x30_x_x10 + u30*x10_x_x20);
    const stk::math::Vector3d negGradW = (1./detJ)*(x20_x_x30 + x30_x_x10 + x10_x_x20); // = -∇w_0, (∇w_0 found by substituting 1,0,0,0 in for u)

    // lumped mass contribution
    sumVolW += 0.25*vol;
    sumContrib += vol*Dot(negGradW, gradU);
  }

  return sumContrib / sumVolW;
}

std::pair<double,double> tri_area_and_divergence_of_node_normals(
    const stk::math::Vector3d & coords0,
    const stk::math::Vector3d & coords1,
    const stk::math::Vector3d & coords2,
    const stk::math::Vector3d & nodeNormals0,
    const stk::math::Vector3d & nodeNormals1,
    const stk::math::Vector3d & nodeNormals2)
{
  const stk::math::Vector3d x10 = coords1 - coords0;
  const stk::math::Vector3d x20 = coords2 - coords0;

  const stk::math::Vector3d n10 = nodeNormals1 - nodeNormals0;
  const stk::math::Vector3d n20 = nodeNormals2 - nodeNormals0;

  const double detJ = (x10[0]*x20[1]-x20[0]*x10[1]);
  const double area = detJ/2.;
  const stk::math::Vector3d gradContrib10(x20[1],-x20[0],0.0);
  const stk::math::Vector3d gradContrib20(-x10[1],x10[0],0.0);
  const double normalDivergence = (1./detJ)*(Dot(n10,gradContrib10) + Dot(n20,gradContrib20));

  return std::make_pair(area, normalDivergence);
}

std::pair<double,double> tet_volume_and_divergence_of_node_normals(
    const stk::math::Vector3d & coords0,
    const stk::math::Vector3d & coords1,
    const stk::math::Vector3d & coords2,
    const stk::math::Vector3d & coords3,
    const stk::math::Vector3d & nodeNormals0,
    const stk::math::Vector3d & nodeNormals1,
    const stk::math::Vector3d & nodeNormals2,
    const stk::math::Vector3d & nodeNormals3)
{
  const stk::math::Vector3d x10 = coords1 - coords0;
  const stk::math::Vector3d x20 = coords2 - coords0;
  const stk::math::Vector3d x30 = coords3 - coords0;

  const stk::math::Vector3d n10 = nodeNormals1 - nodeNormals0;
  const stk::math::Vector3d n20 = nodeNormals2 - nodeNormals0;
  const stk::math::Vector3d n30 = nodeNormals3 - nodeNormals0;

  const stk::math::Vector3d x10_x_x20 = Cross(x10,x20);
  const stk::math::Vector3d x20_x_x30 = Cross(x20,x30);
  const stk::math::Vector3d x30_x_x10 = Cross(x30,x10);

  const double detJ = Dot(x30,x10_x_x20);
  const double vol = detJ/6.;
  const double normalDivergence = (1./detJ)*(Dot(n10,x20_x_x30) + Dot(n20,x30_x_x10) + Dot(n30,x10_x_x20));
  return std::make_pair(vol, normalDivergence);
}

double laplacian_of_field_at_node_of_patch_of_triangles(
    const unsigned node,
    const std::vector<std::array<unsigned,3>> & tris,
    const std::vector<stk::math::Vector3d> & coords,
    const std::vector<double> & u)
{
  double sumAreaW = 0.0;
  double sumContrib = 0.0;

  for(auto & origTri : tris)
  {
    const auto & tri = permute_tri_so_that_node_is_first(origTri, node);

    const double u10 = u[tri[1]] - u[tri[0]];
    const double u20 = u[tri[2]] - u[tri[0]];
    const stk::math::Vector3d x10 = coords[tri[1]] - coords[tri[0]];
    const stk::math::Vector3d x20 = coords[tri[2]] - coords[tri[0]];

    const double detJ = (x10[0]*x20[1]-x20[0]*x10[1]);
    const double area = detJ/2.;
    const stk::math::Vector3d gradContrib10(x20[1],-x20[0],0.0);
    const stk::math::Vector3d gradContrib20(-x10[1],x10[0],0.0);
    const stk::math::Vector3d gradU = (1./detJ)*(u10*gradContrib10 + u20*gradContrib20);
    const stk::math::Vector3d negGradW = (1./detJ)*(gradContrib10 + gradContrib20); // = -∇w_0, (∇w_0 found by substituting 1,0,0 in for u)

    // lumped mass contribution
    sumAreaW += area/3.;
    sumContrib += area*Dot(negGradW, gradU);
  }

  return sumContrib / sumAreaW;
}

double finite_difference_divergence_of_normal_2d(const double delta, const double f0,
    const double fxm, const double fxp,
    const double fym, const double fyp,
    const double fxmym, const double fxmyp, const double fxpym, const double fxpyp)
{
  const double fxx = (fxp - 2*f0 + fxm)/(delta*delta);
  const double fyy = (fyp - 2*f0 + fym)/(delta*delta);
  const double fx = (fxp - fxm)/(2*delta);
  const double fy = (fyp - fym)/(2*delta);
  const double fxy = (fxpyp - fxpym - fxmyp + fxmym)/(4*delta*delta);
  const double normalDivergence = (fx*fx*fyy + fy*fy*fxx - 2*fx*fy*fxy) / std::pow(fx*fx+fy*fy,1.5);
  return normalDivergence;
}

double finite_difference_divergence_of_normal_3d(const double delta, const double f0,
    const double fxm, const double fxp,
    const double fym, const double fyp,
    const double fzm, const double fzp,
    const double fxmym, const double fxmyp, const double fxpym, const double fxpyp,
    const double fxmzm, const double fxmzp, const double fxpzm, const double fxpzp,
    const double fymzm, const double fymzp, const double fypzm, const double fypzp)
{
  const double fxx = (fxp - 2*f0 + fxm)/(delta*delta);
  const double fyy = (fyp - 2*f0 + fym)/(delta*delta);
  const double fzz = (fzp - 2*f0 + fzm)/(delta*delta);
  const double fx = (fxp - fxm)/(2*delta);
  const double fy = (fyp - fym)/(2*delta);
  const double fz = (fzp - fzm)/(2*delta);
  const double fxy = (fxpyp - fxpym - fxmyp + fxmym)/(4*delta*delta);
  const double fxz = (fxpzp - fxpzm - fxmzp + fxmzm)/(4*delta*delta);
  const double fyz = (fypzp - fypzm - fymzp + fymzm)/(4*delta*delta);
  const double normalDivergence = (fx*fx*(fyy+fzz) + fy*fy*(fxx+fzz) + fz*fz*(fxx+fyy) - 2*(fx*fy*fxy + fx*fz*fxz + fy*fz*fyz)) / std::pow(fx*fx+fy*fy+fz*fz,1.5);
  return normalDivergence;
}

template<size_t N>
double divergence_of_node_normals_on_patch_of_elements(
    const std::vector<std::array<unsigned,N>> & elemConn,
    const std::vector<stk::math::Vector3d> & coords,
    const std::vector<stk::math::Vector3d> & nodeNormals)
{
  double sumAreaOrVol = 0.0;
  double sumContrib = 0.0;
  double areaOrVol = 0.;
  double normalDivergence = 0.;

  for(auto & conn : elemConn)
  {
    if constexpr (N == 3)
    {
      std::tie(areaOrVol, normalDivergence) = tri_area_and_divergence_of_node_normals(
          coords[conn[0]], coords[conn[1]], coords[conn[2]],
          nodeNormals[conn[0]], nodeNormals[conn[1]], nodeNormals[conn[2]]);
    }
    else
    {
      static_assert(N == 4);
      std::tie(areaOrVol, normalDivergence) = tet_volume_and_divergence_of_node_normals(
          coords[conn[0]], coords[conn[1]], coords[conn[2]], coords[conn[3]],
          nodeNormals[conn[0]], nodeNormals[conn[1]], nodeNormals[conn[2]], nodeNormals[conn[3]]);
    }
    sumContrib += areaOrVol*normalDivergence;
    sumAreaOrVol += areaOrVol;
  }

  return sumContrib / sumAreaOrVol;
}

template double divergence_of_node_normals_on_patch_of_elements(const std::vector<std::array<unsigned,3>> & triConn,
    const std::vector<stk::math::Vector3d> & coords,
    const std::vector<stk::math::Vector3d> & nodeNormals);
template double divergence_of_node_normals_on_patch_of_elements(const std::vector<std::array<unsigned,4>> & tetConn,
    const std::vector<stk::math::Vector3d> & coords,
    const std::vector<stk::math::Vector3d> & nodeNormals);
}

