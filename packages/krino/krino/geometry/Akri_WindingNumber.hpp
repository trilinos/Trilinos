#ifndef KRINO_KRINO_GEOMETRY_AKRI_WINDINGNUMBER_HPP_
#define KRINO_KRINO_GEOMETRY_AKRI_WINDINGNUMBER_HPP_

#include <vector>
#include <stk_math/StkVector.hpp>

namespace krino {

struct FacetClusterApproximation {
  void clear()
  {
    center = stk::math::Vector3d::ZERO;
    std::fill(areaN.begin(), areaN.end(), 0.);
    std::fill(areaNiXj.begin(), areaNiXj.end(), 0.);
    std::fill(areaDijNk.begin(), areaDijNk.end(), 0.);
  }

  stk::math::Vector3d center;
  std::array<double,3> areaN;
  std::array<double,9> areaNiXj;
  std::array<double,18> areaDijNk;
};

struct ClusterApproximation
{
    void clear()
    {
      center = stk::math::Vector3d::ZERO;
      areaN = stk::math::Vector3d::ZERO;
      std::fill(areaNiXj.begin(), areaNiXj.end(), 0.);

      areaNDdiag = stk::math::Vector3d::ZERO;
      areaNDperm = 0.;
      area2NxDxy_NyDxx = 0.;
      area2NxDxz_NzDxx = 0.;
      area2NyDyz_NzDyy = 0.;
      area2NyDyx_NxDyy = 0.;
      area2NzDzx_NxDzz = 0.;
      area2NzDzy_NyDzz = 0.;
    }

    stk::math::Vector3d center;
    stk::math::Vector3d areaN;

    // data for 1st order term and correction of 2nd order term
    std::array<double,9> areaNiXj;

    // data for 2nd order term
    stk::math::Vector3d areaNDdiag; // NxDxx, NyDyy, NzDzz
    double areaNDperm;              // (NxDyz+NxDzy+NyDzx+NyDxz+NzDxy+NzDyx) = 2*(NxDyz+NyDzx+NzDxy)
    double area2NxDxy_NyDxx; // NxDxy+NxDyx+NyDxx = 2NxDxy+NyDxx
    double area2NxDxz_NzDxx; // NxDxz+NxDzx+NzDxx = 2NxDxz+NzDxx
    double area2NyDyz_NzDyy; // NyDyz+NyDzy+NzDyy = 2NyDyz+NzDyy
    double area2NyDyx_NxDyy; // NyDyx+NyDxy+NxDyy = 2NyDyx+NxDyy
    double area2NzDzx_NxDzz; // NzDzx+NzDxz+NxDzz = 2NzDzx+NxDzz
    double area2NzDzy_NyDzz; // NzDzy+NzDyz+NyDzz = 2NzDzy+NyDzz
};

double compute_facet_winding_number(const stk::math::Vector3d & x0, const stk::math::Vector3d & x1, const stk::math::Vector3d & x2, const stk::math::Vector3d & queryLoc);
double compute_faceted_surface_winding_number(const std::vector<std::array<stk::math::Vector3d,3>> & surfFacets, const stk::math::Vector3d & queryLoc);

void compute_cluster_approximation(const std::vector<std::array<stk::math::Vector3d,3>> & clusterFacets, const stk::math::Vector3d & clusterPt, FacetClusterApproximation & cluster);
double compute_approximate_winding_number(const FacetClusterApproximation & cluster, const stk::math::Vector3d & queryPt);

void compute_cluster_approximation(const std::vector<std::array<stk::math::Vector3d,3>> & clusterFacets, const stk::math::Vector3d & clusterPt, ClusterApproximation & approx);
double compute_approximate_winding_number(const ClusterApproximation & cluster, const stk::math::Vector3d & queryPt);

}

#endif /* KRINO_KRINO_GEOMETRY_AKRI_WINDINGNUMBER_HPP_ */
