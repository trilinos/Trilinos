// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_GEOMETRY_AKRI_WINDINGNUMBERCALCS_HPP_
#define KRINO_KRINO_GEOMETRY_AKRI_WINDINGNUMBERCALCS_HPP_

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

void accumulate_cluster_approximation(const stk::math::Vector3d & facetX0,
    const stk::math::Vector3d & facetX1,
    const stk::math::Vector3d & facetX2,
    const stk::math::Vector3d & clusterPt,
    FacetClusterApproximation & approx);

void compute_cluster_approximation(const std::vector<std::array<stk::math::Vector3d,3>> & clusterFacets, const stk::math::Vector3d & clusterPt, FacetClusterApproximation & cluster);
double compute_approximate_winding_number(const FacetClusterApproximation & cluster, const stk::math::Vector3d & queryPt);

void compute_cluster_approximation(const std::vector<std::array<stk::math::Vector3d,3>> & clusterFacets, const stk::math::Vector3d & clusterPt, ClusterApproximation & approx);
double compute_approximate_winding_number(const ClusterApproximation & cluster, const stk::math::Vector3d & queryPt);

// Build a single-facet approximation about the facet centroid. This is the leaf
// primitive for a bottom-up (post-order) cluster build.
void compute_facet_approximation(const std::array<stk::math::Vector3d,3> & facetCoords, ClusterApproximation & approx);

// Merge a child cluster's moments into a parent cluster about the parent's
// (different) center by applying parallel-axis shift terms. The parent's center
// must already be set; the child carries its own center. This is the O(1) merge
// that lets a parent cluster be assembled from its children instead of being
// re-summed from all of its facets.
void accumulate_cluster_approximation(const ClusterApproximation & childApprox, ClusterApproximation & approx);

double compute_exact_facet_winding_number(const stk::math::Vector3d & x0, const stk::math::Vector3d & x1, const stk::math::Vector3d & x2, const stk::math::Vector3d & queryLoc);

//
// 2D winding number.
//
// The 2D analogs of the 3D types/functions above. A "facet" in 2D is an oriented
// line segment (two endpoints); the winding number about a query point is the
// signed subtended angle / 2*pi, and a cluster carries a 2nd-order multipole
// (Taylor) expansion of that field about a center. Only the x,y components of
// each stk::math::Vector3d are used (matching krino Facet2d storage).
//
// Unlike 3D (where the direct and composable forms differ in size and layout),
// the 2D second-moment tensor is small enough that a single type serves both
// build strategies: Cluster2dApproximation can be accumulated directly from
// segments about a fixed center (top-down build) AND built from a single segment
// (compute_facet_approximation) then merged child-into-parent about a different
// center (accumulate_cluster_approximation, via parallel-axis shift terms) for an
// O(N) bottom-up build.
//
// Moment layout (raw, about the cluster center c):
//   areaN[i]        = sum of length-normals N = (e_y, -e_x), e = x1 - x0
//   areaNiXj[2*i+j] = sum of N_i * (segment_midpoint - c)_j
//   areaNiDjk[3*i+jk] = sum of N_i * Q_jk, where Q is the 2nd position moment of
//                       the segment about c and jk in {xx=0, xy=1, yy=2}.
//

struct Cluster2dApproximation
{
  void clear()
  {
    center = stk::math::Vector3d::ZERO;
    std::fill(areaN.begin(), areaN.end(), 0.);
    std::fill(areaNiXj.begin(), areaNiXj.end(), 0.);
    std::fill(areaNiDjk.begin(), areaNiDjk.end(), 0.);
  }

  stk::math::Vector3d center;
  std::array<double,2> areaN;
  std::array<double,4> areaNiXj;
  std::array<double,6> areaNiDjk;
};

// Accumulate one segment's moments directly into a cluster about clusterPt (the
// top-down build primitive).
void accumulate_cluster_approximation(const stk::math::Vector3d & segX0,
    const stk::math::Vector3d & segX1,
    const stk::math::Vector3d & clusterPt,
    Cluster2dApproximation & approx);

void compute_cluster_approximation(const std::vector<std::array<stk::math::Vector3d,2>> & clusterSegments, const stk::math::Vector3d & clusterPt, Cluster2dApproximation & approx);
double compute_approximate_winding_number(const Cluster2dApproximation & cluster, const stk::math::Vector3d & queryPt);

// Build a single-segment approximation about the segment midpoint. This is the
// leaf primitive for a bottom-up (post-order) cluster build.
void compute_facet_approximation(const std::array<stk::math::Vector3d,2> & segmentCoords, Cluster2dApproximation & approx);

// Merge a child cluster's moments into a parent cluster about the parent's
// (different) center by applying parallel-axis shift terms. The parent's center
// must already be set; the child carries its own center.
void accumulate_cluster_approximation(const Cluster2dApproximation & childApprox, Cluster2dApproximation & approx);

double compute_exact_facet2d_winding_number(const stk::math::Vector3d & x0, const stk::math::Vector3d & x1, const stk::math::Vector3d & queryLoc);

}

#endif /* KRINO_KRINO_GEOMETRY_AKRI_WINDINGNUMBERCALCS_HPP_ */
