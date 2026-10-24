// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_WindingNumberCalcs.hpp>

namespace krino {

// There are two cluster-approximation types here, FacetClusterApproximation and
// ClusterApproximation. Both evaluate the SAME 2nd-order multipole expansion of
// the winding number about a cluster center and agree to round-off (see the
// unit test that cross-checks them). They differ only in how the moment data is
// stored and, more importantly, in how a cluster is assembled:
//
//   FacetClusterApproximation (direct form) -- 33 doubles
//     Stores the raw 2nd-order moment tensor (areaDijNk[18]) and is accumulated
//     directly from facet vertices about a fixed center via
//     accumulate_cluster_approximation(x0,x1,x2,center,...). It has no notion of
//     combining two existing clusters, so a cluster can only be built by summing
//     over its own facets.
//
//   ClusterApproximation (composable / shiftable form) -- 25 doubles
//     Stores the 2nd-order data as pre-combined named scalars and additionally
//     provides accumulate_cluster_approximation(childApprox, parentApprox),
//     which merges a child cluster's moments into a parent about a DIFFERENT
//     center by applying parallel-axis shift terms. A cluster can therefore be
//     built either from facets (via compute_facet_approximation) or by combining
//     already-computed sub-clusters.
//
//   Storage: the two agree on the first 15 doubles (center + areaN + areaNiXj);
//   they differ only in the 2nd-order block. FacetClusterApproximation keeps the
//   full raw tensor (areaDijNk[18]), while ClusterApproximation pre-contracts it
//   into the 10 combinations the evaluator actually needs (areaNDdiag[3],
//   areaNDperm, and 6 mixed area2* scalars), making it ~24% smaller.
//
// How the two trees use these types:
//   WindingNumberTopDownTree uses FacetClusterApproximation and builds top-down:
//   every internal node re-sums ALL facets in its range from scratch, so each
//   facet is visited once per level of the tree -> O(N log N) build work.
//
//   WindingNumberBottomUpTree uses ClusterApproximation and builds post-order:
//   each leaf's moments are computed once from its facets (via
//   compute_facet_approximation, O(N) total), then each internal node's moments
//   are formed as the O(1) merge of its two children -> O(N) build overall.
//   Query cost and accuracy are identical (the evaluated expansion is the same);
//   only build time and per-node storage differ (~24% smaller, 25 vs 33
//   doubles). The bottom-up build stores a dense approximation per node, so it
//   has no "exact leaves carry no approximation" policy.
static std::array<double,6> compute_facet3d_edge_contrib(const std::array<stk::math::Vector3d,3> & facetCoords, const stk::math::Vector3d & centroid)
{
  std::array<double,6> edgeContrib = {0,0,0,0,0,0};
  const stk::math::Vector3d edgeDiff0 = 0.5*(facetCoords[0]+facetCoords[1]) - centroid;
  const stk::math::Vector3d edgeDiff1 = 0.5*(facetCoords[1]+facetCoords[2]) - centroid;
  const stk::math::Vector3d edgeDiff2 = 0.5*(facetCoords[2]+facetCoords[0]) - centroid;

  // Note that this only fills in matrix for j>=i to take advantage of symmetry
  unsigned index = 0;
  for (unsigned i=0; i<3; ++i)
    for (unsigned j=i; j<3; ++j)
      edgeContrib[index++] = (1./3.)*(edgeDiff0[i]*edgeDiff0[j] + edgeDiff1[i]*edgeDiff1[j] + edgeDiff2[i]*edgeDiff2[j]);
  return edgeContrib;
}

double compute_exact_facet_winding_number(const stk::math::Vector3d & facetX0, const stk::math::Vector3d & facetX1, const stk::math::Vector3d & facetX2, const stk::math::Vector3d & queryLoc)
{
    const stk::math::Vector3d x0 = facetX0 - queryLoc;
    const stk::math::Vector3d x1 = facetX1 - queryLoc;
    const stk::math::Vector3d x2 = facetX2 - queryLoc;
    const double len0 = x0.length();
    const double len1 = x1.length();
    const double len2 = x2.length();
    const double num = Dot(Cross(x0,x1), x2);
    const double den = len0*len1*len2 + Dot(x0,x1)*len2 + Dot(x0,x2)*len1 + Dot(x1,x2)*len0;
    return std::atan2(num, den) / (2.*M_PI);
}

void accumulate_cluster_approximation(const stk::math::Vector3d & facetX0,
    const stk::math::Vector3d & facetX1,
    const stk::math::Vector3d & facetX2,
    const stk::math::Vector3d & clusterPt,
    FacetClusterApproximation & approx)
{
  const stk::math::Vector3d facetAreaNormal = 0.5*Cross(facetX1-facetX0,facetX2-facetX0);

  for (int i=0; i<3; ++i)
    approx.areaN[i] += facetAreaNormal[i];

  const stk::math::Vector3d centroidContrib = ((1./3.) * (facetX0+facetX1+facetX2) - clusterPt);
  for (unsigned i=0; i<3; ++i)
    for (unsigned j=0; j<3; ++j)
      approx.areaNiXj[3*i+j] += facetAreaNormal[i]*centroidContrib[j];

  const stk::math::Vector3d edgeDiff0 = 0.5*(facetX0+facetX1) - clusterPt;
  const stk::math::Vector3d edgeDiff1 = 0.5*(facetX1+facetX2) - clusterPt;
  const stk::math::Vector3d edgeDiff2 = 0.5*(facetX2+facetX0) - clusterPt;

  // Note that this only fills in matrix for j>=i to take advantage of symmetry
  unsigned index = 0;
  for (unsigned i=0; i<3; ++i)
    for (unsigned j=i; j<3; ++j)
      for (unsigned k=0; k<3; ++k)
        approx.areaDijNk[index++] += (1./3.) * (edgeDiff0[i]*edgeDiff0[j] + edgeDiff1[i]*edgeDiff1[j] + edgeDiff2[i]*edgeDiff2[j]) * facetAreaNormal[k];
}

double compute_approximate_winding_number(const FacetClusterApproximation & cluster, const stk::math::Vector3d & queryPt)
{
  stk::math::Vector3d dx = cluster.center - queryPt;
  const double len = dx.unitize();

  static constexpr double inv4Pi = 0.25/M_PI;
  const double invLen  = 1./len;
  const double invLen2  = invLen*invLen;
  const double mult = inv4Pi * invLen2;
  double windingNumber = 0.;

  unsigned index = 0;
  for (unsigned i=0; i<3; ++i)
  {
    windingNumber += mult * (cluster.areaN[i]*dx[i] + invLen * cluster.areaNiXj[3*i+i]);
    for (unsigned j=0; j<3; ++j)
    {
      windingNumber -= 3. * mult * invLen* (cluster.areaNiXj[3*i+j]*dx[i]*dx[j]);
      if (j>=i)
      {
        const double deltaij = (i==j) ? 1. : 0.;
        const double symm = (i==j) ? 1.0 : 2.0; // use i-j symmetry
        for (unsigned k=0; k<3; ++k)
        {
          const double deltajk = (j==k) ? 1. : 0.;
          const double deltaik = (i==k) ? 1. : 0.;
          double coeff = 7.5 * dx[i]*dx[j]*dx[k] - 1.5*(dx[i]*deltajk + dx[j]*deltaik + dx[k]*deltaij);
          windingNumber += mult * symm * invLen2 * coeff * cluster.areaDijNk[index++];
        }
      }
    }
  }

  return windingNumber;
}

double compute_approximate_winding_number(const ClusterApproximation & approx, const stk::math::Vector3d & queryPt)
{
  stk::math::Vector3d dx = approx.center - queryPt;
  const double len = dx.unitize();

  static constexpr double inv4Pi = 0.25/M_PI;
  const double invLen  = 1./len;
  const double invLen2  = invLen*invLen;
  const double mult = inv4Pi * invLen2;

  double approxWindingNumber0 = mult*Dot(dx, approx.areaN);

  const double mult1 = mult*invLen;
  double approxWindingNumber1 = 0.;
  for (unsigned i=0; i<3; ++i)
  {
    approxWindingNumber1 += mult1 * approx.areaNiXj[3*i+i];
    for (unsigned j=0; j<3; ++j)
      approxWindingNumber1 -= 3. * mult1 * approx.areaNiXj[3*i+j] * dx[i]*dx[j];
  }

  const stk::math::Vector3d dx2(dx[0]*dx[0], dx[1]*dx[1], dx[2]*dx[2]);
  const stk::math::Vector3d dx3(dx[0]*dx2[0], dx[1]*dx2[1], dx[2]*dx2[2]);
  const stk::math::Vector3d deltaFnCoeff(
    3.*approx.areaNDdiag[0] + approx.area2NyDyx_NxDyy + approx.area2NzDzx_NxDzz,
    3.*approx.areaNDdiag[1] + approx.area2NzDzy_NyDzz + approx.area2NxDxy_NyDxx,
    3.*approx.areaNDdiag[2] + approx.area2NxDxz_NzDxx + approx.area2NyDyz_NzDyy
  );
  const stk::math::Vector3d tmp(
    dx[1]*approx.area2NxDxy_NyDxx + dx[2]*approx.area2NxDxz_NzDxx,
    dx[2]*approx.area2NyDyz_NzDyy + dx[0]*approx.area2NyDyx_NxDyy,
    dx[0]*approx.area2NzDzx_NxDzz + dx[1]*approx.area2NzDzy_NyDzz
  );
  const double approxWindingNumber2 = mult*invLen2*
    (7.5*(Dot(dx3, approx.areaNDdiag) + dx[0]*dx[1]*dx[2]*approx.areaNDperm + Dot(dx2, tmp)) - 1.5*Dot(dx, deltaFnCoeff));

  const double approxWindingNumber = approxWindingNumber0 + approxWindingNumber1 + approxWindingNumber2;

  return approxWindingNumber;
}

void compute_facet_approximation(const std::array<stk::math::Vector3d,3> & facetCoords, ClusterApproximation & approx)
{
  approx.center = (1./3.) * (facetCoords[0]+facetCoords[1]+facetCoords[2]);
  approx.areaN = 0.5*Cross(facetCoords[1]-facetCoords[0],facetCoords[2]-facetCoords[0]);

  std::fill(approx.areaNiXj.begin(), approx.areaNiXj.end(), 0.);

  const stk::math::Vector3d & n = approx.areaN;
  std::array<double,6> edgeContrib = compute_facet3d_edge_contrib(facetCoords, approx.center);
  approx.areaNDdiag = stk::math::Vector3d(n[0]*edgeContrib[0], n[1]*edgeContrib[3], n[2]*edgeContrib[5]);
  approx.areaNDperm = 2*(n[0]*edgeContrib[4] + n[1]*edgeContrib[2] + n[2]*edgeContrib[1]);
  approx.area2NxDxy_NyDxx = 2*n[0]*edgeContrib[1] + n[1]*edgeContrib[0];
  approx.area2NxDxz_NzDxx = 2*n[0]*edgeContrib[2] + n[2]*edgeContrib[0];
  approx.area2NyDyz_NzDyy = 2*n[1]*edgeContrib[4] + n[2]*edgeContrib[3];
  approx.area2NyDyx_NxDyy = 2*n[1]*edgeContrib[1] + n[0]*edgeContrib[3];
  approx.area2NzDzx_NxDzz = 2*n[2]*edgeContrib[2] + n[0]*edgeContrib[5];
  approx.area2NzDzy_NyDzz = 2*n[2]*edgeContrib[4] + n[1]*edgeContrib[5];
}

void accumulate_cluster_approximation(const ClusterApproximation & childApprox, ClusterApproximation & approx)
{
  approx.areaN += childApprox.areaN;

  const stk::math::Vector3d dx = childApprox.center - approx.center;

  const stk::math::Vector3d & n = childApprox.areaN;

  for (unsigned i=0; i<3; ++i)
    for (unsigned j=0; j<3; ++j)
      approx.areaNiXj[3*i+j] += childApprox.areaNiXj[3*i+j] + n[i]*dx[j];

  const double childNxXx = childApprox.areaNiXj[3*0+0];
  const double childNxXy = childApprox.areaNiXj[3*0+1];
  const double childNxXz = childApprox.areaNiXj[3*0+2];
  const double childNyXx = childApprox.areaNiXj[3*1+0];
  const double childNyXy = childApprox.areaNiXj[3*1+1];
  const double childNyXz = childApprox.areaNiXj[3*1+2];
  const double childNzXx = childApprox.areaNiXj[3*2+0];
  const double childNzXy = childApprox.areaNiXj[3*2+1];
  const double childNzXz = childApprox.areaNiXj[3*2+2];

  approx.areaNDdiag += childApprox.areaNDdiag +
      stk::math::Vector3d(2*dx[0]*childNxXx, 2*dx[1]*childNyXy, 2*dx[2]*childNzXz) +
      stk::math::Vector3d(dx[0]*dx[0]*n[0], dx[1]*dx[1]*n[1], dx[2]*dx[2]*n[2]);
  approx.areaNDperm += childApprox.areaNDperm +
      (dx[0]*(2*(childNyXz+childNzXy) + n[1]*dx[2]+n[2]*dx[1]) +
       dx[1]*(2*(childNxXz+childNzXx) + n[0]*dx[2]+n[2]*dx[0]) +
       dx[2]*(2*(childNxXy+childNyXx) + n[0]*dx[1]+n[1]*dx[0]));
  approx.area2NxDxy_NyDxx += childApprox.area2NxDxy_NyDxx + 2*(dx[1]*childNxXx + dx[0]*childNxXy + n[0]*dx[0]*dx[1]) + 2*childNyXx*dx[0] + n[1]*dx[0]*dx[0];
  approx.area2NxDxz_NzDxx += childApprox.area2NxDxz_NzDxx + 2*(dx[2]*childNxXx + dx[0]*childNxXz + n[0]*dx[0]*dx[2]) + 2*childNzXx*dx[0] + n[2]*dx[0]*dx[0];
  approx.area2NyDyz_NzDyy += childApprox.area2NyDyz_NzDyy + 2*(dx[2]*childNyXy + dx[1]*childNyXz + n[1]*dx[1]*dx[2]) + 2*childNzXy*dx[1] + n[2]*dx[1]*dx[1];
  approx.area2NyDyx_NxDyy += childApprox.area2NyDyx_NxDyy + 2*(dx[0]*childNyXy + dx[1]*childNyXx + n[1]*dx[1]*dx[0]) + 2*childNxXy*dx[1] + n[0]*dx[1]*dx[1];
  approx.area2NzDzx_NxDzz += childApprox.area2NzDzx_NxDzz + 2*(dx[0]*childNzXz + dx[2]*childNzXx + n[2]*dx[2]*dx[0]) + 2*childNxXz*dx[2] + n[0]*dx[2]*dx[2];
  approx.area2NzDzy_NyDzz += childApprox.area2NzDzy_NyDzz + 2*(dx[1]*childNzXz + dx[2]*childNzXy + n[2]*dx[2]*dx[1]) + 2*childNyXz*dx[2] + n[1]*dx[2]*dx[2];
}

void compute_cluster_approximation(const std::vector<std::array<stk::math::Vector3d,3>> & clusterFacets, const stk::math::Vector3d & clusterPt, ClusterApproximation & approx)
{
  approx.clear();
  approx.center = clusterPt;
  for (const auto & facetCoords : clusterFacets)
  {
    ClusterApproximation facetApprox;
    compute_facet_approximation(facetCoords, facetApprox);
    accumulate_cluster_approximation(facetApprox, approx);
  }
}

void compute_cluster_approximation(const std::vector<std::array<stk::math::Vector3d,3>> & clusterFacets, const stk::math::Vector3d & clusterPt, FacetClusterApproximation & cluster)
{
  cluster.clear();
  cluster.center = clusterPt;
  for (const auto & facetCoords : clusterFacets)
    accumulate_cluster_approximation(facetCoords[0], facetCoords[1], facetCoords[2], cluster.center, cluster);
}

//
// 2D implementations.
//
// A 2D "facet" is an oriented segment x0 -> x1. Only the x,y components are
// used. The winding-number field is expanded to 2nd order about a cluster
// center c using the length-normal N = (e_y, -e_x), e = x1 - x0, and the segment
// position moments about c:
//   P_j    = (midpoint - c)_j                                 (first moment)
//   Q_jk   = (1/3)(a_j a_k + b_j b_k) + (1/6)(a_j b_k + a_k b_j)  (second moment)
// with a = x0 - c, b = x1 - c. The cluster stores areaN = sum N,
// areaNiXj[2i+j] = sum N_i P_j, and areaNiDjk[3i+m] = sum N_i Q_jk (m: xx=0,
// xy=1, yy=2). Both 2D cluster types share this layout.

// Accumulate one segment's raw moments about clusterPt into the shared arrays.
static void accumulate_segment_moments(
    const stk::math::Vector3d & segX0,
    const stk::math::Vector3d & segX1,
    const stk::math::Vector3d & clusterPt,
    std::array<double,2> & areaN,
    std::array<double,4> & areaNiXj,
    std::array<double,6> & areaNiDjk)
{
  const std::array<double,2> e = {segX1[0]-segX0[0], segX1[1]-segX0[1]};
  const std::array<double,2> n = {e[1], -e[0]}; // length-normal (e_y, -e_x)

  areaN[0] += n[0];
  areaN[1] += n[1];

  const std::array<double,2> a = {segX0[0]-clusterPt[0], segX0[1]-clusterPt[1]};
  const std::array<double,2> b = {segX1[0]-clusterPt[0], segX1[1]-clusterPt[1]};

  const std::array<double,2> p = {0.5*(a[0]+b[0]), 0.5*(a[1]+b[1])}; // midpoint - clusterPt

  for (unsigned i=0; i<2; ++i)
    for (unsigned j=0; j<2; ++j)
      areaNiXj[2*i+j] += n[i]*p[j];

  // Symmetric second moment Q, stored as {xx, xy, yy}.
  const double qxx = (1./3.)*(a[0]*a[0] + b[0]*b[0]) + (1./6.)*(2.*a[0]*b[0]);
  const double qxy = (1./3.)*(a[0]*a[1] + b[0]*b[1]) + (1./6.)*(a[0]*b[1] + a[1]*b[0]);
  const double qyy = (1./3.)*(a[1]*a[1] + b[1]*b[1]) + (1./6.)*(2.*a[1]*b[1]);

  for (unsigned i=0; i<2; ++i)
  {
    areaNiDjk[3*i+0] += n[i]*qxx;
    areaNiDjk[3*i+1] += n[i]*qxy;
    areaNiDjk[3*i+2] += n[i]*qyy;
  }
}

// Evaluate the 2nd-order expansion given the shared moment arrays. Shared by
// both 2D cluster types since they store identical moment data.
static double evaluate_2d_approximate_winding_number(
    const stk::math::Vector3d & center,
    const std::array<double,2> & areaN,
    const std::array<double,4> & areaNiXj,
    const std::array<double,6> & areaNiDjk,
    const stk::math::Vector3d & queryPt)
{
  const double dx = center[0] - queryPt[0];
  const double dy = center[1] - queryPt[1];
  const double rho2 = dx*dx + dy*dy;

  if (rho2 <= 0.)
    return 0.;

  const double invRho2 = 1./rho2;

  // Order 0: (areaN . d) / rho^2
  double twoPiW = invRho2 * (dx*areaN[0] + dy*areaN[1]);

  // Order 1: [ tr(A1) - 2 d^T A1 d / rho^2 ] / rho^2
  const double a00 = areaNiXj[0], a01 = areaNiXj[1], a10 = areaNiXj[2], a11 = areaNiXj[3];
  const double trA1 = a00 + a11;
  const double dA1d = dx*dx*a00 + dx*dy*(a01+a10) + dy*dy*a11;
  twoPiW += invRho2 * (trA1 - 2.*invRho2*dA1d);

  // Order 2: (1/rho^2)[ -2 T1/rho^2 - T2/rho^2 + 4 T3/rho^4 ]
  const double m0xx = areaNiDjk[0], m0xy = areaNiDjk[1], m0yy = areaNiDjk[2];
  const double m1xx = areaNiDjk[3], m1xy = areaNiDjk[4], m1yy = areaNiDjk[5];

  const double t1 = dx*m0xx + dy*m0xy + dx*m1xy + dy*m1yy;
  const double t2 = dx*(m0xx + m0yy) + dy*(m1xx + m1yy);
  const double t3 = dx*(dx*dx*m0xx + 2.*dx*dy*m0xy + dy*dy*m0yy)
                  + dy*(dx*dx*m1xx + 2.*dx*dy*m1xy + dy*dy*m1yy);

  twoPiW += invRho2 * (-2.*t1*invRho2 - t2*invRho2 + 4.*t3*invRho2*invRho2);

  static constexpr double inv2Pi = 0.5/M_PI;
  return inv2Pi * twoPiW;
}

double compute_exact_facet2d_winding_number(const stk::math::Vector3d & x0, const stk::math::Vector3d & x1, const stk::math::Vector3d & queryLoc)
{
  const double ax = x0[0]-queryLoc[0], ay = x0[1]-queryLoc[1];
  const double bx = x1[0]-queryLoc[0], by = x1[1]-queryLoc[1];
  const double cross = ax*by - ay*bx;
  const double dot   = ax*bx + ay*by;
  return std::atan2(cross, dot) / (2.*M_PI);
}

void accumulate_cluster_approximation(const stk::math::Vector3d & segX0,
    const stk::math::Vector3d & segX1,
    const stk::math::Vector3d & clusterPt,
    Cluster2dApproximation & approx)
{
  accumulate_segment_moments(segX0, segX1, clusterPt, approx.areaN, approx.areaNiXj, approx.areaNiDjk);
}

double compute_approximate_winding_number(const Cluster2dApproximation & cluster, const stk::math::Vector3d & queryPt)
{
  return evaluate_2d_approximate_winding_number(cluster.center, cluster.areaN, cluster.areaNiXj, cluster.areaNiDjk, queryPt);
}

void compute_facet_approximation(const std::array<stk::math::Vector3d,2> & segmentCoords, Cluster2dApproximation & approx)
{
  approx.clear();
  approx.center = 0.5*(segmentCoords[0] + segmentCoords[1]);
  accumulate_segment_moments(segmentCoords[0], segmentCoords[1], approx.center, approx.areaN, approx.areaNiXj, approx.areaNiDjk);
}

void accumulate_cluster_approximation(const Cluster2dApproximation & childApprox, Cluster2dApproximation & approx)
{
  // Parallel-axis shift of the child's moments from its center to the parent's.
  const double deltaX = childApprox.center[0] - approx.center[0];
  const double deltaY = childApprox.center[1] - approx.center[1];
  const std::array<double,2> delta = {deltaX, deltaY};

  const std::array<double,2> & childN = childApprox.areaN;
  const std::array<double,4> & childA = childApprox.areaNiXj;

  // First moment: A'_ij = A_ij + delta_j * areaN_i
  for (unsigned i=0; i<2; ++i)
    for (unsigned j=0; j<2; ++j)
      approx.areaNiXj[2*i+j] += childA[2*i+j] + delta[j]*childN[i];

  // Second moment: M'_ijk = M_ijk + delta_j A_ik + delta_k A_ij + delta_j delta_k areaN_i
  // m index: xx->(j,k)=(0,0), xy->(0,1), yy->(1,1)
  for (unsigned i=0; i<2; ++i)
  {
    approx.areaNiDjk[3*i+0] += childApprox.areaNiDjk[3*i+0]
        + 2.*delta[0]*childA[2*i+0] + delta[0]*delta[0]*childN[i];
    approx.areaNiDjk[3*i+1] += childApprox.areaNiDjk[3*i+1]
        + delta[0]*childA[2*i+1] + delta[1]*childA[2*i+0] + delta[0]*delta[1]*childN[i];
    approx.areaNiDjk[3*i+2] += childApprox.areaNiDjk[3*i+2]
        + 2.*delta[1]*childA[2*i+1] + delta[1]*delta[1]*childN[i];
  }

  approx.areaN[0] += childN[0];
  approx.areaN[1] += childN[1];
}

void compute_cluster_approximation(const std::vector<std::array<stk::math::Vector3d,2>> & clusterSegments, const stk::math::Vector3d & clusterPt, Cluster2dApproximation & approx)
{
  approx.clear();
  approx.center = clusterPt;
  for (const auto & segmentCoords : clusterSegments)
    accumulate_cluster_approximation(segmentCoords[0], segmentCoords[1], approx.center, approx);
}

}


