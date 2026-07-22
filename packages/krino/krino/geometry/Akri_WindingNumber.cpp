#include <Akri_Triangle.hpp>
#include <stk_math/StkVector.hpp>
#include <Akri_WindingNumber.hpp>
#include <vector>

namespace krino {

double compute_facet_winding_number(const stk::math::Vector3d & facetX0, const stk::math::Vector3d & facetX1, const stk::math::Vector3d & facetX2, const stk::math::Vector3d & queryLoc)
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

std::array<double,6> compute_facet3d_edge_contrib(const std::array<stk::math::Vector3d,3> & facetCoords, const stk::math::Vector3d & centroid)
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
      (dx[0]*(childNyXz+childNzXy + n[1]*dx[2]+n[2]*dx[1]) +
       dx[1]*(childNxXz+childNzXx + n[0]*dx[2]+n[2]*dx[0]) +
       dx[2]*(childNxXy+childNyXx + n[0]*dx[1]+n[1]*dx[0]));
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

double compute_faceted_surface_winding_number(const std::vector<std::array<stk::math::Vector3d,3>> & surfFacets, const stk::math::Vector3d & queryLoc)
{
  double windingNumber = 0;
  for (const auto & facetCoords : surfFacets)
      windingNumber += compute_facet_winding_number(facetCoords[0], facetCoords[1], facetCoords[2], queryLoc);
  return windingNumber;
}

}
