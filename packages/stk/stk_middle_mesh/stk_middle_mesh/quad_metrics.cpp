#include "quad_metrics.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

void QuadMetrics::compute_det_jacobian(mesh::MeshEntityPtr quad, const std::vector<utils::Point>& xiPts, std::vector<double>& detJ)
{
  assert(quad->get_type() == mesh::MeshEntityType::Quad);
  detJ.resize(xiPts.size());

  std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
  mesh::get_downward(quad, 0, verts.data());
  std::array<utils::Point, mesh::MAX_DOWN> vertCoords;
  for (int i=0; i < 4; ++i)
    vertCoords[i] = verts[i]->get_point_orig(0);

  compute_det_jacobian(vertCoords, xiPts, detJ);
}


void QuadMetrics::compute_det_jacobian(std::array<utils::Point, mesh::MAX_DOWN>& vertCoords, const std::vector<utils::Point>& xiPts,
                          std::vector<double>& detJ)
{
  std::array<double, 2> xiVals, xiDerivs, etaVals, etaDerivs;
  std::array<double, 4> basisDeriv, derivBasis;
  for (size_t i=0; i < xiPts.size(); ++i)
  {
    double xi = xiPts[i].x, eta = xiPts[i].y;

    mesh::compute_lagrange_vals(xi, xiVals.data());
    mesh::compute_lagrange_derivs(xi, xiDerivs.data());
    mesh::compute_lagrange_vals(eta, etaVals.data());
    mesh::compute_lagrange_derivs(eta, etaDerivs.data());

    for (int j=0; j < 4; ++j)
    {
      basisDeriv[j] = 0;
      derivBasis[j] = 0;
    }

    // compute l(xi) * d l(eta)/deta
    basisDeriv[0] = xiVals[0] * etaDerivs[0];
    basisDeriv[1] = xiVals[1] * etaDerivs[0];
    basisDeriv[2] = xiVals[1] * etaDerivs[1];
    basisDeriv[3] = xiVals[0] * etaDerivs[1];

    // compute d l(xi)/dxi * l(eta)
    derivBasis[0] = xiDerivs[0] * etaVals[0];
    derivBasis[1] = xiDerivs[1] * etaVals[0];
    derivBasis[2] = xiDerivs[1] * etaVals[1];
    derivBasis[3] = xiDerivs[0] * etaVals[1];

    utils::Point dxdxi(0, 0, 0), dxdeta(0, 0, 0);
    for (int j=0; j < 4; ++j)
    {
      dxdxi  += derivBasis[j] * vertCoords[j];
      dxdeta += basisDeriv[j] * vertCoords[j];
    }

    utils::Point normal = cross(dxdxi, dxdeta);
    detJ[i] = std::sqrt(dot(normal, normal));
  }
}

}
}
}
}
