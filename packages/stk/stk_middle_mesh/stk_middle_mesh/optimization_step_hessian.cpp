#include "optimization_step_hessian.hpp"

namespace stk {
namespace middle_mesh {
namespace opt {
namespace impl {

// returns updated position of the active vert in parameteric space
utils::Point OptimizationStepHessian::compute_updated_point(ActiveVertData& active)
{
  m_obj->set_active_patch(active);
  std::array<double, 2> deltaX, dfDx;
  utils::impl::Mat2x2<double> hessian;

  mesh::MeshEntityPtr vert = active.get_current_vert();
  auto ptPrime             = m_obj->compute_parameterization(vert->get_point_orig(0));

  // take Newton step
  auto dfDxPt = m_obj->compute_quality_rev(ptPrime);
  dfDx[0]     = dfDxPt.x;
  dfDx[1]     = dfDxPt.y;


  hessian = m_obj->compute_hessian(ptPrime);
  modify_eigenvalues(hessian, 1e-3); 
  matsolve2x2(hessian, deltaX.data(), dfDx.data());

  // inverse2x2(hessian);
  // matvec2x2(hessian, df_dx.data(), delta_x.data());

  deltaX[0] = -deltaX[0];
  deltaX[1] = -deltaX[1];

  auto g = [&](const double alpha) {
    utils::Point pt = utils::Point(ptPrime.x + alpha * deltaX[0], ptPrime.y + alpha * deltaX[1]);
    return m_obj->compute_quality(pt);
  };

  double gGrad0 = dfDx[0] * deltaX[0] + dfDx[1] * deltaX[1];

  double alpha = m_linesearch.search(g, gGrad0);

  ptPrime.x += alpha * deltaX[0];
  ptPrime.y += alpha * deltaX[1];

  return ptPrime;
}

void OptimizationStepHessian::modify_eigenvalues(utils::impl::Mat2x2<double>& a, const double delta)
{
  // compute the modification with the minimum Euclidian norm
  double b = -(a(0, 0) + a(1, 1));
  double c = det2x2(a);

  // the matrix is symmetric, so the eigenvalues have to be real, but
  // sometimes numerical precision pushes them to be slightly imaginary
  double val       = std::sqrt(std::max(b * b - 4 * c, 0.0));
  double lambda1   = (-b + val) / 2;
  double lambda2   = (-b - val) / 2;
  double lambdaMin = std::min(lambda1, lambda2);

  double tau       = std::max(0.0, delta - lambdaMin);
  a(0, 0) += tau;
  a(1, 1) += tau;
}

} // namespace impl

} // namespace opt
} // namespace middle_mesh
} // namespace stk
