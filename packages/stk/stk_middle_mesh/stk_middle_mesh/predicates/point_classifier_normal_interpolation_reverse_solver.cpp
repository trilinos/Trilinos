#include "stk_middle_mesh/predicates/point_classifier_normal_interpolation_reverse_solver.hpp"

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

std::array<double, 3> ReverseClassificationGammaSolver::solve(bool& didConverge)
{
  auto f        = [&](double gamma) { return compute_f_of_gamma(gamma); };
  auto dfdx     = [&](double gamma) { return compute_df_dgamma(gamma); };
  double gamma  = 0.0;
  double retVal = m_newton.solve(f, dfdx, gamma);

  // sometimes f has a large spike at some gamma value and the
  // solution lies beyond that value (when starting from gamma = 0).
  // Try again using a large initial guess for gamma to get past the
  // spike
  if (retVal != 0)
  {
    gamma  = gamma * 100000;
    retVal = m_newton.solve(f, dfdx, gamma);
  }

  didConverge = true;
  if (retVal != 0)
    didConverge = false;
  //  throw std::runtime_error("Newton did not converge");

  std::array<double, 2> lambdas = solve_for_lambda(gamma);

  return {lambdas[0], lambdas[1], gamma};
}

std::array<double, 2> ReverseClassificationGammaSolver::solve_for_lambda(double gamma)
{
  utils::impl::Mat2x2<double> mat;
  mat(0, 0) = m_dx1[0] - gamma * m_dn1[0];
  mat(0, 1) = m_dx2[0] - gamma * m_dn2[0];
  mat(1, 0) = m_dx1[1] - gamma * m_dn1[1];
  mat(1, 1) = m_dx2[1] - gamma * m_dn2[1];

  std::array<double, 2> rhs = {m_pt[0] - m_x0[0] + gamma * m_d0[0], m_pt[1] - m_x0[1] + gamma * m_d0[1]};

  std::array<double, 2> lambdas;
  matsolve2x2(mat, lambdas.data(), rhs.data());

  return lambdas;
}

void ReverseClassificationGammaSolver::solve_for_lambdad_gamma(double gamma, std::array<double, 2>& lambdas,
                                                               std::array<double, 2>& lambdasDot)
{
  utils::impl::Mat2x2<double> mat, matDot;
  mat(0, 0) = m_dx1[0] - gamma * m_dn1[0];
  mat(0, 1) = m_dx2[0] - gamma * m_dn2[0];
  mat(1, 0) = m_dx1[1] - gamma * m_dn1[1];
  mat(1, 1) = m_dx2[1] - gamma * m_dn2[1];

  matDot(0, 0) = -m_dn1[0];
  matDot(0, 1) = -m_dn2[0];
  matDot(1, 0) = -m_dn1[1];
  matDot(1, 1) = -m_dn2[1];

  std::array<double, 2> rhs = {m_pt.x - m_x0.x + gamma * m_d0.x, m_pt.y - m_x0.y + gamma * m_d0.y};

  std::array<double, 2> rhsDot = {m_d0.x, m_d0.y};

  matsolve2x2_dot(mat, matDot, lambdas.data(), lambdasDot.data(), rhs.data(), rhsDot.data());
}

double ReverseClassificationGammaSolver::compute_f_of_gamma(double gamma)
{
  std::array<double, 2> lambdas = solve_for_lambda(gamma);
  double lambda1                = lambdas[0];
  double lambda2                = lambdas[1];

  return -gamma * (m_d0.z + m_dn1.z * lambda1 + m_dn2.z * lambda2) - m_pt.z + m_x0.z + m_dx1.z * lambda1 +
         m_dx2.z * lambda2;
}

double ReverseClassificationGammaSolver::compute_df_dgamma(double gamma)
{
  std::array<double, 2> lambdas, dlambdasDgamma;
  solve_for_lambdad_gamma(gamma, lambdas, dlambdasDgamma);
  const double& lambda1        = lambdas[0];
  const double& lambda2        = lambdas[1];
  const double& dlambda1Dgamma = dlambdasDgamma[0];
  const double& dlambda2Dgamma = dlambdasDgamma[1];

  double partialfPartialgamma = -(m_d0.z + m_dn1.z * lambda1 + m_dn2.z * lambda2);
  double dfDlambda1           = -gamma * m_dn1.z + m_dx1.z;
  double dfDlambda2           = -gamma * m_dn2.z + m_dx2.z;

  return partialfPartialgamma + dfDlambda1 * dlambda1Dgamma + dfDlambda2 * dlambda2Dgamma;
}

} // namespace impl
} // namespace predicates
} // namespace middle_mesh
} // namespace stk
