#ifndef NEWTON_SCALAR_H
#define NEWTON_SCALAR_H

#include <cmath>
#include <iostream>
#include <limits>

#include "backtracking_line_search.hpp"

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

class NewtonScalar
{
  public:
    explicit NewtonScalar(const double tol = 1e-13, const int itermax = 100)
      : m_tol(tol)
      , m_itermax(itermax)
    {}

    // double rhs_func(double): computes f(x). Newton's method attempts
    // solve f(x) = 0
    // double df/dx = jac_func(double x) computes the derivative of f wrt x
    // x0 is updated with the computed solution.
    // returns 0 if Newton's method converged.
    template <typename Trhs, typename Tjac>
    int solve(Trhs rhsFunc, Tjac jacFunc, double& x0)
    {
      // std::cout << "\nEntered Newton2::solve" << std::endl;
      double rhs, dfdx;

      rhs            = rhsFunc(x0);
      double rhsNorm = std::abs(rhs);
      // std::cout << "initial rhs_norm = " << rhs_norm << std::endl;

      int iter = 0;
      while (rhsNorm > m_tol && iter < m_itermax)
      {
        dfdx = jacFunc(x0);
        // x0 -= rhs / dfdx;
        double deltaX = rhs / dfdx;
        globalize(deltaX, rhsFunc, dfdx, x0);

        rhs     = rhsFunc(x0);
        rhsNorm = std::abs(rhs);
        //std::cout << "iteration " << iter << " rhs_norm = " << rhsNorm << std::endl;
        iter++;
      }

      return rhsNorm > m_tol ? 1 : 0;
    }

    template <typename Trhs>
    void globalize(double deltaX, Trhs rhsFunc, double dfdx, double& x)
    {
      auto g = [&](double alpha) {
        double xI     = x - alpha * deltaX;
        double rhsVal = rhsFunc(xI);
        return rhsVal * rhsVal;
      };

      // TODO: we already compute rhs(x), reuse it;
      double rhsI  = rhsFunc(x);
      double grad0 = 2 * rhsI * dfdx * deltaX;

      const double minAlpha        = 0.1;
      const double backtrackingFac = 0.5;
      const double gradientCoeff   = 0.1;
      opt::impl::BacktrackingLineSearch linesearch(backtrackingFac, gradientCoeff);
      double alpha = linesearch.search(g, grad0);
      alpha        = std::max(alpha, minAlpha);

      // std::cout << "alpha = " << alpha << std::endl;

      x -= alpha * deltaX;
    }

  private:
    double m_tol;
    int m_itermax;
};

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
#endif
