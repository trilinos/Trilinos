#ifndef NEWTON_H
#define NEWTON_H

#include <cmath>
#include <iostream>
#include <limits>

#include "backtracking_line_search.hpp"
#include "matrix.hpp"

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

class Newton
{
  public:
    explicit Newton(int n, const double tol = 1e-13, const int itermax = 100)
      : m_size(n)
      , m_jac(n, n)
      , m_rhs(n)
      , m_deltaX(n)
      , m_ipiv(n)
      , m_tol(tol)
      , m_itermax(itermax)
    {}

    // rhs_func(const std::vector<double>& x, std::vector<double>& rhs) overwrites rhs with the function
    // f(x)
    // jac_func(const std::vector<double>& x, Matrix<double>& A) overwrites A with the
    // jacobian
    // on entry, x0 contains the initial guess for the solution x, on
    // exit it will contain the solution (or the most recent iterate if Newton
    // does not converge)
    // all vectors of lengh n, matrices are n x n
    // returns 0 if method converged, 1 otherwise
    template <typename Trhs, typename Tjac>
    int solve(Trhs rhsFunc, Tjac jacFunc, std::vector<double>& x0)
    {
      assert(x0.size() == size_t(m_size));

      // std::cout << "\nEntered Newton2::solve" << std::endl;
      rhsFunc(x0, m_rhs);
      double rhsNorm = norm(m_rhs);

      // double damp_factor = 0.5;

      int iter = 0;
      while (rhsNorm > m_tol && iter < m_itermax)
      {
        jacFunc(x0, m_jac);

        solve_linear_system(m_jac, m_ipiv.data(), m_rhs.data());
        // rhs was overwritten with delta_x;

        globalize(m_rhs, rhsFunc, x0);

        rhsFunc(x0, m_rhs);
        rhsNorm = norm(m_rhs);

        iter++;
      }

      int retVal = rhsNorm > m_tol ? 1 : 0;
      return retVal;
    }

  private:
    double norm(const std::vector<double>& vals)
    {
      double val = 0;
      for (int i = 0; i < m_size; ++i)
        val += vals[i] * vals[i];
      return std::sqrt(val);
    }

    template <typename Trhs>
    void globalize(const std::vector<double>& deltaX, Trhs rhsFunc, std::vector<double>& x)
    {
      std::vector<double> xI(m_size), rhsI(m_size);
      auto g = [&](double alpha) {
        for (int i = 0; i < m_size; ++i)
          xI[i] = x[i] - alpha * deltaX[i];

        rhsFunc(xI, rhsI);

        double rhsSquared = 0;
        for (int i = 0; i < m_size; ++i)
          rhsSquared += rhsI[i] * rhsI[i];

        return rhsSquared;
      };

      // TODO: we already computes rhs(x), reuse it
      rhsFunc(x, rhsI);
      double grad0 = 0;
      for (int i = 0; i < m_size; ++i)
        for (int j = 0; j < m_size; ++j)
          grad0 = 2 * rhsI[i] * m_jac(i, j) * deltaX[j];

      const double minAlpha        = 0.1;
      const double backtrackingFac = 0.5;
      const double gradientCoeff   = 0.1;
      opt::impl::BacktrackingLineSearch linesearch(backtrackingFac, gradientCoeff);
      double alpha = linesearch.search(g, grad0);
      alpha        = std::max(alpha, minAlpha);

      for (int i = 0; i < m_size; ++i)
        x[i] -= alpha * deltaX[i];
    }

    int m_size;
    Matrix<double> m_jac;
    std::vector<double> m_rhs, m_deltaX;
    std::vector<int> m_ipiv;
    double m_tol;
    int m_itermax;
};

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
#endif
