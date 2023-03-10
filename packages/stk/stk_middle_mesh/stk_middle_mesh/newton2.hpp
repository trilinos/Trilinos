#ifndef NEWTON2
#define NEWTON2

#include <cmath>
#include <iostream>
#include <limits>

#include "mat2x2.hpp"

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

class Newton2
{
  public:
    explicit Newton2(const double tol = 1e-13, const int itermax = 100)
      : m_tol(tol)
      , m_itermax(itermax)
    {}

    // jac_func(const double x[2], utils::impl::Mat2x2<double>& A[4]) overwrites A with the
    // jacobian
    // rhs_func(const double x[2], rhs[2]) overwrites rhs with the function
    // f(x)
    template <typename Trhs, typename Tjac>
    void solve(Trhs rhsFunc, Tjac jacFunc, double x0[2])
    {
      // std::cout << "\nEntered Newton2::solve" << std::endl;
      double rhs[2], deltaX[2];
      utils::impl::Mat2x2<double> jac;

      rhsFunc(x0, rhs);
      double rhsNorm = norm(rhs);

      int iter = 0;
      while (rhsNorm > m_tol && iter < m_itermax)
      {
        jacFunc(x0, jac);

        matsolve2x2(jac, deltaX, rhs);

        x0[0] -= deltaX[0];
        x0[1] -= deltaX[1];

        rhsFunc(x0, rhs);
        rhsNorm = norm(rhs);
        iter++;
      }

      if (rhsNorm > m_tol)
        std::cout << "Newton2 did not converge" << std::endl;
    }

  private:
    double norm(const double vals[2]) { return std::sqrt(vals[0] * vals[0] + vals[1] * vals[1]); }

    double m_tol;
    int m_itermax;
};

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
#endif
