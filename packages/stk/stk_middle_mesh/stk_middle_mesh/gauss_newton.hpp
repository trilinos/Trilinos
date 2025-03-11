#ifndef STK_MIDDLE_MESH_GAUSS_NEWTON
#define STK_MIDDLE_MESH_GAUSS_NEWTON

#include <vector>
#include "matrix.hpp"

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

// Solves a nonlinear least squares problems min_x \sum_i r(x)^2
class GaussNewton
{
  public:
    GaussNewton(int m, int n, double tol, int itermax) :
      m_jac(m, n),
      m_residuals(m),
      m_residualsForLapack(m, 1),
      m_lapackWork(m, n),
      m_jacTransposeTimesResiduals(n),
      m_tol(tol),
      m_itermax(itermax)
    {
      assert(m >= n);
    }

    // func must be a function func(const std::vector<double>& x, std::vector<double>& residuals)
    // that computes the residuals r(x)
    // jacFunc must be a function jacFunc(const std::vector<double>& x, Matrix<double>& jac) that
    // computes the rectangular jacobian matrix  J_ij = dr_i(x)/dx_j
    // On entry x0 is the initial guess for the solution.  On exit, it will be the solution computed
    // by Gauss Newton
    // Returns true if the method converged.  Returns false if the method did not converge
    // and allowFailure is true.  Throws and exception otherwise
    template <typename Tfunc, typename Tjac>
    bool solve(Tfunc func, Tjac jacFunc, std::vector<double>& x0, bool allowFailure=false)
    {
      assert(x0.size() == size_t(get_num_variables()));
      func(x0, m_residuals);
      jacFunc(x0, m_jac);

      for (int i=0; i < m_itermax; ++i)
      {
        if (is_converged(m_jac, m_residuals))
          return true;

        for (int j=0; j < get_num_residuals(); ++j)
          m_residualsForLapack(j, 0) = m_residuals[j];

        solve_least_squares(m_jac, m_residualsForLapack, m_lapackWork);
        // now m_residualsForLapack contains delta_x

        for (int j=0; j < get_num_variables(); ++j)
          x0[j] -= m_residualsForLapack(j, 0);

        func(x0, m_residuals);
        jacFunc(x0, m_jac);
      }

      if (is_converged(m_jac, m_residuals))
        return true;

      if (allowFailure)
        return false;
      else
        throw std::runtime_error("Gauss Newton did not converge");
    }

  private:

    bool is_converged(const Matrix<double>& jac, const std::vector<double>& /*residuals*/)
    {        
      matvec(1, jac, m_residuals.data(), 0, m_jacTransposeTimesResiduals.data(), BlasTrans::Trans);

      double norm = 0.0;
      for (int i=0; i < get_num_variables(); ++i)
        norm += m_jacTransposeTimesResiduals[i]*m_jacTransposeTimesResiduals[i];

      return norm < m_tol*m_tol;
    }

    int get_num_residuals() const { return m_jac.extent0(); }

    int get_num_variables() const { return m_jac.extent1(); }

    Matrix<double> m_jac;
    std::vector<double> m_residuals;
    Matrix<double> m_residualsForLapack;
    Matrix<double> m_lapackWork;
    std::vector<double> m_jacTransposeTimesResiduals;
    double m_tol;
    int m_itermax;
};

}
}
}
}

#endif