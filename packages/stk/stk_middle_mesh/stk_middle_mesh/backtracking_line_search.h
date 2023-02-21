#ifndef BACKTRACKING_LINE_SEARCH_H
#define BACKTRACKING_LINE_SEARCH_H

#include <stdexcept>

namespace stk {
namespace middle_mesh {
namespace opt {
namespace impl {

class BacktrackingLineSearch
{
  public:
    BacktrackingLineSearch(const double rho, const double c, const double alpha0 = 1)
      : m_alpha0(alpha0)
      , m_rho(rho)
      , m_c(c)
    {}

    // g_func is a function g(alpha) = f(x + alpha p), where f is the
    // objective function and p is the search direction
    // grad_0 = dg(alpha)/dalpha evaluated at alpha = 0,
    // which is equal to grad f dot pk, evaluated at alpha = 0
    template <typename Tfunc>
    double search(Tfunc gFunc, double grad0)
    {
      double alpha = m_alpha0;
      double g0    = gFunc(0);
      for (int i = 0; i < m_itermax; ++i)
      {
        // std::cout << "\ni = " << i << std::endl;
        double valI = gFunc(alpha);
        // std::cout << "val_i = " << val_i << ", requirement = " << g0 + m_c*alpha*grad_0 << std::endl;
        // std::cout << "diff = " << val_i - (g0 + m_c*alpha*grad_0) << std::endl;
        if (valI <= g0 + m_c * alpha * grad0)
          return alpha;

        alpha *= m_rho;
      }

      throw std::runtime_error("could not find acceptable step length");
    }

  private:
    double m_alpha0;        // starting alpha value
    double m_rho;           // backtracking factor
    double m_c;             // gradient factor
    double m_itermax = 100; // maximum number of iterations, used for
                            // throwing error messages
};

} // namespace impl

} // namespace opt
} // namespace middle_mesh
} // namespace stk
#endif
