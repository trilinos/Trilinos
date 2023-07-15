#ifndef TAYLOR_PATCH_H
#define TAYLOR_PATCH_H

#include "matrix.hpp"
#include "surface_patch.hpp"

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

// TODO: DO NOT USE CUBIC OR HIGHER

// uses a two dimensional Taylor series approximation of z = f(x - x0, y - y0),
// where (x0, y0) is the point at which the approximation is centered
// The degree of the Taylor series is determined by the number of points
// supplied to constructPatch(): it is the minimam degree that will use
// all the points.
class TaylorPatch : public SurfacePatch
{
  public:
    using SurfacePatch::construct_patch;

    // constructs taylor series around pt0, using pts to solve for the function
    // values and derivatives
    void construct_patch(const std::vector<Point>& pts, const Point& pt0) override;

    Point eval_point(const double x, const double y) override;

    // evaluates jacobian of patch dx/du.  derivs must be a 3 x 2 matrix
    void eval_deriv(const double x, const double y, double derivs[2]) override;

  private:
    // returns the u,v coordinates of the point
    Point compute_uv(const Point& pt);

    // returns x,y coordinates of a point in u,v space
    Point compute_xy(const Point& ptUv);

    // computes the pair (degree, number of terms)
    std::pair<int, int> compute_degree(const int npts);

    // evaluates the Taylor polynomial coefficient for a given term
    // of a given degree (0 <= term <= degree)
    double eval_coefficient(const int degree, const int term, const Point& ptUv);

    // dervative of coefficient wrt u
    double eval_coefficient_du(const int degree, const int term, const Point& ptUv);

    // dervative of coefficient wrt v
    double eval_coefficient_dv(const int degree, const int term, const Point& ptUv);

    // computes the coefficient matrix used to solve for f and ints derivatives
    void compute_taylor_matrix(const std::vector<Point>& pts, utils::impl::Matrix<double>& a);

    void compute_taylor_derivs_rank_revealing_qr(const std::vector<Point>& pts);

    int compute_number_singular_values(const Matrix<double>& a, double tol);

    using IntF = int; // integer to use for factorial

    static IntF get_factorial(const unsigned int n);

    int m_degree   = 0;
    int m_numTerms = 0;
    Point m_pt0;
    std::vector<double> m_fvals;
    const static std::vector<IntF> M_FACTORIALS;
};

template <typename T = int>
T factorial(int n)
{
  assert(n >= 0);
  T val = 1;
  for (int i = 1; i <= n; ++i)
    val = val * i;

  return val;
}

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
#endif
