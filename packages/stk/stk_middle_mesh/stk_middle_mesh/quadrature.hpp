#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <cmath>
#include <stdexcept>
#include <vector>

namespace stk {
namespace middle_mesh {
namespace disc1d {
namespace impl {

class LegendrePolynomials
{
  public:
    // evaluate the nth Legendre polynomial at point x
    // n starts at 0 (constant)
    double compute(int n, double x);

    double compute_antideriv(int n, double x);
};

class Quadrature
{
  public:
    explicit Quadrature(int degree) { set_points_and_weights(degree); }

    const std::vector<double>& get_points() const { return m_points; }

    const std::vector<double>& get_weights() const { return m_weights; }

    int get_num_points() const { return m_points.size(); }

    double get_xi_min() const { return -1; }

    double get_xi_max() const { return 1; }

  private:
    // Note: Gauss-Legendre quadrature can exactly integrate 2n-1 degree
    // polynomials, where n is the number of points.  The error is
    // therefore 2n.  We refer to the degree of the quadrature
    // as 2n-1, giving the order as degree + 1 (this is consistent
    // with how finite element basis functions are described:
    // a degree 1 basis can exactly represent degree 1 polynomials,
    // and converges at a rate of 2)
    void set_points_and_weights(int degree);

    std::vector<double> m_points;
    std::vector<double> m_weights;
};

} // namespace impl

} // namespace disc1d
} // namespace middle_mesh
} // namespace stk
#endif