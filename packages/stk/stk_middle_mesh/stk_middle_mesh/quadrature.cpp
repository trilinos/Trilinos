#include "quadrature.hpp"

namespace stk {
namespace middle_mesh {
namespace disc1d {
namespace impl {

double LegendrePolynomials::compute(int n, double x)
{
  switch (n)
  {
    case 0:
      return 1;
    case 1:
      return x;
    case 2:
      return (3 * x * x - 1) / 2;
    case 3:
      return (5 * std::pow(x, 3) - 3 * x) / 2;
    case 4:
      return (25 * std::pow(x, 4) - 30 * x * x + 3) / 8;
    case 5:
      return (63 * std::pow(x, 5) - 70 * std::pow(x, 3) + 15 * x) / 8;
    default:
      throw std::runtime_error("unsupported degree");
  }
}

double LegendrePolynomials::compute_antideriv(int n, double x)
{
  switch (n)
  {
    case 0:
      return x;
    case 1:
      return x * x / 2;
    case 2:
      return (std::pow(x, 3) - x) / 2;
    case 3:
      return (5 * std::pow(x, 4) / 4 - 3 * x * x / 2) / 2;
    case 4:
      return (5 * std::pow(x, 5) - 10 * std::pow(x, 3) + 3 * x) / 8;
    case 5:
      return (63 * std::pow(x, 6) / 6 - 70 * std::pow(x, 4) / 4 + 15 * x * x / 2) / 8;
    default:
      throw std::runtime_error("unsupported degree");
  }
}

void Quadrature::set_points_and_weights(int degree)
{
  switch (degree)
  {
    case 1: {
      m_points  = {0};
      m_weights = {2};
      break;
    }

    case 3: {
      double rad35 = std::sqrt(3.0 / 5.0);
      m_points     = {-rad35, 0, rad35};
      m_weights    = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
      break;
    }

    case 5: {
      double t1 = 3.0 / 7.0;
      double t2 = 2 * std::sqrt(6.0 / 5.0) / 7;

      m_points = {-std::sqrt(t1 + t2), -std::sqrt(t1 - t2), std::sqrt(t1 - t2), std::sqrt(t1 + t2)};

      double w1 = (18 + std::sqrt(30)) / 36;
      double w2 = (18 - std::sqrt(30)) / 36;
      m_weights = {w2, w1, w1, w2};
      break;
    }

    default:
      throw std::runtime_error("unsupported quadrature degree");
  }
}

} // namespace impl
} // namespace disc1d
} // namespace middle_mesh
} // namespace stk
