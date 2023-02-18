#ifndef UTILS_B_SPLINES_H
#define UTILS_B_SPLINES_H

#include <array>
#include <cmath>

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

// evaluates cubic B Splines with a uniform knot vector with parameter
// t ranging from 0 to 1
class CubicBSplines
{
  public:
    void eval(double t, std::array<double, 4>& vals)
    {
      // TODO: use Horners rule
      double t3 = std::pow(t, 3);
      vals[0]   = std::pow(1 - t, 3) / 6;
      vals[1]   = (3 * t3 - 6 * t * t + 4) / 6;
      vals[2]   = (-3 * t3 + 3 * t * t + 3 * t + 1) / 6;
      vals[3]   = t3 / 6;
    }
};

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
#endif