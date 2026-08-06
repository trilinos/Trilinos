#include "stk_search_util/ChangeOfBasis.hpp"
#include <cmath>

namespace stk::search::impl {

std::array<double, 9> compute_rotation_matrix(const std::array<double, 3>& new_z)
{
  double mag = std::sqrt(new_z[0]*new_z[0] + new_z[1]*new_z[1] + new_z[2]*new_z[2]);
  std::array<double, 3> new_z_unit = {new_z[0]/mag, new_z[1]/mag, new_z[2]/mag};
  double alpha = std::asin(-new_z_unit[1]);
  double beta = std::atan2(new_z_unit[0], new_z_unit[2]);

  double cosa = std::cos(alpha);
  double sina = -new_z_unit[1];
  std::array<double, 9> Rx = {1,    0,     0,
                              0, cosa, -sina,
                              0, sina,  cosa};

  double cosb = std::cos(beta);
  double sinb = std::sin(beta);
  std::array<double, 9> Ry = {cosb, 0, sinb,
                                0,  1,    0,
                             -sinb, 0, cosb};

  return matmat(Ry, Rx);
}


}