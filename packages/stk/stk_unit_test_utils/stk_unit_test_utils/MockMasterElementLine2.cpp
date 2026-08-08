#include "stk_unit_test_utils/MockMasterElementLine2.hpp"

namespace stk::unit_test_util {

double Line2::parametric_distance(const std::array<double, 2>& x)
{
  const double ELEMENT_THICKNESS = 0.01;
  double dist = std::fabs(x[0]);
  if(ELEMENT_THICKNESS < x[1] && dist < 1 + x[1]) dist = 1 + x[1];
  return dist;
}

double Line2::is_in_element(const double* elem_nodal_coor, // (3,3)
                            const double* point_coor,      // (3)
                            double* par_coor)
{
  // elem_nodal_coor has the endpoints of the line
  // segment defining this element.  Set the first
  // endpoint to zero.  This means subtrace the
  // first endpoint from the second.
  const double X1 = elem_nodal_coor[1] - elem_nodal_coor[0];
  const double X2 = elem_nodal_coor[3] - elem_nodal_coor[2];
  const double X3 = elem_nodal_coor[5] - elem_nodal_coor[4];

  // Now subtract the first endpoint from the target point
  const double P1 = point_coor[0] - elem_nodal_coor[0];
  const double P2 = point_coor[1] - elem_nodal_coor[2];
  const double P3 = point_coor[2] - elem_nodal_coor[4];

  // Now find the projection along the line of the point
  // This is the parametric coordinate in range (0,1)
  const double norm2 = X1 * X1 + X2 * X2 + X3 * X3;
  STK_ThrowAssert(norm2 != 0.0);

  const double xi = (P1 * X1 + P2 * X2 + P3 * X3) / norm2;
  // rescale to (-1,1)
  par_coor[0] = 2.0 * xi - 1.0;

  // d = ||P x X||/||X||
  // gives the normalized distance from the point to the element.

  const double alpha = std::sqrt((P2 * X1 - P1 * X2) * (P2 * X1 - P1 * X2) + (-P3 * X1 + P1 * X3) * (-P3 * X1 + P1 * X3) +
                                  (P3 * X2 - P2 * X3) * (P3 * X2 - P2 * X3)) / norm2;

  std::array x = { par_coor[0], alpha };
  const double dist = parametric_distance(x);

  return dist;
}

void Line2::interpolate_point(const double* par_coord, // (1)
                              const int& ncomp_field,
                              const double* field,     // (2,ncomp_field)
                              double* result)          // (ncomp_field)
{
  // 'field' is a flat array of dimension (2,ncomp_field) (Fortran ordering);
  double one = 1.0;

  double s = par_coord[0];

  for(int i = 0; i < ncomp_field; i++) {
    int b = 2 * i; // Base 'field array' index for ith component

    result[i] = 0.5 * (one - s) * field[b + 0] + 0.5 * (one + s) * field[b + 1];
  }
}

const std::vector<double>& Line2::coordinate_center()
{
  static const std::vector<double> C(1, 0.);
  return C;
}
}