#include "stk_unit_test_utils/MockMasterElementQuad4.hpp"

namespace stk::unit_test_util {

bool Quad4::within_tol(const double& value, const double& tolerance)
{
  return ( std::fabs(value) < tolerance );
}

template<typename T, std::size_t N>
T Quad4::vector_norm(std::array<T,N> &x)
{
  T norm_sq = 0.0;
  for (std::size_t i = 0; i < N; ++i ) {
    norm_sq += x[i] * x[i];
  }
  return norm_sq;
}

double Quad4::vector_norm(const double* vect, int N)
{
  double norm_sq = 0.0;
  for (int i = 0; i < N; ++i ) {
    norm_sq += vect[i] * vect[i];
  }
  return norm_sq;
}

double Quad4::linear_quad_parametric_distance(const std::array<double, 3>& x)
{
  const double ELEM_THICK = 0.01;
  std::array<double, 3> y = { { std::fabs(x[0]), std::fabs(x[1]), std::fabs(x[2]) } };
  double d = y[0];
  if(d < y[1]) d = y[1];
  if(ELEM_THICK < y[2] && d < 1 + y[2]) d = 1 + y[2];
  return d;
}

double Quad4::compute_linear_quad_cpp_distance(const double* normal, const double* solcur, const double* deltasol,
                                               double* par_coor)
{
  // Rescale the distance vector by the length of the (non-unit) normal vector,
  // which was used above in the NR iteration.
  const double area = std::sqrt(vector_norm(normal, 3));
  const double length = std::sqrt(area);
  const double par_coor_2 = (solcur[2] + deltasol[2]) * length;

  std::array<double, 3> xtmp = { { par_coor[0], par_coor[1], par_coor_2 } };
  return linear_quad_parametric_distance(xtmp);
}

void Quad4::non_unit_face_normal(const double* par_coord,       // (2)
                                 const double* elem_nodal_coor, // (4,3)
                                 double* normal_vector)         // (3)
{
  double xi  = par_coord[0];
  double eta = par_coord[1];

  // Translate element so that node 0 is at (x,y,z) = (0,0,0)
  double x[3] = { elem_nodal_coor[1] - elem_nodal_coor[0],
                  elem_nodal_coor[2] - elem_nodal_coor[0],
                  elem_nodal_coor[3] - elem_nodal_coor[0] };

  double y[3] = { elem_nodal_coor[5] - elem_nodal_coor[4],
                  elem_nodal_coor[6] - elem_nodal_coor[4],
                  elem_nodal_coor[7] - elem_nodal_coor[4] };

  double z[3] = { elem_nodal_coor[9] - elem_nodal_coor[8],
                  elem_nodal_coor[10] - elem_nodal_coor[8],
                  elem_nodal_coor[11] - elem_nodal_coor[8] };

  // Mathematica-generated and simplified code for the normal vector

  double n0 = 0.125* (xi * y[2] * z[0] + y[0] * z[1] + xi * y[0] * z[1] - y[2] * z[1] - xi * y[0] * z[2] +
                      y[1] * (-((1.0 + xi) * z[0]) + (1.0 + eta) * z[2]) +
                      eta * (y[2] * z[0] - y[2] * z[1] - y[0] * z[2]));

  double n1 = 0.125* (-(xi * x[2] * z[0]) - x[0] * z[1] - xi * x[0] * z[1] + x[2] * z[1] + xi * x[0] * z[2] +
                      x[1] * ((1.0 + xi) * z[0] - (1.0 + eta) * z[2]) +
                      eta * (-(x[2] * z[0]) + x[2] * z[1] + x[0] * z[2]));

  double n2 = 0.125* (xi * x[2] * y[0] + x[0] * y[1] + xi * x[0] * y[1] - x[2] * y[1] - xi * x[0] * y[2] +
                      x[1] * (-((1.0 + xi) * y[0]) + (1.0 + eta) * y[2]) +
                      eta * (x[2] * y[0] - x[2] * y[1] - x[0] * y[2]));

  normal_vector[0] = n0;
  normal_vector[1] = n1;
  normal_vector[2] = n2;
}

void Quad4::unit_face_normal(const double* par_coord,
                             const double* elem_nodal_coor, // (4,3)
                             double* normal_vector)         // (3)
{
  non_unit_face_normal(par_coord, elem_nodal_coor, normal_vector);

  double nlen = std::sqrt(vector_norm(normal_vector, 3));

  normal_vector[0] /= nlen;
  normal_vector[1] /= nlen;
  normal_vector[2] /= nlen;
}

void Quad4::interpolate_point(const double* par_coord, // (2)
                              const int& ncomp_field,
                              const double* field,  // (4,ncomp_field)
                              double* result)       // (ncomp_field)
{
  // 'field' is a flat array of dimension (4,ncomp_field) (Fortran ordering);

  double xi = par_coord[0];
  double eta = par_coord[1];

  for(int i = 0; i < ncomp_field; i++) {
    // Base 'field array' index for ith component
    int b = 4 * i;

    result[i] = 0.25 * ((1.0 - eta) * (1.0 - xi) * field[b + 0] +
                        (1.0 - eta) * (1.0 + xi) * field[b + 1] +
                        (1.0 + eta) * (1.0 + xi) * field[b + 2] +
                        (1.0 + eta) * (1.0 - xi) * field[b + 3]);
  }
}

double Quad4::is_in_element(const double* elem_nodal_coor, // (4,3)
                            const double* point_coor,      // (3)
                            double* par_coor)
{
  const double isInElemConverged = 1.0e-16;
  // Translate element so that (x,y,z) coordinates of the first node are (0,0,0)

  double x[3] = { elem_nodal_coor[1] - elem_nodal_coor[0],
                  elem_nodal_coor[2] - elem_nodal_coor[0],
                  elem_nodal_coor[3] - elem_nodal_coor[0] };

  double y[3] = { elem_nodal_coor[5] - elem_nodal_coor[4],
                  elem_nodal_coor[6] - elem_nodal_coor[4],
                  elem_nodal_coor[7] - elem_nodal_coor[4] };

  double z[3] = { elem_nodal_coor[9] - elem_nodal_coor[8],
                  elem_nodal_coor[10] - elem_nodal_coor[8],
                  elem_nodal_coor[11] - elem_nodal_coor[8] };

  // (xp,yp,zp) is the point at which we're searching for (xi,eta,d)
  // (must translate this also)
  // d = (scaled) distance in (x,y,z) space from point (xp,yp,zp) to the
  //     surface defined by the face element (the distance is scaled by
  //     the length of the non-unit normal vector; rescaling of d is done
  //     following the NR iteration below).

  double xp = point_coor[0] - elem_nodal_coor[0];
  double yp = point_coor[1] - elem_nodal_coor[4];
  double zp = point_coor[2] - elem_nodal_coor[8];

  // Newton-Raphson iteration for (xi,eta,d)

  double j[9];
  double gn[3];
  double xcur[3];   // current (x,y,z) point on element surface
  double normal[3]; // (non-unit) normal computed at xcur

  // Solution vector solcur[3] = {xi,eta,d}
  double solcur[3] = { -0.5, -0.5, -0.5 }; // initial guess
  std::array<double,3> deltasol = {{ 1.0, 1.0, 1.0 }};

  unsigned i = 0;
  const unsigned MAX_NR_ITER = 100;
  do {
    // Update guess vector
    solcur[0] += deltasol[0];
    solcur[1] += deltasol[1];
    solcur[2] += deltasol[2];

    interpolate_point(solcur, 3, elem_nodal_coor, xcur);

    // Translate xcur ((x,y,z) point corresponding
    // to current (xi,eta) guess)

    xcur[0] -= elem_nodal_coor[0];
    xcur[1] -= elem_nodal_coor[4];
    xcur[2] -= elem_nodal_coor[8];

    non_unit_face_normal(solcur, elem_nodal_coor, normal);

    gn[0] = xcur[0] - xp + solcur[2] * normal[0];
    gn[1] = xcur[1] - yp + solcur[2] * normal[1];
    gn[2] = xcur[2] - zp + solcur[2] * normal[2];

    // Mathematica-generated code for the jacobian

    j[0] = 0.125* (-2.0* (-1.0 + solcur[1]) * x[0] + (2.0* (1.0 +
            solcur[1]) * (x[1] - x[2]) + solcur[2] * (-(y[1] * z[0]) + y[2] * z[0] + y[0] * z[1] - y[0] * z[2])));

    j[1] = 0.125* (-2.0* (1.0+ solcur[0]) * x[0] + 2.0* (1.0 + solcur[0]) * x[1] - 2.0* (-1.0 + solcur[0]) *
                                    x[2] + (solcur[2] * (y[2] * (z[0] - z[1]) + (-y[0] + y[1]) * z[2])));

    j[2] = normal[0];

    j[3] = 0.125* (-2.0* (-1.0 + solcur[1]) * y[0] + (2.0* (1.0 +
            solcur[1]) * (y[1] - y[2]) + solcur[2] * (x[1] * z[0] - x[2] * z[0] - x[0] * z[1] + x[0] * z[2])));

    j[4] = 0.125* (-2.0* (1.0+ solcur[0]) * y[0] +
                                2.0* (1.0 + solcur[0]) * y[1] -
                                2.0* (-1.0 + solcur[0]) * y[2] +
                                (solcur[2] * (x[2] * (-z[0] + z[1]) + (x[0] - x[1]) * z[2])));

    j[5] = normal[1];

    j[6] = 0.125* ((solcur[2] * (-(x[1] * y[0]) + x[2] * y[0] + x[0] * y[1] - x[0] * y[2])) - 2.0* ((-1.0 + solcur[1]) * z[0] - (1.0 + solcur[1]) * (z[1] - z[2])));

    j[7] = 0.125* ((solcur[2] * (x[2] * (y[0] - y[1]) + (-x[0] + x[1]) * y[2])) - 2.0* (1.0 + solcur[0]) * z[0] + 2.0* (1.0 + solcur[0]) *
                                    z[1] - 2.0* (-1.0 + solcur[0]) * z[2]);

    j[8] = normal[2];

    double jdet = -(j[2] * j[4] * j[6]) + j[1] * j[5] * j[6] + j[2] * j[3] * j[7] - j[0] * j[5] * j[7] - j[1] * j[3] * j[8] +
            j[0] * j[4] * j[8];

    // Solve linear system (j*deltasol = -gn) for deltasol at step n+1

    deltasol[0] = (gn[2] * (j[2] * j[4] - j[1] * j[5]) + gn[1] * (-(j[2] * j[7]) + j[1] * j[8]) +
                    gn[0] * (j[5] * j[7] - j[4] * j[8])) /
                  jdet;
    deltasol[1] = (gn[2] * (-(j[2] * j[3]) + j[0] * j[5]) + gn[1] * (j[2] * j[6] - j[0] * j[8]) +
                    gn[0] * (-(j[5] * j[6]) + j[3] * j[8])) /
                  jdet;
    deltasol[2] = (gn[2] * (j[1] * j[3] - j[0] * j[4]) + gn[1] * (-(j[1] * j[6]) + j[0] * j[7]) +
                    gn[0] * (j[4] * j[6] - j[3] * j[7])) /
                  jdet;

  } while(!within_tol(vector_norm(deltasol), isInElemConverged) && ++i < MAX_NR_ITER);

  // Fill in solution vector; only include the distance (in the third
  // solution vector slot) if npar_coord = 3 (this is how the user
  // requests it)

  par_coor[0] = par_coor[1] = std::numeric_limits<double>::max();

  if(i < MAX_NR_ITER) {
    par_coor[0] = solcur[0] + deltasol[0];
    par_coor[1] = solcur[1] + deltasol[1];
  }

  double dist = compute_linear_quad_cpp_distance(normal, solcur, deltasol.data(), par_coor);

  return dist;
}

const std::vector<double>& Quad4::coordinate_center()
{
  static const std::vector<double> C(2, 0.);
  return C;
}
}