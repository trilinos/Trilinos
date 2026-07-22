#include "stk_unit_test_utils/MockMasterElementQuad4_2D.hpp"

namespace stk::unit_test_util {


bool Quad4_2D::within_tol(const double& value, const double& tolerance)
{
  return ( std::fabs(value) < tolerance );
}

template<typename T, std::size_t N>
T Quad4_2D::vector_norm(std::array<T,N> &x)
{
  T norm_sq = 0.0;
  for (std::size_t i = 0; i < N; ++i ) {
    norm_sq += x[i] * x[i];
  }
  return norm_sq;
}

double Quad4_2D::vector_norm(const double* vect, int N)
{
  double norm_sq = 0.0;
  for (int i = 0; i < N; ++i ) {
    norm_sq += vect[i] * vect[i];
  }
  return norm_sq;
}

double Quad4_2D::linear_quad_parametric_distance(const std::array<double, 2>& x)
{
  std::array<double,2> y = {{ std::fabs(x[0]), std::fabs(x[1]) }};
  const double d = *std::max_element(y.begin(), y.end());
  return d;
}

void Quad4_2D::non_unit_face_normal(const double* par_coord,       // (2)
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

void Quad4_2D::unit_face_normal(const double* par_coord,
                                const double* elem_nodal_coor, // (4,3)
                                double* normal_vector)         // (3)
{
  non_unit_face_normal(par_coord, elem_nodal_coor, normal_vector);

  double nlen = std::sqrt(vector_norm(normal_vector, 3));

  normal_vector[0] /= nlen;
  normal_vector[1] /= nlen;
  normal_vector[2] /= nlen;
}

void Quad4_2D::interpolate_point(const double* par_coord, // (2)
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

double Quad4_2D::is_in_element(const double* elem_nodal_coor, // (4,3)
                               const double* point_coor,      // (3)
                               double* par_coor)
{
const double isInElemConverged = 1.0e-16;

// Translate element so that (x,y) coordinates of the first node are (0,0)
double x[] = { 0., elem_nodal_coor[1] - elem_nodal_coor[0], elem_nodal_coor[2] - elem_nodal_coor[0], elem_nodal_coor[3] - elem_nodal_coor[0] };
double y[] = { 0., elem_nodal_coor[5] - elem_nodal_coor[4], elem_nodal_coor[6] - elem_nodal_coor[4], elem_nodal_coor[7] - elem_nodal_coor[4] };

// (xp,yp) is the point at which we're searching for (xi,eta)
// (must translate this also)

double xp = point_coor[0] - elem_nodal_coor[0];
double yp = point_coor[1] - elem_nodal_coor[4];

// Newton-Raphson iteration for (xi,eta)
double j[4];
double f[2];
double shapefct[4];

double xinew = 0.5; // initial guess
double etanew = 0.5;

double xicur = 0.5;
double etacur = 0.5;

std::array<double,2> xidiff = { 1.0, 1.0 };
unsigned i = 0;

bool converged = false;
const unsigned MAX_NR_ITER = 10; // bilinear so should converge in 1 step

do {
  xicur = xinew;
  etacur = etanew;

  j[0] =  0.25 * (1.0 - etacur) * x[1] + 0.25 * (1.0 + etacur) * x[2] - 0.25 * (1.0 + etacur) * x[3];
  j[1] = -0.25 * (1.0 +  xicur) * x[1] + 0.25 * (1.0 +  xicur) * x[2] + 0.25 * (1.0 -  xicur) * x[3];
  j[2] =  0.25 * (1.0 - etacur) * y[1] + 0.25 * (1.0 + etacur) * y[2] - 0.25 * (1.0 + etacur) * y[3];
  j[3] = -0.25 * (1.0 +  xicur) * y[1] + 0.25 * (1.0 +  xicur) * y[2] + 0.25 * (1.0 -  xicur) * y[3];

  double jdet = j[0] * j[3] - j[1] * j[2];

  shapefct[0] = 0.25 * (1.0 - etacur) * (1.0 - xicur);
  shapefct[1] = 0.25 * (1.0 - etacur) * (1.0 + xicur);
  shapefct[2] = 0.25 * (1.0 + etacur) * (1.0 + xicur);
  shapefct[3] = 0.25 * (1.0 + etacur) * (1.0 - xicur);

  f[0] = (shapefct[1] * x[1] + shapefct[2] * x[2] + shapefct[3] * x[3]) - xp;
  f[1] = (shapefct[1] * y[1] + shapefct[2] * y[2] + shapefct[3] * y[3]) - yp;

  xinew = xicur - (f[0] * j[3] - f[1] * j[1]) / jdet;
  etanew = etacur - (-f[0] * j[2] + f[1] * j[0]) / jdet;

  xidiff[0] = xinew - xicur;
  xidiff[1] = etanew - etacur;

  converged = within_tol(vector_norm(xidiff), isInElemConverged);
} while(!converged && ++i < MAX_NR_ITER);

par_coor[0] = par_coor[1] = std::numeric_limits<double>::max();
double dist = std::numeric_limits<double>::max();

if(i < MAX_NR_ITER) {
  std::array<double,2> xtmp = {{ par_coor[0] = xinew, par_coor[1] = etanew }};
  dist = linear_quad_parametric_distance(xtmp);
}
return dist;
}

const std::vector<double>& Quad4_2D::coordinate_center()
{
  static const std::vector<double> C(2, 0.);
  return C;
}

}