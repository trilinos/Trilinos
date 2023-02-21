#include "change_of_basis.h"
#include "mat2x2.h"
#include "newton.h"
#include "projection.h"
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace stk::middle_mesh;

struct Data
{
    std::array<utils::Point, 3> verts;
    utils::Point pt;
    std::array<utils::Point, 3> normals;
};

Data rotate_data(const Data& data)
{
  std::array<utils::Point, 3> basis;
  int basisType = 0;
  if (basisType >= 0 && basisType <= 2)
  {
    utils::Point dx1 = data.verts[1] - data.verts[0];
    utils::Point dx2 = data.verts[2] - data.verts[0];
    basis[0]         = dx1 / std::sqrt(dot(dx1, dx1));
    basis[1]         = dx2 - dot(dx2, basis[0]) * basis[0];
    basis[2]         = cross(basis[0], basis[1]);

    std::array<utils::Point, 3> basisPerm;
    for (int i = 0; i < 3; ++i)
      basisPerm[i] = basis[(i + basisType) % 3];

    for (int i = 0; i < 3; ++i)
      basis[i] = basisPerm[i];
  } else if (basisType == 3)
  {
    utils::Point dx1    = data.verts[1] - data.verts[0];
    utils::Point normal = (data.normals[0] + data.normals[1] + data.normals[1]) / 3;
    normal              = normal / std::sqrt(dot(normal, normal));
    basis[0]            = normal;
    basis[1]            = dx1 - dot(dx1, basis[0]) * basis[0];
    basis[2]            = cross(basis[0], basis[1]);
  } else
    throw std::runtime_error("unrecognized basis_type");

  std::cout << "basis = " << basis[0] << ", " << basis[1] << ", " << basis[2] << std::endl;

  // the matrix for computing the lambdas can be singular
  // if all the points are in the xz plane.  Do a change
  // of basis so it is full rank.
  utils::impl::ChangeOfBasis q(basis);
  Data data2;
  for (int i = 0; i < 3; ++i)
  {
    data2.verts[i]   = q.project_forward(data.verts[i]);
    data2.normals[i] = q.project_forward(data.normals[i]);
  }
  data2.pt = q.project_forward(data.pt);

  return data2;
}

void newton3x3_func(const Data& data, const std::vector<double>& x, std::vector<double>& f);

void print_inputs(const Data& data)
{
  std::cout << "verts = " << data.verts[0] << ", " << data.verts[1] << ", " << data.verts[2] << std::endl;
  std::cout << "pt = " << data.pt << std::endl;
  std::cout << "normals = " << data.normals[0] << ", " << data.normals[1] << ", " << data.normals[2] << std::endl;
}

void print_inputs_for_visualization(const Data& data)
{
  std::ofstream file;
  file.open("tmp_verts.txt");
  file << std::setprecision(16);

  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
      file << data.verts[i][j] << " ";
    file << std::endl;
  }

  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
      file << (data.verts[i][j] + data.normals[i][j]) << " ";
    file << std::endl;
  }

  for (int j = 0; j < 3; ++j)
    file << data.pt[j] << " ";

  file.close();

  file.open("tmp_edges.txt");
  file << "0 1\n1 2\n0 2\n";
  file << "0 3\n1 4\n2 5\n";
  file.close();
}

void check_solution_reverse(const Data& data, const std::array<double, 3>& solution)
{
  double lambda1 = solution[0];
  double lambda2 = solution[1];
  double lambda0 = 1 - lambda1 - lambda2;
  double gamma   = solution[2];

  utils::Point p1 =
      data.pt + gamma * (lambda0 * data.normals[0] + lambda1 * data.normals[1] + lambda2 * data.normals[2]);
  utils::Point p2   = lambda0 * data.verts[0] + lambda1 * data.verts[1] + lambda2 * data.verts[2];
  utils::Point disp = p2 - p1;
  double dist       = std::sqrt(dot(disp, disp));
  std::cout << "dist = " << dist << std::endl;
  // if (dist > 1e-13)
  //   throw std::runtime_error("Newton produced incorrect solution");
}

std::vector<double> get_initial_guess(const Data& data)
{
  // Idea: use the normal vector at the centroid of the element, and do the linear
  //       projection.  Hopefully this gives an initial guess that is in the
  //       neighborhood of the solution

  auto& verts   = data.verts;
  auto& pt      = data.pt;
  auto& normals = data.normals;

  utils::Point rhs = pt - verts[0];

  utils::Point dx1 = verts[1] - verts[0];
  utils::Point dx2 = verts[2] - verts[0];
  utils::impl::Matrix<double> mat(3, 3);
  utils::Point normal = (normals[0] + normals[1] + normals[2]) / 3;
  mat(0, 0)           = dx1.x;
  mat(0, 1)           = dx2.x;
  mat(0, 2)           = -normal.x;
  mat(1, 0)           = dx1.y;
  mat(1, 1)           = dx2.y;
  mat(1, 2)           = -normal.y;
  mat(2, 0)           = dx1.z;
  mat(2, 1)           = dx2.z;
  mat(2, 2)           = -normal.z;

  std::array<int, 3> ipiv;
  solve_linear_system(mat, ipiv.data(), &(rhs.x));

  // std::cout << "initial guess for gamma = " << rhs.z << std::endl;

  return {rhs.x, rhs.y, rhs.z};
}

std::array<double, 2> solve_for_lambda(const Data& data, double gamma)
{
  auto& verts   = data.verts;
  auto& pt      = data.pt;
  auto& normals = data.normals;

  utils::Point dx1 = verts[1] - verts[0];
  utils::Point dx2 = verts[2] - verts[0];
  utils::Point dn1 = normals[1] - normals[0];
  utils::Point dn2 = normals[2] - normals[0];

  utils::impl::Mat2x2<double> mat;
  mat(0, 0) = dx1[0] - gamma * dn1[0];
  mat(0, 1) = dx2[0] - gamma * dn2[0];
  mat(1, 0) = dx1[1] - gamma * dn1[1];
  mat(1, 1) = dx2[1] - gamma * dn2[1];

  // std::cout << "mat = \n" << mat << std::endl;

  std::array<double, 2> rhs = {pt.x - verts[0].x + gamma * normals[0].x, pt.y - verts[0].y + gamma * normals[0].y};
  // std::cout << "rhs = " << rhs[0] << ", " << rhs[1] << std::endl;

  // std::cout << "dx1 = " << dx1 << ", dx2 = " << dx2 << std::endl;

  std::array<double, 2> lambdas;

  matsolve2x2(mat, lambdas.data(), rhs.data());

  // utils::Point normal = (1 - lambdas[0] - lambdas[1]) * data.normals[0] + lambdas[1] * data.normals[1] + lambdas[2] *
  // data.normals[2]; std::cout << "gamma = " << gamma << ", lambdas = " << lambdas[0] << ", " << lambdas[1] << ",
  // intersection pt = " << (pt + gamma * normal) << std::endl;

  // std::array<double, 2> rhs2;
  // matvec2x2(mat, lambdas.data(), rhs2.data());
  // std::cout << "solve error = " << (rhs[0] - rhs2[0]) << ", " << (rhs[1] - rhs2[1]) << std::endl;

  return lambdas;
}

std::array<double, 2> solve_for_lambdad_gamma(const Data& data, double gamma)
{
  auto& verts   = data.verts;
  auto& pt      = data.pt;
  auto& normals = data.normals;

  utils::Point dx1 = verts[1] - verts[0];
  utils::Point dx2 = verts[2] - verts[0];
  utils::Point dn1 = normals[1] - normals[0];
  utils::Point dn2 = normals[2] - normals[0];

  utils::impl::Mat2x2<double> mat, matDot;
  mat(0, 0) = dx1[0] - gamma * dn1[0];
  mat(0, 1) = dx2[0] - gamma * dn2[0];
  mat(1, 0) = dx1[1] - gamma * dn1[1];
  mat(1, 1) = dx2[1] - gamma * dn2[1];

  matDot(0, 0) = -dn1[0];
  matDot(0, 1) = -dn2[0];
  matDot(1, 0) = -dn1[1];
  matDot(1, 1) = -dn2[1];

  std::array<double, 2> rhs = {pt.x - verts[0].x + gamma * normals[0].x, pt.y - verts[0].y + gamma * normals[0].y};

  std::array<double, 2> rhsDot = {normals[0].x, normals[0].y};

  std::array<double, 2> lambdas, lambdasDot;

  matsolve2x2_dot(mat, matDot, lambdas.data(), lambdasDot.data(), rhs.data(), rhsDot.data());

  return lambdasDot;
}

double compute_fof_gamma(const Data& data, double gamma)
{
  std::array<double, 2> lambdas = solve_for_lambda(data, gamma);
  double lambda1                = lambdas[0];
  double lambda2                = lambdas[1];

  auto& verts   = data.verts;
  auto& pt      = data.pt;
  auto& normals = data.normals;

  utils::Point dx1 = verts[1] - verts[0];
  utils::Point dx2 = verts[2] - verts[0];
  utils::Point dn1 = normals[1] - normals[0];
  utils::Point dn2 = normals[2] - normals[0];

  // std::cout << "x0 = " << data.verts[0]<< std::endl;
  // std::cout << "dx1 = " << dx1 << std::endl;
  // std::cout << "dx2 = " << dx2 << std::endl;
  // std::cout << "d0  = " << data.normals[0] << std::endl;
  // std::cout << "dn1 = " << dn1 << std::endl;
  // std::cout << "dn2 = " << dn2 << std::endl;
  // std::cout << "pt = " <<  pt << std::endl;

  // double term1 = -gamma*(normals[0].z + dn1.z*lambda1 + dn2.z*lambda2);
  // double term2 = - pt.z + verts[0].z + dx1.z*lambda1 + dx2.z*lambda2;
  // std::cout << "gamma = " << gamma << ", term1 = " << term1 << ", term2 = " << term2 << std::endl;
  // std::cout << "  lambda1 = " << lambda1 << ", lambda2 = " << lambda2 << ", normal[0].z = " << normals[0].z << ",
  // dn1.z = " << dn1.z << ", dn2.z = " << dn2.z << std::endl;
  return -gamma * (normals[0].z + dn1.z * lambda1 + dn2.z * lambda2) - pt.z + verts[0].z + dx1.z * lambda1 +
         dx2.z * lambda2;
}

double computed_fd_gamma(const Data& data, double gamma)
{
  std::array<double, 2> lambdas = solve_for_lambda(data, gamma);
  double lambda1                = lambdas[0];
  double lambda2                = lambdas[1];

  std::array<double, 2> dlambdasDgamma = solve_for_lambdad_gamma(data, gamma);
  double dlambda1Dgamma                = dlambdasDgamma[0];
  double dlambda2Dgamma                = dlambdasDgamma[1];

  auto& verts   = data.verts;
  auto& normals = data.normals;

  utils::Point dx1 = verts[1] - verts[0];
  utils::Point dx2 = verts[2] - verts[0];
  utils::Point dn1 = normals[1] - normals[0];
  utils::Point dn2 = normals[2] - normals[0];

  double partialfPartialgamma = -(normals[0].z + dn1.z * lambda1 + dn2.z * lambda2);
  double dfDlambda1           = -gamma * dn1.z + dx1.z;
  double dfDlambda2           = -gamma * dn2.z + dx2.z;

  return partialfPartialgamma + dfDlambda1 * dlambda1Dgamma + dfDlambda2 * dlambda2Dgamma;
}

void output_f_values(double gammaMin, double gammaMax, int npts, const Data& data)
{
  assert(gammaMax > gammaMin);
  assert(npts >= 2);

  double deltaGamma = (gammaMax - gammaMin) / npts;

  std::ofstream of;
  of.open("f_vals.txt");
  of << std::setprecision(16);
  for (int i = 0; i < npts; ++i)
  {
    double gamma = gammaMin + i * deltaGamma;
    double val   = compute_fof_gamma(data, gamma);
    of << gamma << " " << val << std::endl;
    // std::cout << "gamma = " << gamma << ", f(gamma) = " << val << std::endl;
  }

  of.close();
}

// returns a pair (gamma, f(gamma)) where f(gamma) closest to zero
std::pair<double, double> solve_gamma_brute_force(double gammaMin, double gammaMax, int npts, const Data& data)
{
  std::cout << "Finding gamma by brute force" << std::endl;
  assert(gammaMax > gammaMin);
  assert(npts >= 2);

  double deltaGamma = (gammaMax - gammaMin) / npts;

  double fClosest     = std::numeric_limits<double>::max();
  double gammaClosest = 0;
  for (int i = 0; i < npts; ++i)
  {
    double gamma = gammaMin + i * deltaGamma;
    double val   = compute_fof_gamma(data, gamma);
    if (std::abs(val) < std::abs(fClosest))
    {
      fClosest     = val;
      gammaClosest = gamma;
    }
  }

  std::cout << "gamma = " << gammaClosest << " gives f = " << fClosest << std::endl;

  std::array<double, 2> lambdas = solve_for_lambda(data, gammaClosest);
  std::cout << "lambdas = " << utils::Point(1 - lambdas[0] - lambdas[1], lambdas[0], lambdas[1]) << std::endl;
  std::vector<double> newtonRhs(3);
  newton3x3_func(data, {lambdas[0], lambdas[1], gammaClosest}, newtonRhs);
  std::cout << "newton3x3 rhs = " << newtonRhs[0] << ", " << newtonRhs[1] << ", " << newtonRhs[2] << std::endl;

  double rhsNorm = 0;
  for (int i = 0; i < 3; ++i)
    rhsNorm += newtonRhs[i] * newtonRhs[i];
  rhsNorm = std::sqrt(rhsNorm);
  std::cout << "newton 3x3 rhs_norm = " << rhsNorm << std::endl;
  check_solution_reverse(data, {lambdas[0], lambdas[1], gammaClosest});

  return std::make_pair(gammaClosest, fClosest);
}

void solve_newton_gamma(const Data& data)
{
  auto func = [&](const std::vector<double>& x, std::vector<double>& f) {
    // std::cout << "gamma = " << x[0] << std::endl;
    f[0] = compute_fof_gamma(data, x[0]);
  };

  auto jac = [&](const std::vector<double>& x, utils::impl::Matrix<double>& mat) {
    mat(0, 0) = computed_fd_gamma(data, x[0]);
  };

  std::vector<double> x0 = {1e6 /*0*/};
  utils::impl::Newton newton(1, 1e-13, 100000);
  int retVal = newton.solve(func, jac, x0);

  std::array<double, 2> lambdasTwo = solve_for_lambda(data, x0[0]);

  std::array<double, 3> lambdas = {1 - lambdasTwo[0] - lambdasTwo[1], lambdasTwo[0], lambdasTwo[1]};
  utils::Point projectedPt(0, 0, 0);
  for (int i = 0; i < 3; ++i)
    projectedPt += lambdas[i] * data.verts[i];

  std::cout << "NewtonGamma found lambdas " << lambdas[0] << ", " << lambdas[1] << ", " << lambdas[2]
            << ", gamma = " << x0[0] << std::endl;
  std::cout << "found point " << projectedPt << std::endl;
  std::cout << "f(gamma) = " << compute_fof_gamma(data, x0[0]) << std::endl;

  std::vector<double> newtonRhs(3);
  newton3x3_func(data, {lambdas[1], lambdas[2], x0[0]}, newtonRhs);
  std::cout << "newton3x3 rhs = " << newtonRhs[0] << ", " << newtonRhs[1] << ", " << newtonRhs[2] << std::endl;

  double rhsNorm = 0;
  for (int i = 0; i < 3; ++i)
    rhsNorm += newtonRhs[i] * newtonRhs[i];
  rhsNorm = std::sqrt(rhsNorm);
  std::cout << "newton 3x3 rhs_norm = " << rhsNorm << std::endl;
  check_solution_reverse(data, {lambdas[1], lambdas[2], x0[0]});

  if (retVal > 0)
    std::cout << "NEWTON SOLVER DID NOT CONVERGE" << std::endl;
  else
    std::cout << "newton solver converged" << std::endl;
}

void newton3x3_func(const Data& data, const std::vector<double>& x, std::vector<double>& f)
{
  auto& verts   = data.verts;
  auto& pt      = data.pt;
  auto& normals = data.normals;

  utils::Point rhs = -pt + verts[0];
  utils::Point dx1 = verts[1] - verts[0];
  utils::Point dx2 = verts[2] - verts[0];
  utils::Point dn1 = normals[1] - normals[0];
  utils::Point dn2 = normals[2] - normals[0];

  double lambda1 = x[0], lambda2 = x[1], gamma = x[2];
  for (int i = 0; i < 3; ++i)
    f[i] = dx1[i] * lambda1 + dx2[i] * lambda2 - gamma * (normals[0][i] + dn1[i] * lambda1 + dn2[i] * lambda2) + rhs[i];
}

void newton3x3_jac(const Data& data, const std::vector<double>& x, utils::impl::Matrix<double>& mat)
{
  auto& verts   = data.verts;
  auto& normals = data.normals;

  utils::Point dx1 = verts[1] - verts[0];
  utils::Point dx2 = verts[2] - verts[0];
  utils::Point dn1 = normals[1] - normals[0];
  utils::Point dn2 = normals[2] - normals[0];

  double lambda1 = x[0], lambda2 = x[1], gamma = x[2];
  for (int i = 0; i < 3; ++i)
  {
    mat(i, 0) = dx1[i] - gamma * dn1[i];
    mat(i, 1) = dx2[i] - gamma * dn2[i];
    mat(i, 2) = -(normals[i][i] + dn1[i] * lambda1 + dn2[i] * lambda2);
  }
}

void solve_newton3x3(const Data& data)
{
  auto func = [&](const std::vector<double>& x, std::vector<double>& f) { newton3x3_func(data, x, f); };

  auto jac = [&](const std::vector<double>& x, utils::impl::Matrix<double>& mat) { newton3x3_jac(data, x, mat); };

  std::vector<double> x0 = get_initial_guess(data);
  // std::vector<double> x0(3);
  // x0[0] = 0.5;
  // x0[1] = 0.5;
  // x0[2] = 0;  //TODO: maybe there is a better way to do this:
  //       compute vector from pt to centroid of
  //       element, project onto normal vector
  //       at the centroid, compute gamma

  utils::impl::Newton newton(3, 1e-13, 100000);
  int retVal = newton.solve(func, jac, x0);

  std::array<double, 3> lambdas = {1 - x0[0] - x0[1], x0[0], x0[1]};
  utils::Point projectedPt(0, 0, 0);
  for (int i = 0; i < 3; ++i)
    projectedPt += lambdas[i] * data.verts[i];

  std::cout << "found lambdas = " << utils::Point(lambdas[0], lambdas[1], lambdas[2]) << ", gamma = " << x0[2]
            << std::endl;

  std::cout << "found point " << projectedPt << std::endl;
  std::cout << "input point = " << data.pt << std::endl;
  for (int i = 0; i < 3; ++i)
  {
    std::cout << "vert = " << data.verts[i] << ", normal = " << data.normals[i] << std::endl;
    std::cout << "second point on line is " << data.verts[i] + data.normals[i] << std::endl;
  }

  if (retVal > 0)
    std::cout << "NEWTON SOLVER DID NOT CONVERGE" << std::endl;
  else
    std::cout << "newton solver converged" << std::endl;
}

void solve_gamma_fixed_point_iteration(const Data& data, double gamma)
{
  int itermax       = 100000000;
  double gammaTol   = 1e-13;
  double deltaGamma = std::numeric_limits<double>::max();

  utils::Point dx1 = data.verts[1] - data.verts[0];
  utils::Point dx2 = data.verts[2] - data.verts[0];
  utils::Point dn1 = data.normals[1] - data.normals[0];
  utils::Point dn2 = data.normals[2] - data.normals[0];
  std::array<double, 2> lambdas;
  int iter = 0;
  while (deltaGamma > gammaTol)
  {
    std::cout << "iter = " << iter << ", gamma = " << gamma << std::endl;
    lambdas    = solve_for_lambdad_gamma(data, gamma);
    double d   = (1 - lambdas[0] - lambdas[1]) * data.normals[0].z + lambdas[0] * dn1.z + lambdas[1] * dn2.z;
    double rhs = data.pt.z - data.verts[0].z;
    rhs -= dx1.z * lambdas[0] + dx2.z * lambdas[1];

    double gammaNew = -rhs / d;
    deltaGamma      = std::abs(gammaNew - gamma);
    gamma           = 0.5 * (gammaNew + gamma);

    if (iter > itermax)
      break;

    iter++;
  }

  if (deltaGamma > gammaTol)
    std::cout << "FIXED POINT ITERATION DID NOT CONVERGE" << std::endl;
  else
    std::cout << "fixed point ieration converged" << std::endl;

  std::vector<double> newtonRhs(3);
  newton3x3_func(data, {lambdas[0], lambdas[1], gamma}, newtonRhs);
  std::cout << "newton3x3 rhs = " << newtonRhs[0] << ", " << newtonRhs[1] << ", " << newtonRhs[2] << std::endl;

  double rhsNorm = 0;
  for (int i = 0; i < 3; ++i)
    rhsNorm += newtonRhs[i] * newtonRhs[i];
  rhsNorm = std::sqrt(rhsNorm);
  std::cout << "newton 3x3 rhs_norm = " << rhsNorm << std::endl;
}

std::array<utils::impl::Matrix<double>, 3> make_rotation_matrices(double alpha, double beta, double gamma)
{
  std::array<utils::impl::Matrix<double>, 3> mats = {
      utils::impl::Matrix<double>(3, 3), utils::impl::Matrix<double>(3, 3), utils::impl::Matrix<double>(3, 3)};

  auto& rx = mats[0];
  rx.fill(0);
  rx(0, 0) = 1;
  rx(1, 1) = std::cos(alpha);
  rx(1, 2) = std::sin(alpha);
  rx(2, 1) = -std::sin(alpha);
  rx(2, 2) = std::cos(alpha);

  auto& ry = mats[1];
  ry.fill(0);
  ry(0, 1) = 1;
  ry(1, 0) = std::cos(beta);
  ry(0, 2) = -std::sin(beta);
  ry(2, 2) = std::sin(beta);
  ry(2, 2) = std::cos(beta);

  auto& rz = mats[2];
  rz.fill(0);
  rz(2, 2) = 1;
  rz(0, 0) = std::cos(gamma);
  rz(0, 1) = std::sin(gamma);
  rz(1, 0) = -std::sin(gamma);
  rz(1, 1) = std::cos(gamma);

  // std::cout << "Rx = \n" << mats[0] << "\nRy = " << mats[1] << "\nRz = " << mats[2] << std::endl;

  return mats;
}

utils::Point apply_rotation(const std::array<utils::impl::Matrix<double>, 3>& mats, utils::Point ptIn)
{
  for (int j = 0; j < 3; ++j)
  {
    utils::Point ptTmp;
    matvec(1, mats[2 - j], &(ptIn[0]), 0, &(ptTmp[0]));
    ptIn = ptTmp;
  }

  return ptIn;
}

void check_singular(const Data& data, double gamma)
{
  double pi    = std::atan(1) * 4;
  double alpha = 45 * pi / 180;
  double beta  = 45 * pi / 180;
  double delta = 45 * pi / 180;

  Data data2;

  auto rotMats = make_rotation_matrices(alpha, beta, delta);

  for (int i = 0; i < 3; ++i)
  {
    data2.verts[i] = apply_rotation(rotMats, data.verts[i]);
    std::cout << "vert " << data.verts[i] << " rotated to " << data2.verts[i] << std::endl;
    data2.normals[i] = apply_rotation(rotMats, data.normals[i]);
  }
  data2.pt = apply_rotation(rotMats, data.pt);

  utils::Point dx1 = data2.verts[1] - data2.verts[0];
  utils::Point dx2 = data2.verts[2] - data2.verts[0];
  utils::Point dn1 = data2.normals[1] - data2.normals[0];
  utils::Point dn2 = data2.normals[2] - data2.normals[0];

  utils::impl::Mat2x2<double> mat;
  mat(0, 0) = dx1[0] - gamma * dn1[0];
  mat(0, 1) = dx2[0] - gamma * dn2[0];
  mat(1, 0) = dx1[1] - gamma * dn1[1];
  mat(1, 1) = dx2[1] - gamma * dn2[1];

  std::cout << "lambda matrix = " << mat << std::endl;

  double dfdgamma = computed_fd_gamma(data2, gamma);
  std::cout << "dfdgamma = " << dfdgamma << std::endl;
}

Data read_data()
{
  Data data;
  std::ifstream file;
  file.open("newton_data.dat", std::ios_base::binary);

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    {
      double val;
      file.read(reinterpret_cast<char*>(&val), sizeof(double));
      data.verts[i][j] = val;
    }

  for (int i = 0; i < 3; ++i)
    file.read(reinterpret_cast<char*>(&(data.pt[i])), sizeof(double));

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    {
      double val;
      file.read(reinterpret_cast<char*>(&val), sizeof(double));
      data.normals[i][j] = val;
    }

  return data;
}

int main(int argc, char* argv[])
{
  Data data = read_data();
  print_inputs(data);

  Data data2 = rotate_data(data);

  print_inputs_for_visualization(data2);

  std::cout << std::setprecision(16);
  // std::cout << "\nChecking singular" << std::endl;
  check_singular(data, 0);

  std::cout << "\noutputting f values" << std::endl;
  output_f_values(0, 200, 10000, data2);

  std::cout << "\nSolving Gamma brute force" << std::endl;
  solve_gamma_brute_force(60, 70, 1000, data2);

  std::cout << "\nSolving NewtonGamma" << std::endl;
  solve_newton_gamma(data2);

  std::cout << "\nSolving Newton 3x3" << std::endl;
  solve_newton3x3(data2);

  std::cout << "\nSolving Newton 3x3 original" << std::endl;
  solve_newton3x3(data);

  // std::cout << "\nSolving fixed point iteration" << std::endl;
  // solveGammaFixedPointIteration(data, 0);

  return 0;
}
