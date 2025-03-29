#include <Akri_Optimize.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <functional>
#include <tuple>
#include <vector>

namespace krino {

std::tuple<double,double> line_search_armijo_1d(const std::function<double(double)> & fn, const double f0, const double dirDeriv0, const double xmin)
{
  // Uses the interpolation algorithm (Armijo backtracking)

  static constexpr double xmax = 1.;
  static constexpr double c1 = 1.e-4;

  double xk = xmax;
  double fk = fn(xk);

  if (fk <= f0 + c1*xk*dirDeriv0)
    return {xk, fk};

  double xkp1 = -dirDeriv0 * xk*xk / (2.*(fk - f0 - dirDeriv0 * xk));
  double fkp1 = fn(xkp1);

  if (fkp1 <= f0 + c1*xkp1*dirDeriv0)
    return {xkp1, fkp1};

  while (xkp1 > xmin)
  {
    const double factor = 1. / (xk*xk * xkp1*xkp1 * (xkp1-xk));
    const double a = factor * (xk*xk * (fkp1 - f0 - dirDeriv0*xkp1) - xkp1*xkp1 * (fk - f0 - dirDeriv0*xk));
    const double b = factor * (-std::pow(xk,3) * (fkp1 - f0 - dirDeriv0*xkp1) + std::pow(xkp1,3) * (fk - f0 - dirDeriv0*xk));
    double xkp2 = (-b + std::sqrt(std::abs(b*b - 3 * a * dirDeriv0))) / (3.0*a);
    const double fkp2 = fn(xkp2);

    if (fkp2 <= f0 + c1*xkp2*dirDeriv0)
      return {xkp2, fkp2};

    if ((xkp1 - xkp2) > xkp1 / 2.0 || (1 - xkp2/xkp1) < 0.96)
      xkp2 = xkp1 / 2.0;

    xk = xkp1;
    xkp1 = xkp2;
    fk = fkp1;
    fkp1 = fkp2;
  }

  return {0., f0};
}

std::vector<double> xpby(const std::vector<double> & x, const double b, const std::vector<double> & y)
{
  std::vector<double> result;
  result.reserve(x.size());
  for (size_t i=0; i<x.size(); ++i)
    result.push_back(x[i]+b*y[i]);
  return result;
}

stk::math::Vector3d xpby(const stk::math::Vector3d & x, const double b, const stk::math::Vector3d & y)
{
  return x + b*y;
}

double Dot(const std::vector<double> & x, const std::vector<double> & y)
{
  STK_ThrowAssert(x.size() == y.size());
  double dot = 0.;
  for (size_t i=0; i<x.size(); ++i)
    dot += x[i]*y[i];
  return dot;
}

template<typename VEC>
std::tuple<VEC,double> line_search_armijo_vector(const std::function<double(const VEC&)> & fn,
    const VEC& x0,
    const VEC& dir,
    const double f0,
    const VEC& gradf0,
    const double xmin)
{
  const auto fn1d = [&](const double x1d) { return fn(xpby(x0, x1d, dir)); };
  const double dirDeriv = Dot(dir, gradf0);
  const auto & [x1d, obj] = line_search_armijo_1d(fn1d, f0, dirDeriv, xmin);
  return { xpby(x0, x1d, dir), obj };
}

template<typename VEC>
VEC zero_vector(const size_t);

template<>
stk::math::Vector3d zero_vector<stk::math::Vector3d>(const size_t) { return stk::math::Vector3d::ZERO; }

template<>
std::vector<double> zero_vector<std::vector<double>>(const size_t n) { return std::vector<double>(n, 0.); }


template<typename VEC>
VEC bfgs(const std::function<double(const VEC&)> & fn,
    const std::function<VEC(const VEC&)> & gradient,
    const VEC& x0,
    const double tol,
    const int maxIter)
{
  size_t n = x0.size();
  VEC x = x0;
  std::vector<VEC> H(n, zero_vector<VEC>(n));

  // Initialize H as the identity matrix
  for (size_t i = 0; i < n; ++i)
    H[i][i] = 1.0;

  VEC g = gradient(x);
  int iter = 0;

  while (std::sqrt(Dot(g,g)) > tol && iter < maxIter)
  {
    // Compute the search direction
    VEC p = zero_vector<VEC>(n);
    for (size_t i = 0; i < n; ++i)
    {
        p[i] = -H[i][i] * g[i]; // Simplified for identity matrix
    }

    const auto [x_new, f_new] = line_search_armijo_vector(fn, x, p, fn(x), g);

    VEC s = zero_vector<VEC>(n);
    for (size_t i = 0; i < n; ++i)
        s[i] = x_new[i] - x[i];

    if (Dot(s,s) == 0.)
    {
      return x;
    }

    VEC g_new = gradient(x_new);

    VEC y = zero_vector<VEC>(n);
    for (size_t i = 0; i < n; ++i)
        y[i] = g_new[i] - g[i];

    if (Dot(y,s) == 0.)
    {
      return x;
    }

    // Update H using the BFGS formula
    const double rho = 1.0 / Dot(y,s);
    std::vector<VEC> Hs(n, zero_vector<VEC>(n));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            Hs[i][j] = H[i][j] + rho * (s[i] * y[j] + y[i] * s[j]) - (rho * Dot(y,H[i]) * s[i] * s[j]);
        }
    }
    H = Hs;

    // Update x and g
    x = x_new;
    g = g_new;
    iter++;
  }

  return x;
}

template stk::math::Vector3d bfgs(const std::function<double(const stk::math::Vector3d&)> & fn,
    const std::function<stk::math::Vector3d(const stk::math::Vector3d&)> & gradient,
    const stk::math::Vector3d& x0,
    const double tol,
    const int maxIter);

template std::vector<double> bfgs(const std::function<double(const std::vector<double>&)> & fn,
    const std::function<std::vector<double>(const std::vector<double>&)> & gradient,
    const std::vector<double>& x0,
    const double tol,
    const int maxIter);


}

