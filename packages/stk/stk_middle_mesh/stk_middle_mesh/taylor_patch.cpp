#include "taylor_patch.hpp"

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

const std::vector<TaylorPatch::IntF> TaylorPatch::M_FACTORIALS = {
    factorial<TaylorPatch::IntF>(0), factorial<TaylorPatch::IntF>(1), factorial<TaylorPatch::IntF>(2),
    factorial<TaylorPatch::IntF>(3), factorial<TaylorPatch::IntF>(4), factorial<TaylorPatch::IntF>(5),
    factorial<TaylorPatch::IntF>(6), factorial<TaylorPatch::IntF>(7), factorial<TaylorPatch::IntF>(8),
    factorial<TaylorPatch::IntF>(9), factorial<TaylorPatch::IntF>(10)};

void TaylorPatch::construct_patch(const std::vector<Point>& pts, const Point& pt0)
{
  m_pt0 = pt0;
  compute_taylor_derivs(pts);
}

Point TaylorPatch::eval_point(const double x, const double y)
{
  // TODO: I think using TaylorPatch for mesh quality optimization
  //       is numerically unstable.  The z coordinate slowly drifts
  //       away from the value (even for a plane)
  // return Point(x, y, m_pt0.z);

  Point ptUv = compute_uv(Point(x, y));
  double val = 0;
  int jidx   = 0;

  for (int degree = 0; degree <= m_degree; ++degree)
    for (int term = 0; term <= degree; ++term)
    {
      val += eval_coefficient(degree, term, ptUv) * m_fvals[jidx];
      ++jidx;
    }

  return compute_xy(ptUv) + Point(0, 0, val + m_pt0.z);
}

void TaylorPatch::eval_deriv(const double x, const double y, double derivs[2])
{
  Point ptUv  = compute_uv(Point(x, y));
  double dzDu = 0;
  double dzDv = 0;

  int jidx = 0;
  for (int degree = 0; degree <= m_degree; ++degree)
    for (int term = 0; term <= degree; ++term)
    {
      dzDu += eval_coefficient_du(degree, term, ptUv) * m_fvals[jidx];
      dzDv += eval_coefficient_dv(degree, term, ptUv) * m_fvals[jidx];
      ++jidx;
    }

  // u = x - pt0.x -> dz/du = dz/dx
  // z = f(u, v) + pt0.z -> dz/du = df/du, dz/dv = df/dz
  derivs[0] = dzDu;
  derivs[1] = dzDv;
}

Point TaylorPatch::compute_uv(const Point& pt)
{
  return Point(pt.x - m_pt0.x, pt.y - m_pt0.y, 0);
}

Point TaylorPatch::compute_xy(const Point& ptUv)
{
  return Point(ptUv.x + m_pt0.x, ptUv.y + m_pt0.y, 0);
}

std::pair<int, int> TaylorPatch::compute_degree(const int npts)
{
  int nterms = 0;
  int degree = -1;
  while (nterms < npts)
  {
    degree++;
    nterms += degree + 1;
  }

  return std::make_pair(degree, nterms);
}

double TaylorPatch::eval_coefficient(const int degree, const int term, const Point& ptUv)
{
  assert(degree >= 0);
  assert(term >= 0 && term <= degree);

  int ypow = term;
  int xpow = degree - term;
  // protect against case where delta_x = 0 and exponent is zero
  double xval  = xpow == 0 ? 1 : std::pow(ptUv.x, xpow);
  double yval  = ypow == 0 ? 1 : std::pow(ptUv.y, ypow);
  double coeff = (term != 0 && term != degree) ? degree : 1;

  return coeff * xval * yval / get_factorial(degree);
}

double TaylorPatch::eval_coefficient_du(const int degree, const int term, const Point& ptUv)
{
  assert(degree >= 0);
  assert(term >= 0 && term <= degree);

  int ypow = term;
  int xpow = degree - term;
  // protect against case where delta_x = 0 and exponent is zero
  double xval;
  if (xpow == 0)
    xval = 0;
  else if (xpow == 1)
    xval = 1;
  else
    xval = xpow * std::pow(ptUv.x, xpow - 1);
  // double xval = xpow == 0 ? 0 : xpow*std::pow(u, xpow - 1);
  double yval  = ypow == 0 ? 1 : std::pow(ptUv.y, ypow);
  double coeff = (term != 0 && term != degree) ? degree : 1;

  return coeff * xval * yval / get_factorial(degree);
}

double TaylorPatch::eval_coefficient_dv(const int degree, const int term, const Point& ptUv)
{
  assert(degree >= 0);
  assert(term >= 0 && term <= degree);

  int ypow = term;
  int xpow = degree - term;

  // protect against case where delta_x = 0 and exponent is zero
  double yval;
  if (ypow == 0)
    yval = 0;
  else if (ypow == 1)
    yval = 1;
  else
    yval = ypow * std::pow(ptUv.y, ypow - 1);
  // double xval = xpow == 0 ? 0 : xpow*std::pow(u, xpow - 1);
  double xval  = xpow == 0 ? 1 : std::pow(ptUv.x, xpow);
  double coeff = (term != 0 && term != degree) ? degree : 1;

  return coeff * xval * yval / get_factorial(degree);
}

void TaylorPatch::compute_taylor_matrix(const std::vector<Point>& pts, utils::impl::Matrix<double>& a)
{
  a.resize(pts.size(), m_numTerms);
  // A.fill(0);
  for (unsigned int i = 0; i < pts.size(); ++i)
  {
    Point ptUv = compute_uv(pts[i]);
    int jidx   = 0;
    for (int degree = 0; degree <= m_degree; ++degree)
      for (int term = 0; term <= degree; ++term)
      {
        a(i, jidx) = eval_coefficient(degree, term, ptUv);
        ++jidx;
      }

    assert(jidx == m_numTerms);
  }
}

void TaylorPatch::compute_taylor_derivs(const std::vector<Point>& pts)
{
  // TODO: make all dynamic memory objects into class members
  auto p     = compute_degree(pts.size());
  m_degree   = p.first;
  m_numTerms = p.second;

  // get coefficient matrix
  utils::impl::Matrix<double> a(pts.size(), m_numTerms);
  utils::impl::Matrix<double> work(a.extent0() + 1, a.extent1() + 1);
  compute_taylor_matrix(pts, a);

  // setup rhs
  utils::impl::Matrix<double> rhs(std::max(static_cast<int>(pts.size()), m_numTerms), 1);
  for (unsigned int i = 0; i < pts.size(); ++i)
    rhs(i, 0) = pts[i].z - m_pt0.z;

  // check for rank deficiency
  // checks if the diagonal entries of L are zero, replaces the corresponding
  // row in A with the identity, and zeros the rhs
  utils::impl::Matrix<double> b(a);
  std::vector<double> tau(std::max(a.extent0(), a.extent1()));
  compute_qr_factorization(b, work, tau.data());
  const double tol = 1e-13;
  for (int i = 0; i < a.extent0(); ++i)
  {
    if (std::abs(b(i, i)) < tol)
    {
      for (int j = 0; j < a.extent1(); ++j)
        a(i, j) = 0;

      a(i, i)   = 1;
      rhs(i, 0) = 0;
    }
  }
  /*
    // check
    std::cout << "checking" << std::endl;
    for (int i=0; i < A.extent0(); ++i)
      for (int j=0; j < A.extent1(); ++j)
        B(i, j) = A(i, j);

    compute_qr_factorization(B, work, tau.data());
    for (int i=0; i < A.extent0(); ++i)
      std::cout << "i = " << i << ", pivot = " << B(i, i) << std::endl;

    std::cout << "A = " << A << std::endl;
  */

  // solve least squares
  solve_least_squares(a, rhs, work);

  // save results
  m_fvals.resize(m_numTerms);
  for (int i = 0; i < m_numTerms; ++i)
  {
    m_fvals[i] = rhs(i, 0);
  }
}

TaylorPatch::IntF TaylorPatch::get_factorial(const unsigned int n)
{
  assert(n < M_FACTORIALS.size());
  return M_FACTORIALS[n];
}

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
