#include "stk_middle_mesh/predicates/edge_intersection_primitive.hpp"
#include "stk_middle_mesh/mat2x2.hpp"

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

EdgeIntersectionResult::EdgeIntersectionResult(double alpha, double beta)
  : m_intersection1Found(true)
  , m_alpha1(alpha)
  , m_beta1(beta)
{}

EdgeIntersectionResult::EdgeIntersectionResult(double alpha1, double beta1, double alpha2, double beta2)
  : m_intersection1Found(true)
  , m_alpha1(alpha1)
  , m_beta1(beta1)
  , m_intersection2Found(true)
  , m_alpha2(alpha2)
  , m_beta2(beta2)
{}

EdgeIntersectionResult::EdgeIntersectionResult()
  : m_intersection1Found(false)
{}

EdgeIntersectionData::EdgeIntersectionData(const utils::Point& b0, const utils::Point& b1, const utils::Point& g0,
                                           const utils::Point& g1, const utils::Point& d0, const utils::Point& d1,
                                           double tol, double errorTol)
  : m_b0(b0)
  , m_b1(b1)
  , m_g0(g0)
  , m_g1(g1)
  , m_d0(d0)
  , m_d1(d1)
  , m_tol(tol)
  , m_errorTol(errorTol)
{}

void EdgeIntersectionData::compute_solution_in_weird_way()
{
  bool alpha1Valid = false, alpha2Valid = false, beta1Valid = false, beta2Valid = false;
  double alpha1 = 666, alpha2 = 666, beta1 = 666, beta2 = 666;

  utils::Point deltaB = m_b1 - m_b0; // vb in other code
  utils::Point deltaG = m_g1 - m_g0; // vg in other code
  utils::Point deltaD = m_d1 - m_d0; // vn in other code
  utils::Point b0ToG0 = m_g0 - m_b0; // vd in other code

  double c = dot(cross(deltaB, m_d0), b0ToG0);
  double d = dot(cross(deltaB, m_d1), m_g1 - m_b0);
  double a = dot(cross(deltaB, deltaD), deltaG);
  double b = d - c - a;

  double normalizationFac = std::max(std::abs(a), std::abs(b));
  normalizationFac        = std::max(normalizationFac, std::abs(c));
  normalizationFac        = std::max(normalizationFac, std::abs(d));

  if (normalizationFac > 0)
  {
    a /= normalizationFac;
    b /= normalizationFac;
    c /= normalizationFac;
    d /= normalizationFac;
  }

  // std::cout << "a = " << a << std::endl;
  // std::cout << "b = " << b << std::endl;
  // std::cout << "c = " << c << std::endl;
  // std::cout << "d = " << d << std::endl;

  double eps    = 0.1;
  double tolFea = 1;

  if (a == 0)
  {
    if (b == 0)
    {
      beta1Valid = false;
      beta2Valid = false;
      return;
    }

    beta1Valid = true;
    beta2Valid = false;
    beta1      = -c / b;

    if (beta1 >= -tolFea && beta1 <= 1 + tolFea)
    {
      double denominator = (1 - beta1) * dot(cross(deltaG, m_g0), deltaB) + beta1 * dot(cross(deltaG, m_g0), deltaB);
      double numerator   = (1 - beta1) * dot(cross(deltaG, m_d0), deltaD) + beta1 * dot(cross(deltaG, m_d1), deltaD);
      alpha1             = numerator / denominator;
      alpha1Valid        = true;
      alpha2Valid        = false;
    } else
    {
      alpha1Valid = false;
      alpha2Valid = false;
    }
  } else // quadratic case
  {
    // std::cout << "quadratic case" << std::endl;
    double det = b * b - 4 * a * c;
    // std::cout << "det = " << det << std::endl;
    if (det < 0)
    {
      beta1Valid  = false;
      beta2Valid  = false;
      alpha1Valid = false;
      alpha2Valid = false;
      return;
    }

    double sqrtD = std::sqrt(det);

    if (b > 0)
    {
      // std::cout << "b > 0 case " << std::endl;
      double temp = -b - sqrtD;
      beta1       = 2 * c / temp; // Mullers method
      beta2       = temp / 2 / a;

      // std::cout << "if using standard quadratic formulata, betas would be: ";
      // std::cout << (-b + sqrt_d)/(2*a) << ", " << (-b - sqrt_d)/(2*a) << std::endl;
    } else
    {
      // std::cout << "b < 0 case " << std::endl;
      double temp = -b + sqrtD;
      beta1       = temp / 2 / a;
      beta2       = 2 * c / temp;
    }

    beta1Valid = true;
    beta2Valid = true;

    if (std::abs(beta1 - 1) < eps || std::abs(beta2 - 1) < eps)
    {
      // std::cout << "using alternative formula" << std::endl;
      //  use alternative formula to compute intersection
      b     = -(c - d - a);
      sqrtD = std::sqrt(b * b - 4 * a * d);
      if (det < 0) // TODO: det was not updated after redefining b, so this condition has been checked before
      {
        // std::cout << "det < 0, returning" << std::endl;
        alpha1Valid = false;
        alpha2Valid = false;
        beta1Valid  = false;
        beta2Valid  = false;
        return;
      }

      if (b > 0)
      {
        double temp = -b - sqrtD;
        if (std::abs(beta1 - 1) < eps)
          beta1 = 1 + 2 * d / temp;
        if (std::abs(beta2 - 1) < eps)
          beta2 = 1 + temp / 2 / a;
      } else
      {
        double temp = -b + sqrtD;
        if (std::abs(beta1 - 1) < eps)
          beta1 = 1 + temp / 2 / a;
        if (std::abs(beta2 - 1) < eps)
          beta2 = 1 + 2 * d / temp;
      }
    }

    if (beta1 >= -tolFea && beta1 <= 1 + tolFea)
    {
      // std::cout << "beta1 within tol_fea" << std::endl;
      double denominator = (1 - beta1) * dot(cross(deltaG, m_d0), deltaB) + beta1 * dot(cross(deltaG, m_d1), deltaB);
      double numerator   = (1 - beta1) * dot(cross(deltaG, m_d0), b0ToG0) + beta1 * dot(cross(deltaG, m_d1), b0ToG0);

      alpha1      = numerator / denominator;
      alpha1Valid = true;
    } else
    {
      // std::cout << "beta1 not within tol_fea" << std::endl;
      alpha1Valid = false;
    }

    if (beta2 >= -tolFea && beta2 <= 1 + tolFea)
    {
      // std::cout << "beta2 within tol_fea" << std::endl;

      double denominator = (1 - beta2) * dot(cross(deltaG, m_d0), deltaB) + beta2 * dot(cross(deltaG, m_d1), deltaB);
      double numerator   = (1 - beta2) * dot(cross(deltaG, m_d0), b0ToG0) + beta2 * dot(cross(deltaG, m_d1), b0ToG0);
      alpha2             = numerator / denominator;
      alpha2Valid        = true;
    } else
    {
      // std::cout << "beta2 not within tol_fea" << std::endl;
      alpha2Valid = false;
    }

    std::cout << std::boolalpha;
    std::cout << "alphas valid: " << alpha1Valid << ", " << alpha2Valid << std::endl;
    std::cout << "betas valid: " << beta1Valid << ", " << beta2Valid << std::endl;
    std::cout << "first solution: " << alpha1 << ", " << beta1 << std::endl;
    std::cout << "second solution: " << alpha2 << ", " << beta2 << std::endl;
  }
}

void EdgeIntersectionData::compute_beta()
{
  // std::cout << "b0 = " << m_b0 << ", b1 = " << m_b1 << std::endl;
  // std::cout << "g0 = " << m_g0 << ", g1 = " << m_g1 << std::endl;
  // std::cout << "d0 = " << m_d0 << ", d1 = " << m_d1 << std::endl;
  // std::cout << "delta_b = " << m_b1 - m_b0 << ", delta_d = " << m_d1 - m_d0 << std::endl;
  // std::cout << "cross = " << cross(m_b1 - m_b0, m_d1 - m_d0) << std::endl;
  // std::cout << "delta_g = " << m_g1 - m_g0 << std::endl;
  double c2 = dot(cross((m_b1 - m_b0), (m_d1 - m_d0)), (m_g1 - m_g0));
  double c1 = dot(cross((m_b1 - m_b0), m_d0), (m_g1 - m_g0)) + dot(cross((m_b1 - m_b0), (m_d1 - m_d0)), (m_g0 - m_b0));
  double c0 = dot(cross((m_b1 - m_b0), m_d0), (m_g0 - m_b0));

  // std::cout << "initially, c2 = " << c2 << ", c1 = " << c1 << ", c0 = " << c0 << std::endl;
  if (std::abs(c2) < m_tol && std::abs(c1) < m_tol && std::abs(c0) < m_tol)
  {
    // in this case the lines are parallel.  Report no solution, even
    // though it is possible the lines overlap
    m_beta1Valid = false;
    m_beta2Valid = false;
    return;
  }

  double normalizeFac = std::max(std::abs(c2), std::abs(c1));
  normalizeFac        = std::max(normalizeFac, std::abs(c0));
  // std::cout << "normalize_fac = " << normalize_fac << std::endl;
  c2 = c2 / normalizeFac;
  c1 = c1 / normalizeFac;
  c0 = c0 / normalizeFac;

  // check the value of the discriminant to determine if there is a valid (single) solution
  // std::cout << "d0 = " << m_d0 << ", d1 = " << m_d1 << ", diff = " << m_d1 - m_d0 << std::endl;
  // std::cout << "c2 = " << c2 << ", c1 = " << c1 << ", c0 = " << c0 << std::endl;
  // std::cout << "if linear, beta would be " << -c0/c1 << std::endl;

  if (std::abs(c2) <= m_tol)
  {
    m_beta2Valid = false;
    if (std::abs(c1) <= m_tol)
    {
      m_beta1Valid = false;
    } else
    {
      m_beta1 = -c0 / c1;
      // std::cout << "beta1 = " << m_beta1 << std::endl;
      m_beta1Valid = true;
    }
  } else
  {
    double discriminant = c1 * c1 - 4 * c2 * c0;
    // std::cout << "discriminant = " << discriminant << std::endl;

    if (discriminant < 0)
    {
      m_beta1Valid = false;
      m_beta2Valid = false;
    }

    else
    {
      m_beta1      = (-c1 + std::sqrt(discriminant)) / (2 * c2);
      m_beta2      = (-c1 - std::sqrt(discriminant)) / (2 * c2);
      m_beta1Valid = true;
      m_beta2Valid = true;
      // std::cout << "beta1 = " << m_beta1 << std::endl;
      // std::cout << "beta2 = " << m_beta2 << std::endl;
    }
  }
}

namespace {
constexpr double clamp(double val, double minVal, double maxVal)
{
  double val1 = std::max(val, minVal);
  return std::min(val1, maxVal);
}
} // namespace

// solves a linear system for alpha, the distance along the blue edge,
// and gamma, the distance along the direction vector that passes
// through the point on the green edge defined by beta
// Only alpha is returned because that is all we need right now.
// Uses an angle check to see if the solution is an intersection point
double EdgeIntersectionData::solve_alpha_linear_system(double beta, bool& alphaValid)
{
  if (std::abs(beta) > 10)
  {
    alphaValid = false;
    return -1;
  }
  utils::Point c        = m_g0 + beta * (m_g1 - m_g0);
  utils::Point d        = m_d0 + beta * (m_d1 - m_d0);
  utils::Point deltaB   = m_b1 - m_b0;
  double dDeltaB        = dot(d, deltaB);
  utils::Point cMinusB0 = c - m_b0;

  // std::cout << "\nbeta = " << beta << std::endl;
  // std::cout << "c = " << c << std::endl;
  // std::cout << "d = " << d << std::endl;

  double dLen = std::sqrt(dot(d, d));
  // std::cout << "d_len = " << d_len << std::endl;
  if (dLen < m_tol)
  {
    // std::cout << "d = 0, solving single equation" << std::endl;
    double alpha = dot(cMinusB0, deltaB) / dot(deltaB, deltaB);
    alphaValid   = true;
    return alpha;
  }

  utils::impl::Mat2x2<double> a;
  a(0, 0) = dot(deltaB, deltaB);
  a(0, 1) = -dDeltaB;
  a(1, 0) = dDeltaB;
  a(1, 1) = -dot(d, d);

  if (std::abs(det2x2(a)) < m_tol)
  {
    alphaValid = false;
    return -1;
  }

  std::array<double, 2> rhs = {dot(cMinusB0, deltaB), dot(cMinusB0, d)};

  std::array<double, 2> x;
  utils::impl::matsolve2x2(a, x.data(), rhs.data());

  // std::cout << "A = " << A << std::endl;
  // std::cout << "rhs = " << rhs[0] << ", " << rhs[1] << std::endl;

  // check angle
  double alpha = x[0];
  // std::cout << "solved for alpha = " << alpha << ", gamma = " << x[1] << std::endl;

  utils::Point b2 = m_b0 + alpha * deltaB;
  // std::cout << "b2 = " << b2 << std::endl;
  utils::Point bDirection = b2 - c;
  double bDirectionLen    = std::sqrt(dot(bDirection, bDirection));
  // std::cout << "b_direction = " << b_direction << std::endl;
  if (bDirectionLen < m_tol)
  {
    alphaValid = true;
    return alpha;
  }
  bDirection = bDirection / bDirectionLen;

  utils::Point dUnit = d / std::sqrt(dot(d, d));

  // std::cout << "d_unit = " << d_unit << std::endl;
  double cosTheta = std::abs(dot(bDirection, dUnit));
  cosTheta        = clamp(cosTheta, 0, 1);

  // std::cout << "point from gamma = " << c + x[1]*d << std::endl;
  // std::cout << "point from alpha = " << m_b0 + alpha*(m_b1 - m_b0) << std::endl;

  assert(!std::isinf(cosTheta) && !std::isnan(cosTheta));

  alphaValid = std::abs(cosTheta - 1) < m_errorTol;

  return alpha;
}

void EdgeIntersectionData::compute_alpha()
{
  if (m_beta1Valid == false)
  {
    m_alpha1Valid = false;
  } else
  {
    // std::cout << "\nsolving for alpha1" << std::endl;

    m_alpha1 = solve_alpha_linear_system(m_beta1, m_alpha1Valid);
    // std::cout << "error angle = " << compute_error_angle(m_alpha1, m_beta1) * 180/3.14159265 << std::endl;
  }

  if (m_beta2Valid == false)
  {
    m_alpha2Valid = false;
  } else
  {
    // std::cout << "\nsolving for alpha2" << std::endl;

    m_alpha2 = solve_alpha_linear_system(m_beta2, m_alpha2Valid);
    // std::cout << "error angle = " << compute_error_angle(m_alpha2, m_beta2) * 180/3.14159265 << std::endl;
  }
}

/*
void EdgeIntersectionData::compute_alpha()
{
  //std::cout << "beta1 = " << m_beta1 << ", beta2 = " << m_beta2 << std::endl;
  if (m_beta1Valid == false) {
    m_alpha1Valid = false;
  }
  else {
    utils::Point l_1 = cross((m_g1 - m_g0), m_d0 + m_beta1*(m_d1 - m_d0));
    double denom_1 = dot(l_1, (m_b1 - m_b0));

    //std::cout << "denom_1 = " << denom_1 << std::endl;
    //std::cout << "m_tol = " << m_tol << std::endl;
    //std::cout << "possible alpha value = " << (dot(l_1, (m_g0 - m_b0))) / (dot(l_1, (m_b1 - m_b0))) << std::endl;

    if (std::abs(denom_1) > m_tol) {
      m_alpha1Valid = true;
      m_alpha1 = (dot(l_1, (m_g0 - m_b0))) / (dot(l_1, (m_b1 - m_b0)));
    }
    else {
      m_alpha1Valid = false;
    }
  }

  if (m_beta2Valid == false) {
    m_alpha2Valid = false;
  }
  else {
    utils::Point l_2 = cross((m_g1 - m_g0), m_d0 + m_beta2*(m_d1 - m_d0));
    double denom_2 = dot(l_2, (m_b1 - m_b0));

    //std::cout << "num_2   = " << (dot(l_2, (m_g0 - m_b0))) << std::endl;
    //std::cout << "denom_2 = " << denom_2 << std::endl;
    //std::cout << "m_tol = " << m_tol << std::endl;
    //std::cout << "possible alpha value = " << (dot(l_2, (m_g0 - m_b0))) / (dot(l_2, (m_b1 - m_b0))) << std::endl;

    if (std::abs(denom_2) > m_tol) {
      m_alpha2Valid = true;
      m_alpha2 = (dot(l_2, (m_g0 - m_b0))) / (dot(l_2, (m_b1 - m_b0)));
    }
    else {
      m_alpha2Valid = false;
    }
  }

  //std::cout << "alpha1 = " << m_alpha1 << std::endl;
  //std::cout << "alpha2 = " << m_alpha2 << std::endl;
}
*/

double EdgeIntersectionData::compute_error_angle(double alpha, double beta)
{
  // std::cout << "first pt = " << (m_g0 + beta*(m_g1 - m_g0)) << ", second pt = " << (m_b0 + alpha*(m_b1 - m_b0)) <<
  // std::endl; utils::Point displacement_unit = (m_b0 + alpha*(m_b1 - m_b0)) - (m_g0 + beta*(m_g1 - m_g0));
  utils::Point displacementUnit = (m_g0 + beta * (m_g1 - m_g0)) - (m_b0 + alpha * (m_b1 - m_b0));
  double displacementLength     = std::sqrt(dot(displacementUnit, displacementUnit));

  if (displacementLength < 1e-12)
    return 1;

  displacementUnit = displacementUnit / displacementLength;
  // std::cout << "displacement_unit = " << displacement_unit << std::endl;

  utils::Point directionUnit = m_d0 + beta * (m_d1 - m_d0);
  double directionLen        = std::sqrt(dot(directionUnit, directionUnit));
  if (directionLen < 1e-12)
    return 1;

  directionUnit = directionUnit / directionLen;

  // std::cout << "direction_unit = " << direction_unit << std::endl;

  double cosTheta = std::abs(dot(directionUnit, displacementUnit));
  // std::cout << "cos_theta = " << cos_theta << std::endl;
  cosTheta = clamp(cosTheta, -1, 1);

  return cosTheta;
}

EdgeIntersectionResult EdgeIntersectionData::compute_edge_intersection_result()
{
  bool solution1Valid = (m_beta1Valid && m_alpha1Valid);
  bool solution2Valid = (m_beta2Valid && m_alpha2Valid);

  // std::cout << std::boolalpha;
  // std::cout << "solution1 valid: " << m_alpha1Valid << ", " << m_beta1Valid << std::endl;
  // std::cout << "solution2 valid: " << m_alpha2Valid << ", " << m_beta2Valid << std::endl;
  // std::cout << "solution 1: " << m_alpha1 << ", " << m_beta1 << std::endl;
  // std::cout << "solution 2: " << m_alpha2 << ", " << m_beta2 << std::endl;
  // double pi = std::atan(1)*4;

  if (solution1Valid == false && solution2Valid == false)
  {
    return EdgeIntersectionResult();
  }

  else if (solution1Valid == true and solution2Valid == false)
  {
    double cosTheta = compute_error_angle(m_alpha1, m_beta1);

    if (std::abs(cosTheta - 1) < m_errorTol)
    {
      return EdgeIntersectionResult(m_alpha1, m_beta1);
    } else
    {
      return EdgeIntersectionResult();
    }
  }

  else if (solution1Valid == false and solution2Valid == true)
  {
    double cosTheta = compute_error_angle(m_alpha2, m_beta2);
    // std::cout << "solution 2 error = " << error*180/pi << std::endl;

    if (std::abs(cosTheta - 1) < m_errorTol)
    {
      return EdgeIntersectionResult(m_alpha2, m_beta2);
    } else
    {
      return EdgeIntersectionResult();
    }
  }

  else
  { // both solutions valid

    double cosTheta1 = compute_error_angle(m_alpha1, m_beta1);
    double cosTheta2 = compute_error_angle(m_alpha2, m_beta2);
    // std::cout << "solution 1 error = " << error_1*180/pi << std::endl;
    // std::cout << "solution 2 error = " << error_2*180/pi << std::endl;

    if (std::abs(cosTheta1 - 1) < m_errorTol && std::abs(cosTheta2 - 1) < m_errorTol)
    {
      return EdgeIntersectionResult(m_alpha1, m_beta1, m_alpha2, m_beta2);
    } else if (std::abs(cosTheta1 - 1) < m_errorTol)
    {
      return EdgeIntersectionResult(m_alpha1, m_beta1);
    } else if (std::abs(cosTheta2 - 1) < m_errorTol)
    {
      return EdgeIntersectionResult(m_alpha2, m_beta2);
    } else
    {
      return EdgeIntersectionResult();
    }
  }
}

EdgeIntersectionResult compute_edge_intersection(const utils::Point& b0, const utils::Point& b1, const utils::Point& g0,
                                                 const utils::Point& g1, const utils::Point& d0, const utils::Point& d1,
                                                 double tol, double errorTol)
{
  // std::cout << "\nComputing solution in normal way" << std::endl;
  // std::cout << "tol = " << tol << ", error_tol = " << error_tol << std::endl;
  EdgeIntersectionData edgeIntersectionData(b0, b1, g0, g1, d0 / std::sqrt(dot(d0, d0)), d1 / std::sqrt(dot(d1, d1)),
                                            tol, errorTol);
  edgeIntersectionData.compute_beta();
  edgeIntersectionData.compute_alpha();

  auto result = edgeIntersectionData.compute_edge_intersection_result();
  // std::cout << "\nComputing solution in weird way" << std::endl;
  // edgeIntersectionData.compute_solution_in_weird_way();

  return result;
}

} // namespace impl
} // namespace predicates
} // namespace middle_mesh
} // namespace stk
