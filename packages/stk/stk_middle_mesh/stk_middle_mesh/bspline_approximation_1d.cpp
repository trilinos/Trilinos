#include "bspline_approximation_1d.hpp"

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

void BSplineApproximation1D::approximate(const std::vector<double>& xValues, const std::vector<double>& fValues)
{
  assert(xValues.size() > 0);
  assert(xValues.size() == fValues.size());

  setup_local_coordinate_system(xValues[0], xValues[xValues.size() - 1]);
  // computeControlPointValues(x_values, f_values);
  compute_control_points_least_squares(xValues, fValues);
}

double BSplineApproximation1D::eval(double x)
{
  int controlPointStart = get_left_control_point(x) - 3;
  double t              = compute_t_value(x);

  std::array<double, 4> bsplineVals;
  m_bsplines.eval(t, bsplineVals);

  double val = 0;
  for (int i = 0; i < 4; ++i)
    val += m_controlPointVals[i + controlPointStart] * bsplineVals[i];

  return val;
}

void BSplineApproximation1D::setup_local_coordinate_system(double xmin, double xmax)
{
  m_deltaXControlPoints = (xmax - xmin) / (m_numControlPointsInterior - 1);
  // we need 3 extra control points to the left of the domain, so offset m_xmin
  // by that amount
  m_xmin = xmin - 3 * m_deltaXControlPoints;
}

int BSplineApproximation1D::get_left_control_point(double x)
{
  int leftControlPoint = std::floor((x - m_xmin) / m_deltaXControlPoints);
  assert(leftControlPoint >= 3);
  return leftControlPoint;
}

double BSplineApproximation1D::compute_t_value(double x)
{
  int leftControlPoint         = get_left_control_point(x);
  double leftControlPointCoord = m_xmin + leftControlPoint * m_deltaXControlPoints;

  double t = (x - leftControlPointCoord) / m_deltaXControlPoints;
  assert(t >= -1e-13 && t < 1 + 1e-13);

  return t;
}

/*
void BSplineApproximation1D::computeControlPointValues(const std::vector<double>& x_values, const std::vector<double>&
f_values)
{
  std::vector<double> numerator(m_control_point_vals.size(), 0);
  std::vector<double> denominator(m_control_point_vals.size(), 0);
  std::array<double, 4> bspline_vals;
  for (size_t i=0; i < f_values.size(); ++i)
  {
    double t = computeTValue(x_values[i]);
    int starting_control_point = getLeftControlPoint(x_values[i]) - 3;
    std::cout << "\ninput point " << i << "(" << x_values[i] << ", " << f_values[i] << "), t = " << t  << ", f_value = "
<< f_values[i] << std::endl; m_bsplines.eval(t, bspline_vals);

    double w_squared = 0;
    for (int j=0; j < 4; ++j)
    {
      w_squared += bspline_vals[j] * bspline_vals[j];
    }


    for (int j=0; j < 4; ++j)
    {
      double phi = f_values[i] * bspline_vals[j]/w_squared;
      std::cout << "phi " << starting_control_point + j << " = " << phi << std::endl;
      numerator[starting_control_point + j]   += bspline_vals[j]*bspline_vals[j]*phi;
      denominator[starting_control_point + j] += bspline_vals[j]*bspline_vals[j];
    }
  }

  for (size_t i=0; i < m_control_point_vals.size(); ++i)
  {
    if (denominator[i] != 0)
      m_control_point_vals[i] = numerator[i] / denominator[i];
    else
      m_control_point_vals[i] = 0;
    std::cout << "for control point " << i << ", numerator, denominator = " << numerator[i] << ", " << denominator[i] <<
", control point = " << m_control_point_vals[i] << std::endl;

    //m_control_point_vals[i] = 1; //TODO: DEBUGGING

  }
}
*/

void BSplineApproximation1D::compute_control_points_least_squares(const std::vector<double>& xValues,
                                                                  const std::vector<double>& fValues)
{
  Matrix<double> a(xValues.size(), m_controlPointVals.size() - 1);
  Matrix<double> work(xValues.size(), m_controlPointVals.size() - 1);
  Matrix<double> rhs(std::max(xValues.size(), m_controlPointVals.size() - 1), 1);
  a.fill(0);
  for (size_t i = 0; i < xValues.size(); ++i)
  {
    double t                 = compute_t_value(xValues[i]);
    int startingControlPoint = get_left_control_point(xValues[i]) - 3;

    std::array<double, 4> bsplineVals;
    m_bsplines.eval(t, bsplineVals);
    for (int j = 0; j < 4; ++j)
      if (startingControlPoint + j < a.extent1()) // skip the endpoint spline because it is
                                                  // rank deficient
        a(i, startingControlPoint + j) = bsplineVals[j];

    rhs(i, 0) = fValues[i];
  }

  solve_least_squares(a, rhs, work);

  for (size_t i = 0; i < m_controlPointVals.size() - 1; ++i)
    m_controlPointVals[i] = rhs(i, 0);
  m_controlPointVals[m_controlPointVals.size() - 1] = 0;
}

} // namespace impl
} // namespace utils
} // namespace middle_mesh
} // namespace stk
