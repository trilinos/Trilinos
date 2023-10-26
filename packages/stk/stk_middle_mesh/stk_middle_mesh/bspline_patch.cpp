#include "bspline_patch.hpp"
#include <iomanip>

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

void BSplinePatch::construct_patch(const std::vector<Point>& pts, const Point& pt0)
{
  BoundingBox box = compute_bounding_box(pts, pt0);
  setup_local_coordinate_system(box);
  compute_control_points_least_squares(pts);
}

Point BSplinePatch::eval_point(const double x, const double y)
{
  LocalCoords stVals         = compute_st_values(x, y);
  Indices cornerControlPoint = get_corner_control_point(x, y);
  Indices startingControlPoint{cornerControlPoint.i - 3, cornerControlPoint.j - 3};

  std::array<double, 4> bsplineValsS, bsplineValsT;
  m_bsplines.eval(stVals.s, bsplineValsS);
  m_bsplines.eval(stVals.t, bsplineValsT);

  double val = 0;
  for (int m = 0; m < 4; ++m)
    for (int n = 0; n < 4; ++n)
    {
      double basisVals = bsplineValsS[m] * bsplineValsT[n];
      val += m_controlPointVals(startingControlPoint.i + m, startingControlPoint.j + n) * basisVals;
    }

  return Point(x, y, val);
}

BSplinePatch::BoundingBox BSplinePatch::compute_bounding_box(const std::vector<Point>& pts, const Point& pt0)
{
  BoundingBox box{pt0.x, pt0.x, pt0.y, pt0.y};
  for (auto& pt : pts)
  {
    box.xmin = std::min(pt.x, box.xmin);
    box.xmax = std::max(pt.x, box.xmax);
    box.ymin = std::min(pt.y, box.ymin);
    box.ymax = std::max(pt.y, box.ymax);
  }

  return box;
}

void BSplinePatch::setup_local_coordinate_system(const BoundingBox& box)
{
  m_deltaXControlPoints = (box.xmax - box.xmin) / (m_numControlPointsPerDirectionInterior - 1);
  m_deltaYControlPoints = (box.ymax - box.ymin) / (m_numControlPointsPerDirectionInterior - 1);

  // we need 3 extra control points to the left of the domain, so offset m_xmin
  // by that amount
  m_xmin              = box.xmin - 3 * m_deltaXControlPoints;
  m_ymin              = box.ymin - 3 * m_deltaYControlPoints;
  m_xInteriorMidpoint = (box.xmin + box.xmax) / 2;
  m_yInteriorMidpoint = (box.ymin + box.ymax) / 2;
}

BSplinePatch::Indices BSplinePatch::get_corner_control_point(double x, double y)
{
  // these offsets avoid floating point problems where the control point will come out
  // as 2 rather than 3
  double xoffset = x > m_xInteriorMidpoint ? -m_eps : m_eps;
  double yoffset = y > m_yInteriorMidpoint ? -m_eps : m_eps;

  int leftControlPoint   = std::floor((x + xoffset - m_xmin) / m_deltaXControlPoints);
  int bottomControlPoint = std::floor((y + yoffset - m_ymin) / m_deltaYControlPoints);

  assert(leftControlPoint >= 3);
  assert(bottomControlPoint >= 3);
  return {leftControlPoint, bottomControlPoint};
}

BSplinePatch::LocalCoords BSplinePatch::compute_st_values(double x, double y)
{
  Indices idxs                   = get_corner_control_point(x, y);
  double leftControlPointCoord   = m_xmin + idxs.i * m_deltaXControlPoints;
  double bottomControlPointCoord = m_ymin + idxs.j * m_deltaYControlPoints;

  double s = (x - leftControlPointCoord) / m_deltaXControlPoints;
  double t = (y - bottomControlPointCoord) / m_deltaYControlPoints;

  assert(s >= -m_eps && s < 1 + m_eps);
  assert(t >= -m_eps && t < 1 + m_eps);

  return {s, t};
}

int BSplinePatch::compute_linear_index(const Indices& idxs)
{
  int numNonzeroControlPointsPerDirection = (m_numControlPointsPerDirectionInterior + 3 - 1);
  return idxs.i * numNonzeroControlPointsPerDirection + idxs.j;
}

void BSplinePatch::compute_control_points_least_squares(const std::vector<Point>& inputValues)
{
  int numNonzeroControlPointsPerDirection = (m_numControlPointsPerDirectionInterior + 3 - 1);
  int numNonzeroControlPoints             = numNonzeroControlPointsPerDirection * numNonzeroControlPointsPerDirection;

  // TODO: the matrix has tensor-product structure.  Is there a way to use this?
  Matrix<double> a(inputValues.size(), numNonzeroControlPoints);
  Matrix<double> work(inputValues.size(), numNonzeroControlPoints);
  Matrix<double> rhs(std::max(int(inputValues.size()), numNonzeroControlPoints), 1);

  a.fill(0);
  for (size_t i = 0; i < inputValues.size(); ++i)
  {
    LocalCoords stVals         = compute_st_values(inputValues[i].x, inputValues[i].y);
    Indices cornerControlPoint = get_corner_control_point(inputValues[i].x, inputValues[i].y);
    Indices startingControlPoint{cornerControlPoint.i - 3, cornerControlPoint.j - 3};

    std::array<double, 4> bsplineValsS, bsplineValsT;
    m_bsplines.eval(stVals.s, bsplineValsS);
    m_bsplines.eval(stVals.t, bsplineValsT);
    // setup matrix for least squares problem.
    // Need to skip the last spline in the s and t directions because it is always zero
    // inside the domain, and would cause A to be rank-deficient
    for (int m = 0; m < 4; ++m)
      if (startingControlPoint.i + m < numNonzeroControlPointsPerDirection)
        for (int n = 0; n < 4; ++n)
          if (startingControlPoint.j + n < numNonzeroControlPointsPerDirection)
            a(i, compute_linear_index({startingControlPoint.i + m, startingControlPoint.j + n})) =
                bsplineValsS[m] * bsplineValsT[n];

    rhs(i, 0) = inputValues[i].z;
  }

  solve_least_squares(a, rhs, work);

  for (int i = 0; i < numNonzeroControlPointsPerDirection; ++i)
    for (int j = 0; j < numNonzeroControlPointsPerDirection; ++j)
      m_controlPointVals(i, j) = rhs(compute_linear_index({i, j}), 0);

  for (int i = 0; i < m_controlPointVals.extent(0); ++i)
    m_controlPointVals(i, m_controlPointVals.extent(1) - 1) = 0;

  for (int j = 0; j < m_controlPointVals.extent(1); ++j)
    m_controlPointVals(m_controlPointVals.extent(0) - 1, j) = 0;
}

} // namespace impl
} // namespace utils
} // namespace middle_mesh
} // namespace stk
