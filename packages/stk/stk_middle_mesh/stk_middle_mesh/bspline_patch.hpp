#ifndef BSPLINE_PATCH_H
#define BSPLINE_PATCH_H

#include "b_splines.hpp"
#include "matrix.hpp"
#include "projection.hpp"
#include "surface_patch.hpp"

#include <array>
#include <cassert>
#include <vector>

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

// uses cubic B Splines to approximate z = f(x, y)
class BSplinePatch : public SurfacePatch
{
  public:
    // Changing the number of control points per direction is for expert users only!
    // The problem that can occur is there can be no points in the support of one or
    // more of the splines, leading to a singular linear system when we solve for the
    // control point values.
    explicit BSplinePatch(int numControlPointsPerDirectionInterior = 2)
      : m_numControlPointsPerDirectionInterior(numControlPointsPerDirectionInterior)
      , m_controlPointVals(numControlPointsPerDirectionInterior + 3, numControlPointsPerDirectionInterior + 3)
    {}

    using SurfacePatch::construct_patch;

    // constructs patch such that pt0 is included in the patch domain, as well as
    // the points in pts
    void construct_patch(const std::vector<Point>& pts, const Point& pt0) override;

    Point eval_point(const double x, const double y) override;

    void eval_deriv([[maybe_unused]] const double x, [[maybe_unused]] const double y, [[maybe_unused]] double derivs[2]) override
    {
      throw std::runtime_error("BSplinePatch derivative not implemented yet");
    }

  private:
    struct Indices
    {
        int i;
        int j;
    };

    struct LocalCoords
    {
        double s;
        double t;
    };

    struct BoundingBox
    {
        double xmin;
        double xmax;
        double ymin;
        double ymax;
    };

    BoundingBox compute_bounding_box(const std::vector<Point>& pts, const Point& pt0);

    void setup_local_coordinate_system(const BoundingBox& box);

    Indices get_corner_control_point(double x, double y);

    LocalCoords compute_st_values(double x, double y);

    int compute_linear_index(const Indices& idxs);

    void compute_control_points_least_squares(const std::vector<Point>& inputValues);

    int m_numControlPointsPerDirectionInterior;
    double m_deltaXControlPoints = 0;
    double m_deltaYControlPoints = 0;
    double m_xmin                = 0;
    double m_ymin                = 0;
    double m_xInteriorMidpoint   = 0;
    double m_yInteriorMidpoint   = 0;
    double m_eps                 = 1e-13;
    Matrix<double> m_controlPointVals;
    CubicBSplines m_bsplines;
};
} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
#endif