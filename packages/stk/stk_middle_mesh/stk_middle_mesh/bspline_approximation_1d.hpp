#include "b_splines.hpp"

#include <cassert>
#include <iostream>
#include <vector>

#include "matrix.hpp"

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

class BSplineApproximation1D
{
  public:
    // num control points is really the number of control points within the domain of the
    // function.  The actual number of control points is num_control_points + 3
    BSplineApproximation1D(int numControlPoints)
      : m_numControlPointsInterior(numControlPoints)
      , m_controlPointVals(numControlPoints + 3)
    {
      assert(numControlPoints > 1);
    }

    void approximate(const std::vector<double>& xValues, const std::vector<double>& fValues);

    double eval(double x);

  private:
    void setup_local_coordinate_system(double xmin, double xmax);

    int get_left_control_point(double x);

    double compute_t_value(double x);

    // this implements this method from "Scattered Data Interpolation with Mutlilevel B-Splines", Section 3,
    // by Lee, Woldberg, and Shin.  It produces bad results, so don't use it, even though it is much
    // faster than solving a linear system

    // void computeControlPointValues(const std::vector<double>& x_values, const std::vector<double>& f_values);

    void compute_control_points_least_squares(const std::vector<double>& xValues, const std::vector<double>& fValues);

    double m_xmin = 0;
    int m_numControlPointsInterior;
    double m_deltaXControlPoints = 0;
    std::vector<double> m_controlPointVals;
    CubicBSplines m_bsplines;
};

} // namespace impl
} // namespace utils
} // namespace middle_mesh
} // namespace stk
