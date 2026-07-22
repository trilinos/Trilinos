#ifndef PREDICATES_POINT_CLASSIFIER_OPTS
#define PREDICATES_POINT_CLASSIFIER_OPTS

namespace stk {
namespace middle_mesh {
namespace impl {

struct PointClassifierNormalInterpolationTolerances
{
    PointClassifierNormalInterpolationTolerances(double eps)
      : pointClassifierTol(1000 * eps)
      , newtonTol(1e-13)
      , newtonItermax(100000)
    {}

    double pointClassifierTol; // tolerance for determining if a point is classified on
                               // a vertex or edge.  This tolerance is applied to the
                               // barycentric coords, which have range [0, 1] inside the element
    double newtonTol;          // tolerance for Newtons method convergence (ie. norm of function
                               // is less than newtonTol)
    double newtonItermax;      // maximum number of iterations for Newtons method
};

struct PointClassifierNormalWrapperTolerances
{
    explicit PointClassifierNormalWrapperTolerances(double eps = 1e-13, [[maybe_unused]] double smallSlopeTol = 1e-13)
      : quadBreakupTolerance(eps)
      , normalInterpolationTolerances(eps)
    {}

    double quadBreakupTolerance; // when breaking up a quad into 2 triangles, use this tolerance to
                                 // determine which of the triangles a given point is inside
    PointClassifierNormalInterpolationTolerances normalInterpolationTolerances;
};

} // namespace impl
} // namespace middle_mesh
} // namespace stk
#endif