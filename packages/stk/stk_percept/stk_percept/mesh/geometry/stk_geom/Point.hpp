#ifndef Point_hpp
#define Point_hpp

#include <opennurbs.h>
#include <vector>

namespace stk {
  namespace geom {

    typedef ON_2dVector Vector2D;
    typedef std::vector<Vector2D> Vectors2D;

    typedef ON_3dVector Vector3D;
    typedef ON_3dPoint Point3D;
    typedef std::vector<Vector3D> Vectors3D;
    typedef std::vector<Point3D> Points3D;

  }
}
#endif
