#ifndef PROJECTION_H
#define PROJECTION_H

#include <cmath>
#include <iostream>

#include "point.hpp"

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

class Projection
{
  public:
    Projection(const Point& pt1, const Point& pt2, const Point& pt3);

    Point project_to_plane(const Point& pt) const;

    Point project_to_plane_rev(const Point& pt, const Point& ptBar) const;

    // compute coordinates on basis [q1, q2]
    // the point *must* lie in the plane
    // The output point has z coordinate = 0
    Point compute_plane_coords(const Point& pt) const;

    Point compute_plane_coords_rev(const Point& pt, const Point& ptBar) const;

    // combines the above into one function call
    Point project_plane_coords(const Point& pt) const;

    Point project_plane_coords_rev(const Point& pt, const Point& pt2Bar) const;

  private:
    // orthonormal basis for plane
    Point m_q1;
    Point m_q2;
    // scaling for each basis vector
    double m_d1;
    double m_d2;
    Point m_v0;   // a point in the plane
    Point m_norm; // normalized unit vector

    friend std::ostream& operator<<(std::ostream& os, const Projection& p);
};

std::ostream& operator<<(std::ostream& os, const Projection& p);

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
#endif
