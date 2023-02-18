#ifndef ELEMENT_XI_CLASSIFIER_H
#define ELEMENT_XI_CLASSIFIER_H

#include "mesh_entity.h"
#include "predicates/intersection_common.h"

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

// determines which topological entity a given point is
// on given the xi coordinates of the point
class ElementXiClassifier
{
  public:
    explicit ElementXiClassifier(double eps)
      : m_eps(eps)
    {}

    PointRecordForTriangle classify_triangle(const utils::Point& ptXi);

    PointRecordForTriangle classify_quad(const utils::Point& ptXi);

  private:
    double m_eps;
};

} // namespace impl
} // namespace predicates

} // namespace middle_mesh
} // namespace stk
#endif