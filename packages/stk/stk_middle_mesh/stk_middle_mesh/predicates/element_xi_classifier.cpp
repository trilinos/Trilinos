#include "stk_middle_mesh/predicates/element_xi_classifier.hpp"
#include "stk_middle_mesh/predicates/intersection_common.hpp"

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

PointRecordForTriangle ElementXiClassifier::classify_triangle(const utils::Point& ptXi)
{
  double xi1 = ptXi.x, xi2 = ptXi.y, xi3 = 1 - ptXi.get_x() - ptXi.get_y();

  PointClassification type;
  int id;
  // verts
  if (std::abs(xi1 - 1) < m_eps && in_range(xi2, 0, 1, m_eps) && in_range(xi3, 0, 1, m_eps))
  {
    type = PointClassification::Vert;
    id   = 1;
  } else if (std::abs(xi2 - 1) < m_eps && in_range(xi1, 0, 1, m_eps) && in_range(xi3, 0, 1, m_eps))
  {
    type = PointClassification::Vert;
    id   = 2;
  } else if (std::abs(xi3 - 1) < m_eps && in_range(xi1, 0, 1, m_eps) && in_range(xi2, 0, 1, m_eps))
  {
    type = PointClassification::Vert;
    id   = 0;
    // edges
  } else if (std::abs(xi1) < m_eps && in_range(xi2, 0, 1, m_eps) && in_range(xi3, 0, 1, m_eps))
  {
    type = PointClassification::Edge;
    id   = 2;
  } else if (std::abs(xi2) < m_eps && in_range(xi1, 0, 1, m_eps) && in_range(xi3, 0, 1, m_eps))
  {
    type = PointClassification::Edge;
    id   = 0;

  } else if (std::abs(xi3) < m_eps && in_range(xi1, 0, 1, m_eps) && in_range(xi2, 0, 1, m_eps))
  {
    type = PointClassification::Edge;
    id   = 1;
    // interior or exterior
  } else if (in_range(xi1, 0, 1, m_eps) && in_range(xi2, 0, 1, m_eps) && in_range(xi3, 0, 1, m_eps))
  {
    type = PointClassification::Interior;
    id = 0;
  } else
  {
    type = PointClassification::Exterior;
    id = -1;
  }

  return PointRecordForTriangle(type, id, nullptr, utils::Point(xi1, xi2));
}

PointRecordForTriangle ElementXiClassifier::classify_quad(const utils::Point& ptXi)
{
  double xi1 = ptXi.x, xi2 = ptXi.y;

  PointClassification type;
  int id;
  // verts
  if (std::abs(xi1) < m_eps && std::abs(xi2) < m_eps)
  {
    type = PointClassification::Vert;
    id   = 0;
  } else if (std::abs(xi1 - 1) < m_eps && std::abs(xi2) < m_eps)
  {
    type = PointClassification::Vert;
    id   = 1;
  } else if (std::abs(xi1 - 1) < m_eps && std::abs(xi2 - 1) < m_eps)
  {
    type = PointClassification::Vert;
    id   = 2;
  } else if (std::abs(xi1) < m_eps && std::abs(xi2 - 1) < m_eps)
  {
    type = PointClassification::Vert;
    id   = 3;
    // edges
  } else if (std::abs(xi1) < m_eps && in_range(xi2, 0, 1, m_eps))
  {
    type = PointClassification::Edge;
    id   = 3;
  } else if (std::abs(xi1 - 1) < m_eps && in_range(xi2, 0, 1, m_eps))
  {
    type = PointClassification::Edge;
    id   = 1;
  } else if (std::abs(xi2) < m_eps && in_range(xi1, 0, 1, m_eps))
  {
    type = PointClassification::Edge;
    id   = 0;
  } else if (std::abs(xi2 - 1) < m_eps && in_range(xi1, 0, 1, m_eps))
  {
    type = PointClassification::Edge;
    id   = 2;
    // interior or exterior
  } else if (in_range(xi1, 0, 1, m_eps) && in_range(xi2, 0, 1, m_eps))
  {
    type = PointClassification::Interior;
    id = 0;
  } else
  {
    type = PointClassification::Exterior;
    id = -1;
  }

  return PointRecordForTriangle(type, id, nullptr, utils::Point(xi1, xi2));
}

} // namespace impl
} // namespace predicates

} // namespace middle_mesh
} // namespace stk
