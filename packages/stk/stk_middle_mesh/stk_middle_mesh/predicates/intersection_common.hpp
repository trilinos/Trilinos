#ifndef INTERSECTION_COMMON_H
#define INTERSECTION_COMMON_H

#include <stk_middle_mesh/mesh.hpp>
#include <ostream>
#include <string>

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

class MeshRelationalDataScatter;

}
}
}
}

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

enum class IntersectionType
{
  NONE = 0,
  POINT,
  OVERLAP,
};

std::string get_enum_string(IntersectionType type);

enum class PointClassification
{
  Vert = 0,
  Edge,
  Interior,
  Exterior,
};

inline std::ostream& operator<<(std::ostream& os, PointClassification type)
{
  switch (type)
  {
    case PointClassification::Vert: {
      os << "Vert";
      break;
    }
    case PointClassification::Edge: {
      os << "Edge";
      break;
    }
    case PointClassification::Interior: {
      os << "Interior";
      break;
    }
    case PointClassification::Exterior: {
      os << "Exterior";
      break;
    }
  }

  return os;
}

// true if val in [vL - eps, vR + eps]
inline bool in_range(const double val, const double vL, const double vR, double eps)
{
  return val >= vL - eps && val <= vR + eps;
}

// true if val in [vL + eps, vR - eps]
inline bool in_range_e(const double val, const double vL, const double vR, double eps)
{
  return val >= vL + eps && val <= vR - eps;
}

class QuadToTriangles;
class PointClassifierNormalWrapper;
class TriangleCoordUtils;
class PointClassifier;
class PointClassifierForTriangle;

struct PointRecordForTriangle
{
  public:
    explicit PointRecordForTriangle(PointClassification type_ = PointClassification::Exterior, const int id_ = -1,
                                    mesh::MeshEntityPtr el_ = nullptr, const utils::Point& ptXi_ = utils::Point(0, 0))
      : type(type_)
      , id(id_)
      , el(el_)
      , m_ptXi(ptXi_)
    {}

    PointClassification type;
    int id;
    mesh::MeshEntityPtr el;

  private:
    utils::Point m_ptXi;

    friend PointClassifierForTriangle;
    friend TriangleCoordUtils;
    friend QuadToTriangles;
    friend class stk::middle_mesh::nonconformal4::impl::MeshRelationalDataScatter;

    friend std::ostream& operator<<(std::ostream& os, const PointRecordForTriangle& record);
};

inline std::ostream& operator<<(std::ostream& os, const PointRecordForTriangle& record)
{
  os << "PointRecord: type = " << record.type << ", pt_xi = " << record.m_ptXi;
  if (record.type == PointClassification::Vert || record.type == PointClassification::Edge)
    os << ", id = " << record.id;

  return os;
}

mesh::MeshEntityPtr get_entity(const PointRecordForTriangle& record);


struct PointRecord
{
  public:
    explicit PointRecord(PointClassification type_ = PointClassification::Exterior, const int id_ = -1,
                         mesh::MeshEntityPtr el_ = nullptr, const PointRecordForTriangle& r1 = PointRecordForTriangle(),
                         const PointRecordForTriangle& r2 = PointRecordForTriangle())
      : type(type_)
      , id(id_)
      , el(el_)
      , m_r1(r1)
      , m_r2(r2)
    {}

    PointClassification type;
    int id;
    mesh::MeshEntityPtr el;

  private:
    PointRecordForTriangle m_r1;
    PointRecordForTriangle m_r2;

    friend PointClassifier; // TODO: maybe find a way to not use friends
    friend QuadToTriangles;
    friend PointClassifierNormalWrapper;
    friend class stk::middle_mesh::nonconformal4::impl::MeshRelationalDataScatter;

    friend std::ostream& operator<<(std::ostream& os, const PointRecord& record);
};

mesh::MeshEntityPtr get_entity(const PointRecord& record);

int get_entity_id(mesh::MeshEntityPtr el, mesh::MeshEntityPtr entity);

inline std::ostream& operator<<(std::ostream& os, const PointRecord& record)
{
  os << "PointRecord: type = " << record.type;
  if (record.type == PointClassification::Vert || record.type == PointClassification::Edge)
    os << ", id = " << record.id;

  return os;
}

} // namespace impl
} // namespace predicates
} // namespace middle_mesh
} // namespace stk
#endif
