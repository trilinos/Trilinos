#ifndef stk_search_util_stk_mesh_CreateBoundingBox_hpp
#define stk_search_util_stk_mesh_CreateBoundingBox_hpp
#include <vector>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/fem/FieldDeclarations.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>

typedef stk::search::ident::IdentProc<stk::mesh::EntityKey, unsigned> IdentProc;
typedef stk::search::box::AxisAlignedBoundingBox<IdentProc, double, 3> AxisAlignedBoundingBox3D;
typedef stk::search::box::PointBoundingBox<IdentProc, double, 3> PointBoundingBox3D;
typedef std::vector<std::pair<IdentProc, IdentProc> > IdentProcRelation;

namespace stk {

namespace mesh {
class BulkData;
}

namespace search_util {

struct Op
{
  typedef double Data;

  virtual void operator()(Data &xmin, Data &ymin, Data &zmin,
			  Data &xmax, Data &ymax, Data &zmax) const = 0;
  virtual void operator()(AxisAlignedBoundingBox3D &box) const = 0;
  virtual ~Op(){}
};

struct NoOp : public Op
{
  typedef double Data;

  NoOp(){}
  virtual ~NoOp(){}
  virtual void operator()(Data &xmin, Data &ymin, Data &zmin,
			  Data &xmax, Data &ymax, Data &zmax) const
  {}
  virtual void operator()(AxisAlignedBoundingBox3D &box) const
  {}
};

void build_axis_aligned_bbox(stk::mesh::BulkData &bulk_data, stk::mesh::EntityType type,
                             stk::mesh::VectorField *coordinates,
                             std::vector<AxisAlignedBoundingBox3D> &box_vector,
			     bool use_universal_part = false,
			     const Op &op = NoOp());


void build_centroid_bbox(stk::mesh::BulkData &bulk_data,  stk::mesh::EntityType type,
                         stk::mesh::VectorField *coordinates,
                         std::vector<PointBoundingBox3D> &box_vector,
                         bool use_universal_part = false);

/**
 * The bounding box is expanded in all directions by the amount
 * "max_d * scale + offset" where "max_d" is the maximum extent
 * of the bounding box.in any coordinate direction (xmax-xmin,
 * ymax-ymin, zmax-zmin).
 */
struct OffsetScaleOp : public Op
{
  typedef double Data;

  OffsetScaleOp(double scale, double offset)
    : m_scale(scale),
      m_offset(offset)
  {}
  virtual ~OffsetScaleOp(){}
  virtual void operator()(Data &xmin, Data &ymin, Data &zmin,
			  Data &xmax, Data &ymax, Data &zmax) const;
  virtual void operator()(AxisAlignedBoundingBox3D &box) const;

  const double          m_scale;
  const double          m_offset;
};

} // namespace search_util
} // namespace stk
#endif // stk_search_util_stk_mesh_CreateBoundingBox_hpp
