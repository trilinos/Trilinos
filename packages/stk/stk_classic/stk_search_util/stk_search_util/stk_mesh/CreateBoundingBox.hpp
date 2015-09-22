/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_search_util_stk_mesh_CreateBoundingBox_hpp
#define stk_search_util_stk_mesh_CreateBoundingBox_hpp
#include <vector>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>

typedef stk_classic::mesh::Field<double, stk_classic::mesh::Cartesian> CartesianField ;
typedef stk_classic::search::ident::IdentProc<stk_classic::mesh::EntityKey, unsigned> IdentProc;
typedef stk_classic::search::box::AxisAlignedBoundingBox<IdentProc, double, 3> AxisAlignedBoundingBox3D;
typedef stk_classic::search::box::PointBoundingBox<IdentProc, double, 3> PointBoundingBox3D;
typedef std::vector<std::pair<IdentProc, IdentProc> > IdentProcRelation;

namespace stk_classic {

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

void build_axis_aligned_bbox(stk_classic::mesh::BulkData &bulk_data, stk_classic::mesh::EntityRank type,
                             CartesianField *coordinates,
                             std::vector<AxisAlignedBoundingBox3D> &box_vector,
			     bool use_universal_part = false,
			     const Op &op = NoOp());


void build_centroid_bbox(stk_classic::mesh::BulkData &bulk_data,  stk_classic::mesh::EntityRank type,
                         CartesianField *coordinates,
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
} // namespace stk_classic
#endif // stk_search_util_stk_mesh_CreateBoundingBox_hpp
