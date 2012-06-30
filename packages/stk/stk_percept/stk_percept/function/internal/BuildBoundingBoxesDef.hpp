#ifndef stk_percept_BuildBoundingBoxesDef_hpp
#define stk_percept_BuildBoundingBoxesDef_hpp

#include <stk_percept/function/internal/BuildBoundingBoxes.hpp>

namespace stk
{
  namespace percept
  {

    // FIXME
    //template<unsigned SpatialDim>
    //std::ostream &operator<<(std::ostream &out, const typename BuildBoundingBoxes<SpatialDim>::BoundingBox &bbox)
    std::ostream &operator<<(std::ostream &out, const BuildBoundingBoxes<3>::AABoundingBox &bbox)
    {
      out << "bbox_min = { " 
          << bbox.lower(0) << ",  "
          << bbox.lower(1) << ",  "
          << bbox.lower(2) << "}  "
          << "bbox_max = { " 
          << bbox.upper(0) << ",  "
          << bbox.upper(1) << ",  "
          << bbox.upper(2) << "}  ";
      return out;
    }

    //template<unsigned SpatialDim>
    //std::ostream &operator<<(std::ostream &out, const typename BuildBoundingBoxes<SpatialDim>::BoundingBox &bbox)
    std::ostream &operator<<(std::ostream &out, const BuildBoundingBoxes<3>::BoundingPoint &bbox)
    {
      out << "bbox_min = { " 
          << bbox.lower(0) << ",  "
          << bbox.lower(1) << ",  "
          << bbox.lower(2) << "}  "
          << "bbox_max = { " 
          << bbox.upper(0) << ",  "
          << bbox.upper(1) << ",  "
          << bbox.upper(2) << "}  ";
      return out;
    }

    std::ostream &operator<<(std::ostream &out, const BuildBoundingBoxes<2>::AABoundingBox &bbox)
    {
      out << "bbox_min = { " 
          << bbox.lower(0) << ",  "
          << bbox.lower(1) << "}  "
          << "bbox_max = { " 
          << bbox.upper(0) << ",  "
          << bbox.upper(1) << "}  ";
      return out;
    }

    //template<unsigned SpatialDim>
    //std::ostream &operator<<(std::ostream &out, const typename BuildBoundingBoxes<SpatialDim>::BoundingBox &bbox)
    std::ostream &operator<<(std::ostream &out, const BuildBoundingBoxes<2>::BoundingPoint &bbox)
    {
      out << "bbox_min = { " 
          << bbox.lower(0) << ",  "
          << bbox.lower(1) << "}  "
          << "bbox_max = { " 
          << bbox.upper(0) << ",  "
          << bbox.upper(1) << "}  ";
      return out;
    }


#if 0
    template<unsigned SpatialDim>
    std::ostream &operator<<(std::ostream &out, const typename BuildBoundingBoxes<SpatialDim>::AABoundingBox &bbox)
    {
      out << "bbox_min = { " 
          << bbox.lower(0) << ",  "
          << bbox.lower(1) << ",  "
          << bbox.lower(2) << "}  "
          << "bbox_max = { " 
          << bbox.upper(0) << ",  "
          << bbox.upper(1) << ",  "
          << bbox.upper(2) << "}  ";
      return out;
    }
#endif

#if 0
    template<typename BBox_loc>
    std::ostream &operator<<(std::ostream &out, const BBox_loc &bbox)
    {
      out << "bbox_min = { " 
          << bbox.lower(0) << ",  "
          << bbox.lower(1) << ",  "
          << bbox.lower(2) << "}  "
          << "bbox_max = { " 
          << bbox.upper(0) << ",  "
          << bbox.upper(1) << ",  "
          << bbox.upper(2) << "}  ";
      return out;
    }
#endif

  }
}
#endif
