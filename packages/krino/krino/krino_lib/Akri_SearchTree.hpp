// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_SearchTree_h
#define Akri_SearchTree_h

#include <Akri_BoundingBox.hpp>

#include <vector>
#include <functional>

namespace krino {

template <class ENTITY>
class SearchTreeNode {
public:
  typedef ENTITY EntityType;
  typedef BoundingBox BoundingBoxType;
  typedef BoundingBoxType::Real Real;
  typedef BoundingBoxType::VecType VecType;

  SearchTreeNode() : m_entity(0), m_right_index(0) {}

  const BoundingBoxType & bounding_box() const { return m_bbox; }
  BoundingBoxType & bounding_box() { return m_bbox; }
  Real boxLowerBnd2( const VecType & x ) const { return ( m_bbox.SqrDistLowerBnd(x) ); }
  Real boxUpperBnd2( const VecType & x ) const { return ( m_bbox.SqrDistUpperBnd(x) ); }
  Real boxSize2() const { return m_bbox.SqrSize(); }

  bool have_children() const { return (0 != m_right_index); }
  void set_right_index(size_t right_index) { m_right_index = right_index; }
  size_t get_right_index() const { return m_right_index; }
  void set_leaf_entity(EntityType leaf_entity) { m_entity = leaf_entity; }
  EntityType get_leaf_entity() const { return m_entity; }

private:
  BoundingBoxType m_bbox;
  EntityType m_entity;
  size_t m_right_index;
};

template <class ENTITY>
class SearchTree {
public:
  typedef ENTITY EntityType;
  typedef BoundingBox BoundingBoxType;
  typedef std::pair<BoundingBoxType, EntityType> BoundingBoxEntity;
  typedef BoundingBoxType::Real Real;
  typedef BoundingBoxType::VecType VecType;

  SearchTree( std::vector<BoundingBoxEntity> & entities ); // vector is modified by sorting
  SearchTree( const std::vector<EntityType> & entities, const std::function<const BoundingBoxType&(const EntityType &)> & get_bbox );
  SearchTree() = delete;

  // All entities within the max_search_radius that may be the closest to the given point will be returned.
  // If max_search_radius = 0, then all entities that may be the closest to the given point will be returned.
  void find_closest_entities( const Vector3d& search_point, std::vector<EntityType>& return_vec, const Real max_search_radius = 0.0 );
  
  void get_intersecting_entities( const BoundingBoxType& bbox, std::vector<EntityType>& return_vec)
  {
    return_vec.clear();
    if (!empty())
    {
      get_intersecting_entities( bbox, return_vec, 0);
    }
  }

  bool empty() const { return my_nodes.empty(); }
  size_t storage_size() const { return my_nodes.size() * sizeof(SearchTreeNode<EntityType>); }

private:
  size_t build( typename std::vector<BoundingBoxEntity>::iterator entity_begin, typename std::vector<BoundingBoxEntity>::iterator entity_end, const size_t index = 0  );
  double compute_and_update_bounds( const VecType & x,
     Real & lowerBnd2,
     Real & upperBnd2,
     const size_t index = 0 );
  void treeBounds( const VecType & x,
     Real & lowerBnd2,
     Real & upperBnd2,
     const size_t index = 0 );
  void find_closest_entities( const VecType & x,
     const Real maxDist2,
     const Real tol,
     std::vector<EntityType> & nearest,
     const size_t index = 0 );
  void get_intersecting_entities( const BoundingBoxType& bbox, std::vector<EntityType>& return_vec, const size_t index );
  void add_descendants( std::vector<EntityType>& return_vec, const size_t index );
private:
  std::vector<SearchTreeNode<EntityType>> my_nodes;
};

template <class ENTITY>
inline double
SearchTree<ENTITY>::compute_and_update_bounds( const VecType & search_point,
    Real & lowerBnd2,
    Real & upperBnd2,
    const size_t index)
{ /* %TRACE% */  /* %TRACE% */
  const SearchTreeNode<EntityType> & node = my_nodes[index];
  const Real node_lowerBnd2 = node.boxLowerBnd2( search_point );
  if ( node_lowerBnd2 < upperBnd2 )
  {
    const Real node_upperBnd2 = node.boxUpperBnd2( search_point );
    if ( node_upperBnd2 < upperBnd2 )
      upperBnd2 = node_upperBnd2;
    if (!node.have_children() && node_lowerBnd2 < lowerBnd2)
      lowerBnd2 = node_lowerBnd2;
  }
  return node_lowerBnd2;
}

} // namespace krino

#endif // Akri_SearchTree_h
