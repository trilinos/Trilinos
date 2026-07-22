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
#include <stk_math/StkVector.hpp>

#include <algorithm>
#include <vector>
#include <functional>

#include <stk_util/util/ReportHandler.hpp>
#include <Akri_BoundingBoxDistance.hpp>

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
  Real boxLowerBnd2( const VecType & x ) const { return min_possible_closest_squared_distance(m_bbox, x); }
  Real boxUpperBnd2( const VecType & x ) const { return max_possible_closest_squared_distance(m_bbox, x); }
  Real boxSize2() const { return m_bbox.size_squared(); }

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
  typedef BoundingBoxType::Real Real;
  typedef BoundingBoxType::VecType VecType;
  typedef std::pair<EntityType, VecType> EntityCentroid;

  SearchTree( std::vector<EntityCentroid> & entities ); // vector is modified by sorting
  SearchTree( const std::vector<EntityType> & entities,
      const std::function<stk::math::Vector3d(const EntityType)> & get_centroid,
      const std::function<void(const EntityType, BoundingBox & bbox)> & insert_into_bbox);
  SearchTree() = delete;

  // All entities within the max_search_radius that may be the closest to the given point will be returned.
  // If max_search_radius = 0, then all entities that may be the closest to the given point will be returned.
  void find_closest_entities( const stk::math::Vector3d& search_point, const std::function<double(const EntityType, const stk::math::Vector3d &)> & get_entity_distance_squared_fn, std::vector<EntityType>& return_vec, const Real max_search_radius = 0.0 );
  void find_closest_entities( const stk::math::Vector3d& search_point, std::vector<EntityType>& return_vec, const Real max_search_radius = 0.0 );
  
  void recursively_get_entities_within_distance_bound( const stk::math::Vector3d& search_point,
     const double maxSearchSqrDist,
     std::vector<EntityType> & nearest)
  {
    // Algorithm requires tol > std::numeric_limits<double>::epsilon(), here just pick multiplier
    const Real tol = 100.*std::numeric_limits<Real>::epsilon();
    const Real maxDist2 = maxSearchSqrDist;
    const VecType search_pt(search_point.data());
    recursively_get_entities_within_distance_bound( search_pt, maxDist2, tol, nearest, 0 );
  }

  void get_intersecting_entities( const BoundingBoxType& bbox, std::vector<EntityType>& return_vec)
  {
    return_vec.clear();
    if (!empty())
    {
      recursively_get_intersecting_entities( bbox, return_vec, 0);
    }
  }

  bool empty() const { return my_nodes.empty(); }
  size_t size() const { return my_nodes.size(); }
  size_t storage_size() const { return my_nodes.size() * sizeof(SearchTreeNode<EntityType>); }

private:
  size_t build( typename std::vector<EntityCentroid>::iterator entity_begin,
      typename std::vector<EntityCentroid>::iterator entity_end,
      const std::function<void(const EntityType, BoundingBox & bbox)> & insert_into_bbox,
      const size_t index = 0  );
  std::pair<Real,Real> compute_tree_bounds(const VecType & x, const std::function<Real(const EntityType)> & get_entity_distance_squared_fn, const Real max_search_radius2);
  double compute_node_bounds_and_update_tree_bounds( const VecType & x,
     const std::function<Real(const EntityType)> & get_entity_distance_squared_fn,
     Real & lowerBnd2,
     Real & upperBnd2,
     const size_t index = 0 );
  void recursively_update_tree_bounds( const VecType & x,
     const std::function<Real(const EntityType)> & get_entity_distance_squared_fn,
     Real & lowerBnd2,
     Real & upperBnd2,
     const size_t index = 0 );
  void recursively_get_entities_within_distance_bound( const VecType & x,
     const Real maxDist2,
     const Real tol,
     std::vector<EntityType> & nearest,
     const size_t index);
  void recursively_get_intersecting_entities( const BoundingBoxType& bbox, std::vector<EntityType>& return_vec, const size_t index );
  void add_descendants( std::vector<EntityType>& return_vec, const size_t index );
private:
  std::vector<SearchTreeNode<EntityType>> my_nodes;
};

template <class ENTITY>
inline std::pair<typename SearchTree<ENTITY>::Real,typename SearchTree<ENTITY>::Real>
SearchTree<ENTITY>::compute_tree_bounds(const VecType & x,
    const std::function<Real(const EntityType)> & get_entity_distance_squared_fn,
    const Real max_search_radius2)
{
  // first find upper bound for distance
  // the very first estimate for this is the upper bound for the base of the tree
  Real upperBnd2 = my_nodes[0].boxUpperBnd2( x );

  if ( upperBnd2 > max_search_radius2 && max_search_radius2 > 0.0 )
  {
    upperBnd2 = max_search_radius2;
  }

  // now descend relevant parts of tree updating the bounds as we go
  Real lowerBnd2 = upperBnd2;
  recursively_update_tree_bounds(x, get_entity_distance_squared_fn, lowerBnd2, upperBnd2);

  return {lowerBnd2, upperBnd2};
}

template <class ENTITY>
inline double
SearchTree<ENTITY>::compute_node_bounds_and_update_tree_bounds( const VecType & search_point,
    const std::function<Real(const EntityType)> & compute_leaf_entity_distance_squared_fn,
    Real & lowerBnd2,
    Real & upperBnd2,
    const size_t index)
{
  const SearchTreeNode<EntityType> & node = my_nodes[index];
  const Real node_lowerBnd2 = node.boxLowerBnd2( search_point );
  if ( node_lowerBnd2 < upperBnd2 )
  {
    // Make sure that, even if the precision of the leaf query is different than that of the bounding box, we
    // don't get the impossible state where upperBnd2 is less than lowerBnd2.
    const Real node_upperBnd2 = (!node.have_children() && compute_leaf_entity_distance_squared_fn) ?
        std::max(node_lowerBnd2, compute_leaf_entity_distance_squared_fn(node.get_leaf_entity())) :
        node.boxUpperBnd2(search_point);
    if ( node_upperBnd2 < upperBnd2 )
      upperBnd2 = node_upperBnd2;
    if (!node.have_children() && node_lowerBnd2 < lowerBnd2)
      lowerBnd2 = node_lowerBnd2;
  }
  return node_lowerBnd2;
}

template <class ENTITY>
SearchTree<ENTITY>::SearchTree( std::vector<EntityCentroid> & entities )
{
  if (!entities.empty())
  {
    my_nodes.resize(2*entities.size() - 1);
    const size_t tree_size = build( entities.begin(), entities.end() );
    STK_ThrowRequire(tree_size == my_nodes.size());
  }
}

template <class ENTITY>
SearchTree<ENTITY>::SearchTree( const std::vector<EntityType> & entities,
    const std::function<stk::math::Vector3d(const EntityType)> & get_centroid,
    const std::function<void(const EntityType, BoundingBox & bbox)> & insert_into_bbox)
{
  if (!entities.empty())
  {
    std::vector<EntityCentroid> bbox_entities;
    bbox_entities.reserve(entities.size());
    for (auto && entity : entities)
    {
      const stk::math::Vector3d centroid = get_centroid(entity);
      bbox_entities.emplace_back(entity, VecType(centroid[0], centroid[1], centroid[2]));
    }
    my_nodes.resize(2*entities.size() - 1);
    const size_t tree_size = build( bbox_entities.begin(), bbox_entities.end(), insert_into_bbox );
    STK_ThrowRequire(tree_size == my_nodes.size());
  }
}

template <class ENTITY>
void
SearchTree<ENTITY>::find_closest_entities( const stk::math::Vector3d & searchPoint, std::vector<EntityType> & closestEntities, const Real maxSearchRadius )
{
  const std::function<double(const EntityType, const stk::math::Vector3d &)> emptyFunction;
  find_closest_entities(searchPoint, emptyFunction, closestEntities, maxSearchRadius);
}

template <class ENTITY>
void
SearchTree<ENTITY>::find_closest_entities( const stk::math::Vector3d & search_point, const std::function<double(const EntityType, const stk::math::Vector3d &)> & get_entity_distance_squared_fn, std::vector<EntityType> & return_vec, const Real max_search_radius )
{
  return_vec.clear();

  // if there are no entities at all return empty list
  if ( empty() ) return;

  const VecType search_pt(search_point.data());
  const Real max_search_radius2 = max_search_radius * max_search_radius;

  std::function<Real(const EntityType)> get_real_entity_distance_squared_fn;
  if (get_entity_distance_squared_fn)
  {
    get_real_entity_distance_squared_fn = [&search_point, &get_entity_distance_squared_fn](const EntityType entity)->Real { return get_entity_distance_squared_fn(entity, search_point); };
  }

  const auto & [lowerBnd2, upperBnd2] = compute_tree_bounds(search_pt, get_real_entity_distance_squared_fn, max_search_radius2);

  // exit now if no entities within narrow_band
  if ( lowerBnd2 > max_search_radius2 && max_search_radius2 > 0.0 )
  {
    return;
  }

  // with the upper bound known, find all entities within the tolerance
  recursively_get_entities_within_distance_bound( search_point, upperBnd2, return_vec );

  return;
}

template <class ENTITY>
void
SearchTree<ENTITY>::recursively_get_entities_within_distance_bound(
    const VecType & search_point,
    const Real maxDist2,
    const Real tol,
    std::vector<EntityType> & closeEntities,
    const size_t index)
{
  const SearchTreeNode<EntityType> & node = my_nodes[index];
  if ( node.have_children() )
  {
    const size_t Lindex = index + 1;
    const Real L_lowerBnd2 = my_nodes[Lindex].boxLowerBnd2( search_point );
    const Real L_eps = tol*(maxDist2+my_nodes[Lindex].boxSize2());
    if ( L_lowerBnd2 < maxDist2+L_eps )
      recursively_get_entities_within_distance_bound( search_point, maxDist2, tol, closeEntities, Lindex );

    const size_t Rindex = node.get_right_index();
    const Real R_lowerBnd2 = my_nodes[Rindex].boxLowerBnd2( search_point );
    const Real R_eps = tol*(maxDist2+my_nodes[Rindex].boxSize2());
    if ( R_lowerBnd2 < maxDist2+R_eps )
      recursively_get_entities_within_distance_bound( search_point, maxDist2, tol, closeEntities, Rindex );
  }
  else
  {
    EntityType entity = node.get_leaf_entity();
    closeEntities.push_back( entity );
  }
}

template <class ENTITY>
void
SearchTree<ENTITY>::recursively_update_tree_bounds( const VecType & search_point,
    const std::function<Real(const EntityType)> & get_entity_distance_squared_fn,
    Real & lowerBnd2,
    Real & upperBnd2,
    const size_t index)
{
  const SearchTreeNode<EntityType> & node = my_nodes[index];
  if (!node.have_children())
  {
    if (0 == index)
    {
      // Special case of a single node.  Otherwise, this node will already have been handled by parent.
      compute_node_bounds_and_update_tree_bounds(search_point, get_entity_distance_squared_fn,lowerBnd2, upperBnd2, index);
    }
    return;
  }

  const size_t Lindex = index + 1;
  const Real L_lowerBnd2 = compute_node_bounds_and_update_tree_bounds(search_point, get_entity_distance_squared_fn, lowerBnd2, upperBnd2, Lindex);

  const size_t Rindex = node.get_right_index();
  const Real R_lowerBnd2 = compute_node_bounds_and_update_tree_bounds(search_point, get_entity_distance_squared_fn, lowerBnd2, upperBnd2, Rindex);

  if ( L_lowerBnd2 < R_lowerBnd2 )
  {
    if ( L_lowerBnd2 < upperBnd2 )
      recursively_update_tree_bounds( search_point, get_entity_distance_squared_fn, lowerBnd2, upperBnd2, Lindex );
    if ( R_lowerBnd2 < upperBnd2 )
      recursively_update_tree_bounds( search_point, get_entity_distance_squared_fn, lowerBnd2, upperBnd2, Rindex );
  }
  else
  {
    if ( R_lowerBnd2 < upperBnd2 )
      recursively_update_tree_bounds( search_point, get_entity_distance_squared_fn, lowerBnd2, upperBnd2, Rindex );
    if ( L_lowerBnd2 < upperBnd2 )
      recursively_update_tree_bounds( search_point, get_entity_distance_squared_fn, lowerBnd2, upperBnd2, Lindex );
  }
}

template <class ENTITY>
size_t
SearchTree<ENTITY>::build( typename std::vector<EntityCentroid>::iterator entity_begin,
    typename std::vector<EntityCentroid>::iterator entity_end,
    const std::function<void(const EntityType, BoundingBox & bbox)> & insert_into_bbox,
    const size_t index )
{
  const size_t entity_size = std::distance(entity_begin, entity_end);

  // determine bounding box
  BoundingBox & node_bbox = my_nodes[index].bounding_box();
  for ( auto entityIt = entity_begin; entityIt != entity_end; ++entityIt )
    insert_into_bbox(entityIt->first, node_bbox);

  if ( entity_size > 1 )
  {
    BoundingBox centroidBbox;
    for ( auto entityIt = entity_begin; entityIt != entity_end; ++entityIt )
      centroidBbox.accommodate(entityIt->second);

    // determine axis to split tree on
    const unsigned max_spread_dim = centroidBbox.max_span_direction();

    auto compare = [max_spread_dim](const EntityCentroid & lhs, const EntityCentroid & rhs) { return lhs.second[max_spread_dim] < rhs.second[max_spread_dim]; };
    const size_t mid_pt = entity_size / 2;

    // This works around a bug in some versions of nth_element
#ifdef CRAY_LWK
      std::partial_sort( entity_begin, entity_begin + mid_pt, entity_end, compare );
#else
      std::nth_element( entity_begin, entity_begin + mid_pt, entity_end, compare );
#endif

    const size_t Lindex = index + 1;
    const size_t Rindex = build( entity_begin, entity_begin + mid_pt, insert_into_bbox, Lindex );
    my_nodes[index].set_right_index(Rindex);
    const size_t next_index = build( entity_begin + mid_pt, entity_end, insert_into_bbox, Rindex );
    return next_index;
  }
  else
  {
    STK_ThrowAssert(1 == entity_size);
    my_nodes[index].set_leaf_entity(entity_begin->first);
    return index + 1;
  }
}

template <class ENTITY>
void
SearchTree<ENTITY>::recursively_get_intersecting_entities( const BoundingBoxType& bbox, std::vector<EntityType>& return_vec, const size_t index )
{
  const SearchTreeNode<EntityType> & node = my_nodes[index];
  if (!node.bounding_box().intersects(bbox)) return;

  if ( node.have_children() )
  {
    if ( bbox.contains(node.bounding_box()) )
    {
      add_descendants(return_vec, index);
    }
    else
    {
      const size_t Lindex = index + 1;
      recursively_get_intersecting_entities( bbox, return_vec, Lindex );

      const size_t Rindex = node.get_right_index();
      recursively_get_intersecting_entities( bbox, return_vec, Rindex );
    }
  }
  else
  {
    EntityType entity = node.get_leaf_entity();
    return_vec.push_back( entity );
  }
}

template <class ENTITY>
void
SearchTree<ENTITY>::add_descendants( std::vector<EntityType>& return_vec, const size_t index )
{
  const SearchTreeNode<EntityType> & node = my_nodes[index];

  if ( node.have_children() )
  {
    const size_t Lindex = index + 1;
    add_descendants( return_vec, Lindex );

    const size_t Rindex = node.get_right_index();
    add_descendants( return_vec, Rindex );
  }
  else
  {
    EntityType entity = node.get_leaf_entity();
    return_vec.push_back( entity );
  }
}

} // namespace krino

#endif // Akri_SearchTree_h
