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
#include <tuple>
#include <type_traits>

#include <stk_util/util/ReportHandler.hpp>
#include <Akri_BoundingBoxDistance.hpp>

namespace krino {

// Tag passed to the closest-entity search when no per-entity (leaf) distance
// refinement is available. The bounding-box-only bounds are selected at compile
// time via if constexpr, so the no-callback query path carries no indirection
// and no runtime "is a function set?" test.
struct SearchTreeNoLeafDistance {};

template <class ENTITY>
class SearchTreeNode {
public:
  typedef ENTITY EntityType;
  typedef BoundingBox BoundingBoxType;
  typedef BoundingBoxType::Real Real;
  typedef BoundingBoxType::VecType VecType;

  SearchTreeNode() : myEntity{}, myRightIndex(0) {}

  const BoundingBoxType & bounding_box() const { return myBbox; }
  BoundingBoxType & bounding_box() { return myBbox; }
  Real boxLowerBnd2( const VecType & x ) const { return min_possible_closest_squared_distance(myBbox, x); }
  Real boxUpperBnd2( const VecType & x ) const { return max_possible_closest_squared_distance(myBbox, x); }
  Real boxSize2() const { return myBbox.size_squared(); }

  bool have_children() const { return (0 != myRightIndex); }
  void set_right_index(size_t right_index) { myRightIndex = right_index; }
  size_t get_right_index() const { return myRightIndex; }
  void set_leaf_entity(EntityType leaf_entity) { myEntity = leaf_entity; }
  EntityType get_leaf_entity() const { return myEntity; }

private:
  BoundingBoxType myBbox;
  EntityType myEntity;
  // Index of the right child within SearchTree::my_nodes; the left child is
  // always at (this_index + 1). Zero means "no children" (leaf): the root lives
  // at index 0 and can never be a right child, so 0 is an unambiguous sentinel.
  size_t myRightIndex;
};

template <class ENTITY>
class SearchTree {
public:
  typedef ENTITY EntityType;
  typedef BoundingBox BoundingBoxType;
  typedef BoundingBoxType::Real Real;
  typedef BoundingBoxType::VecType VecType;
  typedef std::pair<EntityType, VecType> EntityCentroid;

  // get_centroid and insert_into_bbox are invoked only during construction; they
  // are template parameters so the (statically known) callables inline into the
  // build with no std::function type erasure. Any callable with the signatures
  //   stk::math::Vector3d(const EntityType)
  //   void(const EntityType, BoundingBox & bbox)
  // is accepted.
  template <class GetCentroid, class InsertIntoBbox>
  SearchTree( const std::vector<EntityType> & entities,
      const GetCentroid & get_centroid,
      const InsertIntoBbox & insert_into_bbox);
  SearchTree() = delete;

  // All entities within the max_search_radius that may be the closest to the given point will be returned.
  // If max_search_radius = 0, then all entities that may be the closest to the given point will be returned.
  //
  // get_entity_distance_squared_fn is any callable double(const EntityType, const stk::math::Vector3d &)
  // giving the exact squared distance from a leaf entity to the search point; it
  // is a template parameter so it inlines into the hot query path. The overload
  // without it uses bounding-box-only bounds (SearchTreeNoLeafDistance).
  template <class GetEntityDistanceSquaredFn>
  void find_closest_entities( const stk::math::Vector3d& search_point, const GetEntityDistanceSquaredFn & get_entity_distance_squared_fn, std::vector<EntityType>& return_vec, const Real max_search_radius = 0.0 ) const;
  void find_closest_entities( const stk::math::Vector3d& search_point, std::vector<EntityType>& return_vec, const Real max_search_radius = 0.0 ) const
  {
    find_closest_entities( search_point, SearchTreeNoLeafDistance{}, return_vec, max_search_radius );
  }

  void recursively_get_entities_within_distance_bound( const stk::math::Vector3d& search_point,
     const double maxSearchSqrDist,
     std::vector<EntityType> & nearest) const
  {
    // Algorithm requires tol > std::numeric_limits<double>::epsilon(), here just pick multiplier
    const Real tol = 100.*std::numeric_limits<Real>::epsilon();
    const Real maxDist2 = maxSearchSqrDist;
    const VecType search_pt(search_point.data());
    recursively_get_entities_within_distance_bound( search_pt, maxDist2, tol, nearest, 0 );
  }

  void get_intersecting_entities( const BoundingBoxType& bbox, std::vector<EntityType>& return_vec) const
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
  template <class InsertIntoBbox>
  size_t build( typename std::vector<EntityCentroid>::iterator entity_begin,
      typename std::vector<EntityCentroid>::iterator entity_end,
      const InsertIntoBbox & insert_into_bbox,
      const size_t index = 0  );
  // get_leaf_entity_distance_squared_fn is either a callable Real(const EntityType)
  // or SearchTreeNoLeafDistance; the leaf-refinement branch is selected with
  // if constexpr so the box-only path emits no dead call.
  template <class GetLeafEntityDistanceSquaredFn>
  std::pair<Real,Real> compute_tree_bounds(const VecType & x, const GetLeafEntityDistanceSquaredFn & get_leaf_entity_distance_squared_fn, const Real max_search_radius2) const;
  template <class GetLeafEntityDistanceSquaredFn>
  double compute_node_bounds_and_update_tree_bounds( const VecType & x,
     const GetLeafEntityDistanceSquaredFn & get_leaf_entity_distance_squared_fn,
     Real & lowerBnd2,
     Real & upperBnd2,
     const size_t index = 0 ) const;
  template <class GetLeafEntityDistanceSquaredFn>
  void recursively_update_tree_bounds( const VecType & x,
     const GetLeafEntityDistanceSquaredFn & get_leaf_entity_distance_squared_fn,
     Real & lowerBnd2,
     Real & upperBnd2,
     const size_t index = 0 ) const;
  void recursively_get_entities_within_distance_bound( const VecType & x,
     const Real maxDist2,
     const Real tol,
     std::vector<EntityType> & nearest,
     const size_t index) const;
  void recursively_get_intersecting_entities( const BoundingBoxType& bbox, std::vector<EntityType>& return_vec, const size_t index ) const;
  void add_descendants( std::vector<EntityType>& return_vec, const size_t index ) const;
private:
  std::vector<SearchTreeNode<EntityType>> my_nodes;
};

template <class ENTITY>
template <class GetLeafEntityDistanceSquaredFn>
inline std::pair<typename SearchTree<ENTITY>::Real,typename SearchTree<ENTITY>::Real>
SearchTree<ENTITY>::compute_tree_bounds(const VecType & x,
    const GetLeafEntityDistanceSquaredFn & get_leaf_entity_distance_squared_fn,
    const Real max_search_radius2) const
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
  recursively_update_tree_bounds(x, get_leaf_entity_distance_squared_fn, lowerBnd2, upperBnd2);

  return {lowerBnd2, upperBnd2};
}

template <class ENTITY>
template <class GetLeafEntityDistanceSquaredFn>
inline double
SearchTree<ENTITY>::compute_node_bounds_and_update_tree_bounds( const VecType & search_point,
    const GetLeafEntityDistanceSquaredFn & compute_leaf_entity_distance_squared_fn,
    Real & lowerBnd2,
    Real & upperBnd2,
    const size_t index) const
{
  const SearchTreeNode<EntityType> & node = my_nodes[index];
  const Real node_lowerBnd2 = node.boxLowerBnd2( search_point );
  if ( node_lowerBnd2 < upperBnd2 )
  {
    // Make sure that, even if the precision of the leaf query is different than that of the bounding box, we
    // don't get the impossible state where upperBnd2 is less than lowerBnd2.
    Real node_upperBnd2;
    if constexpr (std::is_same_v<GetLeafEntityDistanceSquaredFn, SearchTreeNoLeafDistance>)
    {
      (void) compute_leaf_entity_distance_squared_fn;
      node_upperBnd2 = node.boxUpperBnd2(search_point);
    }
    else
    {
      node_upperBnd2 = !node.have_children() ?
          std::max(node_lowerBnd2, compute_leaf_entity_distance_squared_fn(node.get_leaf_entity())) :
          node.boxUpperBnd2(search_point);
    }
    if ( node_upperBnd2 < upperBnd2 )
      upperBnd2 = node_upperBnd2;
    if (!node.have_children() && node_lowerBnd2 < lowerBnd2)
      lowerBnd2 = node_lowerBnd2;
  }
  return node_lowerBnd2;
}

template <class ENTITY>
template <class GetCentroid, class InsertIntoBbox>
SearchTree<ENTITY>::SearchTree( const std::vector<EntityType> & entities,
    const GetCentroid & get_centroid,
    const InsertIntoBbox & insert_into_bbox)
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
template <class GetEntityDistanceSquaredFn>
void
SearchTree<ENTITY>::find_closest_entities( const stk::math::Vector3d & search_point, const GetEntityDistanceSquaredFn & get_entity_distance_squared_fn, std::vector<EntityType> & return_vec, const Real max_search_radius ) const
{
  return_vec.clear();

  // if there are no entities at all return empty list
  if ( empty() ) return;

  const VecType search_pt(search_point.data());
  const Real max_search_radius2 = max_search_radius * max_search_radius;

  Real lowerBnd2;
  Real upperBnd2;
  if constexpr (std::is_same_v<GetEntityDistanceSquaredFn, SearchTreeNoLeafDistance>)
  {
    // No leaf-distance refinement: bounds come from bounding boxes only.
    std::tie(lowerBnd2, upperBnd2) = compute_tree_bounds(search_pt, get_entity_distance_squared_fn, max_search_radius2);
  }
  else
  {
    // Bind the search point so the leaf callback matches the Real(EntityType)
    // shape expected by the bounds computation; the lambda type is known here
    // so it inlines with no std::function.
    const auto get_leaf_entity_distance_squared_fn = [&search_point, &get_entity_distance_squared_fn](const EntityType entity)->Real { return get_entity_distance_squared_fn(entity, search_point); };
    std::tie(lowerBnd2, upperBnd2) = compute_tree_bounds(search_pt, get_leaf_entity_distance_squared_fn, max_search_radius2);
  }

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
    const size_t index) const
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
template <class GetLeafEntityDistanceSquaredFn>
void
SearchTree<ENTITY>::recursively_update_tree_bounds( const VecType & search_point,
    const GetLeafEntityDistanceSquaredFn & get_leaf_entity_distance_squared_fn,
    Real & lowerBnd2,
    Real & upperBnd2,
    const size_t index) const
{
  const SearchTreeNode<EntityType> & node = my_nodes[index];
  if (!node.have_children())
  {
    if (0 == index)
    {
      // Special case of a single node.  Otherwise, this node will already have been handled by parent.
      compute_node_bounds_and_update_tree_bounds(search_point, get_leaf_entity_distance_squared_fn,lowerBnd2, upperBnd2, index);
    }
    return;
  }

  const size_t Lindex = index + 1;
  const Real L_lowerBnd2 = compute_node_bounds_and_update_tree_bounds(search_point, get_leaf_entity_distance_squared_fn, lowerBnd2, upperBnd2, Lindex);

  const size_t Rindex = node.get_right_index();
  const Real R_lowerBnd2 = compute_node_bounds_and_update_tree_bounds(search_point, get_leaf_entity_distance_squared_fn, lowerBnd2, upperBnd2, Rindex);

  if ( L_lowerBnd2 < R_lowerBnd2 )
  {
    if ( L_lowerBnd2 < upperBnd2 )
      recursively_update_tree_bounds( search_point, get_leaf_entity_distance_squared_fn, lowerBnd2, upperBnd2, Lindex );
    if ( R_lowerBnd2 < upperBnd2 )
      recursively_update_tree_bounds( search_point, get_leaf_entity_distance_squared_fn, lowerBnd2, upperBnd2, Rindex );
  }
  else
  {
    if ( R_lowerBnd2 < upperBnd2 )
      recursively_update_tree_bounds( search_point, get_leaf_entity_distance_squared_fn, lowerBnd2, upperBnd2, Rindex );
    if ( L_lowerBnd2 < upperBnd2 )
      recursively_update_tree_bounds( search_point, get_leaf_entity_distance_squared_fn, lowerBnd2, upperBnd2, Lindex );
  }
}

template <class ENTITY>
template <class InsertIntoBbox>
size_t
SearchTree<ENTITY>::build( typename std::vector<EntityCentroid>::iterator entity_begin,
    typename std::vector<EntityCentroid>::iterator entity_end,
    const InsertIntoBbox & insert_into_bbox,
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
SearchTree<ENTITY>::recursively_get_intersecting_entities( const BoundingBoxType& bbox, std::vector<EntityType>& return_vec, const size_t index ) const
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
SearchTree<ENTITY>::add_descendants( std::vector<EntityType>& return_vec, const size_t index ) const
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
