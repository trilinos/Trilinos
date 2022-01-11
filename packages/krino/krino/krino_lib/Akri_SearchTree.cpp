// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_SearchTree.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Facet.hpp>
#include <Akri_ProlongationData.hpp>
#include <Akri_Vec.hpp>

#include <algorithm>

namespace krino{

template <class ENTITY>
SearchTree<ENTITY>::SearchTree( std::vector<BoundingBoxEntity> & entities )
{ /* %TRACE[ON]% */ Trace trace__("krino::SearchTree::SearchTree( std::vector<BoundingBoxEntity> & entities )"); /* %TRACE% */

  if (!entities.empty())
  {
    my_nodes.resize(2*entities.size() - 1);
    const size_t tree_size = build( entities.begin(), entities.end() );
    ThrowRequire(tree_size == my_nodes.size());
  }
}

template <class ENTITY>
SearchTree<ENTITY>::SearchTree( const std::vector<EntityType> & entities, const std::function<const BoundingBox&(const EntityType &)> & get_bbox )
{ /* %TRACE[ON]% */ Trace trace__("krino::SearchTree::SearchTree( EntityVec & entities, const std::function<BoundingBox&(EntityType &)> & get_bbox )"); /* %TRACE% */

  if (!entities.empty())
  {
    std::vector<BoundingBoxEntity> bbox_entities;
    bbox_entities.reserve(entities.size());
    for (auto && entity : entities)
    {
      bbox_entities.emplace_back(get_bbox(entity), entity);
    }
    my_nodes.resize(2*entities.size() - 1);
    const size_t tree_size = build( bbox_entities.begin(), bbox_entities.end() );
    ThrowRequire(tree_size == my_nodes.size());
  }
}

template <class ENTITY>
void
SearchTree<ENTITY>::find_closest_entities( const Vector3d & search_point, std::vector<EntityType> & return_vec, const Real max_search_radius )
{
  return_vec.clear();

  // Algorithm requires tol > std::numeric_limits<double>::epsilon(), here just pick multiplier
  const Real tol = 100.*std::numeric_limits<Real>::epsilon();

  // if there are no entities at all return empty list
  if ( empty() ) return;

  const VecType search_pt(search_point.data());

  // first find upper bound for distance
  // the very first estimate for this is the upper bound for the base of the tree
  Real upperBnd2 = my_nodes[0].boxUpperBnd2( search_pt );
  Real max_search_radius2 = max_search_radius * max_search_radius;

  if ( upperBnd2 > max_search_radius2 && max_search_radius2 > 0.0 )
  {
    upperBnd2 = max_search_radius2;
  }

  // now descend relevant parts of tree updating the bounds as we go
  Real lowerBnd2 = upperBnd2;
  treeBounds(search_pt, lowerBnd2, upperBnd2);

  // exit now if no entities within narrow_band
  if ( lowerBnd2 > max_search_radius2 && max_search_radius2 > 0.0 )
  {
    return;
  }

  // with the upper bound known, find all entities within the tolerance
  find_closest_entities( search_pt, upperBnd2, tol, return_vec );

  return;
}

template <class ENTITY>
void
SearchTree<ENTITY>::find_closest_entities(
    const VecType & search_point,
    const Real maxDist2,
    const Real tol,
    std::vector<EntityType> & nearest,
    const size_t index)
{ /* %TRACE% */  /* %TRACE% */
  const SearchTreeNode<EntityType> & node = my_nodes[index];
  if ( node.have_children() )
  {
    const size_t Lindex = index + 1;
    const Real L_lowerBnd2 = my_nodes[Lindex].boxLowerBnd2( search_point );
    const Real L_eps = tol*(maxDist2+my_nodes[Lindex].boxSize2());
    if ( L_lowerBnd2 < maxDist2+L_eps )
      find_closest_entities( search_point, maxDist2, tol, nearest, Lindex );

    const size_t Rindex = node.get_right_index();
    const Real R_lowerBnd2 = my_nodes[Rindex].boxLowerBnd2( search_point );
    const Real R_eps = tol*(maxDist2+my_nodes[Rindex].boxSize2());
    if ( R_lowerBnd2 < maxDist2+R_eps )
      find_closest_entities( search_point, maxDist2, tol, nearest, Rindex );
  }
  else
  {
    EntityType entity = node.get_leaf_entity();
    nearest.push_back( entity );
  }
}

template <class ENTITY>
void
SearchTree<ENTITY>::treeBounds( const VecType & search_point,
    Real & lowerBnd2,
    Real & upperBnd2,
    const size_t index)
{ /* %TRACE% */  /* %TRACE% */
  const SearchTreeNode<EntityType> & node = my_nodes[index];
  if (!node.have_children())
  {
    if (0 == index)
    {
      // Special case of a single node.  Otherwise, this node will already have been handled by parent.
      compute_and_update_bounds(search_point, lowerBnd2, upperBnd2, index);
    }
    return;
  }

  const size_t Lindex = index + 1;
  const Real L_lowerBnd2 = compute_and_update_bounds(search_point, lowerBnd2, upperBnd2, Lindex);

  const size_t Rindex = node.get_right_index();
  const Real R_lowerBnd2 = compute_and_update_bounds(search_point, lowerBnd2, upperBnd2, Rindex);

  if ( L_lowerBnd2 < R_lowerBnd2 )
  {
    if ( L_lowerBnd2 < upperBnd2 )
      treeBounds( search_point, lowerBnd2, upperBnd2, Lindex );
    if ( R_lowerBnd2 < upperBnd2 )
      treeBounds( search_point, lowerBnd2, upperBnd2, Rindex );
  }
  else
  {
    if ( R_lowerBnd2 < upperBnd2 )
      treeBounds( search_point, lowerBnd2, upperBnd2, Rindex );
    if ( L_lowerBnd2 < upperBnd2 )
      treeBounds( search_point, lowerBnd2, upperBnd2, Lindex );
  }
}

template <class ENTITY>
size_t
SearchTree<ENTITY>::build( typename std::vector<BoundingBoxEntity>::iterator entity_begin, typename std::vector<BoundingBoxEntity>::iterator entity_end, const size_t index )
{
  const size_t entity_size = std::distance(entity_begin, entity_end);

  // determine bounding box
  BoundingBox & node_bbox = my_nodes[index].bounding_box();
  for ( auto entity_it = entity_begin; entity_it != entity_end; ++entity_it )
  {
    const BoundingBox & entity_bbox = entity_it->first;
    node_bbox.accommodate( entity_bbox );
  }

  if ( entity_size > 1 )
  {
    // determine axis to split tree on
    int max_spread_dim = node_bbox.max_span_direction();

    auto compare = [max_spread_dim](const BoundingBoxEntity & lhs, const BoundingBoxEntity & rhs) { return lhs.first.center()[max_spread_dim] < rhs.first.center()[max_spread_dim]; };
    const size_t mid_pt = entity_size / 2;

    // This works around a bug in some versions of nth_element
#ifdef CRAY_LWK
      std::partial_sort( entity_begin, entity_begin + mid_pt, entity_end, compare );
#else
      std::nth_element( entity_begin, entity_begin + mid_pt, entity_end, compare );
#endif

    const size_t Lindex = index + 1;
    const size_t Rindex = build( entity_begin, entity_begin + mid_pt, Lindex );
    my_nodes[index].set_right_index(Rindex);
    const size_t next_index = build( entity_begin + mid_pt, entity_end, Rindex );
    return next_index;
  }
  else
  {
    ThrowAssert(1 == entity_size);
    my_nodes[index].set_leaf_entity(entity_begin->second);
    return index + 1;
  }
}

template <class ENTITY>
void
SearchTree<ENTITY>::get_intersecting_entities( const BoundingBoxType& bbox, std::vector<EntityType>& return_vec, const size_t index )
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
      get_intersecting_entities( bbox, return_vec, Lindex );

      const size_t Rindex = node.get_right_index();
      get_intersecting_entities( bbox, return_vec, Rindex );
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

// Explicit template instantiation
template class SearchTree<Facet*>;
template class SearchTree<const ProlongationFacet*>;

} // namespace krino
