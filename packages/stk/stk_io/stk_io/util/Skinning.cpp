/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <assert.h>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <stk_io/util/Skinning.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <Shards_CellTopology.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>

namespace use_case {
#if 0
  std::ostream &pout();               ///< Per-processor output stream
#endif
}

// #include <stk_util/diag/Writer.hpp>
// #include <stk_util/diag/WriterExt.hpp>
// stk::diag::Writer &dw();

namespace {

  static const stk::mesh::EntityRank NODE_RANK = stk::mesh::fem::FEMMetaData::NODE_RANK;

  stk::mesh::Entity *
  get_side_neighbor(const shards::CellTopology &  elem_top ,
		    const stk::mesh::Entity &     elem ,
		    unsigned                      side_id )
  {
    const unsigned side_dimension = elem_top.getDimension() - 1;

    const shards::CellTopology side_top(elem_top.getCellTopologyData(elem_top.getDimension() - 1, side_id));

    const stk::mesh::PairIterRelation elem_nodes = elem.relations( NODE_RANK);

    // Find other element that shares this side...
    stk::mesh::Entity & node = * elem_nodes[ elem_top.getNodeMap(side_dimension, side_id, 0) ].entity();

    const stk::mesh::PairIterRelation node_elems = node.relations( elem.entity_rank());

    stk::mesh::Entity * neighbor = NULL ;

    for ( unsigned i = 0 ; neighbor == NULL && i < node_elems.size() ; ++i ) {

      neighbor = node_elems[i].entity();

      const stk::mesh::PairIterRelation neighbor_nodes = neighbor->relations( NODE_RANK);

      if ( & elem == neighbor ) { neighbor = NULL ; }

      for ( unsigned j = 1 ;
	    neighbor != NULL && j < side_top.getNodeCount() ; ++j ) {

	stk::mesh::Entity * const next_node = elem_nodes[ elem_top.getNodeMap(side_dimension, side_id, j) ].entity();

	// If neighbor does not have node then not this element ...

	bool found = false ;
	for ( unsigned k = 0 ; ! found && k < neighbor_nodes.size() ; ++k ) {
	  found = next_node == neighbor_nodes[k].entity();
	}
	if ( ! found ) { neighbor = NULL ; }
      }

#if 0
      if ( NULL != neighbor ) {
	use_case::pout() << "neighbors( " ;
	use_case::pout() << " Element[ " ;
	use_case::pout() << elem.identifier();
	use_case::pout() << " ]{" ;
	for ( size_t i = 0 ; i < elem_nodes.size() ; ++i ) {
	  use_case::pout() << " " << elem_nodes[i].entity()->identifier();
	}
	use_case::pout() << " } , Element[ " ;
	use_case::pout() << neighbor->identifier();
	use_case::pout() << " ]{" ;
	for ( size_t i = 0 ; i < neighbor_nodes.size() ; ++i ) {
	  use_case::pout() << " " << neighbor_nodes[i].entity()->identifier();
	}
	use_case::pout() << " } , Share { " ;
	for ( unsigned j = 0 ; j < side_top.getNodeCount() ; ++j ) {
	  Entity * const next_node = elem_nodes[ elem_top.getNodeMap(side_dimension, side_id, j) ].entity();
	  use_case::pout() << " " << next_node->identifier();
	}
	use_case::pout() << " } )" ;
	use_case::pout() << std::endl ;
      }
#endif
    }

    return neighbor ;
  }


  unsigned
  determine_local_side_id(const stk::mesh::Entity &     elem,
			  stk::mesh::Entity &           side )
  {
    const shards::CellTopology elem_top(stk::mesh::fem::get_cell_topology( elem ).getCellTopologyData());
    const stk::mesh::PairIterRelation elem_nodes = elem.relations( NODE_RANK);
    const stk::mesh::PairIterRelation side_nodes = side.relations( NODE_RANK);

    const unsigned side_dimension = elem_top.getDimension() - 1;

    int side_id = -1 ;

    for ( unsigned i = 0 ; side_id == -1 && i < elem_top.getSideCount() ; ++i ) {
      const shards::CellTopology side_top(elem_top.getCellTopologyData(2, side_id));

      if ( side_nodes.size() == side_top.getNodeCount() ) {

	side_id = i ;

	for ( unsigned j = 0 ; side_id == static_cast<int>(i) && j < side_top.getNodeCount() ; ++j ) {

	  stk::mesh::Entity * const elem_node = elem_nodes[ elem_top.getNodeMap(side_dimension, i, j) ].entity();

	  bool found = false ;

	  for ( unsigned k = 0 ; ! found && k < side_top.getNodeCount() ; ++k ) {
	    found = elem_node == side_nodes[k].entity();
	  }

	  if ( ! found ) { side_id = -1 ; }
	}
      }
    }

    if ( side_id < 0 ) {
      std::ostringstream msg ;
      msg << "determine_local_side_id( " ;
      msg << elem_top.getName() ;
      msg << " , Element[ " ;
      msg << elem.identifier();
      msg << " ]{" ;
      for ( unsigned i = 0 ; i < elem_nodes.size() ; ++i ) {
	msg << " " << elem_nodes[i].entity()->identifier();
      }
      msg << " } , Side[ " ;
      msg << side.identifier();
      msg << " ]{" ;
      for ( unsigned i = 0 ; i < side_nodes.size() ; ++i ) {
	msg << " " << side_nodes[i].entity()->identifier();
      }
      msg << " } ) FAILED" ;
      throw std::runtime_error( msg.str() );
    }

    return static_cast<unsigned>(side_id) ;
  }

  void generate_element_sides(stk::mesh::BulkData & mesh,
			      stk::mesh::Entity & element,
			      stk::mesh::Part &     side_part,
			      const bool            skin_only )
  {
    // Generate if either the element or its neighbor is owned...
    // Only generate if the element has a smaller identifier
    // If 'skin_only' is true, then only generate if the side is
    // 'exposed' i.e., there is no neighbor sharing the face with this
    // element.

    // For each element
    //   for each element face
    //     get the neighbor element
    //     if no neighbor and owned then generate the face
    //     else if element has smaller id and
    //             either element or neighbor is local
    //       then generate face and attach to neighbor
    //         if either element or neighbor is not local
    //           then add to sharing

    stk::mesh::fem::FEMMetaData * fem_meta = const_cast<stk::mesh::fem::FEMMetaData *>(stk::mesh::MetaData::get(mesh).get_attribute<stk::mesh::fem::FEMMetaData>());

    shards::CellTopology elem_topology = stk::mesh::fem::get_cell_topology(element);
    const shards::CellTopology *elem_top = &elem_topology;
    if (fem_meta && !elem_top->getCellTopologyData()) {
      elem_topology = stk::mesh::fem::get_cell_topology(element);
      elem_top = &elem_topology;
    }

    const unsigned p_rank = mesh.parallel_rank();
    const bool element_owned  = p_rank == element.owner_rank();

    for ( unsigned i = 0 ; i < elem_top->getSideCount() ; ++i ) {

      stk::mesh::Entity * const elem_neighbor = get_side_neighbor(*elem_top, element, i);
      const bool neighbor_owned = elem_neighbor && (elem_neighbor->owner_rank() == p_rank);

      const bool create_side =
	( element_owned || neighbor_owned ) &&
	( ! elem_neighbor ||
	  ( ! skin_only &&
	    element.identifier() < elem_neighbor->identifier() ) );

      if ( create_side ) {
	const shards::CellTopology side_top(elem_top->getCellTopologyData(2, i));

	const stk::mesh::EntityRank side_type = static_cast<stk::mesh::EntityRank>(element.entity_rank() - 1);

	const unsigned side_id = element.identifier() * 10 + i + 1;
        const stk::mesh::EntityId global_side_id(side_id);

	stk::mesh::PartVector parts ;

	parts.push_back(&side_part);

        if (fem_meta) {
	  stk::mesh::Entity & side = stk::mesh::fem::declare_element_side( mesh, global_side_id, element, i, &side_part) ;
	if ( elem_neighbor ) {

	  const unsigned other_side_id = determine_local_side_id( *elem_neighbor , side );

	  mesh.declare_relation( *elem_neighbor , side , other_side_id );
	}
        } else {

	stk::mesh::Entity & side = fem_meta ? stk::mesh::fem::declare_element_side( mesh, global_side_id, element, i, &side_part) :
                                                               mesh.declare_entity( side_type, side_id , parts );

	stk::mesh::PairIterRelation rel = element.relations( NODE_RANK);

	for ( unsigned k = 0 ; k < side_top.getNodeCount() ; ++k ) {
	  stk::mesh::Entity & node = * rel[ elem_top->getNodeMap(elem_top->getDimension() - 1, i, k) ].entity();
	  mesh.declare_relation( side , node , k );
	}

	assert(side.relations(NODE_RANK).size() == side_top.getNodeCount());

	/** \todo REFACTOR Eliminate const_cast... */
	mesh.declare_relation( const_cast<stk::mesh::Entity&>(element), side , i );

	if ( elem_neighbor ) {

	  const unsigned other_side_id = determine_local_side_id( *elem_neighbor , side );

	  mesh.declare_relation( *elem_neighbor , side , other_side_id );
	}
        }
      }
    }
  }

} // namespace <empty>


namespace stk {
  namespace io {
    namespace util {

      void
      generate_sides(stk::mesh::BulkData & mesh,
		     stk::mesh::Part     & side_part,
		     const bool            skin_only )
      {
	// Generate one layer of ghost mesh
        mesh.modification_begin();

	const stk::mesh::fem::FEMMetaData& meta_data = stk::mesh::fem::FEMMetaData::get(mesh);
	const stk::mesh::PartVector & all_parts = meta_data.get_parts();
	for (stk::mesh::PartVector::const_iterator ip = all_parts.begin(); ip != all_parts.end(); ++ip) {
	  stk::mesh::Part *part = *ip;

	  // Filter out parts with "non-solid" (hexes and tets) topology...
	  const CellTopologyData * cell_topo = NULL;
	  cell_topo = meta_data.get_cell_topology(*part).getCellTopologyData();
          stk::mesh::fem::FEMMetaData * fem_meta = const_cast<stk::mesh::fem::FEMMetaData *>(stk::mesh::MetaData::get(mesh).get_attribute<stk::mesh::fem::FEMMetaData>());
          if (fem_meta && !cell_topo) cell_topo = fem_meta->get_cell_topology(*part).getCellTopologyData();
	  if (cell_topo == NULL || cell_topo->dimension != 3)
	    continue;

          stk::mesh::Selector selector = *part & meta_data.locally_owned_part();
          const std::vector<stk::mesh::Bucket*>& all_element_buckets = mesh.buckets(fem_meta->element_rank());
	  std::vector<stk::mesh::Bucket *> elem_buckets;
	  stk::mesh::get_buckets(selector, all_element_buckets, elem_buckets);

          // For each bucket ...
          for(size_t i=0; i<elem_buckets.size(); ++i) {
            // For each element ...
            stk::mesh::Bucket& bucket = *elem_buckets[i];
            stk::mesh::Bucket::iterator
              i_elem = bucket.begin(),
              i_end = bucket.end();

            for(; i_elem != i_end ; ++i_elem ) {
              stk::mesh::Entity & element = *i_elem;
              generate_element_sides(mesh, element, side_part, skin_only);
            }
	  }
	}

        mesh.modification_end();
      }

    } // namespace util
  } // namespace io
} // namespace stk
