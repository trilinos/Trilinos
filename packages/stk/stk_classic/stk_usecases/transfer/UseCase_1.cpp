/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <vector>
#include <fstream>

#include <assert.h>

#include <stk_io/MeshReadWriteUtils.hpp>
#include <init/Ionit_Initializer.h>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>
#include <stk_util/diag/Writer.hpp>
#include <stk_util/diag/WriterExt.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

#include <stk_mesh/diag/EntityKey.hpp>

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/diag/IdentProc.hpp>

#include <stk_search_util/stk_mesh/CreateBoundingBox.hpp>
#include <stk_search_util/stk_mesh/PrintBoundingBox.hpp>
#include <stk_search_util/stk_mesh/PrintEntityProc.hpp>

#include <transfer/UseCaseIsInElement.hpp>

using namespace stk::diag;
using namespace use_case;

typedef stk::mesh::Field<double>                        ScalarField ;
typedef stk::mesh::Field<double, stk::mesh::Cartesian>  VectorField ;

void
use_case_1_driver(
  stk::ParallelMachine  comm,
  const std::string &working_directory,
  const std::string &range_mesh_filename,
  const std::string &range_mesh_type,
  const std::string &range_entity,
  const std::string &domain_mesh_filename,
  const std::string &domain_mesh_type,
  const std::string &domain_entity)
{
  stk::diag::WriterThrowSafe _write_throw_safe(dw());

  stk::diag::Timer timer("Transfer", use_case::TIMER_TRANSFER, use_case::timer());
  stk::diag::Timer timer_transfer(domain_mesh_filename + " to " + range_mesh_filename, timer);
  stk::diag::Timer timer_range_bb("Range bounding box", timer_transfer);
  stk::diag::Timer timer_domain_bb("Domain bounding box", timer_transfer);
  stk::diag::Timer timer_range_search("Range search", timer_transfer);
  stk::diag::Timer timer_ghosting("Ghosting", timer_transfer);

  stk::diag::TimeBlock __timer_transfer(timer_transfer);

  stk::CommAll comm_all( comm );
  const unsigned my_rank = comm_all.parallel_rank();

  dw().m(LOG_TRANSFER) << "Use case 1: Point (range) in Box (domain) Search" << stk::diag::dendl;
  dw().m(LOG_TRANSFER) << "Range  Entity Type = " << range_entity  << stk::diag::dendl;
  dw().m(LOG_TRANSFER) << "Domain Entity Type = " << domain_entity << stk::diag::dendl;

  // Initialize IO system.  Registers all element types and storage
  // types and the exodusII default database type.
  Ioss::Init::Initializer init_db;

  static const size_t spatial_dimension = 3;
  stk::mesh::fem::FEMMetaData range_meta_data( spatial_dimension );
  stk::io::MeshData range_mesh_data;
  std::string filename = working_directory + range_mesh_filename;
  stk::io::create_input_mesh(range_mesh_type, filename, comm,
			     range_meta_data, range_mesh_data);
  range_meta_data.commit();

  stk::mesh::BulkData range_bulk_data(range_meta_data.get_meta_data(range_meta_data) , comm);
  stk::io::populate_bulk_data(range_bulk_data, range_mesh_data);

  stk::mesh::fem::FEMMetaData domain_meta_data( spatial_dimension );
  stk::io::MeshData domain_mesh_data;
  filename = working_directory + domain_mesh_filename;
  stk::io::create_input_mesh(domain_mesh_type, domain_mesh_filename, comm,
			     domain_meta_data, domain_mesh_data);
  const stk::mesh::EntityRank element_rank = domain_meta_data.element_rank();
  stk::mesh::Part & block_hex        = domain_meta_data.declare_part("block_1", element_rank);
  stk::mesh::fem::CellTopology hex_top(shards::getCellTopologyData<shards::Hexahedron<> >());
  stk::mesh::fem::set_cell_topology( block_hex, hex_top );
  domain_meta_data.commit();

  stk::mesh::BulkData domain_bulk_data(domain_meta_data.get_meta_data(domain_meta_data) , comm);
  stk::io::populate_bulk_data(domain_bulk_data, domain_mesh_data);

  // For this use case, the domain consists of an axis-aligned
  // bounding box for each 'domain_entity' in the mesh.  The range is a
  // PointBoundingBox3D at the centroid of each 'range_entity'.  The id of the point
  // will be the same as the id of the containing entity.  If the
  // mesh contains solid elements only, and the range_mesh matches the
  // domain_mesh, then the search should return a single box for each
  // point and the id of the box should match the id of the point.

  VectorField *range_coord_field = range_meta_data.get_field<VectorField>("coordinates");
  std::vector<PointBoundingBox3D> range_vector;

  {
    stk::diag::TimeBlock __timer_range_bb(timer_range_bb);

    // For this use case, if 'node' is the range_entity_rank, then we
    // want to use the universal set nodes, not a nodeset.
    bool node_use_universal_set = true;
    stk::search_util::build_centroid_bbox(range_bulk_data,
                                          range_meta_data.entity_rank(range_entity),
                                          range_coord_field, range_vector,
                                          node_use_universal_set);
  }

  VectorField *domain_coord_field = domain_meta_data.get_field<VectorField>("coordinates");
  std::vector<AxisAlignedBoundingBox3D> domain_vector;

  {
    stk::diag::TimeBlock __timer_domain_bb(timer_domain_bb);

    stk::search_util::build_axis_aligned_bbox(domain_bulk_data,
                                              domain_meta_data.entity_rank(domain_entity),
                                              domain_coord_field, domain_vector,
					      false,
					      stk::search_util::NoOp());
  }

  dw().m(LOG_TRANSFER) << "range  " << range_vector  << stk::diag::dendl;
  dw().m(LOG_TRANSFER) << "domain " << domain_vector << stk::diag::dendl;

  stk::search::FactoryOrder order;
  order.m_communicator = comm;


  dw().m(LOG_TRANSFER) << "Search algorithm " << order.m_algorithm << dendl;

  IdentProcRelation relation;

  {
    stk::diag::TimeBlock __timer_range_search(timer_range_search);
    stk::search::coarse_search(relation, range_vector, domain_vector,order);
  }

  dw().m(LOG_TRANSFER) << "relation :" << relation;

  dw().m(LOG_TRANSFER) << "Pre-Ghosting Info:" << stk::diag::push << stk::diag::dendl;
  stk::search_util::print_entity_proc_map(dw().m(LOG_TRANSFER), range_bulk_data );
  dw().m(LOG_TRANSFER) << stk::diag::pop;

  // Ghosting of range-mesh entity to bring it to the
  // domain-mesh entity processor

  std::vector<stk::mesh::EntityProc> range_to_ghost ;

  IdentProcRelation::const_iterator i=relation.begin(), end=relation.end();

  for (;i!=end; ++i) {
    const std::size_t domain_owning_proc = i->first.proc;
    const std::size_t range_owning_rank  = i->second.proc;
    if (domain_owning_proc != my_rank && range_owning_rank == my_rank) {
      stk::mesh::EntityKey entity_key(i->second.ident);

      stk::mesh::Entity *r_entity = range_bulk_data.get_entity(stk::mesh::entity_rank(entity_key), stk::mesh::entity_id(entity_key));
      if (r_entity->owner_rank() == my_rank) {
        stk::mesh::EntityProc ep(r_entity, domain_owning_proc);
        range_to_ghost.push_back(ep);
      }
    }
  }

  dw().m(LOG_TRANSFER) << "Change ghosts to send:";
  stk::search_util::print_entity_proc_map(dw().m(LOG_TRANSFER), range_to_ghost, "Is ghosting ", " to ");

  {
    stk::diag::TimeBlock __timer_ghosting(timer_ghosting);

    //  Reset the transfer range ghosting by adding the
    //  new to-ghost list and removing the old is-ghosted list.
    //  The 'change_ghosting' method will sort out if anything
    //  is removed and then added.

    range_bulk_data.modification_begin();

    stk::mesh::Ghosting & transfer_range_ghosting =
      range_bulk_data.create_ghosting( std::string("transter_test") );

    {
      std::vector<stk::mesh::Entity*> receive;
      transfer_range_ghosting.receive_list( receive );
      range_bulk_data.change_ghosting( transfer_range_ghosting ,
                                       range_to_ghost ,
                                       receive );
    }

    range_bulk_data.modification_end();

    dw().m(LOG_TRANSFER) << "Post-Ghosting Info" << stk::diag::push << stk::diag::dendl;
    stk::search_util::print_entity_proc_map(dw().m(LOG_TRANSFER), range_bulk_data );
    dw().m(LOG_TRANSFER) << stk::diag::pop;

    // Copy coordinates to the newly ghosted nodes

    std::vector<const stk::mesh::FieldBase *> fields ;
    fields.push_back(range_coord_field);

    stk::mesh::communicate_field_data( transfer_range_ghosting , fields);
  }

  // Create the on-processor domain-to-range (range possibly ghosted) map.
  std::vector<std::pair<stk::mesh::Entity*, stk::mesh::Entity*> > entity_map;
  {
    IdentProcRelation::const_iterator I=relation.begin(), rend=relation.end();
    for ( ; I!=rend; ++I) {
      const std::size_t domain_owning_proc = I->first.proc;
      if (domain_owning_proc == my_rank) {
        stk::mesh::EntityKey domain_entity_key(I->first.ident);
        stk::mesh::EntityKey range_entity_key(I->second.ident);

        stk::mesh::Entity * const d_entity = domain_bulk_data.get_entity(stk::mesh::entity_rank(domain_entity_key), stk::mesh::entity_id(domain_entity_key));
        stk::mesh::Entity * const r_entity = range_bulk_data.get_entity (stk::mesh::entity_rank(range_entity_key), stk::mesh::entity_id(range_entity_key));
        if ( d_entity == 0 || r_entity == 0 ) {
          std::ostringstream msg ;
          msg << "P" << my_rank << ": Ghosting failed for " ;
          msg << "( " ;
          stk::mesh::print_entity_key( msg , range_meta_data.get_meta_data(range_meta_data) , domain_entity_key );
          if ( d_entity == NULL ) { msg << " -> NULL" ; }
          msg << " , " << I->first.proc ;
          msg << " ) ( " ;
          stk::mesh::print_entity_key( msg , range_meta_data.get_meta_data(range_meta_data) , range_entity_key );
          if ( r_entity == NULL ) { msg << " -> NULL" ; }
          msg << " , " << I->second.proc ;
          msg << " )" ;
          throw std::runtime_error( msg.str() );
        }
        std::pair<stk::mesh::Entity*, stk::mesh::Entity*> e(d_entity, r_entity);
        entity_map.push_back(e);
      }
    }
  }

  std::sort(entity_map.begin(), entity_map.end());
  entity_map.erase(std::unique (entity_map.begin(), entity_map.end()), entity_map.end());

  std::vector<std::size_t> not_in_element;
  stk::usecase::is_in_element(domain_coord_field, range_coord_field, entity_map, not_in_element);

  dw().m(LOG_TRANSFER) << "[" << my_rank << "] Detailed search found " << not_in_element.size()
                       << " nodes that were not in the element following the detailed search.";
  dw().m(LOG_TRANSFER) << not_in_element << stk::diag::dendl;

  std::vector<std::size_t>::reverse_iterator I=not_in_element.rbegin();
  for (; I!=not_in_element.rend(); ++I)
  {
    std::vector<std::pair<stk::mesh::Entity*, stk::mesh::Entity*> >::iterator del = entity_map.begin()+*I;
    entity_map.erase(del);
  }

  dw().m(LOG_TRANSFER) << "[" << my_rank << "]  Detailed search found " << entity_map.size()
                       << " range nodes in the " << domain_vector.size() << " domain elements.\n";
  stk::search_util::print_entity_map(dw().m(LOG_TRANSFER), entity_map, " contains ");
}

