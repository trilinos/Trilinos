/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <vector>
#include <fstream>

#include <assert.h>

#include <stk_io/StkMeshIoBroker.hpp>
#include <init/Ionit_Initializer.h>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>
#include <stk_util/diag/Writer.hpp>
#include <stk_util/diag/WriterExt.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/diag/IdentProc.hpp>

#include <stk_search_util/stk_mesh/CreateBoundingBox.hpp>
#include <stk_search_util/stk_mesh/PrintBoundingBox.hpp>
#include <stk_search_util/stk_mesh/PrintEntityProc.hpp>

#include <transfer/UseCaseIsInElement.hpp>

#include <Shards_BasicTopologies.hpp>

using namespace stk::diag;
using namespace stk::search;
using namespace use_case;

typedef stk::mesh::Field<double>                        ScalarField ;

void
use_case_3_driver(stk::ParallelMachine  comm,
                  const std::string &working_directory,
                  const std::string &range_mesh_filename,
                  const std::string &range_mesh_type,
                  const std::string &range_entity,
                  const std::string &domain_mesh_filename,
                  const std::string &domain_mesh_type,
                  const std::string &domain_entity,
                  double offset,
                  double scale)
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
  const int my_rank = comm_all.parallel_rank();

  dw().m(LOG_TRANSFER) << "Use case 3: Point (range) in Box (domain) Search" << stk::diag::dendl;
  dw().m(LOG_TRANSFER) << "Range  Entity Type = " << range_entity  << stk::diag::dendl;
  dw().m(LOG_TRANSFER) << "Domain Entity Type = " << domain_entity << stk::diag::dendl;

  stk::io::StkMeshIoBroker range_mesh_data(comm);
  std::string filename = working_directory + range_mesh_filename;
  range_mesh_data.open_mesh_database(filename, range_mesh_type);
  range_mesh_data.create_input_mesh();
  range_mesh_data.populate_bulk_data();

  stk::mesh::MetaData &range_meta_data = range_mesh_data.meta_data();
  stk::mesh::BulkData &range_bulk_data = range_mesh_data.bulk_data();

  stk::io::StkMeshIoBroker domain_mesh_data(comm);
  filename = working_directory + domain_mesh_filename;
  domain_mesh_data.open_mesh_database(filename, domain_mesh_type);
  domain_mesh_data.create_input_mesh();

  stk::mesh::MetaData &domain_meta_data = domain_mesh_data.meta_data();
  const stk::mesh::EntityRank element_rank = stk::mesh::MetaData::ELEMENT_RANK;
  const stk::mesh::EntityRank side_rank    = domain_meta_data.side_rank();
  stk::mesh::Part & block_hex        = domain_meta_data.declare_part("block_1", element_rank);
  stk::mesh::Part & block_quad       = domain_meta_data.declare_part("block_2", side_rank);
  stk::mesh::CellTopology hex_top (shards::getCellTopologyData<shards::Hexahedron<> >());
  stk::mesh::CellTopology quad_top(shards::getCellTopologyData<shards::Quadrilateral<> >());
  stk::mesh::set_cell_topology( block_hex,  hex_top );
  stk::mesh::set_cell_topology( block_quad, quad_top );
  stk::mesh::Part & block_skin       = domain_meta_data.declare_part("skin", side_rank);
  stk::mesh::set_cell_topology( block_skin, quad_top );
  domain_meta_data.commit();

  domain_mesh_data.populate_bulk_data();
  stk::mesh::BulkData &domain_bulk_data = domain_mesh_data.bulk_data();
  {
    stk::mesh::PartVector add_parts(1,&block_skin);
    stk::mesh::skin_mesh(domain_bulk_data, add_parts);
  }

  // For this use case, the domain consists of an axis-aligned
  // bounding box for each 'domain_entity' in the mesh.  The range is a
  // PointBoundingBox3D at the centroid of each 'range_entity'.  The id of the point
  // will be the same as the id of the containing entity.  If the
  // mesh contains solid elements only, and the range_mesh matches the
  // domain_mesh, then the search should return a single box for each
  // point and the id of the box should match the id of the point.

  CartesianField const& range_coord_field = static_cast<CartesianField const&>(range_mesh_data.get_coordinate_field());
  std::vector<PointBoundingBox3D> range_vector;

  {
    stk::diag::TimeBlock __timer_range_bb(timer_range_bb);

    // For this use case, if 'node' is the range_entity_rank, then we
    // want to use the universal set nodes, not a nodeset.
    bool use_universal_set = true;
    stk::search_util::build_centroid_bbox(range_bulk_data,
                                          range_meta_data.entity_rank(range_entity),
                                          range_coord_field, range_vector,
                                          use_universal_set);
  }

  CartesianField const& domain_coord_field = static_cast<CartesianField const&>(domain_mesh_data.get_coordinate_field());
  std::vector<AxisAlignedBoundingBox3D> domain_vector;

  {
    stk::diag::TimeBlock __timer_domain_bb(timer_domain_bb);

    bool use_universal_set = true; // Use the 'skin' part
    stk::search_util::build_axis_aligned_bbox(domain_bulk_data,
                                              domain_meta_data.entity_rank(domain_entity),
                                              domain_coord_field, domain_vector,
					      use_universal_set,
					      stk::search_util::OffsetScaleOp(scale, offset));
  }

  dw().m(LOG_TRANSFER) << "range  " << range_vector  << stk::diag::dendl;
  dw().m(LOG_TRANSFER) << "domain " << domain_vector << stk::diag::dendl;

  FactoryOrder order;
  order.m_communicator = comm;



  dw().m(LOG_TRANSFER) << "Search algorithm " << order.m_algorithm << dendl;

  IdentProcRelation relation;

  {
    stk::diag::TimeBlock __timer_range_search(timer_range_search);
    stk::search::coarse_search(relation, range_vector, domain_vector, order);
  }

  dw().m(LOG_TRANSFER) << "relation :" << relation;

  dw().m(LOG_TRANSFER) << "Pre-Ghosting Info:" << stk::diag::push << stk::diag::dendl;
  stk::search_util::print_entity_proc_map(dw().m(LOG_TRANSFER), range_bulk_data );
  dw().m(LOG_TRANSFER) << stk::diag::pop;

  std::vector<stk::mesh::EntityProc> range_to_ghost ;

  IdentProcRelation::const_iterator i=relation.begin(), end=relation.end();
  for (;i!=end; ++i) {
    stk::mesh::EntityKey range_entity_key(i->second.ident);

    const int domain_owning_proc = i->first.proc;
    const int range_owning_rank  = i->second.proc;
    if (domain_owning_proc != my_rank && range_owning_rank == my_rank) {
      stk::mesh::Entity r_entity = range_bulk_data.get_entity(range_entity_key.rank(), range_entity_key.id());
      if (range_bulk_data.parallel_owner_rank(r_entity) == my_rank) {
        stk::mesh::EntityProc ep(r_entity, domain_owning_proc);
        range_to_ghost.push_back(ep);
      }
    }
  }

  dw().m(LOG_TRANSFER) << "Change ghosts to send:";
  stk::search_util::print_entity_proc_map(dw().m(LOG_TRANSFER), range_bulk_data, range_to_ghost, "Is ghosting ", " to ");


  {
    stk::diag::TimeBlock __timer_ghosting(timer_ghosting);

    range_bulk_data.modification_begin();

    stk::mesh::Ghosting & transfer_range_ghosting =
      range_bulk_data.create_ghosting( std::string("transter_test") );

    {
      std::vector<stk::mesh::EntityKey> receive ;
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
    fields.push_back(&range_coord_field);

    stk::mesh::communicate_field_data( transfer_range_ghosting , fields);
  }

  std::vector<std::pair<stk::mesh::Entity , stk::mesh::Entity> > entity_map;
  {
    IdentProcRelation::const_iterator I=relation.begin(), rend=relation.end();
    for ( ; I!=rend; ++I) {
      const int domain_owning_proc = I->first.proc;
      if (domain_owning_proc == my_rank) {
        stk::mesh::EntityKey domain_entity_key(I->first.ident);
        stk::mesh::EntityKey range_entity_key(I->second.ident);

        stk::mesh::Entity d_entity = domain_bulk_data.get_entity(domain_entity_key.rank(), domain_entity_key.id());
        stk::mesh::Entity r_entity = range_bulk_data.get_entity (range_entity_key.rank(), range_entity_key.id());
        assert(domain_bulk_data.is_valid(d_entity));
        assert(range_bulk_data.is_valid(r_entity));
        std::pair<stk::mesh::Entity , stk::mesh::Entity> e(d_entity, r_entity);
        entity_map.push_back(e);
      }
    }
  }

  std::sort(entity_map.begin(), entity_map.end());
  entity_map.erase(std::unique (entity_map.begin(), entity_map.end()), entity_map.end());

  if (dw().shouldPrint(LOG_TRANSFER)) {
    dw() << "[" << my_rank << "]  Detailed search found " << entity_map.size()
          << " range nodes in the " << domain_vector.size() << " domain faces." << stk::diag::push << stk::diag::dendl;
    stk::search_util::print_entity_map(dw().m(LOG_TRANSFER), domain_bulk_data, range_bulk_data,
                                       entity_map, " may contain ");
    dw() << stk::diag::pop << stk::diag::dendl;
  }

  std::vector<std::size_t> not_in_element;
  stk::usecase::is_in_element(domain_bulk_data, range_bulk_data,
                              domain_coord_field, range_coord_field, entity_map, not_in_element);

  dw().m(LOG_TRANSFER) << "[" << my_rank << "] Detailed search found " << not_in_element.size()
                       << " nodes that were not in the face following the detailed search.";
  dw().m(LOG_TRANSFER) << not_in_element << stk::diag::dendl;

  std::vector<std::size_t>::reverse_iterator I=not_in_element.rbegin();
  for (; I!=not_in_element.rend(); ++I)
  {
    std::vector<std::pair<stk::mesh::Entity , stk::mesh::Entity> >::iterator del = entity_map.begin()+*I;
    entity_map.erase(del);
  }

  dw().m(LOG_TRANSFER) << "[" << my_rank << "]  Detailed search found " << entity_map.size()
                       << " range nodes in the " << domain_vector.size() << " domain faces.\n";
  stk::search_util::print_entity_map(dw().m(LOG_TRANSFER), domain_bulk_data, range_bulk_data,
                                     entity_map, " contains ");
}
