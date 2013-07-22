/*------------------------------------------------------------------------*/
/*                 Copyright 2013 Sandia Corporation.                     */      
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */      
/*  United States Government.                                             */      
/*------------------------------------------------------------------------*/

#include <Intrepid_FieldContainer.hpp>
#include <boost/shared_ptr.hpp>


#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_io/MeshReadWriteUtils.hpp>
#include <init/Ionit_Initializer.h>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>
#include <stk_util/diag/PrintTimer.hpp>

#include <stk_transfer/STKNode.hpp>
#include <stk_transfer/LinearInterpolate.hpp>
#include <stk_transfer/Transfer.hpp>

namespace bopt = boost::program_options;

typedef stk::mesh::Field<double>                       ScalarField ;
typedef stk::mesh::Field<double, stk::mesh::Cartesian> CartesianField ;


bool use_case_6_driver(stk::ParallelMachine  comm,
                      const std::string &working_directory,
                      const std::string &range_mesh,
                      const std::string &range_filetype,
                      const std::string &domain_mesh,
                      const std::string &domain_filetype)
{
  stk::diag::Timer timer("Transfer Use Case 6", 
                          use_case::TIMER_TRANSFER, 
                          use_case::timer());
  stk::diag::Timer timer_node_to_node(" Node To Node", timer);
  use_case::timerSet().setEnabledTimerMask(use_case::TIMER_ALL);

  const double TOLERANCE = 0.000001;

  bool status = true;

  //const double TOLERANCE = 0.000001;
  //const double  rand_max = RAND_MAX;
  enum {           DIM = 3  };

  stk::io::MeshData range_mesh_data(comm);
  stk::io::MeshData domain_mesh_data(comm);

  std::string filename = working_directory + range_mesh;
  range_mesh_data.open_mesh_database(filename, range_filetype);
  range_mesh_data.create_input_mesh();

  stk::mesh::MetaData &range_meta_data = range_mesh_data.meta_data();
  const stk::mesh::EntityRank node_rank = stk::mesh::MetaData::NODE_RANK;
  stk::mesh::Part & range_block         = range_meta_data.declare_part("nodes", node_rank);
  stk::mesh::CellTopology node_top (shards::getCellTopologyData<shards::Node>());
  stk::mesh::set_cell_topology( range_block,  node_top );


  const std::string data_field_name = "Sum_Of_Coordinates";
  ScalarField &range_coord_sum_field = stk::mesh::put_field( 
                        range_meta_data.declare_field<ScalarField>(data_field_name), 
                        stk::mesh::MetaData::NODE_RANK , 
                        range_meta_data.universal_part() );

  range_meta_data.commit();
  range_mesh_data.populate_bulk_data();
  stk::mesh::BulkData &range_bulk_data = range_mesh_data.bulk_data();

  filename = working_directory + domain_mesh;
  domain_mesh_data.open_mesh_database(filename, domain_filetype);
  domain_mesh_data.create_input_mesh();

  stk::mesh::MetaData &domain_meta_data = domain_mesh_data.meta_data();
  stk::mesh::Part & domain_block        = domain_meta_data.declare_part("nodes", node_rank);
  stk::mesh::CellTopology hex_top (shards::getCellTopologyData<shards::Hexahedron<> >());
  stk::mesh::CellTopology quad_top(shards::getCellTopologyData<shards::Quadrilateral<> >());
  stk::mesh::set_cell_topology( domain_block,      hex_top );
  stk::mesh::set_cell_topology( domain_block,      quad_top );
  const stk::mesh::EntityRank side_rank    = domain_meta_data.side_rank();
  stk::mesh::Part & block_skin       = domain_meta_data.declare_part("skin", side_rank);
  stk::mesh::set_cell_topology( block_skin, quad_top );
  
  ScalarField &domain_coord_sum_field = stk::mesh::put_field( 
                        domain_meta_data.declare_field<ScalarField>(data_field_name), 
                        stk::mesh::MetaData::NODE_RANK , 
                        domain_meta_data.universal_part() );
  domain_meta_data.commit();

  domain_mesh_data.populate_bulk_data();
  stk::mesh::BulkData &domain_bulk_data = domain_mesh_data.bulk_data();
  stk::mesh::skin_mesh(domain_bulk_data, stk::mesh::MetaData::ELEMENT_RANK, &block_skin);
  // For this use case, the domain consists of an axis-aligned
  // bounding box for each 'domain_entity' in the mesh.  The range is a
  // PointBoundingBox3D at the centroid of each 'range_entity'.  The id of the point
  // will be the same as the id of the containing entity.  If the
  // mesh contains solid elements only, and the range_mesh matches the
  // domain_mesh, then the search should return a single box for each
  // point and the id of the box should match the id of the point.
  
  CartesianField &range_coord_field  = static_cast<CartesianField&>( range_mesh_data.get_coordinate_field());
  CartesianField &domain_coord_field = static_cast<CartesianField&>(domain_mesh_data.get_coordinate_field());

  stk::mesh::Selector range_nodes = range_meta_data.locally_owned_part();
  stk::mesh::Selector domain_nodes= domain_meta_data.locally_owned_part();
  
  std::vector<stk::mesh::Entity> domain_entities;
  {
    stk::mesh::get_selected_entities(domain_nodes, domain_bulk_data.buckets(stk::mesh::MetaData::NODE_RANK), domain_entities);
    const size_t num_entities = domain_entities.size();
    for (size_t i = 0; i < num_entities; ++i) {
      const stk::mesh::Entity entity = domain_entities[i];
      double *entity_coordinates = domain_bulk_data.field_data(domain_coord_field, entity);
      double *entity_coord_sum   = domain_bulk_data.field_data(domain_coord_sum_field, entity);
      *entity_coord_sum = entity_coordinates[0] + entity_coordinates[1] + entity_coordinates[2];
    }
  }
  std::vector<stk::mesh::Entity> range_entities;
  {
    stk::mesh::get_selected_entities(range_nodes, range_bulk_data.buckets(stk::mesh::MetaData::NODE_RANK), range_entities);
    const size_t num_entities = range_entities.size();
    const double rand_max = RAND_MAX;
    for (size_t i = 0; i < num_entities; ++i) {
      const stk::mesh::Entity entity = range_entities[i];
      double *entity_coord_sum   = range_bulk_data.field_data(range_coord_sum_field, entity);
      *entity_coord_sum = rand()/rand_max;
    }
  }

  const double radius=.25;
  const std::vector<stk::mesh::FieldBase*> from_fields(1, &domain_coord_sum_field);
  boost::shared_ptr<stk::transfer::STKNode<3> >
    transfer_domain_mesh (new stk::transfer::STKNode<3>(domain_entities, domain_coord_field, from_fields, radius));

  const std::vector<stk::mesh::FieldBase*> to_fields  (1, &range_coord_sum_field);
  boost::shared_ptr<stk::transfer::STKNode<3> >
    transfer_range_mesh (new stk::transfer::STKNode<3>(range_entities, range_coord_field, to_fields, radius));

  
  stk::transfer::GeometricTransfer<
    class stk::transfer::LinearInterpolate<
      class stk::transfer::STKNode<3>, 
      class stk::transfer::STKNode<3>
    >
  >
  transfer(transfer_domain_mesh, transfer_range_mesh, "STK Transfer test Use case 6");
  
  {
    stk::diag::TimeBlock __timer_node_to_node(timer_node_to_node);
    try {
      transfer.initialize();
      transfer.apply();
    } catch (std::exception &e) {
      std::cout <<__FILE__<<":"<<__LINE__
                <<" Caught an std::exception with what string:"
                <<e.what()
                <<"      rethrowing....."
                <<std::endl;
      status = status && false;
    } catch (...) {
      std::cout <<__FILE__<<":"<<__LINE__
                <<" Caught an exception, rethrowing..."
                <<std::endl;
      status = status && false;
    } 
  }

  if (status) {
    std::vector<stk::mesh::Entity> range_points;
    stk::mesh::get_selected_entities(range_nodes,
                                     range_bulk_data.buckets(stk::mesh::MetaData::NODE_RANK),
                                     range_points);

    const unsigned TONUMPOINTS = range_points.size();

    bool success = true;
    for (unsigned i=0 ; i<TONUMPOINTS; ++i) {
      const  stk::mesh::Entity  entity = range_points[i];

      const double *ToValues = static_cast<double*>(range_bulk_data.field_data(range_coord_sum_field, entity));
      const double *ToPoints = static_cast<double*>(range_bulk_data.field_data(range_coord_field, entity));

      double check_l = 0;
      for (unsigned j=0 ; j<DIM; ++j) check_l += ToPoints[j];   
      if (TOLERANCE < fabs(check_l-ToValues[0]) ) { 
        const stk::mesh::EntityKey k = range_bulk_data.entity_key(entity);
        std::cout <<__FILE__<<":"<<__LINE__
                  <<" EntityKey:"<<k
                  <<" ToPoints:"<<ToPoints[0]<<" "<<ToPoints[1]<<" "<<ToPoints[2]
                  <<" ToValues:"<<ToValues[0]
                  <<" check:"<<check_l
                  <<" error:"<<fabs(check_l-ToValues[0])
                  <<std::endl;
        success = false;
      }   
    }
    status = status && success;
  }
  timer.stop();
//stk::diag::printTimersTable(std::cout, timer, 
//      stk::diag::METRICS_CPU_TIME | stk::diag::METRICS_WALL_TIME, false, comm);


  const bool collective_result = use_case::print_status(comm, status);
  return collective_result;
}
