// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <memory>

#include <Intrepid_FieldContainer.hpp>

#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/FieldBase.hpp>

#include <stk_io/StkMeshIoBroker.hpp>
#include <init/Ionit_Initializer.h>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>
#include <stk_util/diag/PrintTimer.hpp>

#include "STKNode.hpp"
#include "LinearInterpolate.hpp"
#include <stk_transfer/GeometricTransfer.hpp>

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

  stk::io::StkMeshIoBroker range_mesh_data(comm);
  stk::io::StkMeshIoBroker domain_mesh_data(comm);

  std::string filename = working_directory + range_mesh;
  range_mesh_data.add_mesh_database(filename, range_filetype, stk::io::READ_MESH);
  range_mesh_data.create_input_mesh();

  stk::mesh::MetaData &range_meta_data = range_mesh_data.meta_data();
  const stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;
  stk::mesh::Part & range_block         = range_meta_data.declare_part("nodes", node_rank);
  stk::mesh::CellTopology node_top (shards::getCellTopologyData<shards::Node>());
  stk::mesh::set_cell_topology( range_block,  node_top );


  const std::string data_field_name = "Sum_Of_Coordinates";
  ScalarField &range_coord_sum_field = stk::mesh::put_field(
                        range_meta_data.declare_field<ScalarField>(stk::topology::NODE_RANK, data_field_name),
                        range_meta_data.universal_part() );

  range_meta_data.commit();
  range_mesh_data.populate_bulk_data();
  stk::mesh::BulkData &range_bulk_data = range_mesh_data.bulk_data();

  filename = working_directory + domain_mesh;
  domain_mesh_data.add_mesh_database(filename, domain_filetype, stk::io::READ_MESH);
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
                        domain_meta_data.declare_field<ScalarField>(stk::topology::NODE_RANK, data_field_name),
                        domain_meta_data.universal_part() );
  domain_meta_data.commit();

  domain_mesh_data.populate_bulk_data();
  stk::mesh::BulkData &domain_bulk_data = domain_mesh_data.bulk_data();
  stk::mesh::PartVector add_parts(1,&block_skin);
  stk::mesh::skin_mesh(domain_bulk_data, add_parts);
  // For this use case, the domain consists of an axis-aligned
  // bounding box for each 'domain_entity' in the mesh.  The range is a
  // PointBoundingBox3D at the centroid of each 'range_entity'.  The id of the point
  // will be the same as the id of the containing entity.  If the
  // mesh contains solid elements only, and the range_mesh matches the
  // domain_mesh, then the search should return a single box for each
  // point and the id of the box should match the id of the point.

  CartesianField const& range_coord_field  = static_cast<CartesianField const&>( range_mesh_data.get_coordinate_field());
  CartesianField const& domain_coord_field = static_cast<CartesianField const&>(domain_mesh_data.get_coordinate_field());

  stk::mesh::Selector range_nodes = range_meta_data.locally_owned_part();
  stk::mesh::Selector domain_nodes= domain_meta_data.locally_owned_part();

  std::vector<stk::mesh::Entity> domain_entities;
  {
    stk::mesh::get_selected_entities(domain_nodes, domain_bulk_data.buckets(stk::topology::NODE_RANK), domain_entities);
    const size_t num_entities = domain_entities.size();
    for (size_t i = 0; i < num_entities; ++i) {
      const stk::mesh::Entity entity = domain_entities[i];
      double *entity_coordinates = stk::mesh::field_data(domain_coord_field, entity);
      double *entity_coord_sum   = stk::mesh::field_data(domain_coord_sum_field, entity);
      *entity_coord_sum = entity_coordinates[0] + entity_coordinates[1] + entity_coordinates[2];
    }
  }
  std::vector<stk::mesh::Entity> range_entities;
  {
    stk::mesh::get_selected_entities(range_nodes, range_bulk_data.buckets(stk::topology::NODE_RANK), range_entities);
    const size_t num_entities = range_entities.size();
    const double rand_max = RAND_MAX;
    for (size_t i = 0; i < num_entities; ++i) {
      const stk::mesh::Entity entity = range_entities[i];
      double *entity_coord_sum   = stk::mesh::field_data(range_coord_sum_field, entity);
      *entity_coord_sum = rand()/rand_max;
    }
  }

  const double radius=.25;
  const std::vector<stk::mesh::FieldBase*> from_fields(1, &domain_coord_sum_field);
  std::shared_ptr<stk::transfer::STKNode >
    transfer_domain_mesh (new stk::transfer::STKNode(domain_entities, domain_coord_field, from_fields, radius));

  const std::vector<stk::mesh::FieldBase*> to_fields  (1, &range_coord_sum_field);
  std::shared_ptr<stk::transfer::STKNode >
    transfer_range_mesh (new stk::transfer::STKNode(range_entities, range_coord_field, to_fields, radius));


  stk::transfer::GeometricTransfer<
    class stk::transfer::LinearInterpolate<
      class stk::transfer::STKNode,
      class stk::transfer::STKNode
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
                                     range_bulk_data.buckets(stk::topology::NODE_RANK),
                                     range_points);

    const unsigned TONUMPOINTS = range_points.size();

    bool success = true;
    for (unsigned i=0 ; i<TONUMPOINTS; ++i) {
      const  stk::mesh::Entity  entity = range_points[i];

      const double *ToValues = static_cast<double*>(stk::mesh::field_data(range_coord_sum_field, entity));
      const double *ToPoints = static_cast<double*>(stk::mesh::field_data(range_coord_field, entity));

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
