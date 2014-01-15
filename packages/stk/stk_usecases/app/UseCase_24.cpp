/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stdexcept>
#include <sstream>
#include <vector>
#include <iostream>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/TopologyDimensions.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <init/Ionit_Initializer.h>

namespace stk {
namespace app {

//----------------------------------------------------------------------
// This file contains the implementation of use-case 13
// The function 'use_case_24_driver' below is the equivalent of 'main'.
// The use case is to define some fields here and there and iterate
// elements to populate this data.
//----------------------------------------------------------------------

enum { SpatialDim = 3 };

//----------------------------------------------------------------------

typedef mesh::Field<double>                    ScalarFieldType ;
typedef mesh::Field<double,mesh::Cartesian>    VectorFieldType ;
typedef mesh::Field<double,mesh::QuadratureTag>   ScalarIPFieldType ;
typedef mesh::Field<double,mesh::Cartesian, mesh::QuadratureTag> VectorIPFieldType ;

typedef mesh::EntityArray<ScalarFieldType> ScalarArray;
typedef mesh::EntityArray<VectorFieldType> VectorArray;
typedef mesh::EntityArray<ScalarIPFieldType> ScalarIPArray;
typedef mesh::EntityArray<VectorIPFieldType> VectorIPArray;

//--------------------------------
// prototype for the function that will generate the use-case mesh.

void use_case_24_boundary_algorithm(
  mesh::BulkData & bulkData ,
  mesh::Part & quad_block,
  const ScalarIPFieldType & side_pressure_ip,
  const ScalarIPFieldType & element_pressure_ip,
  const VectorIPFieldType & side_mass_flux_ip,
  const VectorIPFieldType & element_mass_flux_ip,
  const VectorIPFieldType & momentum_flux_bip,
  const ScalarFieldType   & pressure,
  const VectorFieldType   & velocity,
  const VectorFieldType   & nodal_momentum_flux );

//--------------------------------------------------------------------
//
// main driver for use-case 24: face:element looping
//

bool use_case_24_driver(
  MPI_Comm comm,
  const std::string &working_directory,
  const std::string &meshName,
  const std::string &outputName )
{
  const unsigned p_rank = parallel_machine_rank( comm );

  //--------------------------------------------------------------------

  if ( ! p_rank ) {
    std::cout << "stk_mesh Use Case #24, begin" << std::endl
              << "  Number Processes = " << parallel_machine_size( comm )
              << std::endl ;
    std::cout.flush();
  }

  //-------------------------------- Read the exodusII mesh file.
  //Create parts and fields for all io parts, attributes, and fields
  //found in the input mesh file.
  const std::string dbtype("exodusii");
  std::string filename = working_directory + meshName;
  stk::io::StkMeshIoBroker mesh_data(comm);
  mesh_data.open_mesh_database(filename, dbtype, stk::io::READ_MESH);
  mesh_data.create_input_mesh();

  //------------------------------------------------------------------
  // Declare the mesh meta data: element blocks and associated fields
  stk::mesh::MetaData &fem_meta = mesh_data.meta_data();
  const mesh::EntityRank element_rank = stk::mesh::MetaData::ELEMENT_RANK;
  const mesh::EntityRank side_rank    = fem_meta.side_rank();
  const mesh::EntityRank edge_rank    = mesh::MetaData::EDGE_RANK;

  //--------------------------------
  // Element-block declarations typically occur when reading the
  // mesh-file meta-data, and thus won't usually appear in application code.
  // Declaring the element blocks and associating an element traits
  // with each element block.

  mesh::Part & allParts = fem_meta.universal_part();
  mesh::Part & hex_io1 = fem_meta.declare_part("block_1", element_rank);

  stk::mesh::CellTopology hex_top(shards::getCellTopologyData<shards::Hexahedron<> >());
  stk::mesh::set_cell_topology( hex_io1, hex_top );

    //--------------------------------
    // Surface-block declarations
  mesh::Part & quad_io1 = fem_meta.declare_part("surface_1", side_rank);
  mesh::Part & quad_io2 = fem_meta.declare_part("surface_2", side_rank);
  mesh::Part & quad_io3 = fem_meta.declare_part("surface_3", side_rank);

  //========================================
  // Define where the field lives and its' size;
  //========================================

  //-----------------------------------------
  // coordinates_field, velocity and pressure
  // live on all nodes in the domain
  //-----------------------------------------

  VectorFieldType & coordinates_field =
    fem_meta.declare_field< VectorFieldType >(stk::topology::NODE_RANK, "coordinates" );

  VectorFieldType & velocity =
    fem_meta.declare_field< VectorFieldType >(stk::topology::NODE_RANK, "velocity" );

  ScalarFieldType & pressure =
    fem_meta.declare_field< ScalarFieldType >(stk::topology::NODE_RANK, "pressure" );

  // add these nodal fields to allParts that have nodes in the universe

  // NOTE SPD 10/24/2008
  // how about put_field_on_part() rather than put_field()
  mesh::put_field( coordinates_field , allParts , SpatialDim );
  mesh::put_field( velocity , allParts , SpatialDim );
  mesh::put_field( pressure , allParts );

  //-----------------------------------------
  // nodal_face_momentum_flux lives on all
  // sideset nodes in the domain
  //-----------------------------------------

  VectorFieldType & nodal_momentum_flux =
    fem_meta.declare_field< VectorFieldType >(stk::topology::NODE_RANK, "nodal_momentum_flux" );

  mesh::put_field( nodal_momentum_flux , quad_io1 );
  mesh::put_field( nodal_momentum_flux , quad_io2 );
  mesh::put_field( nodal_momentum_flux , quad_io3 );

  //========================================
  // hard code dimensions for now for ip
  //========================================
  const int numBoundaryIps = 4;    // valid for a linear Quad element
  const int numElementSCSIps = 12; // valid for a linear hex CVFEM element

  // momentum_flux_bip lives only on the boundary points

  stk::topology::rank_t sideRank = static_cast<stk::topology::rank_t>(side_rank);
  VectorIPFieldType & momentum_flux_bip =
    fem_meta.declare_field< VectorIPFieldType >(sideRank, "momentum_flux_bip" );

  // NOTE SPD 10/22/2008
  // I am not sure that I like all of these specialty put_field() methods
  // Perhaps it would be better to define a set of attributes that a field
  // requires like dimension, and where it lives. Then, we can have some
  // sort of method on the field that sets it...
  // so, below would look like:

  /*

  mom_flux_bip.spatial_dimension(SpatialDim);
  mom_flux_bip.size_of_where_it_lives(numBoundaryIps);

  */

  // This would also help in extracting the size of things rather than remembering
  // what index should be passed to dimension()

  mesh::put_field(momentum_flux_bip , quad_io1 , SpatialDim, numBoundaryIps );
  mesh::put_field(momentum_flux_bip , quad_io2 , SpatialDim, numBoundaryIps );
  mesh::put_field(momentum_flux_bip , quad_io3 , SpatialDim, numBoundaryIps );

  //-----------------------------------------
  // pressure_ip lives both on boundary
  // and interior integration points
  //-----------------------------------------

  ScalarIPFieldType & side_pressure_ip =
    fem_meta.declare_field< ScalarIPFieldType >(sideRank, "side_pressure_ip" );

  // first on boundary ips
  mesh::put_field( side_pressure_ip , quad_io1 , numBoundaryIps );
  mesh::put_field( side_pressure_ip , quad_io2 , numBoundaryIps );
  mesh::put_field( side_pressure_ip , quad_io3 , numBoundaryIps );

  ScalarIPFieldType & element_pressure_ip =
    fem_meta.declare_field< ScalarIPFieldType >(stk::topology::ELEMENT_RANK, "element_pressure_ip" );
  // second on all interior element ips
  mesh::put_field( element_pressure_ip , hex_io1 , numElementSCSIps );

  //-----------------------------------------
  // mass_flux_ip lives both on boundary
  // and interior integration points
  //-----------------------------------------

  VectorIPFieldType & side_mass_flux_ip =
    fem_meta.declare_field< VectorIPFieldType >(sideRank, "side_mass_flux_ip" );

  // first on boundary ips
  mesh::put_field(side_mass_flux_ip , quad_io1 , SpatialDim, numBoundaryIps );
  mesh::put_field(side_mass_flux_ip , quad_io2 , SpatialDim, numBoundaryIps );
  mesh::put_field(side_mass_flux_ip , quad_io3 , SpatialDim, numBoundaryIps );

  VectorIPFieldType & element_mass_flux_ip =
    fem_meta.declare_field< VectorIPFieldType >(stk::topology::ELEMENT_RANK, "element_mass_flux_ip" );
  // second on all interior element ips
  mesh::put_field(element_mass_flux_ip , hex_io1 , SpatialDim, numElementSCSIps );

  // Commit (finalize) the meta data (part and field definitions).
  // Is now ready to be used in the creation and management of mesh bulk data.
  fem_meta.commit();

  //------------------------------------------------------------------
  // mesh::BulkData bulk data (nodes, faces, node date, etc)
  // conforming to the meta data.


  // Read bulk data from mesh file
  mesh_data.populate_bulk_data();
  mesh::BulkData &mesh_bulk_data = mesh_data.bulk_data();

  {
    std::vector<unsigned> count ;
    mesh::Selector selector = fem_meta.locally_owned_part() |
                              fem_meta.globally_shared_part();
    count_entities( selector, mesh_bulk_data, count );

    std::cout << "  P" << p_rank << ": Uses {" ;
    std::cout << " Node = " << count[ mesh::MetaData::NODE_RANK] ;
    std::cout << " Edge = " << count[ edge_rank] ;
    std::cout << " Face = " << count[ side_rank] ;
    std::cout << " Elem = " << count[ element_rank] ;
    std::cout << " }" << std::endl ;
    std::cout.flush();
  }

  //------------------------------------------------------------------
  // process only two of the sidesets
  use_case_24_boundary_algorithm( mesh_bulk_data , quad_io1,
                                  side_pressure_ip, element_pressure_ip, side_mass_flux_ip, element_mass_flux_ip,
                                  momentum_flux_bip, pressure, velocity,
                                  nodal_momentum_flux );

  use_case_24_boundary_algorithm( mesh_bulk_data , quad_io3,
                                  side_pressure_ip, element_pressure_ip, side_mass_flux_ip, element_mass_flux_ip,
                                  momentum_flux_bip, pressure, velocity,
                                  nodal_momentum_flux );

  size_t result_output_index = mesh_data.create_output_mesh(outputName, stk::io::WRITE_RESULTS);

  // Output the pressure and velocity fields...
  mesh_data.add_field(result_output_index, pressure);
  mesh_data.add_field(result_output_index, velocity);

  mesh_data.process_output_request(result_output_index,0.0);

  return true;
}

//--------------------------------------------------------------------

void use_case_24_boundary_algorithm(
  mesh::BulkData & bulkData ,
  mesh::Part & quad_block,
  const ScalarIPFieldType & side_pressure_ip,
  const ScalarIPFieldType & element_pressure_ip,
  const VectorIPFieldType & side_mass_flux_ip,
  const VectorIPFieldType & element_mass_flux_ip,
  const VectorIPFieldType & momentum_flux_bip,
  const ScalarFieldType   & pressure,
  const VectorFieldType   & velocity,
  const VectorFieldType   & nodal_momentum_flux )
{
  mesh::MetaData &fem_meta = mesh::MetaData::get(bulkData);
  const mesh::EntityRank side_rank    = fem_meta.side_rank();

  // create intersection of faces and locally owned
  mesh::Selector select_owned_block = quad_block & fem_meta.locally_owned_part();

  // Get vector of buckets ( entities and field data)
  // for which the sides (quad_block) are all locally owned.

  const std::vector<mesh::Bucket*> & buckets = bulkData.buckets( side_rank );

  unsigned count = 0 ;

  for ( std::vector<mesh::Bucket *>::const_iterator
          ik = buckets.begin() ; ik != buckets.end() ; ++ik ) if ( select_owned_block( **ik ) ) {

      const mesh::Bucket & faceBucket = **ik ;

      // Number of sides in this bucket of sides and side field data

      const int number = faceBucket.size();

      count += number ;

      for ( int i = 0 ; i < number ; ++i ) {

        mesh::Entity const *face_elem = faceBucket.begin_elements(i);
        int num_face_elems = faceBucket.num_elements(i);
        mesh::ConnectivityOrdinal const *face_elem_ordinal = faceBucket.begin_element_ordinals(i);
        mesh::Entity const *face_nodes = faceBucket.begin_nodes(i);
        int num_face_nodes = faceBucket.num_nodes(i);

        if ( num_face_elems != 1 ) {
          std::cout << " There appears to be more than one element connected to this boundary face.." << std::endl;
          throw std::exception();
        }
        else {

          mesh::Entity elem = face_elem[0];
          const mesh::Bucket &elem_bucket = bulkData.bucket(elem);
          const mesh::Ordinal elem_bordinal = bulkData.bucket_ordinal(elem);

          //=======================================
          // populate face bip data with values
          //=======================================

          // define some arrays for bip data
          ScalarIPArray p_bip_array(side_pressure_ip, faceBucket, i);
          VectorIPArray mf_bip_array(side_mass_flux_ip, faceBucket, i);
          VectorIPArray mom_flux_bip(momentum_flux_bip, faceBucket, i);

          // FIXME SPD 10/22/2008
          // I do not like to remember specific indexes, e.g. what dimension(0),
          // dimension(1), etc. means
          // For example dimension(o) for ScalarIPFieldType provides numBoundaryIps
          // whereas dimension(0) for a VectorIPFieldType, provides the SpatialDim..
          // Again, I might like to see some specific method off of the field that returns
          // things such as SpatialDim, and HowManyThereAre, e.g., numBoundaryIps
          const int numBoundaryIps = p_bip_array.dimension(0);
          const int nDim = mom_flux_bip.dimension(0);
          for ( int k = 0; k < numBoundaryIps; ++k )
          {
            p_bip_array(k) = 1.0;
            for (int j = 0; j < nDim; ++j){
              mf_bip_array(j,k) = 1.0;
              mom_flux_bip(j,k) = 1.0;
            }

          }

          //=======================================
          // populate face nodal data with values
          //=======================================
          for ( int k = 0; k < num_face_nodes; ++k )
          {
            mesh::Entity node = face_nodes[k];
            const mesh::Bucket &node_bucket = bulkData.bucket(node);
            const mesh::Ordinal node_bordinal = bulkData.bucket_ordinal(node);

            double * const pNode = bulkData.field_data( pressure ,node );
            *pNode = 1.0;

            VectorArray v_face_node_array(velocity, node_bucket, node_bordinal);
            VectorArray mf_face_node_array(nodal_momentum_flux, node_bucket, node_bordinal);
            for ( int j = 0; j < nDim; ++j ){
              v_face_node_array(j) = 1.0;
              mf_face_node_array(j) = 1.0;
            }

          }


          //=======================================
          // populate element ip values
          //=======================================

          // define some arrays for bip data
          ScalarIPArray p_ip_array(element_pressure_ip, elem_bucket, elem_bordinal);
          VectorIPArray mf_ip_array(element_mass_flux_ip, elem_bucket, elem_bordinal);

          // FIXME SPD 10/22/2008
          // I do not like to remember specific indexes, e.g. dimension(0) means,
          // for a ScalarIPFieldType, numBoundaryIps whereas for a VectorIPFieldType,
          // dimension(0) provides the SpatialDim..
          const int numInteriorIps = p_ip_array.dimension(0);
          const int nnDim = mf_ip_array.dimension(0);
          for ( int k = 0; k < numInteriorIps; ++k )
          {
            p_ip_array(k) = 1.0;
            for (int j = 0; j < nnDim; ++j){
              mf_ip_array(j,k) = 2.0;
            }

          }

          //================================================
          // populate off elem nodes data with values
          //================================================


          mesh::Entity const *elem_nodes = bulkData.begin_nodes(elem);
          int num_elem_nodes = bulkData.num_nodes(elem);
          const int localFaceNumber = face_elem_ordinal[0];
          const CellTopologyData *topo = mesh::get_cell_topology(elem_bucket).getCellTopologyData();
          std::vector<int> skipNode(num_elem_nodes, 0);

          for ( int j = 0; j < numBoundaryIps; ++j )
          {
            const int faceNode = topo->side[localFaceNumber].node[j];
            skipNode[faceNode] = 1;
          }

          // define some arrays for bip data
          for (int k = 0; k < num_elem_nodes; ++k )
          {
            mesh::Entity node = elem_nodes[k];
            const mesh::Bucket &node_bucket = bulkData.bucket(node);
            const mesh::Ordinal node_bordinal = bulkData.bucket_ordinal(node);

            double * const pNode = bulkData.field_data( pressure, node );
            VectorArray v_face_node_array(velocity, node_bucket, node_bordinal);

            if ( !skipNode[k]  )
            {
              *pNode = 2.0;

              for ( int j = 0; j < nDim; ++j )
                v_face_node_array(j) = 2.0;
            }
          }
        }
      }
    }

  std::cout << "use_case_24_boundary_algorithm applied to '"
            << quad_block.name()
            << "' processed " << count << " faces" << std::endl ;
}

} // end namespace stk
} // end namespace app

