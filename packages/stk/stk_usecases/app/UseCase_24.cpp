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

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/MeshReadWriteUtils.hpp>
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

// Specification for the aggressive gather pointer-field for elements.

typedef mesh::Field<double*,mesh::ElementNode> ElementNodePointerFieldType ;

//--------------------------------
// prototype for the function that will generate the use-case mesh.

void use_case_24_boundary_algorithm(
  mesh::BulkData & bulkData ,
  mesh::Part & quad_block,
  const ScalarIPFieldType & pressure_ip,
  const VectorIPFieldType & mass_flux_ip,
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

  //--------------------------------------------------------------------

  //------------------------------------------------------------------
  // Declare the mesh meta data: element blocks and associated fields
  stk::mesh::fem::FEMMetaData fem_meta;
  fem_meta.FEM_initialize(SpatialDim);
  const mesh::EntityRank element_rank = fem_meta.element_rank();
  const mesh::EntityRank side_rank    = fem_meta.side_rank();
  const mesh::EntityRank edge_rank    = fem_meta.edge_rank();

  //--------------------------------
  // Element-block declarations typically occur when reading the
  // mesh-file meta-data, and thus won't usually appear in application code.
  // Declaring the element blocks and associating an element traits
  // with each element block.

  mesh::Part & allParts = fem_meta.universal_part();
  mesh::Part & hex_io1 = fem_meta.declare_part("block_1", element_rank);

  stk::mesh::fem::CellTopology hex_top(shards::getCellTopologyData<shards::Hexahedron<> >());
  stk::mesh::fem::set_cell_topology( hex_io1, hex_top );

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
    fem_meta.declare_field< VectorFieldType >( "coordinates" );

  VectorFieldType & velocity =
    fem_meta.declare_field< VectorFieldType >( "velocity" );

  ScalarFieldType & pressure =
    fem_meta.declare_field< ScalarFieldType >( "pressure" );

  // add these nodal fields to allParts that have nodes in the universe

  // NOTE SPD 10/24/2008
  // how about put_field_on_part() rather than put_field()
  mesh::put_field( coordinates_field , mesh::fem::FEMMetaData::NODE_RANK , allParts , SpatialDim );
  mesh::put_field( velocity , mesh::fem::FEMMetaData::NODE_RANK , allParts , SpatialDim );
  mesh::put_field( pressure , mesh::fem::FEMMetaData::NODE_RANK , allParts );

  //-----------------------------------------
  // nodal_face_momentum_flux lives on all
  // sideset nodes in the domain
  //-----------------------------------------

  VectorFieldType & nodal_momentum_flux =
    fem_meta.declare_field< VectorFieldType >( "nodal_momentum_flux" );

  mesh::put_field( nodal_momentum_flux , mesh::fem::FEMMetaData::NODE_RANK , quad_io1 );
  mesh::put_field( nodal_momentum_flux , mesh::fem::FEMMetaData::NODE_RANK , quad_io2 );
  mesh::put_field( nodal_momentum_flux , mesh::fem::FEMMetaData::NODE_RANK , quad_io3 );

  //========================================
  // hard code dimensions for now for ip
  //========================================
  const int numBoundaryIps = 4;    // valid for a linear Quad element
  const int numElementSCSIps = 12; // valid for a linear hex CVFEM element

  // momentum_flux_bip lives only on the boundary points

  VectorIPFieldType & momentum_flux_bip =
    fem_meta.declare_field< VectorIPFieldType >( "momentum_flux_bip" );

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

  mesh::put_field(momentum_flux_bip , side_rank , quad_io1 , SpatialDim, numBoundaryIps );
  mesh::put_field(momentum_flux_bip , side_rank , quad_io2 , SpatialDim, numBoundaryIps );
  mesh::put_field(momentum_flux_bip , side_rank , quad_io3 , SpatialDim, numBoundaryIps );

  //-----------------------------------------
  // pressure_ip lives both on boundary
  // and interior integration points
  //-----------------------------------------

  ScalarIPFieldType & pressure_ip =
    fem_meta.declare_field< ScalarIPFieldType >( "pressure_ip" );

  // first on boundary ips
  mesh::put_field( pressure_ip , side_rank , quad_io1 , numBoundaryIps );
  mesh::put_field( pressure_ip , side_rank , quad_io2 , numBoundaryIps );
  mesh::put_field( pressure_ip , side_rank , quad_io3 , numBoundaryIps );

  // second on all interior element ips
  mesh::put_field( pressure_ip , element_rank , hex_io1 , numElementSCSIps );

  //-----------------------------------------
  // mass_flux_ip lives both on boundary
  // and interior integration points
  //-----------------------------------------

  VectorIPFieldType & mass_flux_ip =
    fem_meta.declare_field< VectorIPFieldType >( "mass_flux_ip" );

  // first on boundary ips
  mesh::put_field(mass_flux_ip , side_rank , quad_io1 , SpatialDim, numBoundaryIps );
  mesh::put_field(mass_flux_ip , side_rank , quad_io2 , SpatialDim, numBoundaryIps );
  mesh::put_field(mass_flux_ip , side_rank , quad_io3 , SpatialDim, numBoundaryIps );

  // second on all interior element ips
  mesh::put_field(mass_flux_ip , element_rank , hex_io1 , SpatialDim, numElementSCSIps );

  // Initialize IO system.  Registers all element types and storage
  // types and the exodusII default database type.
  Ioss::Init::Initializer init_db;

  //-------------------------------- Read the exodusII mesh file.
  //Create parts and fields for all io parts, attributes, and fields
  //found in the input mesh file.

  const std::string dbtype("exodusii");
  stk::io::MeshData mesh_data;
  stk::mesh::MetaData & meta_data =
        stk::mesh::fem::FEMMetaData::get_meta_data(fem_meta);
  std::string filename = working_directory + meshName;
  stk::io::create_input_mesh(dbtype, filename, comm, fem_meta, mesh_data);

  // Commit (finalize) the meta data (part and field definitions).
  // Is now ready to be used in the creation and management of mesh bulk data.

  // Output the pressure and velocity fields...
  stk::io::set_field_role(pressure, Ioss::Field::TRANSIENT);
  stk::io::set_field_role(velocity, Ioss::Field::TRANSIENT);

  fem_meta.commit();

  //------------------------------------------------------------------
  // mesh::BulkData bulk data (nodes, faces, node date, etc)
  // conforming to the meta data.

  mesh::BulkData mesh_bulk_data( meta_data , comm );

  // Read bulk data from mesh file
  stk::io::populate_bulk_data(mesh_bulk_data, mesh_data);

  {
    std::vector<unsigned> count ;
    mesh::Selector selector = fem_meta.locally_owned_part() |
                              fem_meta.globally_shared_part();
    count_entities( selector, mesh_bulk_data, count );

    std::cout << "  P" << p_rank << ": Uses {" ;
    std::cout << " Node = " << count[ mesh::fem::FEMMetaData::NODE_RANK] ;
    std::cout << " Edge = " << count[ edge_rank] ;
    std::cout << " Face = " << count[ side_rank] ;
    std::cout << " Elem = " << count[ element_rank] ;
    std::cout << " }" << std::endl ;
    std::cout.flush();
  }

  //------------------------------------------------------------------
  // process only two of the sidesets
  use_case_24_boundary_algorithm( mesh_bulk_data , quad_io1,
                                  pressure_ip, mass_flux_ip,
                                  momentum_flux_bip, pressure, velocity,
                                  nodal_momentum_flux );

  use_case_24_boundary_algorithm( mesh_bulk_data , quad_io3,
                                  pressure_ip, mass_flux_ip,
                                  momentum_flux_bip, pressure, velocity,
                                  nodal_momentum_flux );

  stk::io::create_output_mesh(outputName, comm, mesh_bulk_data, mesh_data);
  stk::io::define_output_fields(mesh_data, fem_meta);
  stk::io::process_output_request(mesh_data, mesh_bulk_data, 0.0);

  return true;
}

//--------------------------------------------------------------------

void use_case_24_boundary_algorithm(
  mesh::BulkData & bulkData ,
  mesh::Part & quad_block,
  const ScalarIPFieldType & pressure_ip,
  const VectorIPFieldType & mass_flux_ip,
  const VectorIPFieldType & momentum_flux_bip,
  const ScalarFieldType   & pressure,
  const VectorFieldType   & velocity,
  const VectorFieldType   & nodal_momentum_flux )
{
  mesh::fem::FEMMetaData &fem_meta = mesh::fem::FEMMetaData::get(bulkData);
  const mesh::EntityRank element_rank = fem_meta.element_rank();
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

      mesh::Entity & face = faceBucket[i] ;

      const mesh::PairIterRelation face_elem = face.relations( element_rank );
      const mesh::PairIterRelation face_nodes = face.relations( mesh::fem::FEMMetaData::NODE_RANK);

      if ( face_elem.size() != 1 ) {
	std::cout << " There appears to be more than one element connected to this boundary face.." << std::endl;
	throw std::exception();
      }
      else {

	mesh::Entity &elem = *face_elem[0].entity();

	//=======================================
	// populate face bip data with values
	//=======================================

	// define some arrays for bip data
	ScalarIPArray p_bip_array(pressure_ip, face);
	VectorIPArray mf_bip_array(mass_flux_ip, face);
	VectorIPArray mom_flux_bip(momentum_flux_bip, face);

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

	for ( unsigned int k = 0; k < face_nodes.size(); ++k )
	{
	  mesh::Entity &node = *face_nodes[k].entity();

	  double * const pNode = field_data( pressure ,node );
	  *pNode = 1.0;

	  VectorArray v_face_node_array(velocity, node);
	  VectorArray mf_face_node_array(nodal_momentum_flux, node);
	  for ( int j = 0; j < nDim; ++j ){
	    v_face_node_array(j) = 1.0;
	    mf_face_node_array(j) = 1.0;
	  }

	}


	//=======================================
	// populate element ip values
	//=======================================

	// define some arrays for bip data
	ScalarIPArray p_ip_array(pressure_ip, elem);
	VectorIPArray mf_ip_array(mass_flux_ip, elem);

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

        const mesh::PairIterRelation elem_nodes = elem.relations( mesh::fem::FEMMetaData::NODE_RANK);
	const int localFaceNumber = face_elem[0].identifier();
        const CellTopologyData *topo = mesh::fem::get_cell_topology(elem).getCellTopologyData();
	std::vector<int> skipNode(elem_nodes.size(),0);

	for ( int j = 0; j < numBoundaryIps; ++j )
	{
	  const int faceNode = topo->side[localFaceNumber].node[j];
	  skipNode[faceNode] = 1;
	}

	// define some arrays for bip data
	for ( unsigned int k = 0; k < elem_nodes.size(); ++k )
	{
	  mesh::Entity &node = *elem_nodes[k].entity();

	  double * const pNode = field_data( pressure ,node );
	  VectorArray v_face_node_array(velocity, node);

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

