/*--------------------------------------------------------------------*/
/*    Copyright 2009, 2010, 2011 Sandia Corporation.                        */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Shards_BasicTopologies.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_mesh/fem/FEMHelpers.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/Comm.hpp>

#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

#include <stk_percept/TopologyVerifier.hpp>
#include <stk_percept/GeometryVerifier.hpp>
#include <stk_percept/mesh/gen/SweepMesher.hpp>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <math.h>

//static  bool isParallel = true; // FIXME
static  bool isParallel()
{
  stk::ParallelMachine parallel_machine = MPI_COMM_WORLD ;
  int mpi_size = stk::parallel_machine_size(parallel_machine);
  return mpi_size > 0;
}

namespace stk_percept_unit
{
void use_encr_case_1_driver( MPI_Comm comm );

int myMain()
{
  stk::ParallelMachine parallel_machine = MPI_COMM_WORLD ;

  if (isParallel() ) return 0;

  use_encr_case_1_driver( parallel_machine );

  // Need to check whether all use-cases actually did run ok.  If they didn't, then don't
  // print the following "OK" string.  (The "OK" string is for the runtest tool to
  // determine whether the test gets a "pass" or "fail" status.)
  std::cout << "unit test cases ran OK" << std::endl;

  return 0;
}

#ifndef REDS
STKUNIT_UNIT_TEST(topo, test1)
{
  /* %TRACE[OFF]% */  /* %TRACE% */

  // totally different
  STKUNIT_EXPECT_FALSE( 0 );


  if (!isParallel())
  {
    myMain();
  }
}
#endif

} // namespace SEncr


namespace stk
{
namespace percept
{
namespace unit_tests
{

#ifndef REDS
STKUNIT_UNIT_TEST(topo, testCrossedElems)
{

  if (isParallel() ) return;
  using namespace stk::percept::util;
  using namespace stk::percept::interface_table;
  /* %TRACE[OFF]% */  /* %TRACE% */
  stk::ParallelMachine parallel_machine = MPI_COMM_WORLD ;
  bool verbose = false;

  // totally different
  STKUNIT_EXPECT_FALSE( 0 );

  stk::percept::TopologyVerifier topoVerifier;

  unsigned line2Elems[] = {
    0,1,
    1,2,
    2,3,
    3,4,
    0,0  // sentinel
  };
  enum {
    numElemsL2 = sizeof(line2Elems)/(2*sizeof(line2Elems[0])) - 1  // drop sentinel
  };
  if(verbose) std::cout << "numElemsL2= " << numElemsL2 << std::endl;

  boost::array<double,3> coordsLine[] = {
    //{{0,0,0}}, {{1,0,0}}, {{2,0,0}}, {{3,0,0}}, {{4,0,0}}
    {{0,0,0}}, {{1.1,0,0}}, {{2.22,0,0}}, {{3.333,0,0}}, {{4.4444,0,0}}
  };
  unsigned numNodesLine = sizeof(coordsLine)/sizeof(coordsLine[0]);
  if(verbose) std::cout << "numNodesLine= " << numNodesLine << std::endl;

  //------ create a valid heterogeneous mesh
  // line2 mesh
  SweepMesher tp2;
  tp2.dump(verbose);

  tp2.initNodes(coordsLine, numNodesLine);
  tp2.initElems(shards_Line_2, line2Elems, numElemsL2);
  if(verbose) std::cout << "line2 mesh\n";
  tp2.dump();

  // sweep to make a quad mesh from line mesh
  boost::array< double, 3> dir = {{0,1.234,0}};
  TransformDir xf( dir );
  std::vector<Transform *> xforms(3, &xf);
  //xforms.push_back(xf);
  tp2.sweep(shards_Line_2, shards_Quadrilateral_4,  xforms);
  if(verbose) std::cout << "after line to quad sweep\n";
  tp2.dump();

  SweepMesher quadMeshCopy;
  quadMeshCopy.CopyFromBasicMesh(tp2);

  // break one of the quads into tris
  tp2.breakElement<shards_Quadrilateral_4, shards_Triangle_3>(0);
  //tp2.breakElement<shards_Quadrilateral_4, shards_Triangle_3>(1);
  if(verbose) std::cout << "after break\n";
  tp2.dump();
  tp2.stkMeshCreate(parallel_machine);
  tp2.writeSTKMesh("topo-goodQuadTri.e");

  // verify valid topology
  bool isBad = false;
  if(verbose) std::cout << "verify this is a good topology " << std::endl;
  isBad = topoVerifier.isTopologyBad( *tp2.getBulkData() );
  STKUNIT_EXPECT_FALSE(isBad);

  //------- a bad topology with a duplicated node
  tp2.initialize();
  tp2.dump(verbose);
  tp2.CopyFromBasicMesh(quadMeshCopy);
  if(verbose) std::cout << "before creating invalid mesh\n";
  tp2.dump();

  //ShardsInterfaceTable::s_elemInfo[shards_Quadrilateral_4].vertex_count;
  // { 0 1 6 5  }
  // { 1 2 7 6  }
  int whichElem = 1;

  STKUNIT_EXPECT_EQ(tp2.m_elems[shards_Quadrilateral_4][ 4*whichElem + 3 ], 6u);
  // create a duplicated node
  tp2.m_elems[shards_Quadrilateral_4][ 4*whichElem + 3 ] = 7;

  if(verbose) std::cout << "after creating invalid mesh with duplicate node\n";
  tp2.dump();
  tp2.stkMeshCreate(parallel_machine);
  tp2.writeSTKMesh("topo-badQuadDupl.e");

  // verify bad topology
  isBad = topoVerifier.isTopologyBad( *tp2.getBulkData() );
  STKUNIT_EXPECT_TRUE(isBad);

  //------ create a bad topology with crossed elements
  tp2.initialize();
  tp2.dump(verbose);
  tp2.CopyFromBasicMesh(quadMeshCopy);
  if(verbose) std::cout << "before creating invalid mesh\n";
  tp2.dump();

  //ShardsInterfaceTable::s_elemInfo[shards_Quadrilateral_4].vertex_count;
  // { 0 1 6 5  }
  // { 1 2 7 6  }
  whichElem = 1;

  STKUNIT_EXPECT_EQ(tp2.m_elems[shards_Quadrilateral_4][ 4*whichElem + 3 ], 6u);
  // create a crossed mesh
  tp2.m_elems[shards_Quadrilateral_4][ 4*whichElem + 3 ] = 5;

  if(verbose) std::cout << "after creating invalid mesh\n";
  tp2.dump();
  tp2.stkMeshCreate(parallel_machine);
  tp2.writeSTKMesh("topo-badQuadCrossed.e");

  // verify bad topology
  isBad = topoVerifier.isTopologyBad( *tp2.getBulkData() );
  STKUNIT_EXPECT_TRUE(isBad);
}

STKUNIT_UNIT_TEST(geom, geomPrints)
{
  if (isParallel() ) return;
  using namespace stk::percept::util;
  using namespace stk::percept::interface_table;

  /* %TRACE[OFF]% */  /* %TRACE% */
  stk::ParallelMachine parallel_machine = MPI_COMM_WORLD ;
  bool verbose = false;

  // totally different
  STKUNIT_EXPECT_FALSE( 0 );

  stk::percept::GeometryVerifier geomVerifier(false); // true for more dumps

  unsigned line2Elems[] = {
    0,1,
    1,2,
    2,3,
    3,4,
    0,0  // sentinel
  };
  enum {
    numElemsL2 = sizeof(line2Elems)/(2*sizeof(line2Elems[0])) - 1  // drop sentinel
  };
  if(verbose) std::cout << "numElemsL2= " << numElemsL2 << std::endl;

  boost::array<double,3> coordsLine[] = {
    //{{0,0,0}}, {{1,0,0}}, {{2,0,0}}, {{3,0,0}}, {{4,0,0}}
    {{0,0,0}}, {{1.1,0,0}}, {{2.22,0,0}}, {{3.333,0,0}}, {{4.4444,0,0}}
  };
  unsigned numNodesLine = sizeof(coordsLine)/sizeof(coordsLine[0]);
  if(verbose) std::cout << "numNodesLine= " << numNodesLine << std::endl;

  //------ create a valid heterogeneous mesh
  // line2 mesh
  SweepMesher tp2;
  tp2.dump(verbose);

  tp2.initNodes(coordsLine, numNodesLine);
  tp2.initElems(shards_Line_2, line2Elems, numElemsL2);
  if(verbose) std::cout << "line2 mesh\n";
  tp2.dump();

  // sweep to make a quad mesh from line mesh
  boost::array< double,3> dir = {{0,1.234,0}};
  TransformDir xf ( dir );
  std::vector<Transform *> xforms(3, &xf);
  //xforms.push_back(xf);
  tp2.sweep(shards_Line_2, shards_Quadrilateral_4,  xforms);

  //if(true || verbose) std::cout << "after line to quad sweep\n";
  if( verbose) std::cout << "after line to quad sweep\n";
  //tp2.dump(true);
  tp2.dump();

  SweepMesher quadMeshCopy;
  quadMeshCopy.CopyFromBasicMesh(tp2);

  // break one of the quads into tris
  tp2.breakElement<shards_Quadrilateral_4, shards_Triangle_3>(0);
  //tp2.breakElement<shards_Quadrilateral_4, shards_Triangle_3>(1);
  if(verbose) std::cout << "after break\n";
  tp2.dump();
  tp2.stkMeshCreate(parallel_machine);
  tp2.writeSTKMesh("geom-goodQuadTri.e");

  // verify valid geometry
  bool isBad = false;
  if(verbose) std::cout << "verify this is a good geometry " << std::endl;
  isBad = geomVerifier.isGeometryBad( *tp2.getBulkData() );
  STKUNIT_EXPECT_FALSE(isBad);

  /////////////// path test 3
  tp2.initialize();
  tp2.CopyFromBasicMesh(quadMeshCopy);
  double rad = 10.0;
  boost::array<double, 3> dirT = {{0,rad,0}};
  xf = TransformDir( dirT );
  tp2.transform(xf);

  VectorOfCoord path3;
  VectorOfCoord dir3;

  boost::array<double, 3> pt0 =  {{0,rad,0}} ;
  path3.push_back(pt0);
  boost::array<double, 3> dr0 =  {{0,0,1}} ;
  dir3.push_back(dr0);

  unsigned ntheta = 8;
  for (unsigned ith = 1; ith <= ntheta; ith++)
  {
    double th = M_PI*((double)ith)/((double)ntheta);
    boost::array<double, 3> pt1 = {{0, rad*cos(th), rad*sin(th)}};
    boost::array<double, 3> dr1 = {{0, -sin(th), cos(th)}};;
    path3.push_back(pt1);
    dir3.push_back(dr1);
  }
  tp2.sweep(path3, dir3);

  tp2.stkMeshCreate(parallel_machine);
  tp2.writeSTKMesh("geom-all-hex-path3.e");
  geomVerifier = GeometryVerifier(false);
  geomVerifier.isGeometryBad(*tp2.getBulkData() );

  // break path3 of the hexes into tets
  tp2.breakAllElements<shards_Hexahedron_8, shards_Tetrahedron_4>();
  tp2.stkMeshCreate(parallel_machine);
  tp2.writeSTKMesh("geom-all-hex-tet-path3.e");
  geomVerifier.isGeometryBad(*tp2.getBulkData() );

}

STKUNIT_UNIT_TEST(geom, geomEqui)
{
  if (isParallel() ) return;
  using namespace stk::percept::util;
  using namespace stk::percept::interface_table;

  /* %TRACE[OFF]% */  /* %TRACE% */
  stk::ParallelMachine parallel_machine = MPI_COMM_WORLD ;
  bool verbose = false;

  stk::percept::GeometryVerifier geomVerifier(false); // true for more dumps

  unsigned tetElems[] = {
    0,1,2,3
  };
  enum {
    numElems = 1
  };

  boost::array<double, 3> coordsTet[] = {
    {{0,0,0}}, {{1,0,0}}, {{0.5,sqrt(3.)/2.,0}}, {{0.5,sqrt(3.)/6.,sqrt(6.)/3.}}
  };
  unsigned numNodesTet = sizeof(coordsTet)/sizeof(coordsTet[0]);
  if(verbose) std::cout << "numNodesTet= " << numNodesTet << std::endl;

  //------ create a valid mesh of one tet element
  //
  SweepMesher tp2;
  tp2.dump(verbose);

  tp2.initNodes(coordsTet, numNodesTet);
  tp2.initElems(shards_Tetrahedron_4, tetElems, numElems);
  if(verbose) std::cout << "tet4 mesh\n";
  tp2.dump();

  tp2.stkMeshCreate(parallel_machine);
  tp2.writeSTKMesh("equi-tet.e");
  geomVerifier = GeometryVerifier(false);
  geomVerifier.isGeometryBad(*tp2.getBulkData(), true );

  //------ scale the mesh
  double sf= 4.0;
  boost::array<double, 3> coordsTetScaled[] = {
    {{0,0,0}}, {{sf*1,0,0}}, {{sf*0.5,sf*sqrt(3.)/2.,0}}, {{sf*0.5, sf*sqrt(3.)/6., sf*sqrt(6.)/3.}}
  };

  tp2.initialize();
  tp2.initNodes(coordsTetScaled, numNodesTet);
  tp2.initElems(shards_Tetrahedron_4, tetElems, numElems);
  if(verbose) std::cout << "tet4 mesh scaled\n";
  tp2.dump();

  tp2.stkMeshCreate(parallel_machine);
  tp2.writeSTKMesh("equi-tet-scaled.e");
  geomVerifier = GeometryVerifier(false);
  geomVerifier.isGeometryBad(*tp2.getBulkData(), true );

}
#endif


}
}
}



using namespace stk ;


namespace stk_percept_unit {

enum { SpatialDim = 3 };

//----------------------------------------------------------------------

typedef mesh::Field<double>                    ScalarFieldType ;
typedef mesh::Field<double,mesh::Cartesian>    VectorFieldType ;

// Specification for the aggressive gather pointer-field for elements.

typedef mesh::Field<double*,mesh::ElementNode> ElementNodePointerFieldType ;

//--------------------------------
// prototype for the function that will generate the use-case mesh.
// copied from stk::mesh use cases

void use_encr_case_1_generate_mesh(
  mesh::BulkData & mesh ,
  const unsigned N[] ,
  const VectorFieldType & node_coord ,
  const ElementNodePointerFieldType & elem_node_coord ,
  mesh::Part & hex_block );



//--------------------------------------------------------------------
//
// main driver for use-case 13: heterogeneous element mesh.
//

void use_encr_case_1_driver( MPI_Comm comm )
{
  const unsigned p_rank = parallel_machine_rank( comm );

  //reset_malloc_stats();

  const unsigned box_size[3] = { 1 , 1 , 1 };

  //--------------------------------------------------------------------

  if ( ! p_rank ) {
    std::cout << "stk_percept_utest begin" << std::endl
              << "  Number Processes = " << parallel_machine_size( comm )
              << std::endl ;
    std::cout.flush();
  }

  //--------------------------------------------------------------------

  {
    //------------------------------------------------------------------
    // Declare the mesh meta data: element blocks and associated fields

    mesh::fem::FEMMetaData mesh_meta_data(3, mesh::fem::entity_rank_names(3) );

    //--------------------------------
    // Element-block declarations typically occur when reading the
    // mesh-file meta-data, and thus won't usually appear in application code.
    // Declaring the element blocks and associating an element traits
    // with each element block.

    mesh::Part & universal = mesh_meta_data.universal_part();
    mesh::Part & block_hex = mesh_meta_data.declare_part("block_1", mesh_meta_data.element_rank());

    /// set cell topology for the part block_1
    stk::mesh::fem::set_cell_topology< shards::Hexahedron<8>  >( block_hex );

    //--------------------------------
    // Declare coordinates field on all nodes with 3D:

    VectorFieldType & coordinates_field =
      mesh_meta_data.declare_field< VectorFieldType >( "coordinates" );

    stk::mesh::put_field(
      coordinates_field , mesh::fem::FEMMetaData::NODE_RANK , universal , SpatialDim );

    //--------------------------------

    //--------------------------------
    // Declare an aggressive "gather" field which is an
    // array of pointers to the element's nodes' coordinate field data.
    // The declaration specifies:
    //
    //     double * elem_node_coord[number_of_nodes]

    ElementNodePointerFieldType & elem_node_coord =
      mesh_meta_data.
      declare_field< ElementNodePointerFieldType >( "elem_node_coord" );

    // Declare that the 'elem_node_coord' pointer field data
    // points to the 'coordinates_field' data on the nodes.

    mesh_meta_data.declare_field_relation(
      elem_node_coord ,
      stk::mesh::fem::get_element_node_stencil(3) ,
      coordinates_field );

    // Declare the size of the aggressive "gather" field
    //     double * elem_node_coord[ size = number_of_nodes ]
    // is the number of nodes of the elements.
    // This size is different for each element block.

    stk::mesh::put_field(
                         elem_node_coord , mesh_meta_data.element_rank() , block_hex , shards::Hexahedron<8> ::node_count );

    //--------------------------------
    // Commit (finalize) the meta data.  Is now ready to be used
    // in the creation and management of mesh bulk data.

    mesh_meta_data.commit();

    //------------------------------------------------------------------
    // mesh::BulkData bulk data conforming to the meta data.

    mesh::BulkData mesh_bulk_data( stk::mesh::fem::FEMMetaData::get_meta_data(mesh_meta_data) , MPI_COMM_WORLD );

    // In a typical app, the mesh would be read from file at this point.
    // But in this use-case, we generate the mesh and initialize
    // field data to use-case defined values.

    use_encr_case_1_generate_mesh(
      mesh_bulk_data ,
      box_size ,
      coordinates_field ,
      elem_node_coord ,
      block_hex );

    //use_encr_case_1_generate_sides( mesh_bulk_data , false );


    {
      std::vector<unsigned> count ;
      mesh::Selector selector(mesh_meta_data.globally_shared_part());
      count_entities( selector, mesh_bulk_data, count );

      std::cout << "  P" << p_rank << ": Uses {" ;
      std::cout << " Node = " << count[ 0 ] ;
      std::cout << " Edge = " << count[ 1 ] ;
      std::cout << " Face = " << count[ 2 ] ;
      std::cout << " Elem = " << count[ 3 ] ;
      std::cout << " }" << std::endl ;
      std::cout.flush();
    }

    //------------------------------------------------------------------

#ifdef USE_GNU_MALLOC_HOOKS
    if (parallel_machine_rank(comm) == 0) {
      double net_alloc = alloc_MB() - freed_MB();
      std::cout << "Mesh creation:" << "\n   Total allocated: "
                << alloc_MB()<<"MB in "<<alloc_blks() << " blocks."
                << "\n   Total freed: " << freed_MB() << "MB in "
                << freed_blks() << " blocks."
                << "\n   Net allocated: "<<net_alloc << "MB."<<std::endl;
    }
#endif

    //------------------------------------------------------------------



    //------------------------------------------------------------------

  }
}

//--------------------------------------------------------------------


// PairIterRelation is a typedef for:
//
//   PairIter< std::vector<Relation>::const_iterator >
//
// Where PairIter is a wrapper for a pair of iterators.


//--------------------------------------------------------------------
//----------------------------------------------------------------------


} // namespace stk_percept_unit

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#include <stk_util/unit_test_support/GeneratedMesh.hpp>

namespace stk_percept_unit {

void use_encr_case_1_generate_mesh(
  mesh::BulkData & mesh ,
  const unsigned N[] ,
  const VectorFieldType & node_coord_field ,
  const ElementNodePointerFieldType & elem_node_coord_field,
  mesh::Part & hex_block )
{
  //EXCEPTWATCH;
  mesh.modification_begin();

  const unsigned parallel_size = mesh.parallel_size();
  const unsigned parallel_rank = mesh.parallel_rank();

  stk::percept::TopologyVerifier topoVerifier;

  double t = 0 ;
  size_t num_hex = 0 ;
  size_t num_shell = 0 ;
  size_t num_nodes = 0 ;
  size_t num_block = 0 ;
  int error_flag = 0 ;

  try {

    stk::io::util::GeneratedMesh gmesh( N[0], N[1], N[2], parallel_size, parallel_rank );

    num_nodes = gmesh.node_count_proc();
    num_block = gmesh.block_count();

    t = wall_time();

    std::vector<int> node_map( num_nodes , 0 );

    // global nodes
    gmesh.node_map( node_map );

    // create a bad element
    node_map[1] = node_map[0];

    //TODO CPPUNIT_ASSERT_EQUAL( num_nodes , node_map.size() );

    {
      for ( size_t i = 1 ; i <= num_block ; ++i ) {
        const size_t                     num_elem = gmesh.element_count_proc(i);
        const std::pair<std::string,int> top_info = gmesh.topology_type(i);

        std::vector<int> elem_map( num_elem , 0 );
        gmesh.element_map( i, elem_map );

        std::vector<int> elem_conn( num_elem * top_info.second );

        gmesh.connectivity( i , elem_conn );

        //std::cout << "top_info.first= " <<top_info.first << std::endl;

        if ( top_info.second == 8 ) {

          for ( size_t j = 0 ; j < num_elem ; ++j ) {

            const int * const local_node_id = & elem_conn[ j * 8 ] ;

            const stk::mesh::EntityId node_id[8] = {
              node_map[ local_node_id[0] - 1 ] ,
              node_map[ local_node_id[1] - 1 ] ,
              node_map[ local_node_id[2] - 1 ] ,
              // TODO
              node_map[ local_node_id[3] - 1 ] ,
              node_map[ local_node_id[4] - 1 ] ,
              node_map[ local_node_id[5] - 1 ] ,
              node_map[ local_node_id[6] - 1 ] ,
              node_map[ local_node_id[7] - 1 ]
            };

            const stk::mesh::EntityId elem_id = elem_map[ j ];

            stk::mesh::fem::declare_element( mesh , hex_block , elem_id , node_id );

            mesh::Entity * const elem = mesh.get_entity( stk::mesh::fem::FEMMetaData::get(mesh).element_rank() , elem_id );


            if (!topoVerifier.isTopologyBad(*elem))
            {
              std::cout << "no bad element and bad element expected" << "\n";
              std::cout.flush();
#ifndef REDS
              STKUNIT_EXPECT_TRUE( 0 );
#endif
            }


            ++num_hex ;
          }
        }
      }
    }


    std::vector<double> node_coordinates( 3 * node_map.size() );

    gmesh.coordinates( node_coordinates );

    if ( 3 * node_map.size() != node_coordinates.size() ) {
      std::ostringstream msg ;
      msg << "  P" << mesh.parallel_rank()
          << ": ERROR, node_map.size() = "
          << node_map.size()
          << " , node_coordinates.size() / 3 = "
          << ( node_coordinates.size() / 3 );
      throw std::runtime_error( msg.str() );
    }

    for ( unsigned i = 0 ; i < node_map.size() ; ++i ) {
      const unsigned i3 = i * 3 ;

      mesh::Entity * const node = mesh.get_entity( mesh::fem::FEMMetaData::NODE_RANK , node_map[i] );

      if ( NULL == node ) {
        std::ostringstream msg ;
        msg << "  P:" << mesh.parallel_rank()
            << " ERROR, Node not found: "
            << node_map[i] << " = node_map[" << i << "]" ;
        throw std::runtime_error( msg.str() );
      }

      double * const data = field_data( node_coord_field , *node );
      data[0] = node_coordinates[ i3 + 0 ];
      data[1] = node_coordinates[ i3 + 1 ];
      data[2] = node_coordinates[ i3 + 2 ];

      //FIXME ??? elem_node_coord_field
    }

  }
  catch ( const std::exception & X ) {
    std::cout << "  P:" << mesh.parallel_rank() << ": " << X.what()
              << std::endl ;
    std::cout.flush();
    error_flag = 1 ;
  }
  catch( ... ) {
    std::cout << "  P:" << mesh.parallel_rank()
              << " Caught unknown exception"
              << std::endl ;
    std::cout.flush();
    error_flag = 1 ;
  }

  all_reduce( mesh.parallel() , ReduceMax<1>( & error_flag ) );

  if ( error_flag ) {
    std::string msg( "Failed mesh generation" );
    throw std::runtime_error( msg );
  }

  mesh.modification_end();

  double dt = wall_dtime( t );

  all_reduce( mesh.parallel() , ReduceMax<1>( & dt ) );

  std::cout << "  P" << mesh.parallel_rank()
            << ": Meshed Hex = " << num_hex
            << " , Shell = " << num_shell
            << " , Node = " << num_nodes
            << " in " << dt << " sec"
            << std::endl ;
  std::cout.flush();
}

} // namespace stk_percept_unit

