/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/PrintTable.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>
#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/ConstantFunction.hpp>
#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/fixtures/QuadFixture.hpp>
#include <stk_percept/fixtures/WedgeFixture.hpp>

#include <stk_adapt/UniformRefinerPattern.hpp>
#include <stk_adapt/UniformRefiner.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_io/IossBridge.hpp>

#include <boost/lexical_cast.hpp>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>

#include <math.h>
#include <mpi.h>

namespace stk {
namespace adapt {
namespace unit_tests {

#define EXTRA_PRINT 0

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(unit_uniformRefiner, quad4_quad4_4_test_1)
{
  EXCEPTWATCH;
  stk::ParallelMachine pm = MPI_COMM_WORLD ;

  //const unsigned p_rank = stk::parallel_machine_rank( pm );
  const unsigned p_size = stk::parallel_machine_size( pm );
  //if (p_size == 1 || p_size == 3)
  if (p_size <= 3)
  {
    //const unsigned p_rank = stk::parallel_machine_rank( pm );
    //const unsigned p_size = stk::parallel_machine_size( pm );

    const unsigned n = 12;
    //const unsigned nx = n , ny = n , nz = p_size*n ;
    const unsigned nx = n , ny = n;

    percept::QuadFixture<double> fixture( pm , nx , ny, true);
    fixture.meta_data.commit();
    fixture.generate_mesh();

    percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
    eMesh.printInfo("quad fixture");
    eMesh.saveAs("quad_fixture.e");
  }
}

/// uses the Sierra-ported tables from framework/{element,mesh_modification}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(unit_uniformRefiner, break_quad_to_quad_sierra)
{
  EXCEPTWATCH;
  stk::ParallelMachine pm = MPI_COMM_WORLD ;

  //const unsigned p_rank = stk::parallel_machine_rank( pm );
  const unsigned p_size = stk::parallel_machine_size( pm );
  //if (p_size == 1 || p_size == 3)
  if (p_size <= 3)
  {
    //const unsigned p_rank = stk::parallel_machine_rank( pm );
    //const unsigned p_size = stk::parallel_machine_size( pm );

    const unsigned n = 12;
    //const unsigned nx = n , ny = n , nz = p_size*n ;
    const unsigned nx = n , ny = n;

    percept::QuadFixture<double> fixture( pm , nx , ny, true);

    percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, false);

    UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > break_quad_to_quad_4(eMesh);
    // FIXME
    int scalarDimension = 0; // a scalar
    FieldBase* proc_rank_field = eMesh.addField("proc_rank", mesh::Element, scalarDimension);

    //fixture.meta_data.commit();
    eMesh.commit();

    fixture.generate_mesh();

    eMesh.printInfo("quad mesh");

    UniformRefiner breaker(eMesh, break_quad_to_quad_4, proc_rank_field);
    //breaker.setRemoveOldElements(false);
    breaker.doBreak();

    //!!eMesh.saveAs("./square_quad4_ref_sierra_out.e");
    // end_demo
  }

}

STKUNIT_UNIT_TEST(unit_uniformRefiner, break_tri_to_tri_sierra)
{
  EXCEPTWATCH;
  stk::ParallelMachine pm = MPI_COMM_WORLD ;

  //const unsigned p_rank = stk::parallel_machine_rank( pm );
  const unsigned p_size = stk::parallel_machine_size( pm );
  if (p_size <= 3)
  {
    const unsigned n = 12;
    //const unsigned nx = n , ny = n , nz = p_size*n ;
    const unsigned nx = n , ny = n;

    percept::QuadFixture<double, shards::Triangle<3> > fixture( pm , nx , ny, true);

    percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, false);

    //             UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > break_tri_to_tri_4(eMesh);
    //             // FIXME
    //             int scalarDimension = 0; // a scalar
    //             FieldBase* proc_rank_field = eMesh.addField("proc_rank", mesh::Element, scalarDimension);

    //fixture.meta_data.commit();
    eMesh.commit();

    fixture.generate_mesh();

    eMesh.printInfo("tri mesh");

    //             UniformRefiner breaker(eMesh, break_tri_to_tri_4, proc_rank_field);
    //             breaker.doBreak();

    eMesh.saveAs("./quad_fixture_tri3.e");
    // end_demo
  }

}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(unit_uniformRefiner, hex8_hex8_8_1)
{
  EXCEPTWATCH;
  MPI_Barrier( MPI_COMM_WORLD );

  // start_demo_uniformRefiner_hex8_hex8_8_1

  percept::PerceptMesh eMesh;

  unsigned p_size = eMesh.getParallelSize();

  // generate a 4x4x(4*p_size) mesh
  std::string gmesh_spec = std::string("4x4x")+toString(4*p_size)+std::string("|bbox:0,0,0,1,1,1");
  eMesh.newMesh(percept::PerceptMesh::GMeshSpec(gmesh_spec));

  Hex8_Hex8_8 break_hex_to_hex(eMesh);

  int scalarDimension = 0; // a scalar
  FieldBase* proc_rank_field = eMesh.addField("proc_rank", mesh::Element, scalarDimension);

  eMesh.commit();

  UniformRefiner breaker(eMesh, break_hex_to_hex, proc_rank_field);
  breaker.doBreak();
  // end_demo
}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(unit_uniformRefiner, wedge6_1)
{
  EXCEPTWATCH;
  MPI_Barrier( MPI_COMM_WORLD );

  // start_demo_uniformRefiner_wedge6_1

  percept::PerceptMesh eMesh;

  unsigned p_size = eMesh.getParallelSize();
  if (p_size == 1)
  {

    //         void createMesh(stk::ParallelMachine parallel_machine,
    //                         unsigned n_nodes_x, unsigned n_nodes_y, unsigned n_nodes_z,
    //                         double xmin, double xmax,
    //                         double ymin, double ymax,
    //                         double zmin, double zmax,
    //                         std::string output_filename
    //                         )

    percept::WedgeFixture wedgeFixture;
    wedgeFixture.createMesh(MPI_COMM_WORLD,
                            4, 3, 2,
                            0, 1,
                            0, 1,
                            0, 1,
                            std::string("swept-wedge_0.e") );
  }
}


} // namespace regression_tests
} // namespace adapt
} // namespace stk

