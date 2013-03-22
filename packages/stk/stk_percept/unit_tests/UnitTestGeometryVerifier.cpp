/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_percept/function/FunctionOperator.hpp>
#include <stk_percept/function/ConstantFunction.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/TopologyVerifier.hpp>
#include <stk_percept/GeometryVerifier.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <stk_percept/fixtures/Fixture.hpp>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>
#include <math.h>

namespace stk {
namespace percept {
namespace unit_tests {

// on pathscale platform this doesn't work (something to do with static variables)

static int dw_enabled = 1;
static stk::diag::Writer s_diagWriter(std::cout.rdbuf(), dw_enabled);
static stk::diag::Writer &
dw()
{
  //static stk::diag::Writer s_diagWriter(dwout().rdbuf(), 0);

  s_diagWriter.setPrintMask(LOG_NORM+LOG_ALWAYS);

  return s_diagWriter;
}

#define EXTRA_PRINT 0

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(geom, volume)
{
  dw().m(LOG_GEOMETRY_VERIFIER) << "TEST::geom::volume " << stk::diag::dendl;

  const size_t num_x = 3;
  const size_t num_y = 3;
  const size_t num_z = 3;
  std::string config_mesh =
    Ioss::Utils::to_string(num_x) + "x" +
    Ioss::Utils::to_string(num_y) + "x" +
    Ioss::Utils::to_string(num_z) + "|bbox:0,0,0,1,1,1";
	
  PerceptMesh eMesh(3u);
  eMesh.new_mesh(GMeshSpec(config_mesh));
  eMesh.commit();
  // no need for this in create mode: eMesh.readBulkData();

  //FEMMetaData& metaData = *eMesh.get_fem_meta_data();
  mesh::BulkData& bulkData = *eMesh.get_bulk_data();

  eMesh.dump();
  GeometryVerifier geomVerifier(false);
  geomVerifier.isGeometryBad(bulkData, true);
  //setDoPause(true);
  //pause();
}

}
}
}
