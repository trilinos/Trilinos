//
//  ESEAS_PyramidBasisAgreementTests.cpp
//  VerificationDriver
//
//  Created by Roberts, Nathan V on 9/14/22.
//

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_TestUtils.hpp"
#include "Intrepid2_Types.hpp"
#include "Shards_CellTopology.hpp"

#include "Intrepid2_ESEAS_Interface.hpp"
#include "VerificationDriverHelpers.hpp"

using namespace Intrepid2;
using DeviceType = Kokkos::DefaultExecutionSpace::device_type;

namespace
{

TEUCHOS_UNIT_TEST( BasisAgreement, Pyramid_HGRAD )
{
  const EFunctionSpace fs = FUNCTION_SPACE_HGRAD;
  shards::CellTopology pyramidTopo = shards::CellTopology(shards::getCellTopologyData<shards::Pyramid<> >() );
  
  success = testBasis(pyramidTopo, fs, out, success);
}

TEUCHOS_UNIT_TEST( BasisAgreement, Pyramid_HCURL )
{
  const EFunctionSpace fs = FUNCTION_SPACE_HCURL;
  shards::CellTopology pyramidTopo = shards::CellTopology(shards::getCellTopologyData<shards::Pyramid<> >() );
  
  success = testBasis(pyramidTopo, fs, out, success);
}

TEUCHOS_UNIT_TEST( BasisAgreement, Pyramid_HDIV )
{
  const EFunctionSpace fs = FUNCTION_SPACE_HDIV;
  shards::CellTopology pyramidTopo = shards::CellTopology(shards::getCellTopologyData<shards::Pyramid<> >() );
  
  success = testBasis(pyramidTopo, fs, out, success);
}

TEUCHOS_UNIT_TEST( BasisAgreement, Pyramid_HVOL )
{
  const EFunctionSpace fs = FUNCTION_SPACE_HVOL;
  shards::CellTopology pyramidTopo = shards::CellTopology(shards::getCellTopologyData<shards::Pyramid<> >() );
  
  success = testBasis(pyramidTopo, fs, out, success);
}

} // end anonymous namespace

