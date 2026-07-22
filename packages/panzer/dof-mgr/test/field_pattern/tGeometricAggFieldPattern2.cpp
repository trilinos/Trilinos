// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <string>
#include <iostream>

#include "Intrepid2_config.h"
#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"

#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_GeometricAggFieldPattern.hpp"

#include "TestFieldPattern.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {

std::string note = "***   NOTE: UNIT TEST BASED ON SEPT 2010   ***\n"
                   "***   INTREPID AND SHARDS Trilinos-dev     ***\n"
                   "***   DOXYGEN WEBSITE                      ***\n";

typedef Kokkos::DynRankView<double,PHX::Device> FieldContainer;

/////////////////////////////////////////////
// 2D tests
/////////////////////////////////////////////

// triangle tests
TEUCHOS_UNIT_TEST(tGeometricFieldPattern, test2d)
{
   out << note << std::endl;

   // basis to build patterns from
   const int maxOrder = Intrepid2::Parameters::MaxOrder; 
   const int order = std::min(6, maxOrder);
   RCP<Intrepid2::Basis<PHX::Device,double,double> > basis = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<PHX::Device,double,double>(order));
   RCP<const FieldPattern> pattern = rcp(new Intrepid2FieldPattern(basis));

   std::vector<std::pair<FieldType,RCP<const FieldPattern>>> patterns;
   std::vector<int> indices;

   {
      // test unbuilt exceptions
      GeometricAggFieldPattern gfp;

      TEST_THROW(gfp.getDimension(),std::logic_error);
      TEST_THROW(gfp.getSubcellCount(0),std::logic_error);
      TEST_THROW(gfp.getSubcellIndices(0,0),std::logic_error);
      TEST_THROW(gfp.getSubcellClosureIndices(0,0,indices),std::logic_error);
   }

   {
      patterns.clear();
      patterns.push_back(std::make_pair(FieldType::CG,pattern));

      GeometricAggFieldPattern gfp;
      gfp.buildPattern(patterns);

      TEST_NOTHROW(gfp.getDimension());
      TEST_NOTHROW(gfp.getSubcellCount(0));
      TEST_NOTHROW(gfp.getSubcellIndices(0,0));
      TEST_THROW(gfp.getSubcellClosureIndices(0,0,indices),std::logic_error);

      TestFieldPattern tfp;
      tfp.cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >());
      tfp.subcellIndices.resize(3); // 3 geometric types: vertex, edge, interior

      // vertex
      tfp.subcellIndices[0].resize(4);
      tfp[0][0].push_back(0); tfp[0][1].push_back(1); tfp[0][2].push_back(2); tfp[0][3].push_back(3);

      // edge and interior
      tfp.subcellIndices[1].resize(4);
      tfp.subcellIndices[2].resize(1);
      if (order > 1) {
        tfp[1][0].push_back(4); tfp[1][1].push_back(5); tfp[1][2].push_back(6); tfp[1][3].push_back(7);
        tfp[2][0].push_back(8); 
      }
      TEST_ASSERT(gfp.equals(tfp));
   }
}

}
