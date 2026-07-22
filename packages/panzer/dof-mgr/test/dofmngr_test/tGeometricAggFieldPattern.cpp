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

#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_GeometricAggFieldPattern.hpp"

// include some intrepid basis functions
// 2D basis 
#include "Intrepid2_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid2_HDIV_TRI_I1_FEM.hpp"

#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid2_HDIV_HEX_I1_FEM.hpp"

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
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisA = rcp(new Intrepid2::Basis_HGRAD_TRI_C1_FEM<PHX::exec_space,double,double>);
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisB = rcp(new Intrepid2::Basis_HGRAD_TRI_C2_FEM<PHX::exec_space,double,double>);
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisC = rcp(new Intrepid2::Basis_HGRAD_TRI_C1_FEM<PHX::exec_space,double,double>);

   RCP<const FieldPattern> patternA = rcp(new Intrepid2FieldPattern(basisA));
   RCP<const FieldPattern> patternB = rcp(new Intrepid2FieldPattern(basisB));
   RCP<const FieldPattern> patternC = rcp(new Intrepid2FieldPattern(basisB));

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
      patterns.push_back(std::make_pair(FieldType::CG,patternA));
      patterns.push_back(std::make_pair(FieldType::CG,patternB));
      patterns.push_back(std::make_pair(FieldType::CG,patternC));

      GeometricAggFieldPattern gfp;
      gfp.buildPattern(patterns);

      TEST_NOTHROW(gfp.getDimension());
      TEST_NOTHROW(gfp.getSubcellCount(0));
      TEST_NOTHROW(gfp.getSubcellIndices(0,0));
      TEST_THROW(gfp.getSubcellClosureIndices(0,0,indices),std::logic_error);

      TestFieldPattern tfp;
      tfp.cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3> >());
      tfp.subcellIndices.resize(3); // two dimensional

      tfp.subcellIndices[0].resize(3);
      tfp[0][0].push_back(0); tfp[0][1].push_back(1); tfp[0][2].push_back(2);

      tfp.subcellIndices[1].resize(3);
      tfp[1][0].push_back(3); tfp[1][1].push_back(4); tfp[1][2].push_back(5);

      tfp.subcellIndices[2].resize(1);
      // empty

      TEST_ASSERT(gfp.equals(tfp));
   }

   {
      RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2::Basis_HDIV_TRI_I1_FEM<PHX::exec_space,double,double>);
      patternB = rcp(new Intrepid2FieldPattern(basis));

      patterns.clear();
      patterns.push_back(std::make_pair(FieldType::CG,patternA));
      patterns.push_back(std::make_pair(FieldType::CG,patternB));

      GeometricAggFieldPattern gfp;
      gfp.buildPattern(patterns);

      TEST_NOTHROW(gfp.getDimension());
      TEST_NOTHROW(gfp.getSubcellCount(0));
      TEST_NOTHROW(gfp.getSubcellIndices(0,0));
      TEST_THROW(gfp.getSubcellClosureIndices(0,0,indices),std::logic_error);

      TestFieldPattern tfp;
      tfp.cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3> >());
      tfp.subcellIndices.resize(3); // two dimensional

      tfp.subcellIndices[0].resize(3);
      tfp[0][0].push_back(0); tfp[0][1].push_back(1); tfp[0][2].push_back(2);

      tfp.subcellIndices[1].resize(3);
      tfp[1][0].push_back(3); tfp[1][1].push_back(4); tfp[1][2].push_back(5);

      tfp.subcellIndices[2].resize(1);
      // empty

      TEST_ASSERT(gfp.equals(tfp));
   }
}

// HEX tests
TEUCHOS_UNIT_TEST(tGeometricFieldPattern, test3d)
{

   out << note << std::endl;

   // basis to build patterns from
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisA = rcp(new Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::exec_space,double,double>);
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisB = rcp(new Intrepid2::Basis_HGRAD_HEX_C2_FEM<PHX::exec_space,double,double>);
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisC = rcp(new Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::exec_space,double,double>);

   RCP<const FieldPattern> patternA = rcp(new Intrepid2FieldPattern(basisA));
   RCP<const FieldPattern> patternB = rcp(new Intrepid2FieldPattern(basisB));
   RCP<const FieldPattern> patternC = rcp(new Intrepid2FieldPattern(basisB));

   std::vector<std::pair<FieldType,RCP<const FieldPattern>>> patterns;
   std::vector<int> indices;

   {
      patterns.clear();
      patterns.push_back(std::make_pair(FieldType::CG,patternA));
      patterns.push_back(std::make_pair(FieldType::CG,patternB));
      patterns.push_back(std::make_pair(FieldType::CG,patternC));

      GeometricAggFieldPattern gfp;
      gfp.buildPattern(patterns);

      TEST_NOTHROW(gfp.getDimension());
      TEST_NOTHROW(gfp.getSubcellCount(0));
      TEST_NOTHROW(gfp.getSubcellIndices(0,0));
      TEST_THROW(gfp.getSubcellClosureIndices(0,0,indices),std::logic_error);

      TestFieldPattern tfp;
      tfp.cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >());
      tfp.subcellIndices.resize(4); // three dimensional

      tfp.subcellIndices[0].resize(8);
      for(std::size_t i=0;i<8;i++)
         tfp[0][i].push_back(i);

      tfp.subcellIndices[1].resize(12);
      for(std::size_t i=0;i<12;i++)
         tfp[1][i].push_back(i+8);

      tfp.subcellIndices[2].resize(6);
      for(std::size_t i=0;i<6;i++)
         tfp[2][i].push_back(i+20);

      tfp.subcellIndices[3].resize(1);
      tfp[3][0].push_back(26);

      TEST_ASSERT(gfp.equals(tfp));
   }

   {
      RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2::Basis_HDIV_HEX_I1_FEM<PHX::exec_space,double,double>);
      patternB = rcp(new Intrepid2FieldPattern(basis));

      patterns.clear();
      patterns.push_back(std::make_pair(FieldType::CG,patternA));
      patterns.push_back(std::make_pair(FieldType::CG,patternB));

      GeometricAggFieldPattern gfp;
      gfp.buildPattern(patterns);

      TEST_NOTHROW(gfp.getDimension());
      TEST_NOTHROW(gfp.getSubcellCount(0));
      TEST_NOTHROW(gfp.getSubcellIndices(0,0));
      TEST_THROW(gfp.getSubcellClosureIndices(0,0,indices),std::logic_error);

      TestFieldPattern tfp;
      tfp.cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >());
      tfp.subcellIndices.resize(4); // three dimensional

      tfp.subcellIndices[0].resize(8);
      for(std::size_t i=0;i<8;i++)
         tfp[0][i].push_back(i);

      tfp.subcellIndices[1].resize(12);
      // empty

      tfp.subcellIndices[2].resize(6);
      for(std::size_t i=0;i<6;i++)
         tfp[2][i].push_back(i+8);

      tfp.subcellIndices[3].resize(1);

      tfp.print(out);
      gfp.print(out);

      TEST_ASSERT(gfp.equals(tfp));
   }
}

}
