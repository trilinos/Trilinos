#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <string>
#include <iostream>

#include "dofmngr/Panzer_IntrepidFieldPattern.hpp"
#include "dofmngr/Panzer_GeometricAggFieldPattern.hpp"

// include some intrepid basis functions
// 2D basis 
#include "Intrepid_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid_HDIV_TRI_I1_FEM.hpp"

#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid_HDIV_HEX_I1_FEM.hpp"

#include "TestFieldPattern.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {

std::string note = "***   NOTE: UNIT TEST BASED ON SEPT 2010   ***\n"
                   "***   INTREPID AND SHARDS Trilinos-dev     ***\n"
                   "***   DOXYGEN WEBSITE                      ***\n";

typedef Intrepid::FieldContainer<double> FieldContainer;

/////////////////////////////////////////////
// 2D tests
/////////////////////////////////////////////

// triangle tests
TEUCHOS_UNIT_TEST(tGeometricFieldPattern, test2d)
{
   out << note << std::endl;

   // basis to build patterns from
   RCP<Intrepid::Basis<double,FieldContainer> > basisA = rcp(new Intrepid::Basis_HGRAD_TRI_C1_FEM<double,FieldContainer>);
   RCP<Intrepid::Basis<double,FieldContainer> > basisB = rcp(new Intrepid::Basis_HGRAD_TRI_C2_FEM<double,FieldContainer>);
   RCP<Intrepid::Basis<double,FieldContainer> > basisC = rcp(new Intrepid::Basis_HGRAD_TRI_C1_FEM<double,FieldContainer>);

   RCP<const FieldPattern> patternA = rcp(new IntrepidFieldPattern(basisA));
   RCP<const FieldPattern> patternB = rcp(new IntrepidFieldPattern(basisB));
   RCP<const FieldPattern> patternC = rcp(new IntrepidFieldPattern(basisB));

   std::vector<RCP<const FieldPattern> > patterns;
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
      patterns.push_back(patternA);
      patterns.push_back(patternB);
      patterns.push_back(patternC);

      GeometricAggFieldPattern gfp;
      gfp.buildPattern(patterns);

      TEST_NOTHROW(gfp.getDimension());
      TEST_NOTHROW(gfp.getSubcellCount(0));
      TEST_NOTHROW(gfp.getSubcellIndices(0,0));
      TEST_THROW(gfp.getSubcellClosureIndices(0,0,indices),std::logic_error);

      TestFieldPattern tfp;
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
      RCP<Intrepid::Basis<double,FieldContainer> > basis = rcp(new Intrepid::Basis_HDIV_TRI_I1_FEM<double,FieldContainer>);
      patternB = rcp(new IntrepidFieldPattern(basis));

      patterns.clear();
      patterns.push_back(patternA);
      patterns.push_back(patternB);

      GeometricAggFieldPattern gfp;
      gfp.buildPattern(patterns);

      TEST_NOTHROW(gfp.getDimension());
      TEST_NOTHROW(gfp.getSubcellCount(0));
      TEST_NOTHROW(gfp.getSubcellIndices(0,0));
      TEST_THROW(gfp.getSubcellClosureIndices(0,0,indices),std::logic_error);

      TestFieldPattern tfp;
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
   RCP<Intrepid::Basis<double,FieldContainer> > basisA = rcp(new Intrepid::Basis_HGRAD_HEX_C1_FEM<double,FieldContainer>);
   RCP<Intrepid::Basis<double,FieldContainer> > basisB = rcp(new Intrepid::Basis_HGRAD_HEX_C2_FEM<double,FieldContainer>);
   RCP<Intrepid::Basis<double,FieldContainer> > basisC = rcp(new Intrepid::Basis_HGRAD_HEX_C1_FEM<double,FieldContainer>);

   RCP<const FieldPattern> patternA = rcp(new IntrepidFieldPattern(basisA));
   RCP<const FieldPattern> patternB = rcp(new IntrepidFieldPattern(basisB));
   RCP<const FieldPattern> patternC = rcp(new IntrepidFieldPattern(basisB));

   std::vector<RCP<const FieldPattern> > patterns;
   std::vector<int> indices;

   {
      patterns.clear();
      patterns.push_back(patternA);
      patterns.push_back(patternB);
      patterns.push_back(patternC);

      GeometricAggFieldPattern gfp;
      gfp.buildPattern(patterns);

      TEST_NOTHROW(gfp.getDimension());
      TEST_NOTHROW(gfp.getSubcellCount(0));
      TEST_NOTHROW(gfp.getSubcellIndices(0,0));
      TEST_THROW(gfp.getSubcellClosureIndices(0,0,indices),std::logic_error);

      TestFieldPattern tfp;
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
      RCP<Intrepid::Basis<double,FieldContainer> > basis = rcp(new Intrepid::Basis_HDIV_HEX_I1_FEM<double,FieldContainer>);
      patternB = rcp(new IntrepidFieldPattern(basis));

      patterns.clear();
      patterns.push_back(patternA);
      patterns.push_back(patternB);

      GeometricAggFieldPattern gfp;
      gfp.buildPattern(patterns);

      TEST_NOTHROW(gfp.getDimension());
      TEST_NOTHROW(gfp.getSubcellCount(0));
      TEST_NOTHROW(gfp.getSubcellIndices(0,0));
      TEST_THROW(gfp.getSubcellClosureIndices(0,0,indices),std::logic_error);

      TestFieldPattern tfp;
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
