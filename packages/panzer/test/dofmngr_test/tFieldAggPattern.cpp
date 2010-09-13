#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <string>
#include <iostream>

#include "dofmngr/Panzer_FieldAggPattern.hpp"
#include "dofmngr/Panzer_IntrepidFieldPattern.hpp"
#include "dofmngr/Panzer_GeometricAggFieldPattern.hpp"

// include some intrepid basis functions
// 2D basis 
#include "Intrepid_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"

// 3D basis 
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C2_FEM.hpp"

#include "Intrepid_HDIV_TRI_I1_FEM.hpp"

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

// tests:
//    - getGeometricAggPattern
//    - getFieldPattern
//    - getDimension 
TEUCHOS_UNIT_TEST(tFieldAggPattern, testA)
{
   out << note << std::endl;

   // basis to build patterns from
   RCP<Intrepid::Basis<double,FieldContainer> > basisA = rcp(new Intrepid::Basis_HGRAD_TRI_C1_FEM<double,FieldContainer>);
   RCP<Intrepid::Basis<double,FieldContainer> > basisB = rcp(new Intrepid::Basis_HGRAD_TRI_C2_FEM<double,FieldContainer>);
   RCP<Intrepid::Basis<double,FieldContainer> > basisC = rcp(new Intrepid::Basis_HGRAD_TRI_C1_FEM<double,FieldContainer>);

   RCP<const FieldPattern> patternA = rcp(new IntrepidFieldPattern(basisA));
   RCP<const FieldPattern> patternB = rcp(new IntrepidFieldPattern(basisB));
   RCP<const FieldPattern> patternC = rcp(new IntrepidFieldPattern(basisB));

   std::vector<int> closureIndices;
   std::vector<RCP<const FieldPattern> > patternV;
   std::vector<std::pair<int,RCP<const FieldPattern> > > patternM;

   patternV.push_back(patternA);
   patternV.push_back(patternB);
   patternV.push_back(patternC);

   GeometricAggFieldPattern geom(patternV);

   patternM.push_back(std::make_pair(7,patternA));
   patternM.push_back(std::make_pair(3,patternB));
   patternM.push_back(std::make_pair(4,patternC));

   // test build of geometric field pattern
   {
      FieldAggPattern agg; 
   
      TEST_THROW(agg.getSubcellClosureIndices(2,0,closureIndices),std::logic_error);
      TEST_THROW(*agg.getGeometricAggFieldPattern(),Teuchos::NullReferenceError);

      TEST_EQUALITY(agg.getGeometricAggFieldPattern(),Teuchos::null);

      agg.buildPattern(patternM);

      bool equality = false;
      TEST_NOTHROW(equality = geom.equals(*agg.getGeometricAggFieldPattern()));
      TEST_ASSERT(equality);

      TEST_EQUALITY(geom.getDimension(),agg.getGeometricAggFieldPattern()->getDimension());
   }

   // test build of geometric field pattern
   {
      FieldAggPattern agg(patternM); 

      TEST_THROW(agg.getSubcellClosureIndices(2,0,closureIndices),std::logic_error);
   
      bool equality = false;
      TEST_NOTHROW(equality = geom.equals(*agg.getGeometricAggFieldPattern()));
      TEST_ASSERT(equality);

      TEST_EQUALITY(geom.getDimension(),agg.getGeometricAggFieldPattern()->getDimension());
   }

   {
      FieldAggPattern agg(patternM); 
   
      TEST_THROW(agg.getFieldPattern(5),std::logic_error); // no field 5

      TEST_NOTHROW(agg.getFieldPattern(3));
      TEST_NOTHROW(agg.getFieldPattern(7));
      TEST_NOTHROW(agg.getFieldPattern(4));

      TEST_EQUALITY(agg.getFieldPattern(3),patternB);
      TEST_EQUALITY(agg.getFieldPattern(7),patternA);
      TEST_EQUALITY(agg.getFieldPattern(4),patternC);
   }
}

// tests:
//    - buildPattern --> different and same geometries
TEUCHOS_UNIT_TEST(tFieldAggPattern, testB)
{
   out << note << std::endl;

   // basis to build patterns from
   RCP<Intrepid::Basis<double,FieldContainer> > basisA = rcp(new Intrepid::Basis_HGRAD_TRI_C1_FEM<double,FieldContainer>);
   RCP<Intrepid::Basis<double,FieldContainer> > basisB = rcp(new Intrepid::Basis_HGRAD_TRI_C2_FEM<double,FieldContainer>);
   RCP<Intrepid::Basis<double,FieldContainer> > basisC = rcp(new Intrepid::Basis_HGRAD_TRI_C1_FEM<double,FieldContainer>);

   RCP<const FieldPattern> patternA = rcp(new IntrepidFieldPattern(basisA));
   RCP<const FieldPattern> patternB = rcp(new IntrepidFieldPattern(basisB));
   RCP<const FieldPattern> patternC = rcp(new IntrepidFieldPattern(basisB));

   // test that single construction gives the same pattern
   {
      std::vector<std::pair<int,RCP<const FieldPattern> > > patternM;
      patternM.push_back(std::make_pair(3,patternB));

      FieldAggPattern agg(patternM); 
      agg.print(out); // for debugging purposes

      TEST_EQUALITY(agg.getDimension(),patternB->getDimension()); 
      TEST_EQUALITY((int) agg.fieldIds().size(),patternB->numberIds());
      TEST_EQUALITY((int) agg.numFieldsPerId().size(),patternB->numberIds());

      const std::vector<int> & fieldsPerId = agg.numFieldsPerId();
      const std::vector<int> & fieldIds = agg.fieldIds();

      int numberFields = 0;
      for(std::size_t i=0;i<fieldsPerId.size();i++)
         numberFields += fieldsPerId[i];

      TEST_EQUALITY((int) fieldIds.size(),numberFields);

      // check fields per ID
      {
         bool allOnes = true;
         for(std::size_t i=0;i<fieldsPerId.size();i++)
            allOnes &= (fieldsPerId[i]==1);
         TEST_ASSERT(allOnes);
      }

      // check fields ids
      {
         bool allIds = true;
         for(std::size_t i=0;i<fieldIds.size();i++)
            allIds &= (fieldIds[i]==3);
         TEST_ASSERT(allIds);
      }

      for(int i=0;i<agg.getDimension();i++)
         TEST_EQUALITY(agg.getSubcellCount(i),patternB->getSubcellCount(i)); 
      TEST_ASSERT(agg.equals(*patternB)); 
   }

   // test that single construction gives the same pattern
   {
      std::vector<std::pair<int,RCP<const FieldPattern> > > patternM;
      patternM.push_back(std::make_pair(3,patternB));
      patternM.push_back(std::make_pair(7,patternA));

      FieldAggPattern agg(patternM); 
      agg.print(out); // for debugging purposes

      TEST_ASSERT(patternB->sameGeometry(agg));
      TEST_EQUALITY(agg.getDimension(),patternB->getDimension()); 
      TEST_EQUALITY((int) agg.numFieldsPerId().size(),agg.getGeometricAggFieldPattern()->numberIds());

      const std::vector<int> & fieldsPerId = agg.numFieldsPerId();
      const std::vector<int> & fieldIds = agg.fieldIds();

      int numberFields = 0;
      for(std::size_t i=0;i<fieldsPerId.size();i++)
         numberFields += fieldsPerId[i];

      TEST_EQUALITY((int) fieldIds.size(),numberFields);
      TEST_EQUALITY((int) fieldsPerId.size(),agg.getGeometricAggFieldPattern()->numberIds());

      int offset, dimension = 0;
      RCP<const FieldPattern> geomPattern = agg.getGeometricAggFieldPattern();

      // check out fieldsPerId
      offset = 0;
      dimension = 0;
      for(int nodes=0;nodes<geomPattern->getSubcellCount(dimension);nodes++)
         TEST_EQUALITY(fieldsPerId[offset+nodes],2); 

      offset += geomPattern->getSubcellCount(dimension);
      dimension = 1;
      for(int edges=0;edges<geomPattern->getSubcellCount(dimension);edges++)
         TEST_EQUALITY(fieldsPerId[offset+edges],1); 

      // this tests length of fieldsPerId is not longer than expected
      TEST_EQUALITY(offset+geomPattern->getSubcellCount(dimension),(int) fieldsPerId.size());

      // check out fieldIds
      offset = 0;
      dimension = 0;
      for(int nodes=0;nodes<geomPattern->getSubcellCount(dimension);nodes++) {
         TEST_EQUALITY(fieldIds[offset+2*nodes],3);
         TEST_EQUALITY(fieldIds[offset+2*nodes+1],7);
      }
 
      offset += 2*geomPattern->getSubcellCount(dimension);
      dimension = 1;
      for(int edges=0;edges<geomPattern->getSubcellCount(dimension);edges++)
         TEST_EQUALITY(fieldIds[offset+edges],3);

      // this tests length of fieldIds is not longer than expected
      TEST_EQUALITY(offset+geomPattern->getSubcellCount(dimension),(int) fieldIds.size());

      // check out subcell information
      int index=0;
      dimension = 0;
      for(int nodes=0;nodes<geomPattern->getSubcellCount(dimension);++nodes) {
         const std::vector<int> & indices = agg.getSubcellIndices(dimension,nodes);
         TEST_EQUALITY((int) indices.size(),2);
         for(std::size_t i=0;i<indices.size();i++,index++) 
            TEST_EQUALITY(indices[i],index);
      }

      dimension = 1;
      for(int edges=0;edges<geomPattern->getSubcellCount(dimension);++edges) {
         const std::vector<int> & indices = agg.getSubcellIndices(dimension,edges);
         TEST_EQUALITY((int) indices.size(),1);
         for(std::size_t i=0;i<indices.size();i++,index++) 
            TEST_EQUALITY(indices[i],index);
      }

      dimension = 2;
      for(int cells=0;cells<geomPattern->getSubcellCount(dimension);++cells) {
         const std::vector<int> & indices = agg.getSubcellIndices(dimension,cells);
         TEST_EQUALITY((int) indices.size(),0);
      }
   }
}

}
