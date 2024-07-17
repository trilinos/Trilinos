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

#include "Panzer_FieldAggPattern.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_GeometricAggFieldPattern.hpp"
#include "Panzer_NodalFieldPattern.hpp"

// include some intrepid basis functions
// 2D basis 
#include "Intrepid2_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_HCURL_QUAD_I1_FEM.hpp"

// 3D basis 
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_Cn_FEM.hpp"
#include "Intrepid2_HDIV_HEX_I1_FEM.hpp"
#include "Intrepid2_HCURL_HEX_I1_FEM.hpp"

#include "Intrepid2_HDIV_TRI_I1_FEM.hpp"

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

// tests:
//    - getGeometricAggPattern
//    - getFieldPattern
//    - getDimension 
TEUCHOS_UNIT_TEST(tFieldAggPattern, testA)
{

   out << note << std::endl;

   // basis to build patterns from
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisA = rcp(new Intrepid2::Basis_HGRAD_TRI_C1_FEM<PHX::exec_space,double,double>);
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisB = rcp(new Intrepid2::Basis_HGRAD_TRI_C2_FEM<PHX::exec_space,double,double>);
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisC = rcp(new Intrepid2::Basis_HGRAD_TRI_C1_FEM<PHX::exec_space,double,double>);

   RCP<const FieldPattern> patternA = rcp(new Intrepid2FieldPattern(basisA));
   RCP<const FieldPattern> patternB = rcp(new Intrepid2FieldPattern(basisB));
   RCP<const FieldPattern> patternC = rcp(new Intrepid2FieldPattern(basisB));

   std::vector<int> closureIndices;
   std::vector<std::pair<FieldType,RCP<const FieldPattern>>> patternV;
   std::vector<std::tuple<int,FieldType,RCP<const FieldPattern> > > patternM;

   patternV.push_back(std::make_pair(FieldType::CG,patternA));
   patternV.push_back(std::make_pair(FieldType::CG,patternB));
   patternV.push_back(std::make_pair(FieldType::CG,patternC));

   GeometricAggFieldPattern geom(patternV);

   patternM.push_back(std::make_tuple(7,FieldType::CG,patternA));
   patternM.push_back(std::make_tuple(3,FieldType::CG,patternB));
   patternM.push_back(std::make_tuple(4,FieldType::CG,patternC));

   // test build of geometric field pattern
   {
      FieldAggPattern agg; 
   
      TEST_THROW(agg.getSubcellClosureIndices(2,0,closureIndices),std::logic_error);
      // TEST_THROW(*agg.getGeometricAggFieldPattern(),Teuchos::NullReferenceError);

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
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisA = rcp(new Intrepid2::Basis_HGRAD_TRI_C1_FEM<PHX::exec_space,double,double>);
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisB = rcp(new Intrepid2::Basis_HGRAD_TRI_C2_FEM<PHX::exec_space,double,double>);
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisC = rcp(new Intrepid2::Basis_HGRAD_TRI_C1_FEM<PHX::exec_space,double,double>);

   RCP<const FieldPattern> patternA = rcp(new Intrepid2FieldPattern(basisA));
   RCP<const FieldPattern> patternB = rcp(new Intrepid2FieldPattern(basisB));
   RCP<const FieldPattern> patternC = rcp(new Intrepid2FieldPattern(basisB));

   // test that single construction gives the same pattern
   {
      std::vector<std::tuple<int,FieldType,RCP<const FieldPattern> > > patternM;
      patternM.push_back(std::make_tuple(3,FieldType::CG,patternB));

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
      std::vector<std::tuple<int,FieldType,RCP<const FieldPattern> > > patternM;
      patternM.push_back(std::make_tuple(3,FieldType::CG,patternB));
      patternM.push_back(std::make_tuple(7,FieldType::CG,patternA));

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

// tests:
//    - localOffsets --> different and same geometries
TEUCHOS_UNIT_TEST(tFieldAggPattern, testC)
{
   out << note << std::endl;

   // basis to build patterns from

   {
      RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2::Basis_HGRAD_TRI_C2_FEM<PHX::exec_space,double,double>);
      RCP<const FieldPattern> pattern = rcp(new Intrepid2FieldPattern(basis));
      std::vector<std::tuple<int,FieldType,RCP<const FieldPattern> > > patternM;
      patternM.push_back(std::make_tuple(3,FieldType::CG,pattern));

      FieldAggPattern agg(patternM);
      agg.print(out);
      TEUCHOS_ASSERT(agg.equals(*pattern));

      TEST_THROW(agg.localOffsets(7),std::logic_error); // check for failure if ID does not exist
 
      std::vector<int> v = agg.localOffsets(3);
      for(std::size_t i=0;i<v.size();i++)
         TEST_EQUALITY(v[i],(int) i);
   }

   {
      RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2::Basis_HGRAD_HEX_C2_FEM<PHX::exec_space,double,double>);
      RCP<const FieldPattern> pattern = rcp(new Intrepid2FieldPattern(basis));
      std::vector<std::tuple<int,FieldType,RCP<const FieldPattern> > > patternM;
      patternM.push_back(std::make_tuple(3,FieldType::CG,pattern));

      FieldAggPattern agg(patternM);
      agg.print(out);

      TEST_THROW(agg.localOffsets(7),std::logic_error); // check for failure if ID does not exist

      std::vector<int> v_true(27);
      // nodes
      v_true[0] = 0; v_true[1] = 1; v_true[2] = 2; v_true[3] = 3;
      v_true[4] = 4; v_true[5] = 5; v_true[6] = 6; v_true[7] = 7;
 
      // edges     
      v_true[8]  = 8; v_true[9]  = 9; v_true[10] = 10; v_true[11] = 11;
      v_true[12] = 16; v_true[13] = 17; v_true[14] = 18 ;v_true[15] = 19;
      v_true[16] = 12; v_true[17] = 13; v_true[18] = 14; v_true[19] = 15;

      // volume
      v_true[20] = 26;

      // faces
      v_true[21] = 24; v_true[22] = 25; v_true[23] = 23; 
      v_true[24] = 21; v_true[25] = 20; v_true[26] = 22;
 
      std::vector<int> v = agg.localOffsets(3);
      for(std::size_t i=0;i<v.size();i++)
         TEST_EQUALITY(v[i],v_true[i]);
   }

   {
      RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisC1 = rcp(new Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::exec_space,double,double>);
      RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisC2 = rcp(new Intrepid2::Basis_HGRAD_HEX_C2_FEM<PHX::exec_space,double,double>);
      RCP<const FieldPattern> patternC1 = rcp(new Intrepid2FieldPattern(basisC1));
      RCP<const FieldPattern> patternC2 = rcp(new Intrepid2FieldPattern(basisC2));
      std::vector<std::tuple<int,FieldType,RCP<const FieldPattern> > > patternM;
      patternM.push_back(std::make_tuple(0,FieldType::CG,patternC2));
      patternM.push_back(std::make_tuple(1,FieldType::CG,patternC2));
      patternM.push_back(std::make_tuple(2,FieldType::CG,patternC2));
      patternM.push_back(std::make_tuple(3,FieldType::CG,patternC1));

      FieldAggPattern agg(patternM);
      agg.print(out);

      TEST_THROW(agg.localOffsets(7),std::logic_error); // check for failure if ID does not exist

      {
         out << "\nTesting Q2-Q1 NS: Pressure" << std::endl;
         std::vector<int> v_true(8);
         v_true[0] = 3; v_true[1] = 7; v_true[2] = 11; v_true[3] = 15;
         v_true[4] = 19; v_true[5] = 23; v_true[6] = 27; v_true[7] = 31;
         std::vector<int> v = agg.localOffsets(3);
         for(std::size_t i=0;i<v.size();i++)
            TEST_EQUALITY(v[i],v_true[i]);
      }

      for(int d=0;d<3;d++) {
         out << "\nTesting Q2-Q1 NS: Velocity d=" << d << std::endl;
         std::vector<int> v_true(27);

         // nodes
         v_true[0] = 4*0+0; v_true[1] = 4*1+0; v_true[2] = 4*2+0; v_true[3] = 4*3+0;
         v_true[4] = 4*4+0; v_true[5] = 4*5+0; v_true[6] = 4*6+0; v_true[7] = 4*7+0;
    
         // edges     
         v_true[8]  = 3*0+32; v_true[9]  = 3*1+32; v_true[10] = 3*2+32; v_true[11] = 3*3+32;
         v_true[12] = 3*8+32; v_true[13] = 3*9+32; v_true[14] = 3*10+32; v_true[15] = 3*11+32;
         v_true[16] = 3*4+32; v_true[17] = 3*5+32; v_true[18] = 3*6+32; v_true[19] = 3*7+32;
   
         // volume
         v_true[20] = 86;
   
         // faces
         v_true[21] = 3*4+68; v_true[22] = 3*5+68; v_true[23] = 3*3+68; 
         v_true[24] = 3*1+68; v_true[25] = 3*0+68; v_true[26] = 3*2+68;

         std::vector<int> v = agg.localOffsets(d);
         for(std::size_t i=0;i<v.size();i++)
            TEST_EQUALITY(v[i],v_true[i]+d);
      }
   }

   {
      out << "Crazy HEX basis test" << std::endl;
 
      RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisC1 = rcp(new Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::exec_space,double,double>);
      RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisDivI1 = rcp(new Intrepid2::Basis_HDIV_HEX_I1_FEM<PHX::exec_space,double,double>);
      RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisCurlI1 = rcp(new Intrepid2::Basis_HCURL_HEX_I1_FEM<PHX::exec_space,double,double>);
      RCP<const FieldPattern> patternC1 = rcp(new Intrepid2FieldPattern(basisC1));
      RCP<const FieldPattern> patternDivI1 = rcp(new Intrepid2FieldPattern(basisDivI1));
      RCP<const FieldPattern> patternCurlI1 = rcp(new Intrepid2FieldPattern(basisCurlI1));
      std::vector<std::tuple<int,FieldType,RCP<const FieldPattern> > > patternM;
      patternM.push_back(std::make_tuple(3,FieldType::CG,patternC1));
      patternM.push_back(std::make_tuple(9,FieldType::CG,patternDivI1));
      patternM.push_back(std::make_tuple(1,FieldType::CG,patternCurlI1));

      FieldAggPattern agg(patternM);
      agg.print(out);

      TEST_THROW(agg.localOffsets(7),std::logic_error); // check for failure if ID does not exist

      {
         out << "\nTesting C1" << std::endl;
         std::vector<int> v = agg.localOffsets(3);
         for(std::size_t i=0;i<v.size();i++)
            TEST_EQUALITY(v[i],(int) i);
      }

      {
         out << "\nTesting Div-I1" << std::endl;
         std::vector<int> v = agg.localOffsets(9);
         for(std::size_t i=0;i<v.size();i++)
            TEST_EQUALITY(v[i],(int) i+20);
      }

      {
         out << "\nTesting Curl-I1" << std::endl;
         std::vector<int> v = agg.localOffsets(1);
         for(std::size_t i=0;i<v.size();i++)
            TEST_EQUALITY(v[i],(int) i+8);
      }
   }

   {
      out << "Highorder HEX basis test: FAILS!!!! DISABLED FOR NOW!!!!!" << std::endl;
 
/*
      RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2::Basis_HGRAD_HEX_Cn_FEM<PHX::exec_space,double,double>(4,Intrepid2::POINTTYPE_EQUISPACED));
      RCP<const FieldPattern> pattern = rcp(new Intrepid2FieldPattern(basis));
      std::vector<std::tuple<int,FieldType,RCP<const FieldPattern> > > patternM;
      patternM.push_back(std::make_tuple(3,FieldType::CG,pattern));

      FieldAggPattern agg(patternM);
      out << "4th order hex" << std::endl;
      pattern->print(out);

      out << "Aggregated hex" << std::endl;
      agg.print(out);

      out << "Geometric hex" << std::endl;
      agg.getGeometricAggFieldPattern()->print(out);

      TEST_THROW(agg.localOffsets(7),std::logic_error); // check for failure if ID does not exist
*/
   }
}

// tests:
//    - localOffsets_closure --> different and same geometries
TEUCHOS_UNIT_TEST(tFieldAggPattern, testD)
{
   {
      RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2::Basis_HGRAD_HEX_C2_FEM<PHX::exec_space,double,double>);
      RCP<const FieldPattern> pattern = rcp(new Intrepid2FieldPattern(basis));
      std::vector<std::tuple<int,FieldType,RCP<const FieldPattern> > > patternM;
      patternM.push_back(std::make_tuple(3,FieldType::CG,pattern));

      FieldAggPattern agg(patternM);
      agg.print(out);

      out << "Testing throws of localOffsets_closure" << std::endl;
      TEST_THROW(agg.localOffsets_closure(7,0,1),std::logic_error); // check for failure if ID does not exist
      TEST_THROW(agg.localOffsets_closure(3,2,6),std::logic_error); // check for failure if sub cell doesn't exist

      std::vector<int> v_true;
      std::vector<int> v;

      // test (2,0) closure
      out << "Testing localOffsets_closure(3,2,0) on a C2 HEX" << std::endl;
      v_true.resize(9);
      v = agg.localOffsets_closure(3,2,0).first;

      // nodes
      v_true[0] = 0; v_true[1] = 1; v_true[2] = 4; v_true[3] = 5;
 
      // edges     
      v_true[4]  = 8; v_true[5]  = 17; v_true[6] = 12; v_true[7] = 16;

      // areas
      v_true[8] = 20;

      TEST_EQUALITY(v.size(),v_true.size());

      std::sort(v.begin(),v.end());
      std::sort(v_true.begin(),v_true.end());
      for(std::size_t i=0;i<v.size();i++)
         TEST_EQUALITY(v[i],v_true[i]);

      // test (1,7) closure
      out << "Testing localOffsets_closure(3,1,7) on a C2 HEX" << std::endl;
      v_true.resize(3);
      v = agg.localOffsets_closure(3,1,7).first;

      // nodes
      v_true[0] = 4; v_true[1] = 7;
 
      // edges     
      v_true[2]  = 15; 

      TEST_EQUALITY(v.size(),v_true.size());

      std::sort(v.begin(),v.end());
      std::sort(v_true.begin(),v_true.end());
      for(std::size_t i=0;i<v.size();i++)
         TEST_EQUALITY(v[i],v_true[i]);

      // test all dimension zero closures
      out << "Testing localOffsets_closure(3,0,*) on a C2 HEX" << std::endl;
      v_true.resize(1);
      for(int i=0;i<8;i++) {
         v = agg.localOffsets_closure(3,0,i).first;

         // nodes
         v_true[0] = i;
 
         TEST_EQUALITY(v.size(),v_true.size());
         TEST_EQUALITY(v[0],v_true[0]);
      }
   }

   {
      RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisC1 = rcp(new Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::exec_space,double,double>);
      RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisC2 = rcp(new Intrepid2::Basis_HGRAD_HEX_C2_FEM<PHX::exec_space,double,double>);
      RCP<const FieldPattern> patternP = rcp(new Intrepid2FieldPattern(basisC1));
      RCP<const FieldPattern> patternU = rcp(new Intrepid2FieldPattern(basisC2));
      std::vector<std::tuple<int,FieldType,RCP<const FieldPattern> > > patternM;
      int numP = 8;
      int numU = 4;
      patternM.push_back(std::make_tuple(numP,FieldType::CG,patternP));
      patternM.push_back(std::make_tuple(numU,FieldType::CG,patternU));

      FieldAggPattern agg(patternM);
      agg.print(out);

      out << "Testing throws of localOffsets_closure" << std::endl;
      TEST_THROW(agg.localOffsets_closure(7,0,1),std::logic_error); // check for failure if ID does not exist
      TEST_THROW(agg.localOffsets_closure(3,2,6),std::logic_error); // check for failure if sub cell doesn't exist

      TEST_EQUALITY(agg.localOffsets_closure(numP,1,6).first.size(),2); // pressure basis only has nodes
      TEST_EQUALITY(agg.localOffsets_closure(numP,2,1).first.size(),4); // pressure basis only has nodes

      std::vector<int> v_true;
      std::vector<int> v;

      // test numP first
      ///////////////////////////////////////

      // test (2,2) closure
      out << "Testing localOffsets_closure(numP,2,2) on a C2 HEX" << std::endl;
      v_true.resize(4);
      v = agg.localOffsets_closure(numP,2,2).first;

      // nodes
      v_true[0] = 2*2; v_true[1] = 2*3; v_true[2] = 2*7; v_true[3] = 2*6;
 
      TEST_EQUALITY(v.size(),v_true.size());

      std::sort(v.begin(),v.end());
      std::sort(v_true.begin(),v_true.end());
      for(std::size_t i=0;i<v.size();i++)
         TEST_EQUALITY(v[i],v_true[i]);

      // test (1,11) closure
      out << "Testing localOffsets_closure(numP,1,11) on a C2 HEX" << std::endl;
      v_true.resize(2);
      v = agg.localOffsets_closure(numP,1,11).first;

      // nodes
      v_true[0] = 2*3; v_true[1] = 2*7;
 
      TEST_EQUALITY(v.size(),v_true.size());

      std::sort(v.begin(),v.end());
      std::sort(v_true.begin(),v_true.end());
      for(std::size_t i=0;i<v.size();i++)
         TEST_EQUALITY(v[i],v_true[i]);

      // test all dimension zero closures
      out << "Testing localOffsets_closure(numP,0,*) on a C2 HEX" << std::endl;
      v_true.resize(1);
      for(int i=0;i<8;i++) {
         v = agg.localOffsets_closure(numP,0,i).first;

         // nodes
         v_true[0] = 2*i;
 
         TEST_EQUALITY(v.size(),v_true.size());
         TEST_EQUALITY(v[0],v_true[0]);
      }

      // test numU second
      ///////////////////////////////////////

      // test (2,4) closure
      out << "Testing localOffsets_closure(numU,2,4) on a C2 HEX" << std::endl;
      v_true.resize(9);
      v = agg.localOffsets_closure(numU,2,4).first;

      // nodes
      v_true[0] = 2*0+1; v_true[1] = 2*1+1; v_true[2] = 2*2+1; v_true[3] = 2*3+1;
 
      // edges     
      v_true[4]  = 16; v_true[5]  = 17; v_true[6] = 18; v_true[7] = 19;

      // areas
      v_true[8] = 32;

      TEST_EQUALITY(v.size(),v_true.size());

      std::sort(v.begin(),v.end());
      std::sort(v_true.begin(),v_true.end());
      for(std::size_t i=0;i<v.size();i++)
         TEST_EQUALITY(v[i],v_true[i]);

      // test (1,10) closure
      out << "Testing localOffsets_closure(numU,1,10) on a C2 HEX" << std::endl;
      v_true.resize(3);
      v = agg.localOffsets_closure(numU,1,10).first;

      // nodes
      v_true[0] = 2*2+1; v_true[1] = 2*6+1;
 
      // edges     
      v_true[2]  = 26; 

      TEST_EQUALITY(v.size(),v_true.size());

      std::sort(v.begin(),v.end());
      std::sort(v_true.begin(),v_true.end());
      for(std::size_t i=0;i<v.size();i++)
         TEST_EQUALITY(v[i],v_true[i]);

      // test all dimension zero closures
      out << "Testing localOffsets_closure(numU,0,*) on a C2 HEX" << std::endl;
      v_true.resize(1);
      for(int i=0;i<8;i++) {
         v = agg.localOffsets_closure(numU,0,i).first;

         // nodes
         v_true[0] = 2*i+1;
 
         TEST_EQUALITY(v.size(),v_true.size());
         TEST_EQUALITY(v[0],v_true[0]);
      }
   }
}

// Tests the HCURL case where the geometry includes the nodes, but the pattern does not
TEUCHOS_UNIT_TEST(tFieldAggPattern, testE)
{
   out << note << std::endl;

   // basis to build patterns from
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisA = rcp(new Intrepid2::Basis_HCURL_QUAD_I1_FEM<PHX::exec_space,double,double>);

   RCP<const FieldPattern> patternA = rcp(new Intrepid2FieldPattern(basisA));
   RCP<const FieldPattern> patternNode = rcp(new NodalFieldPattern(basisA->getBaseCellTopology()));

   std::vector<int> closureIndices;
   std::vector<std::pair<FieldType,RCP<const FieldPattern>>> patternV;
   std::vector<std::tuple<int,FieldType,RCP<const FieldPattern> > > patternM;

   patternV.push_back(std::make_pair(FieldType::CG,patternA));
   patternV.push_back(std::make_pair(FieldType::CG,patternNode));

   GeometricAggFieldPattern geom(patternV);

   patternM.push_back(std::make_tuple(7,FieldType::CG,patternA));

   // test build of geometric field pattern
   {
      FieldAggPattern agg; 
   
      agg.buildPattern(patternM,Teuchos::rcpFromRef(geom));

      bool equality = false;
      TEST_NOTHROW(equality = geom.equals(*agg.getGeometricAggFieldPattern()));
      TEST_ASSERT(equality);

      TEST_EQUALITY(geom.getDimension(),agg.getGeometricAggFieldPattern()->getDimension());

      out << agg << std::endl;
   }
}

}
