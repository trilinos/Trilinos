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

#include "Intrepid2_HVOL_C0_FEM.hpp"
#include "Panzer_FieldAggPattern.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_GeometricAggFieldPattern.hpp"
#include "Panzer_GlobalIndexer_Utilities.hpp"

// include some intrepid basis functions
// 2D basis 
#include "Intrepid2_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid2_HCURL_TRI_I1_FEM.hpp"
#include "Intrepid2_HDIV_TRI_I1_FEM.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_HCURL_QUAD_I1_FEM.hpp"
#include "Intrepid2_HDIV_QUAD_I1_FEM.hpp"

#include "Intrepid2_HDIV_TET_I1_FEM.hpp"
#include "Intrepid2_HGRAD_TET_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TET_C2_FEM.hpp"

#include "Intrepid2_HDIV_HEX_I1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"

#include "Shards_BasicTopologies.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {

std::string note = "***   NOTE: UNIT TEST BASED ON MARCH 2012  ***\n"
                   "***   INTREPID AND SHARDS Trilinos-dev     ***\n"
                   "***   DOXYGEN WEBSITE                      ***\n";
typedef Kokkos::DynRankView<double,PHX::Device> FieldContainer;

/////////////////////////////////////////////
// 2D tests
/////////////////////////////////////////////

TEUCHOS_UNIT_TEST(tOrientation, testEdgeBasis_tri)
{

   out << note << std::endl;

   // basis to build patterns from
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisA = rcp(new Intrepid2::Basis_HGRAD_TRI_C1_FEM<PHX::exec_space,double,double>);
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisB = rcp(new Intrepid2::Basis_HCURL_TRI_I1_FEM<PHX::exec_space,double,double>);
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisC = rcp(new Intrepid2::Basis_HGRAD_TRI_C2_FEM<PHX::exec_space,double,double>); // used further down

   RCP<const FieldPattern> patternA = rcp(new Intrepid2FieldPattern(basisA));
   RCP<const FieldPattern> patternB = rcp(new Intrepid2FieldPattern(basisB));
   RCP<const FieldPattern> patternC = rcp(new Intrepid2FieldPattern(basisC)); // used further down

   TEST_EQUALITY(patternA->numberIds(),3);
   TEST_EQUALITY(patternB->numberIds(),3);

   std::vector<std::pair<int,int> > topEdgeIndices;
   orientation_helpers::computePatternEdgeIndices(*patternA,topEdgeIndices);

   TEST_EQUALITY(topEdgeIndices.size(),3);
   TEST_EQUALITY(topEdgeIndices[0].first,0); TEST_EQUALITY(topEdgeIndices[0].second,1);
   TEST_EQUALITY(topEdgeIndices[1].first,1); TEST_EQUALITY(topEdgeIndices[1].second,2);
   TEST_EQUALITY(topEdgeIndices[2].first,2); TEST_EQUALITY(topEdgeIndices[2].second,0);

   std::vector<std::vector<panzer::GlobalOrdinal>> connectivity(4);
   connectivity[0].resize(patternA->numberIds());
   connectivity[1].resize(patternA->numberIds());
   connectivity[2].resize(patternA->numberIds());
   connectivity[3].resize(patternA->numberIds());

   // Topologocally the four elements look like:
   //
   //    3-----1
   //    |\   /|
   //    | \ / |
   //    |  6  |
   //    | / \ |
   //    |/   \|
   //    5-----0
   //
   // all that matters is the global
   // node numbering and the local ordering
   // The local ordering is defined by the following connectivity
   connectivity[0][0] = 0; connectivity[0][1] = 6; connectivity[0][2] = 5; 
   connectivity[1][0] = 6; connectivity[1][1] = 0; connectivity[1][2] = 1; 
   connectivity[2][0] = 1; connectivity[2][1] = 3; connectivity[2][2] = 6; 
   connectivity[3][0] = 3; connectivity[3][1] = 5; connectivity[3][2] = 6; 

   // LOCAL element orientations are always set so that they flow in the positive
   // direction along an edge from node 0 to node 1. As a result if the GID of
   // node 0 is larger then node 1 then the GLOBAL orientation is -1 (and positive
   // otherwise. The local definition of the edge direction is defined by 
   // the shards cell topology.
   {
      std::vector<signed char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity[0], *patternB, orientations);
      TEST_EQUALITY(orientations[0],1);
      TEST_EQUALITY(orientations[1],-1);
      TEST_EQUALITY(orientations[2],-1);
   }

   {
      std::vector<signed char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity[1], *patternB, orientations);
      TEST_EQUALITY(orientations[0],-1);
      TEST_EQUALITY(orientations[1],1);
      TEST_EQUALITY(orientations[2],1);
   }

   {
      std::vector<signed char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity[2], *patternB, orientations);
      TEST_EQUALITY(orientations[0],1);
      TEST_EQUALITY(orientations[1],1);
      TEST_EQUALITY(orientations[2],-1);
   }

   {
      std::vector<signed char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity[3], *patternB, orientations);
      TEST_EQUALITY(orientations[0],1);
      TEST_EQUALITY(orientations[1],1);
      TEST_EQUALITY(orientations[2],-1);
   }

   // lets now test an aggregate basis
   std::vector<std::tuple<int,FieldType,Teuchos::RCP<const FieldPattern> > > patterns;
   patterns.push_back(std::make_tuple(0,FieldType::CG,patternB));
   patterns.push_back(std::make_tuple(1,FieldType::CG,patternC));
   patterns.push_back(std::make_tuple(2,FieldType::CG,patternA));
   Teuchos::RCP<const FieldPattern> aggPattern = Teuchos::rcp(new FieldAggPattern(patterns));

   std::vector<signed char> orientations(aggPattern->numberIds(),0);
   orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity[2], *aggPattern, orientations);
   int nonzeroCount = 0;
   for(std::size_t s=0;s<orientations.size();s++)
      nonzeroCount += orientations[s]*orientations[s]; // should be +1 only if it is an edge
   TEST_EQUALITY(nonzeroCount,6);

   // loop over edges
   signed char tests[] = { 1,1,-1}; // for element 2
   for(std::size_t e=0;e<3;e++) {
      const std::vector<int> & edgeIndices = aggPattern->getSubcellIndices(1,e);
      for(std::size_t s=0;s<edgeIndices.size();s++)
         TEST_EQUALITY(orientations[edgeIndices[s]],tests[e]);
   }       
}

TEUCHOS_UNIT_TEST(tOrientation, testEdgeBasis_quad)
{

   out << note << std::endl;

   // basis to build patterns from
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisA = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double>);
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisB = rcp(new Intrepid2::Basis_HCURL_QUAD_I1_FEM<PHX::exec_space,double,double>);
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisC = rcp(new Intrepid2::Basis_HGRAD_QUAD_C2_FEM<PHX::exec_space,double,double>); // used further down

   RCP<const FieldPattern> patternA = rcp(new Intrepid2FieldPattern(basisA));
   RCP<const FieldPattern> patternB = rcp(new Intrepid2FieldPattern(basisB));
   RCP<const FieldPattern> patternC = rcp(new Intrepid2FieldPattern(basisC)); // used further down

   TEST_EQUALITY(patternA->numberIds(),4);
   TEST_EQUALITY(patternB->numberIds(),4);

   std::vector<std::pair<int,int> > topEdgeIndices;
   orientation_helpers::computePatternEdgeIndices(*patternA,topEdgeIndices);

   TEST_EQUALITY(topEdgeIndices.size(),4);
   TEST_EQUALITY(topEdgeIndices[0].first,0); TEST_EQUALITY(topEdgeIndices[0].second,1);
   TEST_EQUALITY(topEdgeIndices[1].first,1); TEST_EQUALITY(topEdgeIndices[1].second,2);
   TEST_EQUALITY(topEdgeIndices[2].first,2); TEST_EQUALITY(topEdgeIndices[2].second,3);
   TEST_EQUALITY(topEdgeIndices[3].first,3); TEST_EQUALITY(topEdgeIndices[3].second,0);

   std::vector<std::vector<panzer::GlobalOrdinal>> connectivity(4);
   connectivity[0].resize(patternA->numberIds());
   connectivity[1].resize(patternA->numberIds());
   connectivity[2].resize(patternA->numberIds());
   connectivity[3].resize(patternA->numberIds());

   // Topologocally the four elements look like:
   //
   //    7-----6-----8
   //    |     |     |
   //    |     |     |
   //    3-----4-----5
   //    |     |     |
   //    |     |     |
   //    0-----1-----2
   //
   // all that matters is the global
   // node numbering and the local ordering
   // The local ordering is defined by the following connectivity
   connectivity[0][0] = 0; connectivity[0][1] = 1; connectivity[0][2] = 4; connectivity[0][3] = 3; 
   connectivity[1][0] = 1; connectivity[1][1] = 2; connectivity[1][2] = 5; connectivity[1][3] = 4; 
   connectivity[2][0] = 3; connectivity[2][1] = 4; connectivity[2][2] = 6; connectivity[2][3] = 7; 
   connectivity[3][0] = 4; connectivity[3][1] = 5; connectivity[3][2] = 8; connectivity[3][3] = 6; 

   // LOCAL element orientations are always set so that they flow in the positive
   // direction along an edge from node 0 to node 1. As a result if the GID of
   // node 0 is larger then node 1 then the GLOBAL orientation is -1 (and positive
   // otherwise. The local definition of the edge direction is defined by 
   // the shards cell topology.
   {
      std::vector<signed char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity[0], *patternB, orientations);
      TEST_EQUALITY(orientations[0],1);
      TEST_EQUALITY(orientations[1],1);
      TEST_EQUALITY(orientations[2],-1);
      TEST_EQUALITY(orientations[3],-1);
   }

   {
      std::vector<signed char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity[1], *patternB, orientations);
      TEST_EQUALITY(orientations[0],1);
      TEST_EQUALITY(orientations[1],1);
      TEST_EQUALITY(orientations[2],-1);
      TEST_EQUALITY(orientations[3],-1);
   }

   {
      std::vector<signed char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity[2], *patternB, orientations);
      TEST_EQUALITY(orientations[0],1);
      TEST_EQUALITY(orientations[1],1);
      TEST_EQUALITY(orientations[2],1);
      TEST_EQUALITY(orientations[3],-1);
   }

   {
      std::vector<signed char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity[3], *patternB, orientations);
      TEST_EQUALITY(orientations[0],1);
      TEST_EQUALITY(orientations[1],1);
      TEST_EQUALITY(orientations[2],-1);
      TEST_EQUALITY(orientations[3],-1);
   }

   // lets now test an aggregate basis
   std::vector<std::tuple<int,FieldType,Teuchos::RCP<const FieldPattern> > > patterns;
   patterns.push_back(std::make_tuple(0,FieldType::CG,patternB));
   patterns.push_back(std::make_tuple(1,FieldType::CG,patternC));
   patterns.push_back(std::make_tuple(2,FieldType::CG,patternA));
   Teuchos::RCP<const FieldPattern> aggPattern = Teuchos::rcp(new FieldAggPattern(patterns));

   std::vector<signed char> orientations(aggPattern->numberIds(),0);
   orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity[2], *aggPattern, orientations);
   int nonzeroCount = 0;
   for(std::size_t s=0;s<orientations.size();s++)
      nonzeroCount += orientations[s]*orientations[s]; // should be +1 only if it is an edge
   TEST_EQUALITY(nonzeroCount,8);

   // loop over edges
   signed char tests[] = { 1,1,1,-1}; // for element 2
   for(std::size_t e=0;e<4;e++) {
      const std::vector<int> & edgeIndices = aggPattern->getSubcellIndices(1,e);
      for(std::size_t s=0;s<edgeIndices.size();s++)
         TEST_EQUALITY(orientations[edgeIndices[s]],tests[e]);
   }       
}

/////////////////////////////////////////////
// 2D tests - face basis
/////////////////////////////////////////////

TEUCHOS_UNIT_TEST(tOrientation, testFaceBasis_tri)
{

   out << note << std::endl;

   // basis to build patterns from
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisA = rcp(new Intrepid2::Basis_HGRAD_TRI_C1_FEM<PHX::exec_space,double,double>);
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisB = rcp(new Intrepid2::Basis_HDIV_TRI_I1_FEM<PHX::exec_space,double,double>);
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisC = rcp(new Intrepid2::Basis_HGRAD_TRI_C2_FEM<PHX::exec_space,double,double>); // used further down

   RCP<const FieldPattern> patternA = rcp(new Intrepid2FieldPattern(basisA));
   RCP<const FieldPattern> patternB = rcp(new Intrepid2FieldPattern(basisB));
   RCP<const FieldPattern> patternC = rcp(new Intrepid2FieldPattern(basisC)); // used further down

   TEST_EQUALITY(patternA->numberIds(),3);
   TEST_EQUALITY(patternB->numberIds(),3);

   std::vector<std::pair<int,int> > topFaceIndices;
   orientation_helpers::computePatternEdgeIndices(*patternA,topFaceIndices);

   TEST_EQUALITY(topFaceIndices.size(),3);
   TEST_EQUALITY(topFaceIndices[0].first,0); TEST_EQUALITY(topFaceIndices[0].second,1);
   TEST_EQUALITY(topFaceIndices[1].first,1); TEST_EQUALITY(topFaceIndices[1].second,2);
   TEST_EQUALITY(topFaceIndices[2].first,2); TEST_EQUALITY(topFaceIndices[2].second,0);

   std::vector<std::vector<panzer::GlobalOrdinal>> connectivity(4);
   connectivity[0].resize(patternA->numberIds());
   connectivity[1].resize(patternA->numberIds());
   connectivity[2].resize(patternA->numberIds());
   connectivity[3].resize(patternA->numberIds());

   // Topologocally the four elements look like:
   //
   //    3-----1
   //    |\   /|
   //    | \ / |
   //    |  6  |
   //    | / \ |
   //    |/   \|
   //    5-----0
   //
   // all that matters is the global
   // node numbering and the local ordering
   // The local ordering is defined by the following connectivity
   connectivity[0][0] = 0; connectivity[0][1] = 6; connectivity[0][2] = 5; 
   connectivity[1][0] = 6; connectivity[1][1] = 0; connectivity[1][2] = 1; 
   connectivity[2][0] = 1; connectivity[2][1] = 3; connectivity[2][2] = 6; 
   connectivity[3][0] = 3; connectivity[3][1] = 5; connectivity[3][2] = 6; 

   // LOCAL element orientations are always set so that they flow in the positive
   // direction away from the cell center from node 0 to node 1. As a result if the GID of
   // node 0 is larger then node 1 then the GLOBAL orientation is -1 (and positive
   // otherwise. The local definition of the face direction is defined by 
   // the shards cell topology.
   {
      std::vector<signed char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topFaceIndices, connectivity[0], *patternB, orientations);
      TEST_EQUALITY(orientations[0],1);
      TEST_EQUALITY(orientations[1],-1);
      TEST_EQUALITY(orientations[2],-1);
   }

   {
      std::vector<signed char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topFaceIndices, connectivity[1], *patternB, orientations);
      TEST_EQUALITY(orientations[0],-1);
      TEST_EQUALITY(orientations[1],1);
      TEST_EQUALITY(orientations[2],1);
   }

   {
      std::vector<signed char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topFaceIndices, connectivity[2], *patternB, orientations);
      TEST_EQUALITY(orientations[0],1);
      TEST_EQUALITY(orientations[1],1);
      TEST_EQUALITY(orientations[2],-1);
   }

   {
      std::vector<signed char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topFaceIndices, connectivity[3], *patternB, orientations);
      TEST_EQUALITY(orientations[0],1);
      TEST_EQUALITY(orientations[1],1);
      TEST_EQUALITY(orientations[2],-1);
   }

   // lets now test an aggregate basis
   std::vector<std::tuple<int,FieldType,Teuchos::RCP<const FieldPattern> > > patterns;
   patterns.push_back(std::make_tuple(0,FieldType::CG,patternB));
   patterns.push_back(std::make_tuple(1,FieldType::CG,patternC));
   patterns.push_back(std::make_tuple(2,FieldType::CG,patternA));
   Teuchos::RCP<const FieldPattern> aggPattern = Teuchos::rcp(new FieldAggPattern(patterns));

   std::vector<signed char> orientations(aggPattern->numberIds(),0);
   orientation_helpers::computeCellEdgeOrientations(topFaceIndices, connectivity[2], *aggPattern, orientations);
   int nonzeroCount = 0;
   for(std::size_t s=0;s<orientations.size();s++)
      nonzeroCount += orientations[s]*orientations[s]; // should be +1 only if it is an face
   TEST_EQUALITY(nonzeroCount,6);

   // loop over faces
   signed char tests[] = { 1,1,-1}; // for element 2
   for(std::size_t e=0;e<3;e++) {
      const std::vector<int> & faceIndices = aggPattern->getSubcellIndices(1,e);
      for(std::size_t s=0;s<faceIndices.size();s++)
         TEST_EQUALITY(orientations[faceIndices[s]],tests[e]);
   }       
}

TEUCHOS_UNIT_TEST(tOrientation, testFaceBasis_quad)
{

   out << note << std::endl;

   // basis to build patterns from
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisA = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double>);
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisB = rcp(new Intrepid2::Basis_HDIV_QUAD_I1_FEM<PHX::exec_space,double,double>);
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisC = rcp(new Intrepid2::Basis_HGRAD_QUAD_C2_FEM<PHX::exec_space,double,double>); // used further down

   RCP<const FieldPattern> patternA = rcp(new Intrepid2FieldPattern(basisA));
   RCP<const FieldPattern> patternB = rcp(new Intrepid2FieldPattern(basisB));
   RCP<const FieldPattern> patternC = rcp(new Intrepid2FieldPattern(basisC)); // used further down

   TEST_EQUALITY(patternA->numberIds(),4);
   TEST_EQUALITY(patternB->numberIds(),4);

   std::vector<std::pair<int,int> > topFaceIndices;
   orientation_helpers::computePatternEdgeIndices(*patternA,topFaceIndices);

   TEST_EQUALITY(topFaceIndices.size(),4);
   TEST_EQUALITY(topFaceIndices[0].first,0); TEST_EQUALITY(topFaceIndices[0].second,1);
   TEST_EQUALITY(topFaceIndices[1].first,1); TEST_EQUALITY(topFaceIndices[1].second,2);
   TEST_EQUALITY(topFaceIndices[2].first,2); TEST_EQUALITY(topFaceIndices[2].second,3);
   TEST_EQUALITY(topFaceIndices[3].first,3); TEST_EQUALITY(topFaceIndices[3].second,0);

   std::vector<std::vector<panzer::GlobalOrdinal>> connectivity(4);
   connectivity[0].resize(patternA->numberIds());
   connectivity[1].resize(patternA->numberIds());
   connectivity[2].resize(patternA->numberIds());
   connectivity[3].resize(patternA->numberIds());

   // Topologocally the four elements look like:
   //
   //    7-----6-----8
   //    |     |     |
   //    |     |     |
   //    3-----4-----5
   //    |     |     |
   //    |     |     |
   //    0-----1-----2
   //
   // all that matters is the global
   // node numbering and the local ordering
   // The local ordering is defined by the following connectivity
   connectivity[0][0] = 0; connectivity[0][1] = 1; connectivity[0][2] = 4; connectivity[0][3] = 3; 
   connectivity[1][0] = 1; connectivity[1][1] = 2; connectivity[1][2] = 5; connectivity[1][3] = 4; 
   connectivity[2][0] = 3; connectivity[2][1] = 4; connectivity[2][2] = 6; connectivity[2][3] = 7; 
   connectivity[3][0] = 4; connectivity[3][1] = 5; connectivity[3][2] = 8; connectivity[3][3] = 6; 

   // LOCAL element orientations are always set so that they flow in the positive
   // direction along an face from node 0 to node 1. As a result if the GID of
   // node 0 is larger then node 1 then the GLOBAL orientation is -1 (and positive
   // otherwise. The local definition of the face direction is defined by 
   // the shards cell topology.
   {
      std::vector<signed char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topFaceIndices, connectivity[0], *patternB, orientations);
      TEST_EQUALITY(orientations[0],1);
      TEST_EQUALITY(orientations[1],1);
      TEST_EQUALITY(orientations[2],-1);
      TEST_EQUALITY(orientations[3],-1);
   }

   {
      std::vector<signed char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topFaceIndices, connectivity[1], *patternB, orientations);
      TEST_EQUALITY(orientations[0],1);
      TEST_EQUALITY(orientations[1],1);
      TEST_EQUALITY(orientations[2],-1);
      TEST_EQUALITY(orientations[3],-1);
   }

   {
      std::vector<signed char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topFaceIndices, connectivity[2], *patternB, orientations);
      TEST_EQUALITY(orientations[0],1);
      TEST_EQUALITY(orientations[1],1);
      TEST_EQUALITY(orientations[2],1);
      TEST_EQUALITY(orientations[3],-1);
   }

   {
      std::vector<signed char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topFaceIndices, connectivity[3], *patternB, orientations);
      TEST_EQUALITY(orientations[0],1);
      TEST_EQUALITY(orientations[1],1);
      TEST_EQUALITY(orientations[2],-1);
      TEST_EQUALITY(orientations[3],-1);
   }

   // lets now test an aggregate basis
   std::vector<std::tuple<int,FieldType,Teuchos::RCP<const FieldPattern> > > patterns;
   patterns.push_back(std::make_tuple(0,FieldType::CG,patternB));
   patterns.push_back(std::make_tuple(1,FieldType::CG,patternC));
   patterns.push_back(std::make_tuple(2,FieldType::CG,patternA));
   Teuchos::RCP<const FieldPattern> aggPattern = Teuchos::rcp(new FieldAggPattern(patterns));

   std::vector<signed char> orientations(aggPattern->numberIds(),0);
   orientation_helpers::computeCellEdgeOrientations(topFaceIndices, connectivity[2], *aggPattern, orientations);
   int nonzeroCount = 0;
   for(std::size_t s=0;s<orientations.size();s++)
      nonzeroCount += orientations[s]*orientations[s]; // should be +1 only if it is an face
   TEST_EQUALITY(nonzeroCount,8);

   // loop over faces
   signed char tests[] = { 1,1,1,-1}; // for element 2
   for(std::size_t e=0;e<4;e++) {
      const std::vector<int> & faceIndices = aggPattern->getSubcellIndices(1,e);
      for(std::size_t s=0;s<faceIndices.size();s++)
         TEST_EQUALITY(orientations[faceIndices[s]],tests[e]);
   }       
}

TEUCHOS_UNIT_TEST(tOrientation, testFaceBasis_tri2)
{

   out << note << std::endl;

   shards::CellTopology tri(shards::getCellTopologyData<shards::Triangle<3> >());

   // basis to build patterns from
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisA = rcp(new Intrepid2::Basis_HVOL_C0_FEM<PHX::exec_space,double,double>(tri));

   RCP<const FieldPattern> patternA = rcp(new Intrepid2FieldPattern(basisA));

   TEST_EQUALITY(patternA->numberIds(),1);

   std::vector<std::vector<int> > topFaceIndices;
   orientation_helpers::computePatternFaceIndices(*patternA,topFaceIndices);

   TEST_EQUALITY(topFaceIndices.size(),1);
   TEST_EQUALITY(topFaceIndices[0].size(),3); TEST_EQUALITY(topFaceIndices[0][0],0); TEST_EQUALITY(topFaceIndices[0][1],1); TEST_EQUALITY(topFaceIndices[0][2],2);
}

TEUCHOS_UNIT_TEST(tOrientation, testFaceBasis_quad2)
{

   out << note << std::endl;

   shards::CellTopology quad(shards::getCellTopologyData<shards::Quadrilateral<4> >()); 

   // basis to build patterns from
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisA = rcp(new Intrepid2::Basis_HVOL_C0_FEM<PHX::exec_space,double,double>(quad));

   RCP<const FieldPattern> patternA = rcp(new Intrepid2FieldPattern(basisA));

   TEST_EQUALITY(patternA->numberIds(),1);

   std::vector<std::vector<int> > topFaceIndices;
   orientation_helpers::computePatternFaceIndices(*patternA,topFaceIndices);

   TEST_EQUALITY(topFaceIndices.size(),1);
   TEST_EQUALITY(topFaceIndices[0].size(),4); 
   TEST_EQUALITY(topFaceIndices[0][0],0); TEST_EQUALITY(topFaceIndices[0][1],1); TEST_EQUALITY(topFaceIndices[0][2],2); TEST_EQUALITY(topFaceIndices[0][3],3); 
}

/////////////////////////////////////////////
// 3D tests - face basis
/////////////////////////////////////////////

TEUCHOS_UNIT_TEST(tOrientation, testFaceBasis_tet)
{

   out << note << std::endl;

   // basis to build patterns from
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisA = rcp(new Intrepid2::Basis_HGRAD_TET_C1_FEM<PHX::exec_space,double,double>);
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisB = rcp(new Intrepid2::Basis_HDIV_TET_I1_FEM<PHX::exec_space,double,double>);
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisC = rcp(new Intrepid2::Basis_HGRAD_TET_C2_FEM<PHX::exec_space,double,double>); // used further down

   RCP<const FieldPattern> patternA = rcp(new Intrepid2FieldPattern(basisA));
   RCP<const FieldPattern> patternB = rcp(new Intrepid2FieldPattern(basisB));
   RCP<const FieldPattern> patternC = rcp(new Intrepid2FieldPattern(basisC)); // used further down

   TEST_EQUALITY(patternA->numberIds(),4);
   TEST_EQUALITY(patternB->numberIds(),4);

   std::vector<std::vector<int> > topFaceIndices;
   orientation_helpers::computePatternFaceIndices(*patternA,topFaceIndices);

   TEST_EQUALITY(topFaceIndices.size(),4);
   TEST_EQUALITY(topFaceIndices[0].size(),3); TEST_EQUALITY(topFaceIndices[0][0],0); TEST_EQUALITY(topFaceIndices[0][1],1); TEST_EQUALITY(topFaceIndices[0][2],3);
   TEST_EQUALITY(topFaceIndices[1].size(),3); TEST_EQUALITY(topFaceIndices[1][0],1); TEST_EQUALITY(topFaceIndices[1][1],2); TEST_EQUALITY(topFaceIndices[1][2],3);
   TEST_EQUALITY(topFaceIndices[2].size(),3); TEST_EQUALITY(topFaceIndices[2][0],0); TEST_EQUALITY(topFaceIndices[2][1],3); TEST_EQUALITY(topFaceIndices[2][2],2);
   TEST_EQUALITY(topFaceIndices[3].size(),3); TEST_EQUALITY(topFaceIndices[3][0],0); TEST_EQUALITY(topFaceIndices[3][1],2); TEST_EQUALITY(topFaceIndices[3][2],1);

   // Topologically the first elements look like (shown by looking at each face), note
   // that the expected orientation is included as a +/- sign in the element
   //                                              //
   //      7          7          7          2      //
   //     / \        / \        / \        / \     //
   //    / - \      / - \      / + \      / + \    //
   //   /     \    /     \    /     \    /     \   //
   //  6 ----- 2  2 ----- 9  9 ----- 6  6 ----- 9  //
   //                                              //
   // all that matters is the global
   // node numbering and the local ordering

   // The local ordering is defined by the following connectivity
   std::vector<std::vector<panzer::GlobalOrdinal> > connectivity(1);
   connectivity[0].resize(patternA->numberIds());

   connectivity[0][0] = 6; connectivity[0][1] = 2; connectivity[0][2] = 9; connectivity[0][3] = 7; 

   {
      std::vector<signed char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellFaceOrientations(topFaceIndices, connectivity[0], *patternB, orientations);
      TEST_EQUALITY(orientations[0],-1);
      TEST_EQUALITY(orientations[1],-1);
      TEST_EQUALITY(orientations[2],1);
      TEST_EQUALITY(orientations[3],1);
   }
}

TEUCHOS_UNIT_TEST(tOrientation, testFaceBasis_hex)
{
   out << note << std::endl;

   // basis to build patterns from
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisA = rcp(new Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::exec_space,double,double>);
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisB = rcp(new Intrepid2::Basis_HDIV_HEX_I1_FEM<PHX::exec_space,double,double>);
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basisC = rcp(new Intrepid2::Basis_HGRAD_HEX_C2_FEM<PHX::exec_space,double,double>); // used further down

   RCP<const FieldPattern> patternA = rcp(new Intrepid2FieldPattern(basisA));
   RCP<const FieldPattern> patternB = rcp(new Intrepid2FieldPattern(basisB));
   RCP<const FieldPattern> patternC = rcp(new Intrepid2FieldPattern(basisC)); // used further down

   TEST_EQUALITY(patternA->numberIds(),8);
   TEST_EQUALITY(patternB->numberIds(),6);

   std::vector<std::vector<int> > topFaceIndices;
   orientation_helpers::computePatternFaceIndices(*patternA,topFaceIndices);

   // this checks to ensure that each face is oriented in a counter clockwise direction

   TEST_EQUALITY(topFaceIndices.size(),6);
   TEST_EQUALITY(topFaceIndices[0].size(),4); 
     TEST_EQUALITY(topFaceIndices[0][0],0); TEST_EQUALITY(topFaceIndices[0][1],1); TEST_EQUALITY(topFaceIndices[0][2],5); TEST_EQUALITY(topFaceIndices[0][3],4); 
   TEST_EQUALITY(topFaceIndices[1].size(),4);
     TEST_EQUALITY(topFaceIndices[1][0],1); TEST_EQUALITY(topFaceIndices[1][1],2); TEST_EQUALITY(topFaceIndices[1][2],6); TEST_EQUALITY(topFaceIndices[1][3],5); 
   TEST_EQUALITY(topFaceIndices[2].size(),4);
     TEST_EQUALITY(topFaceIndices[2][0],2); TEST_EQUALITY(topFaceIndices[2][1],3); TEST_EQUALITY(topFaceIndices[2][2],7); TEST_EQUALITY(topFaceIndices[2][3],6); 
   TEST_EQUALITY(topFaceIndices[3].size(),4);
     TEST_EQUALITY(topFaceIndices[3][0],0); TEST_EQUALITY(topFaceIndices[3][1],4); TEST_EQUALITY(topFaceIndices[3][2],7); TEST_EQUALITY(topFaceIndices[3][3],3); 
   TEST_EQUALITY(topFaceIndices[4].size(),4);
     TEST_EQUALITY(topFaceIndices[4][0],0); TEST_EQUALITY(topFaceIndices[4][1],3); TEST_EQUALITY(topFaceIndices[4][2],2); TEST_EQUALITY(topFaceIndices[4][3],1); 
   TEST_EQUALITY(topFaceIndices[5].size(),4);
     TEST_EQUALITY(topFaceIndices[5][0],4); TEST_EQUALITY(topFaceIndices[5][1],5); TEST_EQUALITY(topFaceIndices[5][2],6); TEST_EQUALITY(topFaceIndices[5][3],7); 

   // Topologically the first elements look like (shown by looking at each face), note
   // that the expected orientation is included as a +/- sign in the element
   // 
   //  0 ----- 8   8 ----- 9   9 ----- 1
   //  |       |   |       |   |       |
   //  |   +   |   |   +   |   |   -   |
   //  |       |   |       |   |       |
   //  5 ----- 2   2 ----- 6   6 ----- 7
   // 
   //  1 ----- 0   5 ----- 2   1 ----- 9
   //  |       |   |       |   |       |
   //  |   +   |   |   +   |   |   -   |
   //  |       |   |       |   |       |
   //  7 ----- 5   7 ----- 6   0 ----- 8
   //
   // all that matters is the global
   // node numbering and the local ordering

   // The local ordering is defined by the following connectivity
   std::vector<std::vector<panzer::GlobalOrdinal>> connectivity(1);
   connectivity[0].resize(patternA->numberIds());

   connectivity[0][0] = 5; connectivity[0][1] = 2; connectivity[0][2] = 6; connectivity[0][3] = 7; 
   connectivity[0][4] = 0; connectivity[0][5] = 8; connectivity[0][6] = 9; connectivity[0][7] = 1; 

   {
      std::vector<signed char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellFaceOrientations(topFaceIndices, connectivity[0], *patternB, orientations);
      TEST_EQUALITY(orientations[0],1);
      TEST_EQUALITY(orientations[1],1);
      TEST_EQUALITY(orientations[2],-1);
      TEST_EQUALITY(orientations[3],1);
      TEST_EQUALITY(orientations[4],1);
      TEST_EQUALITY(orientations[5],-1);
   }
}

}
