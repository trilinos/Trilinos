// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <string>
#include <iostream>

#include "Panzer_FieldAggPattern.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_GeometricAggFieldPattern.hpp"
#include "Panzer_UniqueGlobalIndexer_Utilities.hpp"

// include some intrepid basis functions
// 2D basis 
#include "Intrepid_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid_HCURL_TRI_I1_FEM.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HCURL_QUAD_I1_FEM.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {

std::string note = "***   NOTE: UNIT TEST BASED ON MARCH 2012  ***\n"
                   "***   INTREPID AND SHARDS Trilinos-dev     ***\n"
                   "***   DOXYGEN WEBSITE                      ***\n";
typedef Intrepid::FieldContainer<double> FieldContainer;

/////////////////////////////////////////////
// 2D tests
/////////////////////////////////////////////

TEUCHOS_UNIT_TEST(tOrientation, testEdgeBasis_tri)
{
   out << note << std::endl;

   // basis to build patterns from
   RCP<Intrepid::Basis<double,FieldContainer> > basisA = rcp(new Intrepid::Basis_HGRAD_TRI_C1_FEM<double,FieldContainer>);
   RCP<Intrepid::Basis<double,FieldContainer> > basisB = rcp(new Intrepid::Basis_HCURL_TRI_I1_FEM<double,FieldContainer>);
   RCP<Intrepid::Basis<double,FieldContainer> > basisC = rcp(new Intrepid::Basis_HGRAD_TRI_C2_FEM<double,FieldContainer>); // used further down

   RCP<const FieldPattern> patternA = rcp(new IntrepidFieldPattern(basisA));
   RCP<const FieldPattern> patternB = rcp(new IntrepidFieldPattern(basisB));
   RCP<const FieldPattern> patternC = rcp(new IntrepidFieldPattern(basisC)); // used further down

   TEST_EQUALITY(patternA->numberIds(),3);
   TEST_EQUALITY(patternB->numberIds(),3);

   std::vector<std::pair<int,int> > topEdgeIndices;
   orientation_helpers::computePatternEdgeIndices(*patternA,topEdgeIndices);

   TEST_EQUALITY(topEdgeIndices.size(),3);
   TEST_EQUALITY(topEdgeIndices[0].first,0); TEST_EQUALITY(topEdgeIndices[0].second,1);
   TEST_EQUALITY(topEdgeIndices[1].first,1); TEST_EQUALITY(topEdgeIndices[1].second,2);
   TEST_EQUALITY(topEdgeIndices[2].first,2); TEST_EQUALITY(topEdgeIndices[2].second,0);

   std::vector<std::vector<long> > connectivity(4);
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
      std::vector<char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity[0], *patternB, orientations);
      TEST_EQUALITY(orientations[0],char(1));
      TEST_EQUALITY(orientations[1],char(-1));
      TEST_EQUALITY(orientations[2],char(-1));
   }

   {
      std::vector<char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity[1], *patternB, orientations);
      TEST_EQUALITY(orientations[0],char(-1));
      TEST_EQUALITY(orientations[1],char(1));
      TEST_EQUALITY(orientations[2],char(1));
   }

   {
      std::vector<char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity[2], *patternB, orientations);
      TEST_EQUALITY(orientations[0],char(1));
      TEST_EQUALITY(orientations[1],char(1));
      TEST_EQUALITY(orientations[2],char(-1));
   }

   {
      std::vector<char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity[3], *patternB, orientations);
      TEST_EQUALITY(orientations[0],char(1));
      TEST_EQUALITY(orientations[1],char(1));
      TEST_EQUALITY(orientations[2],char(-1));
   }

   // lets now test an aggregate basis
   std::vector<std::pair<int,Teuchos::RCP<const FieldPattern> > > patterns;
   patterns.push_back(std::make_pair(0,patternB));
   patterns.push_back(std::make_pair(1,patternC));
   patterns.push_back(std::make_pair(2,patternA));
   Teuchos::RCP<const FieldPattern> aggPattern = Teuchos::rcp(new FieldAggPattern(patterns));

   std::vector<char> orientations(aggPattern->numberIds(),0);
   orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity[2], *aggPattern, orientations);
   int nonzeroCount = 0;
   for(std::size_t s=0;s<orientations.size();s++)
      nonzeroCount += orientations[s]*orientations[s]; // should be +1 only if it is an edge
   TEST_EQUALITY(nonzeroCount,6);

   // loop over edges
   char tests[] = { 1,1,-1}; // for element 2
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
   RCP<Intrepid::Basis<double,FieldContainer> > basisA = rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer>);
   RCP<Intrepid::Basis<double,FieldContainer> > basisB = rcp(new Intrepid::Basis_HCURL_QUAD_I1_FEM<double,FieldContainer>);
   RCP<Intrepid::Basis<double,FieldContainer> > basisC = rcp(new Intrepid::Basis_HGRAD_QUAD_C2_FEM<double,FieldContainer>); // used further down

   RCP<const FieldPattern> patternA = rcp(new IntrepidFieldPattern(basisA));
   RCP<const FieldPattern> patternB = rcp(new IntrepidFieldPattern(basisB));
   RCP<const FieldPattern> patternC = rcp(new IntrepidFieldPattern(basisC)); // used further down

   TEST_EQUALITY(patternA->numberIds(),4);
   TEST_EQUALITY(patternB->numberIds(),4);

   std::vector<std::pair<int,int> > topEdgeIndices;
   orientation_helpers::computePatternEdgeIndices(*patternA,topEdgeIndices);

   TEST_EQUALITY(topEdgeIndices.size(),4);
   TEST_EQUALITY(topEdgeIndices[0].first,0); TEST_EQUALITY(topEdgeIndices[0].second,1);
   TEST_EQUALITY(topEdgeIndices[1].first,1); TEST_EQUALITY(topEdgeIndices[1].second,2);
   TEST_EQUALITY(topEdgeIndices[2].first,2); TEST_EQUALITY(topEdgeIndices[2].second,3);
   TEST_EQUALITY(topEdgeIndices[3].first,3); TEST_EQUALITY(topEdgeIndices[3].second,0);

   std::vector<std::vector<long> > connectivity(4);
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
      std::vector<char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity[0], *patternB, orientations);
      TEST_EQUALITY(orientations[0],char(1));
      TEST_EQUALITY(orientations[1],char(1));
      TEST_EQUALITY(orientations[2],char(-1));
      TEST_EQUALITY(orientations[3],char(-1));
   }

   {
      std::vector<char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity[1], *patternB, orientations);
      TEST_EQUALITY(orientations[0],char(1));
      TEST_EQUALITY(orientations[1],char(1));
      TEST_EQUALITY(orientations[2],char(-1));
      TEST_EQUALITY(orientations[3],char(-1));
   }

   {
      std::vector<char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity[2], *patternB, orientations);
      TEST_EQUALITY(orientations[0],char(1));
      TEST_EQUALITY(orientations[1],char(1));
      TEST_EQUALITY(orientations[2],char(1));
      TEST_EQUALITY(orientations[3],char(-1));
   }

   {
      std::vector<char> orientations(patternB->numberIds(),0);
      orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity[3], *patternB, orientations);
      TEST_EQUALITY(orientations[0],char(1));
      TEST_EQUALITY(orientations[1],char(1));
      TEST_EQUALITY(orientations[2],char(-1));
      TEST_EQUALITY(orientations[3],char(-1));
   }

   // lets now test an aggregate basis
   std::vector<std::pair<int,Teuchos::RCP<const FieldPattern> > > patterns;
   patterns.push_back(std::make_pair(0,patternB));
   patterns.push_back(std::make_pair(1,patternC));
   patterns.push_back(std::make_pair(2,patternA));
   Teuchos::RCP<const FieldPattern> aggPattern = Teuchos::rcp(new FieldAggPattern(patterns));

   std::vector<char> orientations(aggPattern->numberIds(),0);
   orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity[2], *aggPattern, orientations);
   int nonzeroCount = 0;
   for(std::size_t s=0;s<orientations.size();s++)
      nonzeroCount += orientations[s]*orientations[s]; // should be +1 only if it is an edge
   TEST_EQUALITY(nonzeroCount,8);

   // loop over edges
   char tests[] = { 1,1,1,-1}; // for element 2
   for(std::size_t e=0;e<4;e++) {
      const std::vector<int> & edgeIndices = aggPattern->getSubcellIndices(1,e);
      for(std::size_t s=0;s<edgeIndices.size();s++)
         TEST_EQUALITY(orientations[edgeIndices[s]],tests[e]);
   }       
}

}
