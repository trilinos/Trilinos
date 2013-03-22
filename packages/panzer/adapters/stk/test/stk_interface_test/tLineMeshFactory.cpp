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
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_LineMeshFactory.hpp"

#include "Shards_BasicTopologies.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

namespace panzer_stk {

TEUCHOS_UNIT_TEST(tLineMeshFactory, defaults)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;
   
   LineMeshFactory factory; 
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);

   TEST_EQUALITY(mesh->getPeriodicBCVector().size(),0);
 
   if(mesh->isWritable());
      mesh->writeToExodus("Line.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),1);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),2);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),5);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),6);
}

TEUCHOS_UNIT_TEST(tLineMeshFactory, element_counts)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",1);
   pl->set("X Elements",2);
   
   LineMeshFactory factory; 
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);
 
   if(mesh->isWritable());
      mesh->writeToExodus("Line_oddelmt.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),1);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),2);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),2);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),3);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),3);
}

TEUCHOS_UNIT_TEST(tLineMeshFactory, allblock)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   int xe = 4;
   int bx = 4;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",bx);
   pl->set("X Elements",xe);

   xe *= bx;
   
   LineMeshFactory factory; 
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);
 
   if(mesh->isWritable());
      mesh->writeToExodus("Line_allblock.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),1);
   TEST_EQUALITY(mesh->getNumElementBlocks(),(std::size_t) bx);
   TEST_EQUALITY(mesh->getNumSidesets(),2);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),(std::size_t) xe);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),(std::size_t) xe+1);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(std::size_t) xe+1);
}

TEUCHOS_UNIT_TEST(tLineMeshFactory, two_block)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",2);
   pl->set("X Elements",5);
   
   LineMeshFactory factory; 
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);
 
   if(mesh->isWritable());
      mesh->writeToExodus("Line_2block.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),1);
   TEST_EQUALITY(mesh->getNumElementBlocks(),2);
   TEST_EQUALITY(mesh->getNumSidesets(),2);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),10);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),11);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),11);
}

}
