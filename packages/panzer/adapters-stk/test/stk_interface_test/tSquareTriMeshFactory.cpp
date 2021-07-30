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
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareTriMeshFactory.hpp"

#include "Shards_BasicTopologies.hpp"

#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Selector.hpp"

#include "Ioss_DatabaseIO.h"
#include "Ioss_IOFactory.h"
#include "Ioss_Region.h"
#include "Ioss_EdgeBlock.h"

namespace panzer_stk {

void edge_block_test_helper(Teuchos::FancyOStream &out,
                            bool &success,
                            Teuchos::RCP<Teuchos::ParameterList> pl,
                            std::string exodus_filename,
                            uint32_t expected_edge_block_count)
{
   SquareTriMeshFactory factory; 
   factory.setParameterList(pl);
   Teuchos::RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);
 
   if(mesh->isWritable())
      mesh->writeToExodus(exodus_filename.c_str());

   {
   Ioss::DatabaseIO *db_io = Ioss::IOFactory::create("exodus", 
                                                     exodus_filename.c_str(), 
                                                     Ioss::READ_MODEL);
   TEST_ASSERT(db_io);

   Ioss::Region region(db_io);
   TEST_ASSERT(db_io->ok() == true);

   auto all_edge_blocks = region.get_edge_blocks();
   TEST_ASSERT(all_edge_blocks.size() == expected_edge_block_count);
   }
}

TEUCHOS_UNIT_TEST(tSquareTriMeshFactory, defaults)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;
   
   SquareTriMeshFactory factory; 
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
 
   if(mesh->isWritable())
      mesh->writeToExodus("square-tri.exo");

   // minimal requirements
   TEST_ASSERT(mesh!=Teuchos::null);
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),4);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),2*25);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),25+60);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),36);

   int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
   int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);

   int mpi_numprocs = -1;
   MPI_Comm_size(MPI_COMM_WORLD, &mpi_numprocs);
   int mpi_rank = -1;
   MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
   TEST_EQUALITY(numprocs,mpi_numprocs);
   TEST_EQUALITY(rank,mpi_rank);

   // check for nodeset
   std::vector<std::string> nodesets;
   mesh->getNodesetNames(nodesets);
 
   TEST_EQUALITY(nodesets.size(),1);
   TEST_EQUALITY(nodesets[0],"origin");
}

TEUCHOS_UNIT_TEST(tSquareTriMeshFactory, default_edge_face_blocks)
{
   using Teuchos::RCP;

   int xe = 2, ye = 2;
   int bx = 1, by = 1;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",bx);
   pl->set("Y Blocks",by);
   pl->set("X Elements",xe);
   pl->set("Y Elements",ye);

   edge_block_test_helper(out, success, pl, "SquareTri_default_edge_blocks.exo", 0);
}

TEUCHOS_UNIT_TEST(tSquareTriMeshFactory, create_edge_blocks_pl)
{
   using Teuchos::RCP;

   int xe = 2, ye = 2;
   int bx = 1, by = 1;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",bx);
   pl->set("Y Blocks",by);
   pl->set("X Elements",xe);
   pl->set("Y Elements",ye);
   pl->set("Create Edge Blocks",true);

   edge_block_test_helper(out, success, pl, "SquareTri_create_edge_blocks_pl.exo", 1);
}

}
