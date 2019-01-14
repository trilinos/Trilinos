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

#include <stk_mesh/base/GetBuckets.hpp>

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListExceptions.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"

#include "Shards_BasicTopologies.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

#ifdef PANZER_HAVE_IOSS

namespace panzer_stk {

TEUCHOS_UNIT_TEST(tExodusReaderFactory, basic_test)
{
   int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
   int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
   out << "Running numprocs = " << numprocs << " rank = " << rank << std::endl;

   std::vector<Teuchos::RCP<STK_ExodusReaderFactory> > facts;
   facts.push_back(Teuchos::rcp(new STK_ExodusReaderFactory("meshes/basic.gen")));
   facts.push_back(Teuchos::rcp(new STK_ExodusReaderFactory));

   Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
   pl->set("File Name","meshes/basic.gen");
   facts[1]->setParameterList(pl);

   out << "\n***reading from meshes/basic.gen ... writes to meshes/outputcheck.gen" << std::endl;
   for(std::size_t i=0;i<facts.size();i++) {
      {
         // read from file and build mesh
         Teuchos::RCP<STK_Interface> mesh = facts[i]->buildUncommitedMesh(MPI_COMM_WORLD);
         facts[i]->completeMeshConstruction(*mesh,MPI_COMM_WORLD);

         TEST_ASSERT(mesh!=Teuchos::null);
         TEST_ASSERT(mesh->getDimension()==2);
         TEST_ASSERT(mesh->isWritable());
         TEST_ASSERT(not mesh->isModifiable());

         out << "Begin writing to meshes/outputcheck.gen" << std::endl;
         mesh->writeToExodus("meshes/outputcheck.gen");
         out << "Finished writing to meshes/outputcheck.gen" << std::endl;

         // check element blocks
         std::vector<std::string> eBlocks;
         mesh->getElementBlockNames(eBlocks);
         TEST_EQUALITY((int) eBlocks.size(),2);
         out << "E-Blocks: ";
         for(std::size_t j=0;j<eBlocks.size();++j)
            out << "\"" << eBlocks[j] << "\" ";
         out << std::endl;

         // check side sets
         std::vector<std::string> sidesets;
         mesh->getSidesetNames(sidesets);
         TEST_EQUALITY((int) sidesets.size(),7);
         out << "Sides: ";
         for(std::size_t j=0;j<sidesets.size();++j)
            out << "\"" << sidesets[j] << "\" ";
         out << std::endl;

         // check node sets
         std::vector<std::string> nodesets;
         mesh->getNodesetNames(nodesets);
         TEST_EQUALITY((int) nodesets.size(),2);
         out << "Nodesets: ";
         for(std::size_t j=0;j<nodesets.size();++j)
            out << "\"" << nodesets[j] << "\" ";
         out << std::endl;

         TEST_EQUALITY(mesh->getSideRank(),mesh->getEdgeRank());
         TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),8);
         TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),22);
         TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),15);
      }

      // in an effort to be as cerebral as possible I read in the
      // outputed mesh and then re-output it
      out << "\n***reading from meshes/outputcheck.gen ... writes to meshes/outputcheck2.gen" << std::endl;
      {
         // read from file and build mesh
         Teuchos::RCP<STK_Interface> mesh = facts[i]->buildMesh(MPI_COMM_WORLD);

         // check element blocks
         std::vector<std::string> eBlocks;
         mesh->getElementBlockNames(eBlocks);
         TEST_EQUALITY((int) eBlocks.size(),2);
         out << "E-Blocks: ";
         for(std::size_t j=0;j<eBlocks.size();++j)
            out << "\"" << eBlocks[j] << "\" ";
         out << std::endl;

         // check side sets
         std::vector<std::string> sidesets;
         mesh->getSidesetNames(sidesets);
         TEST_EQUALITY((int) sidesets.size(),7);
         out << "Sides: ";
         for(std::size_t j=0;j<sidesets.size();++j)
            out << "\"" << sidesets[j] << "\" ";
         out << std::endl;

         // check node sets
         std::vector<std::string> nodesets;
         mesh->getNodesetNames(nodesets);
         TEST_EQUALITY((int) nodesets.size(),2);
         out << "Nodesets: ";
         for(std::size_t j=0;j<nodesets.size();++j)
            out << "\"" << nodesets[j] << "\" ";
         out << std::endl;

         mesh->writeToExodus("meshes/outputcheck2.gen");
         TEST_EQUALITY(mesh->getSideRank(),mesh->getEdgeRank());
         TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),8);
         TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),22);
         TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),15);
      }
   }
}

TEUCHOS_UNIT_TEST(tExodusReaderFactory, exo_scaling)
{
  {
    // These should correspond to the node coordinates, in order, in the
    // mesh file "basic.gen" as read above.
    double good_node_coords[15][2] = {{0.0, 0.0},
                                      {0.0, 0.5},
                                      {0.5, 0.5},
                                      {0.5, 0.0},
                                      {0.0, 1.0},
                                      {0.5, 1.0},
                                      {1.0, 0.5},
                                      {1.0, 0.0},
                                      {1.0, 1.0},
                                      {1.5, 0.5},
                                      {1.5, 0.0},
                                      {1.5, 1.0},
                                      {2.0, 0.5},
                                      {2.0, 0.0},
                                      {2.0, 1.0}};

    STK_ExodusReaderFactory factory;

    Teuchos::RCP<Teuchos::ParameterList> inp_pl = Teuchos::rcp(new Teuchos::ParameterList);
    inp_pl->set("File Name", "meshes/basic.gen");
    inp_pl->set("Scale Factor", 0.5);

    TEST_NOTHROW(factory.setParameterList(inp_pl));

    Teuchos::RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    TEST_ASSERT(mesh!=Teuchos::null);

    // Make sure the node coordinates as they exist in the data
    // structure have been scaled by the 1/"Scale Factor" from above.
    double sf = 2.0;  // 1/(scale factor)
    for (stk::mesh::EntityId id=1; id <= 15; ++id)
    {

      stk::mesh::Entity node = mesh->getBulkData()->get_entity(mesh->getNodeRank(), id);
      if (mesh->isValid(node))
      {
        double const* coords = mesh->getNodeCoordinates(id);
        TEST_EQUALITY(coords[0], sf*good_node_coords[id-1][0]);
        TEST_EQUALITY(coords[1], sf*good_node_coords[id-1][1]);
      }
    }
  }

}

TEUCHOS_UNIT_TEST(tExodusReaderFactory, periodic_bc)
{
   // correct setting of parameter list
   {
      STK_ExodusReaderFactory factory;

      Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
      pl->set("File Name","meshes/basic.gen");
      TEST_NOTHROW(factory.setParameterList(pl));

      Teuchos::RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
      TEST_ASSERT(mesh!=Teuchos::null);
      TEST_EQUALITY(mesh->getPeriodicBCVector().size(),0);
   }

   {
      STK_ExodusReaderFactory factory;
      Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
      pl->set("File Name","meshes/basic.gen");
      Teuchos::ParameterList & pbcs = pl->sublist("Periodic BCs");
      pbcs.set<int>("Count",1);
      pbcs.set("Periodic Condition 1","x-coord left;right");

      TEST_NOTHROW(factory.setParameterList(pl));

      Teuchos::RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
      TEST_ASSERT(mesh!=Teuchos::null);
      TEST_EQUALITY(mesh->getPeriodicBCVector().size(),1);
   }
}

TEUCHOS_UNIT_TEST(tExodusReaderFactory, parameter_list_construction)
{
   // correct setting of parameter list
   {
      STK_ExodusReaderFactory factory;

      Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
      pl->set("File Name","meshes/basic.gen");
      TEST_NOTHROW(factory.setParameterList(pl));
   }

   // first incorrect paramter list ... extra parameteres
   {
      STK_ExodusReaderFactory factory;

      Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
      pl->set("File Name","meshes/basic.gen");
      pl->set("Foo","Bar");
      TEST_THROW(factory.setParameterList(pl),Teuchos::Exceptions::InvalidParameter);
   }

   // second incorrect paramter list ... no file name
   {
      STK_ExodusReaderFactory factory;

      Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
      pl->set("No File Name","meshes/basic.gen");
      TEST_THROW(factory.setParameterList(pl),Teuchos::Exceptions::InvalidParameter);
   }
}

}

#endif
