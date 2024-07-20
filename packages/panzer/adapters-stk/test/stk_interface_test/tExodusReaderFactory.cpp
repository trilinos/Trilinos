// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#include <stk_mesh/base/Comm.hpp>

#include "Ioss_DatabaseIO.h"
#include "Ioss_IOFactory.h"
#include "Ioss_Region.h"
#include "Ioss_EdgeBlock.h"
#include "Ioss_FaceBlock.h"

#ifdef PANZER_HAVE_IOSS

namespace panzer_stk {

void edge_face_block_test_helper(Teuchos::FancyOStream &out,
                                 bool &success,
                                 std::string exodus_filename,
                                 uint32_t expected_edge_block_count,
                                 uint32_t expected_face_block_count)
{
   Ioss::DatabaseIO *db_io = Ioss::IOFactory::create("exodus",
                                                     exodus_filename.c_str(),
                                                     Ioss::READ_MODEL);
   TEST_ASSERT(db_io);

   Ioss::Region region(db_io);
   TEST_ASSERT(db_io->ok() == true);

   auto all_edge_blocks = region.get_edge_blocks();
   TEST_ASSERT(all_edge_blocks.size() == expected_edge_block_count);
   auto all_face_blocks = region.get_face_blocks();
   TEST_ASSERT(all_face_blocks.size() == expected_face_block_count);
}


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
   pl->set("File Type","Exodus");
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

/*
 * This is a much simplified copy of the "basic_test"
 * which confirms that by default the edge and face
 * blocks are NOT created when reading in an Exodus
 * file that doesn't already have edge or face blocks.
*/
TEUCHOS_UNIT_TEST(tExodusReaderFactory, default_edge_face_block_test)
{
   auto erf = Teuchos::rcp(new STK_ExodusReaderFactory());

   Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
   pl->set("File Name","meshes/basic3d.gen");
   erf->setParameterList(pl);

   // read from file and build mesh
   Teuchos::RCP<STK_Interface> mesh = erf->buildUncommitedMesh(MPI_COMM_WORLD);
   erf->completeMeshConstruction(*mesh,MPI_COMM_WORLD);

   TEST_ASSERT(mesh!=Teuchos::null);
   TEST_ASSERT(mesh->getDimension()==3);
   TEST_ASSERT(mesh->isWritable());
   TEST_ASSERT(not mesh->isModifiable());

   mesh->writeToExodus("meshes/default_edge_face_block_check.gen");

   // check edge blocks
   std::vector<std::string> edgeblocks;
   mesh->getEdgeBlockNames(edgeblocks);
   TEST_EQUALITY((int) edgeblocks.size(),0);

   // check face blocks
   std::vector<std::string> faceblocks;
   mesh->getFaceBlockNames(faceblocks);
   TEST_EQUALITY((int) faceblocks.size(),0);

   edge_face_block_test_helper(out,
                               success,
                               "meshes/default_edge_face_block_check.gen",
                               0,
                               0);
}

/*
 * This is a much simplified copy of the "basic_test"
 * which confirms that the edge block is created in
 * step 1 and copied in step 2.
*/
TEUCHOS_UNIT_TEST(tExodusReaderFactory, edge_block_test)
{
   {
   auto erf = Teuchos::rcp(new STK_ExodusReaderFactory());

   Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
   pl->set("File Name","meshes/basic.gen");
   pl->set("Create Edge Blocks",true);
   erf->setParameterList(pl);

   // read from file and build mesh
   Teuchos::RCP<STK_Interface> mesh = erf->buildUncommitedMesh(MPI_COMM_WORLD);
   erf->completeMeshConstruction(*mesh,MPI_COMM_WORLD);

   TEST_ASSERT(mesh!=Teuchos::null);
   TEST_ASSERT(mesh->getDimension()==2);
   TEST_ASSERT(mesh->isWritable());
   TEST_ASSERT(not mesh->isModifiable());

   mesh->writeToExodus("meshes/edge_block_check.gen");

   // check edge blocks
   std::vector<std::string> edgeblocks;
   mesh->getEdgeBlockNames(edgeblocks);
   TEST_EQUALITY((int) edgeblocks.size(),1);

   edge_face_block_test_helper(out,
                               success,
                               "meshes/edge_block_check.gen",
                               1,
                               0);
   }
   {
   // in an effort to be as cerebral as possible I read in the
   // outputed mesh and then re-output it

   // read from file and build mesh
   auto erf = Teuchos::rcp(new STK_ExodusReaderFactory());

   Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
   pl->set("File Name","meshes/edge_block_check.gen");
   pl->set("Create Edge Blocks",true);
   erf->setParameterList(pl);

   Teuchos::RCP<STK_Interface> mesh = erf->buildMesh(MPI_COMM_WORLD);

   // check edge blocks
   std::vector<std::string> edgeblocks;
   mesh->getEdgeBlockNames(edgeblocks);
   TEST_EQUALITY((int) edgeblocks.size(),1);

   mesh->writeToExodus("meshes/edge_block_check2.gen");

   edge_face_block_test_helper(out,
                               success,
                               "meshes/edge_block_check2.gen",
                               1,
                               0);
   }
}

/*
 * This is a much simplified copy of the "basic_test"
 * which confirms that the face block is created in
 * step 1 and copied in step 2.
*/
TEUCHOS_UNIT_TEST(tExodusReaderFactory, face_block_test)
{
   {
   auto erf = Teuchos::rcp(new STK_ExodusReaderFactory());

   Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
   pl->set("File Name","meshes/basic3d.gen");
   pl->set("Create Face Blocks",true);
   erf->setParameterList(pl);

   // read from file and build mesh
   Teuchos::RCP<STK_Interface> mesh = erf->buildUncommitedMesh(MPI_COMM_WORLD);
   erf->completeMeshConstruction(*mesh,MPI_COMM_WORLD);

   TEST_ASSERT(mesh!=Teuchos::null);
   TEST_ASSERT(mesh->getDimension()==3);
   TEST_ASSERT(mesh->isWritable());
   TEST_ASSERT(not mesh->isModifiable());

   mesh->writeToExodus("meshes/face_block_check.gen");

   // check face blocks
   std::vector<std::string> faceblocks;
   mesh->getFaceBlockNames(faceblocks);
   TEST_EQUALITY((int) faceblocks.size(),1);

   edge_face_block_test_helper(out,
                               success,
                               "meshes/face_block_check.gen",
                               0,
                               1);
   }
   {
   // in an effort to be as cerebral as possible I read in the
   // outputed mesh and then re-output it

   // read from file and build mesh
   auto erf = Teuchos::rcp(new STK_ExodusReaderFactory());

   Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
   pl->set("File Name","meshes/face_block_check.gen");
   pl->set("Create Face Blocks",true);
   erf->setParameterList(pl);

   Teuchos::RCP<STK_Interface> mesh = erf->buildMesh(MPI_COMM_WORLD);

   // check face blocks
   std::vector<std::string> faceblocks;
   mesh->getFaceBlockNames(faceblocks);
   TEST_EQUALITY((int) faceblocks.size(),1);

   mesh->writeToExodus("meshes/face_block_check2.gen");

   edge_face_block_test_helper(out,
                               success,
                               "meshes/face_block_check2.gen",
                               0,
                               1);
   }
}

/*
 * This is a much simplified copy of the "basic_test"
 * which confirms that the face block is created in
 * step 1 and copied in step 2.
*/
TEUCHOS_UNIT_TEST(tExodusReaderFactory, multiblock_edge_block_phase1_test)
{
   {
   auto erf = Teuchos::rcp(new STK_ExodusReaderFactory());

   Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
   pl->set("File Name","meshes/twoblock_cube.gen");
   pl->set("Create Edge Blocks",true);
   erf->setParameterList(pl);

   // read from file and build mesh
   Teuchos::RCP<STK_Interface> mesh = erf->buildUncommitedMesh(MPI_COMM_WORLD);
   erf->completeMeshConstruction(*mesh,MPI_COMM_WORLD);

   TEST_ASSERT(mesh!=Teuchos::null);
   TEST_ASSERT(mesh->getDimension()==3);
   TEST_ASSERT(mesh->isWritable());
   TEST_ASSERT(not mesh->isModifiable());

   mesh->writeToExodus("meshes/multiblock_edge_block_check.gen");

   // check edge blocks
   std::vector<std::string> edgeblocks;
   mesh->getEdgeBlockNames(edgeblocks);
   TEST_EQUALITY((int) edgeblocks.size(),1);

   edge_face_block_test_helper(out,
                               success,
                               "meshes/multiblock_edge_block_check.gen",
                               1,
                               0);
   }
}
TEUCHOS_UNIT_TEST(tExodusReaderFactory, multiblock_edge_block_phase2_test)
{
   {
   // in an effort to be as cerebral as possible I read in the
   // outputed mesh and then re-output it

   // read from file and build mesh
   auto erf = Teuchos::rcp(new STK_ExodusReaderFactory());

   Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
   pl->set("File Name","meshes/multiblock_edge_block_check.gen");
   pl->set("Create Edge Blocks",true);
   erf->setParameterList(pl);

   Teuchos::RCP<STK_Interface> mesh = erf->buildMesh(MPI_COMM_WORLD);

   // check edge blocks
   std::vector<std::string> edgeblocks;
   mesh->getEdgeBlockNames(edgeblocks);
   TEST_EQUALITY((int) edgeblocks.size(),1);

   mesh->writeToExodus("meshes/multiblock_edge_block_check2.gen");

   edge_face_block_test_helper(out,
                               success,
                               "meshes/multiblock_edge_block_check2.gen",
                               1,
                               0);
   }
}

/*
 * This is a much simplified copy of the "basic_test"
 * which confirms that the face block is created in
 * step 1 and copied in step 2.
*/
TEUCHOS_UNIT_TEST(tExodusReaderFactory, multiblock_face_block_phase1_test)
{
   {
   auto erf = Teuchos::rcp(new STK_ExodusReaderFactory());

   Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
   pl->set("File Name","meshes/twoblock_cube.gen");
   pl->set("Create Face Blocks",true);
   erf->setParameterList(pl);

   // read from file and build mesh
   Teuchos::RCP<STK_Interface> mesh = erf->buildUncommitedMesh(MPI_COMM_WORLD);
   erf->completeMeshConstruction(*mesh,MPI_COMM_WORLD);

   TEST_ASSERT(mesh!=Teuchos::null);
   TEST_ASSERT(mesh->getDimension()==3);
   TEST_ASSERT(mesh->isWritable());
   TEST_ASSERT(not mesh->isModifiable());

   mesh->writeToExodus("meshes/multiblock_face_block_check.gen");

   // check face blocks
   std::vector<std::string> faceblocks;
   mesh->getFaceBlockNames(faceblocks);
   TEST_EQUALITY((int) faceblocks.size(),2);

   edge_face_block_test_helper(out,
                               success,
                               "meshes/multiblock_face_block_check.gen",
                               0,
                               2);
   }
}
TEUCHOS_UNIT_TEST(tExodusReaderFactory, multiblock_face_block_phase2_test)
{
   {
   // in an effort to be as cerebral as possible I read in the
   // outputed mesh and then re-output it

   // read from file and build mesh
   auto erf = Teuchos::rcp(new STK_ExodusReaderFactory());

   Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
   pl->set("File Name","meshes/multiblock_face_block_check.gen");
   pl->set("Create Face Blocks",true);
   erf->setParameterList(pl);

   Teuchos::RCP<STK_Interface> mesh = erf->buildMesh(MPI_COMM_WORLD);

   // check face blocks
   std::vector<std::string> faceblocks;
   mesh->getFaceBlockNames(faceblocks);
   TEST_EQUALITY((int) faceblocks.size(),2);

   mesh->writeToExodus("meshes/multiblock_face_block_check2.gen");

   edge_face_block_test_helper(out,
                               success,
                               "meshes/multiblock_face_block_check2.gen",
                               0,
                               2);
   }
}

/*
 * This is a much simplified copy of the "basic_test"
 * which confirms that the face block is created in
 * step 1 and copied in step 2.
*/
TEUCHOS_UNIT_TEST(tExodusReaderFactory, multiblock_edge_face_block_phase1_test)
{
   {
   auto erf = Teuchos::rcp(new STK_ExodusReaderFactory());

   Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
   pl->set("File Name","meshes/twoblock_cube.gen");
   pl->set("Create Edge Blocks",true);
   pl->set("Create Face Blocks",true);
   erf->setParameterList(pl);

   // read from file and build mesh
   Teuchos::RCP<STK_Interface> mesh = erf->buildUncommitedMesh(MPI_COMM_WORLD);
   erf->completeMeshConstruction(*mesh,MPI_COMM_WORLD);

   TEST_ASSERT(mesh!=Teuchos::null);
   TEST_ASSERT(mesh->getDimension()==3);
   TEST_ASSERT(mesh->isWritable());
   TEST_ASSERT(not mesh->isModifiable());

   mesh->writeToExodus("meshes/multiblock_edge_face_block_check.gen");

   // check edge blocks
   std::vector<std::string> edgeblocks;
   mesh->getEdgeBlockNames(edgeblocks);
   TEST_EQUALITY((int) edgeblocks.size(),1);
   // check face blocks
   std::vector<std::string> faceblocks;
   mesh->getFaceBlockNames(faceblocks);
   TEST_EQUALITY((int) faceblocks.size(),2);

   edge_face_block_test_helper(out,
                               success,
                               "meshes/multiblock_edge_face_block_check.gen",
                               1,
                               2);
   }
}
TEUCHOS_UNIT_TEST(tExodusReaderFactory, multiblock_edge_face_block_phase2_test)
{
   {
   // in an effort to be as cerebral as possible I read in the
   // outputed mesh and then re-output it

   // read from file and build mesh
   auto erf = Teuchos::rcp(new STK_ExodusReaderFactory());

   Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
   pl->set("File Name","meshes/multiblock_edge_face_block_check.gen");
   pl->set("Create Edge Blocks",true);
   pl->set("Create Face Blocks",true);
   erf->setParameterList(pl);

   Teuchos::RCP<STK_Interface> mesh = erf->buildMesh(MPI_COMM_WORLD);

   // check edge blocks
   std::vector<std::string> edgeblocks;
   mesh->getEdgeBlockNames(edgeblocks);
   TEST_EQUALITY((int) edgeblocks.size(),1);
   // check face blocks
   std::vector<std::string> faceblocks;
   mesh->getFaceBlockNames(faceblocks);
   TEST_EQUALITY((int) faceblocks.size(),2);

   mesh->writeToExodus("meshes/multiblock_edge_face_block_check2.gen");

   edge_face_block_test_helper(out,
                               success,
                               "meshes/multiblock_edge_face_block_check2.gen",
                               1,
                               2);
   }
}

TEUCHOS_UNIT_TEST(tExodusReaderFactory, getMeshDimension)
{
   TEST_EQUALITY(panzer_stk::getMeshDimension("meshes/basic.gen",MPI_COMM_WORLD,"Exodus"),2);
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

#ifdef PANZER_HAVE_PERCEPT
TEUCHOS_UNIT_TEST(tExodusReaderFactory, percept)
{
  STK_ExodusReaderFactory factory;

  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
  pl->set("File Name","meshes/basic.gen");
  pl->set("Levels of Uniform Refinement",1);
  TEST_NOTHROW(factory.setParameterList(pl));

  Teuchos::RCP<STK_Interface> mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);
  const bool delayCommitForRefinement = true;
  mesh->initialize(MPI_COMM_WORLD,false,delayCommitForRefinement);
  factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);

  // test  number of total elements to make sure refinement works
  TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),32);
}
#endif

#ifdef PANZER_HAVE_UMR
TEUCHOS_UNIT_TEST(tExodusReaderFactory, umr_refine_once_with_geometry)
{
  STK_ExodusReaderFactory factory;

  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
  pl->set("File Type","Exodus Refinement");
  pl->set("File Name","meshes/2block_cylinders_30deg.g");
  pl->set("Geometry File Name","meshes/2block_cylinders_30deg.stp");
  pl->set("Levels of Uniform Refinement",1);
  TEST_NOTHROW(factory.setParameterList(pl));

  Teuchos::RCP<STK_Interface> mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);
  factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);

  TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),12*8);

  std::vector<std::string> sidesets;
  mesh->getSidesetNames(sidesets);
  TEST_EQUALITY((int) sidesets.size(),5);

  std::vector<std::string> nodesets;
  mesh->getNodesetNames(nodesets);
  TEST_EQUALITY((int) nodesets.size(),3);

  std::map<std::string,int> sidesets_counts = {{"top",4},  {"interface",2}, {"surface_3",2}, {"surface_4",4}, {"inner",2}};

  for (auto const& x : sidesets_counts) {
    std::vector<size_t> globalCounts;
    stk::mesh::Part * ss_part = mesh->getSideset(x.first);
    stk::mesh::Selector selector(*ss_part);
    stk::mesh::comm_mesh_counts( *(mesh->getBulkData()), globalCounts, &selector);

    TEST_EQUALITY((int) globalCounts[mesh->getFaceRank()],x.second*4);
  }
  mesh->writeToExodus("meshes/2block_cylinders_30deg.r1.g");
}
#endif

}

#endif
