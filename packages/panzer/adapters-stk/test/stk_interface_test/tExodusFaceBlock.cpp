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
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"

#include "Shards_BasicTopologies.hpp"

#include "Kokkos_DynRankView.hpp"

namespace panzer_stk {

TEUCHOS_UNIT_TEST(tExodusFaceBlock, face_count)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
   int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
   out << "Running numprocs = " << numprocs << " rank = " << rank << std::endl;

   std::size_t xelems=2;
   std::size_t yelems=4;
   std::size_t zelems=5;
   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",1);
   pl->set("Y Blocks",1);
   pl->set("Z Blocks",1);
   pl->set("X Elements",(int)xelems);
   pl->set("Y Elements",(int)yelems);
   pl->set("Z Elements",(int)zelems);
   pl->set("Create Face Blocks",true);

   CubeHexMeshFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);

   if(mesh->isWritable())
      mesh->writeToExodus("FaceBlock1.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),3);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),xelems*yelems*zelems);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),xelems*yelems*(zelems+1)+xelems*zelems*(yelems+1)+yelems*zelems*(xelems+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getEdgeRank()),xelems*(yelems+1)*(zelems+1)+yelems*(xelems+1)*(zelems+1)+zelems*(xelems+1)*(yelems+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(yelems+1)*(xelems+1)*(zelems+1));

   std::vector<stk::mesh::Entity> all_faces;
   mesh->getAllFaces("quad_4_"+panzer_stk::STK_Interface::faceBlockString, all_faces);
   TEST_EQUALITY(all_faces.size(),xelems*yelems*(zelems+1)
                                 +xelems*zelems*(yelems+1)
                                 +yelems*zelems*(xelems+1));

   std::vector<stk::mesh::Entity> my_faces;
   if (numprocs==1) {
     // all faces belong to rank0, so getMyFaces() is equivalent to getAllFaces()
     mesh->getMyFaces("quad_4_"+panzer_stk::STK_Interface::faceBlockString, my_faces);
     TEST_EQUALITY(my_faces.size(),all_faces.size());
     TEST_EQUALITY(my_faces.size(),xelems*yelems*(zelems+1)
                                  +xelems*zelems*(yelems+1)
                                  +yelems*zelems*(xelems+1));
   }
   else if(numprocs==2 && rank==0) {
     // rank0 owns all faces in it's half of the mesh including the
     // faces on the plane shared with rank1
     std::size_t my_xelems=xelems/2;
     std::size_t my_yelems=yelems;
     std::size_t my_zelems=zelems;
     mesh->getMyFaces("quad_4_"+panzer_stk::STK_Interface::faceBlockString, my_faces);
     TEST_EQUALITY(my_faces.size(),my_xelems*my_yelems*(my_zelems+1)
                                  +my_xelems*my_zelems*(my_yelems+1)
                                  +my_yelems*my_zelems*(my_xelems+1));
   }
   else if(numprocs==2 && rank==1) {
     // rank1 doesn't own the faces on the plane shared with rank0
     std::size_t my_xelems=xelems/2;
     std::size_t my_yelems=yelems;
     std::size_t my_zelems=zelems;
     mesh->getMyFaces("quad_4_"+panzer_stk::STK_Interface::faceBlockString, my_faces);
     TEST_EQUALITY(my_faces.size(),my_xelems*my_yelems*(my_zelems+1)
                                  +my_xelems*my_zelems*(my_yelems+1)
                                  +my_yelems*my_zelems*(my_xelems+1)
                                  -my_yelems*my_zelems);             // remove faces owned by rank0
   }
   else {
     // fail!
     TEST_ASSERT(false && "This test must run with either 1 or 2 ranks.");
   }
}

TEUCHOS_UNIT_TEST(tExodusFaceBlock, is_face_local)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
   int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
   out << "Running numprocs = " << numprocs << " rank = " << rank << std::endl;

   std::size_t xelems=2;
   std::size_t yelems=4;
   std::size_t zelems=5;
   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",1);
   pl->set("Y Blocks",1);
   pl->set("Z Blocks",1);
   pl->set("X Elements",(int)xelems);
   pl->set("Y Elements",(int)yelems);
   pl->set("Z Elements",(int)zelems);
   pl->set("Create Face Blocks",true);

   CubeHexMeshFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);

   if(mesh->isWritable())
      mesh->writeToExodus("FaceBlock2.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),3);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),xelems*yelems*zelems);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),xelems*yelems*(zelems+1)+xelems*zelems*(yelems+1)+yelems*zelems*(xelems+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getEdgeRank()),xelems*(yelems+1)*(zelems+1)+yelems*(xelems+1)*(zelems+1)+zelems*(xelems+1)*(yelems+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(yelems+1)*(xelems+1)*(zelems+1));

   std::vector<stk::mesh::Entity> my_faces;
   mesh->getMyFaces("quad_4_"+panzer_stk::STK_Interface::faceBlockString, my_faces);
   for ( auto face : my_faces ) {
      TEST_ASSERT(mesh->isFaceLocal(face));
   }

   std::vector<stk::mesh::Entity> all_faces;
   mesh->getAllFaces("quad_4_"+panzer_stk::STK_Interface::faceBlockString, all_faces);
   for ( auto face : all_faces ) {
      if (mesh->getBulkData()->parallel_owner_rank(face)==rank) {
        TEST_ASSERT(mesh->isFaceLocal(face));
      } else {
        TEST_ASSERT(not mesh->isFaceLocal(face));
      }
   }
}

TEUCHOS_UNIT_TEST(tExodusFaceBlock, add_face_field)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",1);
   pl->set("Y Blocks",1);
   pl->set("Z Blocks",1);
   pl->set("X Elements",2);
   pl->set("Y Elements",4);
   pl->set("Z Elements",5);
   pl->set("Create Face Blocks",true);

   CubeHexMeshFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);

   mesh->addFaceField("face_field_1", "eblock-0_0_0");
   mesh->addFaceField("face_field_2", "eblock-0_0_0");

   stk::mesh::Field<double> * face_field_1 = mesh->getFaceField("face_field_1", "eblock-0_0_0");
   stk::mesh::Field<double> * face_field_2 = mesh->getFaceField("face_field_2", "eblock-0_0_0");

   std::vector<stk::mesh::Entity> faces;
   mesh->getAllFaces("quad_4_"+panzer_stk::STK_Interface::faceBlockString, faces);
   for(auto face : faces) {
     double* data = stk::mesh::field_data(*face_field_1, face);
     // set the face's field value to face's entity ID
     *data = mesh->getBulkData()->identifier(face);
   }
   for(auto face : faces) {
     double* data = stk::mesh::field_data(*face_field_2, face);
     // set the face's field value to face's entity ID * 2
     *data = 2*mesh->getBulkData()->identifier(face);
   }

   if(mesh->isWritable())
      mesh->writeToExodus("FaceBlock3.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),3);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),4*2*5);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),2*4*(5+1)+2*5*(4+1)+4*5*(2+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getEdgeRank()),2*(4+1)*(5+1)+4*(2+1)*(5+1)+5*(2+1)*(4+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(4+1)*(2+1)*(5+1));

   mesh->getAllFaces("quad_4_"+panzer_stk::STK_Interface::faceBlockString, faces);
   TEST_EQUALITY(faces.size(),2*4*(5+1)+4*5*(2+1)+5*2*(4+1));
}

TEUCHOS_UNIT_TEST(tExodusFaceBlock, set_face_field_data)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",1);
   pl->set("Y Blocks",1);
   pl->set("Z Blocks",1);
   pl->set("X Elements",2);
   pl->set("Y Elements",4);
   pl->set("Z Elements",5);
   pl->set("Create Face Blocks",true);

   CubeHexMeshFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);

   mesh->addFaceField("face_field_3", "eblock-0_0_0");
   mesh->addFaceField("face_field_4", "eblock-0_0_0");

   std::vector<stk::mesh::Entity> faces;
   mesh->getMyFaces("quad_4_"+panzer_stk::STK_Interface::faceBlockString, faces);

   Kokkos::DynRankView<double,PHX::Device> faceValues;
   faceValues = Kokkos::createDynRankView(faceValues,"faceValues",faces.size());
   auto faceValues_h = Kokkos::create_mirror_view(faceValues);

   std::vector<std::size_t> faceIds;
   for(auto face : faces) {
     faceIds.push_back(mesh->faceLocalId(face));
   }
   sort(faceIds.begin(),faceIds.end());

   for(std::size_t i=0;i<faceIds.size();i++) {
     faceValues_h(i) = 3*mesh->faceGlobalId(faceIds[i]);
   }
   mesh->setFaceFieldData("face_field_3",
                          "eblock-0_0_0",
                          faceIds,
                          faceValues_h);

   for(std::size_t i=0;i<faceIds.size();i++) {
     faceValues_h(i) = 4*mesh->faceGlobalId(faceIds[i]);
   }
   mesh->setFaceFieldData("face_field_4",
                          "eblock-0_0_0",
                          faceIds,
                          faceValues_h);

   if(mesh->isWritable())
      mesh->writeToExodus("FaceBlock4.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),3);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),4*2*5);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),2*4*(5+1)+2*5*(4+1)+4*5*(2+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getEdgeRank()),2*(4+1)*(5+1)+4*(2+1)*(5+1)+5*(2+1)*(4+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(4+1)*(2+1)*(5+1));

   mesh->getAllFaces("quad_4_"+panzer_stk::STK_Interface::faceBlockString, faces);
   TEST_EQUALITY(faces.size(),2*4*(5+1)+4*5*(2+1)+5*2*(4+1));
}

}
