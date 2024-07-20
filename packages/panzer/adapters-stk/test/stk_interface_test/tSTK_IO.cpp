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
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Kokkos_ViewFactory.hpp"

#include "Kokkos_DynRankView.hpp"

#include "Ioss_DatabaseIO.h"
#include "Ioss_IOFactory.h"
#include "Ioss_Region.h"


#ifdef PANZER_HAVE_IOSS

using Teuchos::RCP;
using Teuchos::rcp;

typedef Kokkos::DynRankView<double,PHX::Device> FieldContainer;

namespace panzer_stk {

RCP<STK_Interface> buildMesh(int xElements,int yElements);
RCP<STK_Interface> buildMesh_cells(int xElements,int yElements);
void buildLocalIds(const STK_Interface & mesh,
                   std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > > & localIds);

void assignBlock(FieldContainer & block,FieldContainer & vertices, double (* func)(double,double));
void assignBlock(FieldContainer & block,FieldContainer & vertices, double val);

double xval(double x, double /* y */) { return x; }
double yval(double /* x */, double y) { return y; }
double block2(double x,double y) { return (x-0.5)*(x-0.5)+y; }

// triangle tests
TEUCHOS_UNIT_TEST(tSTK_IO, fields)
{
   RCP<STK_Interface> mesh = buildMesh(8,8);

   std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > > localIds;
   buildLocalIds(*mesh,localIds);

   FieldContainer vert0, vert1;
   out << "get vertices" << std::endl;
   mesh->getElementVertices(*localIds["eblock-0_0"],vert0);
   mesh->getElementVertices(*localIds["eblock-1_0"],vert1);

   FieldContainer ublock0, tblock0, tblock1;
   ublock0 = Kokkos::createDynRankView(ublock0,"ublock0",localIds["eblock-0_0"]->size(),4);
   tblock0 = Kokkos::createDynRankView(tblock0,"tblock0",localIds["eblock-0_0"]->size(),4);
   tblock1 = Kokkos::createDynRankView(tblock1,"tblock1",localIds["eblock-1_0"]->size(),4);
   out << "assigning" << std::endl;

   assignBlock(ublock0,vert0,xval);
   assignBlock(tblock0,vert0,yval);
   assignBlock(tblock1,vert1,block2);

   mesh->setSolutionFieldData("u","eblock-0_0",*localIds["eblock-0_0"],ublock0);
   mesh->setSolutionFieldData("T","eblock-0_0",*localIds["eblock-0_0"],tblock0);
   mesh->setSolutionFieldData("T","eblock-1_0",*localIds["eblock-1_0"],tblock1);

   out << "write to exodus" << std::endl;

   mesh->writeToExodus("output.exo");
}

TEUCHOS_UNIT_TEST(tSTK_IO, cell_fields)
{
   RCP<STK_Interface> mesh = buildMesh_cells(8,8);

   std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > > localIds;
   buildLocalIds(*mesh,localIds);

   out << "write to exodus" << std::endl;

   FieldContainer vert0, vert1;
   out << "get vertices" << std::endl;
   mesh->getElementVertices(*localIds["eblock-0_0"],vert0);
   mesh->getElementVertices(*localIds["eblock-1_0"],vert1);

   FieldContainer ublock0, tblock0, tblock1;
   ublock0 = Kokkos::createDynRankView(ublock0,"ublock0",localIds["eblock-0_0"]->size(),4);
   tblock0 = Kokkos::createDynRankView(tblock0,"tblock0",localIds["eblock-0_0"]->size(),4);
   tblock1 = Kokkos::createDynRankView(tblock1,"tblock1",localIds["eblock-1_0"]->size(),4);
   out << "assigning" << std::endl;

   assignBlock(ublock0,vert0,xval);
   assignBlock(tblock0,vert0,yval);
   assignBlock(tblock1,vert1,block2);

   mesh->setCellFieldData("u","eblock-0_0",*localIds["eblock-0_0"],ublock0);
   mesh->setCellFieldData("T","eblock-0_0",*localIds["eblock-0_0"],tblock0);
   mesh->setCellFieldData("T","eblock-1_0",*localIds["eblock-1_0"],tblock1);
   mesh->setSolutionFieldData("T","eblock-1_0",*localIds["eblock-1_0"],tblock1);

   mesh->writeToExodus("output-cells.exo");
}

TEUCHOS_UNIT_TEST(tSTK_IO, exodus_factory_transient_fields)
{
   STK_ExodusReaderFactory factory("meshes/basic.gen");
   RCP<STK_Interface> mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);
   mesh->addSolutionField("u","block_1");
   mesh->addSolutionField("T","block_1");
   mesh->addSolutionField("T","block_2");
   factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);

   std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > > localIds;
   buildLocalIds(*mesh,localIds);

   FieldContainer vert0, vert1;
   out << "get vertices" << std::endl;
   mesh->getElementVertices(*localIds["block_1"],vert0);
   mesh->getElementVertices(*localIds["block_2"],vert1);

   FieldContainer ublock0, tblock0, tblock1;
   ublock0 = Kokkos::createDynRankView(ublock0,"ublock0",localIds["block_1"]->size(),4);
   tblock0 = Kokkos::createDynRankView(tblock0,"tblock0",localIds["block_1"]->size(),4);
   tblock1 = Kokkos::createDynRankView(tblock1,"tblock1",localIds["block_2"]->size(),4);

   mesh->setupExodusFile("transient_exo.exo");

   out << "assigning 4.5" << std::endl;
   {
      assignBlock(ublock0,vert0,6.0);
      assignBlock(tblock0,vert0,7.0);
      assignBlock(tblock1,vert1,8.0);

      mesh->setSolutionFieldData("u","block_1",*localIds["block_1"],ublock0);
      mesh->setSolutionFieldData("T","block_1",*localIds["block_1"],tblock0);
      mesh->setSolutionFieldData("T","block_2",*localIds["block_2"],tblock1);
   }

   out << "write to exodus: 4.5" << std::endl;
   mesh->writeToExodus(4.5);
}

TEUCHOS_UNIT_TEST(tSTK_IO, transient_fields)
{
   RCP<STK_Interface> mesh = buildMesh(20,20);

   std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > > localIds;
   buildLocalIds(*mesh,localIds);

   FieldContainer vert0, vert1;
   out << "get vertices" << std::endl;
   mesh->getElementVertices(*localIds["eblock-0_0"],vert0);
   mesh->getElementVertices(*localIds["eblock-1_0"],vert1);

   FieldContainer ublock0, tblock0, tblock1;
   ublock0 = Kokkos::createDynRankView(ublock0,"ublock0",localIds["eblock-0_0"]->size(),4);
   tblock0 = Kokkos::createDynRankView(tblock0,"tblock0",localIds["eblock-0_0"]->size(),4);
   tblock1 = Kokkos::createDynRankView(tblock1,"tblock1",localIds["eblock-1_0"]->size(),4);

   mesh->setupExodusFile("transient.exo");

   out << "assigning 3.0" << std::endl;
   {
      assignBlock(ublock0,vert0,1.0);
      assignBlock(tblock0,vert0,2.0);
      assignBlock(tblock1,vert1,3.0);

      mesh->setSolutionFieldData("u","eblock-0_0",*localIds["eblock-0_0"],ublock0);
      mesh->setSolutionFieldData("T","eblock-0_0",*localIds["eblock-0_0"],tblock0);
      mesh->setSolutionFieldData("T","eblock-1_0",*localIds["eblock-1_0"],tblock1);
   }

   out << "write to exodus: 3.0" << std::endl;
   mesh->writeToExodus(3.0);
   TEST_EQUALITY(mesh->getCurrentStateTime(),3.0);

   out << "assigning 4.5" << std::endl;
   {
      assignBlock(ublock0,vert0,6.0);
      assignBlock(tblock0,vert0,7.0);
      assignBlock(tblock1,vert1,8.0);

      mesh->setSolutionFieldData("u","eblock-0_0",*localIds["eblock-0_0"],ublock0);
      mesh->setSolutionFieldData("T","eblock-0_0",*localIds["eblock-0_0"],tblock0);
      mesh->setSolutionFieldData("T","eblock-1_0",*localIds["eblock-1_0"],tblock1);
   }

   out << "write to exodus: 4.5" << std::endl;
   mesh->writeToExodus(4.5);
   TEST_EQUALITY(mesh->getCurrentStateTime(),4.5);
   // Data can be buffered in writeToExodus() call. Flush to file by closing.
   mesh = Teuchos::null;


   STK_ExodusReaderFactory factory("transient.exo",2);
   RCP<STK_Interface> mesh_read = factory.buildMesh(MPI_COMM_WORLD);
   TEST_EQUALITY(mesh_read->getInitialStateTime(),4.5);
   TEST_EQUALITY(mesh_read->getCurrentStateTime(),0.0); // writeToExodus has not yet been called
}

TEUCHOS_UNIT_TEST(tSTK_IO, addInformationRecords)
{
   using Teuchos::RCP;

   std::vector<std::string> info_records_1 = {
     "DG::eblock-0_0_0::basis::Basis_HGRAD_HEX_C1_FEM",
     "DG::eblock-0_0_0::field::Ex"
   };
   std::vector<std::string> info_records_2 = {
     "DG::eblock-0_0_0::basis::Basis_HGRAD_HEX_C1_FEM",
     "DG::eblock-0_0_0::field::Ey"
   };
   std::vector<std::string> info_records_3 = {
     "DG::eblock-0_0_0::basis::Basis_HGRAD_HEX_C1_FEM",
     "DG::eblock-0_0_0::field::Ez"
   };
   std::vector<std::string> info_records_4 = {
     "DG::eblock-1_1_1::basis::Basis_HGRAD_HEX_C1_FEM",
     "DG::eblock-1_1_1::field::Ex"
   };
   std::vector<std::string> info_records_5 = {
     "DG::eblock-1_1_1::basis::Basis_HGRAD_HEX_C1_FEM",
     "DG::eblock-1_1_1::field::Ey"
   };
   std::vector<std::string> info_records_6 = {
     "DG::eblock-1_1_1::basis::Basis_HGRAD_HEX_C1_FEM",
     "DG::eblock-1_1_1::field::Ez"
   };

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",1);
   pl->set("Y Blocks",1);
   pl->set("Z Blocks",1);
   pl->set("X Elements",2);
   pl->set("Y Elements",4);
   pl->set("Z Elements",5);

   CubeHexMeshFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);
   mesh->addInformationRecords(info_records_1);
   mesh->addInformationRecords(info_records_2);
   mesh->addInformationRecords(info_records_3);
   mesh->addInformationRecords(info_records_4);
   mesh->addInformationRecords(info_records_5);
   mesh->addInformationRecords(info_records_6);
   factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);

   TEST_ASSERT(mesh!=Teuchos::null);

   mesh->writeToExodus("info_records.exo");

   {
   Ioss::DatabaseIO *db_io = Ioss::IOFactory::create("exodus",
                                                     "info_records.exo",
                                                     Ioss::READ_MODEL);
   TEST_ASSERT(db_io);

   Ioss::Region region(db_io);
   TEST_ASSERT(db_io->ok() == true);

   auto info_records_from_ioss = db_io->get_information_records();

   // put all the info records in one vector to make them easier to iterate over
   std::vector<std::string> all_expected_records;
   all_expected_records.insert(all_expected_records.end(), info_records_1.begin(), info_records_1.end());
   all_expected_records.insert(all_expected_records.end(), info_records_2.begin(), info_records_2.end());
   all_expected_records.insert(all_expected_records.end(), info_records_3.begin(), info_records_3.end());
   all_expected_records.insert(all_expected_records.end(), info_records_4.begin(), info_records_4.end());
   all_expected_records.insert(all_expected_records.end(), info_records_5.begin(), info_records_5.end());
   all_expected_records.insert(all_expected_records.end(), info_records_6.begin(), info_records_6.end());
   // make sure all the input records appear in the results returned from IOSS
   for ( auto r : all_expected_records ) {
      auto iter = std::find(info_records_from_ioss.begin(), info_records_from_ioss.end(), r);
      TEST_ASSERT(iter != info_records_from_ioss.end());
   }
   }
}

RCP<STK_Interface> buildMesh(int xElements,int yElements)
{
    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",xElements);  // in each block
    pl->set("Y Elements",yElements);  // in each block

    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<STK_Interface> mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);

    mesh->addSolutionField("u","eblock-0_0");
    mesh->addSolutionField("T","eblock-0_0");
    mesh->addSolutionField("T","eblock-1_0");

    factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);

    return mesh;
}

RCP<STK_Interface> buildMesh_cells(int xElements,int yElements)
{
    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",xElements);  // in each block
    pl->set("Y Elements",yElements);  // in each block

    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<STK_Interface> mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);

    mesh->addCellField("u","eblock-0_0");
    mesh->addCellField("T","eblock-0_0");
    mesh->addCellField("T","eblock-1_0");

    mesh->addSolutionField("T","eblock-1_0");

    factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);

    return mesh;
}

void buildLocalIds(const STK_Interface & mesh,
                   std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > > & localIds)
{
   // defines ordering of blocks
   std::vector<std::string> blockIds;
   mesh.getElementBlockNames(blockIds);

   std::vector<std::string>::const_iterator idItr;
   for(idItr=blockIds.begin();idItr!=blockIds.end();++idItr) {
      std::string blockId = *idItr;

      localIds[blockId] = Teuchos::rcp(new std::vector<std::size_t>);
      std::vector<std::size_t> & localBlockIds = *localIds[blockId];

      // grab elements on this block
      std::vector<stk::mesh::Entity> blockElmts;
      mesh.getMyElements(blockId,blockElmts);

      std::vector<stk::mesh::Entity>::const_iterator itr;
      for(itr=blockElmts.begin();itr!=blockElmts.end();++itr)
         localBlockIds.push_back(mesh.elementLocalId(*itr));

      std::sort(localBlockIds.begin(),localBlockIds.end());
   }
}

void assignBlock(FieldContainer & block,FieldContainer & vertices, double val)
{
   TEUCHOS_ASSERT(block.extent(0)==vertices.extent(0));
   TEUCHOS_ASSERT(block.extent(1)==vertices.extent(1));

   std::size_t cellCnt = block.extent(0);
   std::size_t nodeCnt = block.extent(1);

   Kokkos::parallel_for(
       cellCnt, KOKKOS_LAMBDA(size_t cell) {
          for (std::size_t node = 0; node < nodeCnt; node++) {
             block(cell, node) = val;
          }
       }
   );
}

void assignBlock(FieldContainer & block,FieldContainer & vertices, double (* func)(double,double))
{
   TEUCHOS_ASSERT(block.extent(0)==vertices.extent(0));
   TEUCHOS_ASSERT(block.extent(1)==vertices.extent(1));

   std::size_t cellCnt = block.extent(0);
   std::size_t nodeCnt = block.extent(1);

   auto vertices_h = Kokkos::create_mirror_view(vertices);
   Kokkos::deep_copy(vertices_h, vertices);

   auto block_h = Kokkos::create_mirror_view(block);
   Kokkos::deep_copy(block_h, block);

   for(std::size_t cell = 0; cell < cellCnt; cell++) {
      for(std::size_t node = 0; node < nodeCnt; node++) {
         double x = vertices_h(cell, node, 0);
         double y = vertices_h(cell, node, 1);
         block_h(cell, node) = func(x, y);
      }
   }

   Kokkos::deep_copy(block, block_h);
   Kokkos::deep_copy(vertices, vertices_h);
}

}

#endif
