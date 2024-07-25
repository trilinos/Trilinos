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
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include "Panzer_STK_SquareTriMeshFactory.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_CubeTetMeshFactory.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"

#include "Panzer_FaceToElement.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {
namespace disc_fe {
namespace core_tests {

  void setParamList(RCP<Teuchos::ParameterList> mesh_pl, int nx, int ny, int nz=0, int bx=1, int by=1, int bz=0)
  {
    mesh_pl->set("X Blocks", bx);
    mesh_pl->set("Y Blocks", by);
    mesh_pl->set("X Elements", nx);
    mesh_pl->set("Y Elements", ny);

    mesh_pl->set("X0", 0.0);
    mesh_pl->set("Y0", 0.0);
    mesh_pl->set("Xf", 1.0);
    mesh_pl->set("Yf", 1.0);

    if ( nz > 0 ) {
      mesh_pl->set("Z Blocks", bz);
      mesh_pl->set("Z Elements", nz);
      mesh_pl->set("Z0", 0.0);
      mesh_pl->set("Zf", 1.0);
    }
  }

#ifndef PANZER_HIDE_DEPRECATED_CODE
  TEUCHOS_UNIT_TEST(FaceToElement, test_default_communicator)
  {
    Teuchos::RCP<const Teuchos::MpiComm<int>> comm_world(new Teuchos::MpiComm< int>(MPI_COMM_WORLD));

    RCP<Teuchos::ParameterList> mesh_pl(new Teuchos::ParameterList("Mesh"));
    int nx=4, ny=8, nz=-1, bx=1, by=2, bz=-1;
    setParamList(mesh_pl, nx, ny, nz, bx, by, bz);

    RCP<panzer_stk::STK_MeshFactory> mesh_factory = rcp(new panzer_stk::SquareTriMeshFactory);
    mesh_factory->setParameterList(mesh_pl);

    RCP<panzer_stk::STK_Interface> mesh = mesh_factory->buildUncommitedMesh(*comm_world->getRawMpiComm());
    mesh_factory->completeMeshConstruction(*mesh,*comm_world->getRawMpiComm());
    const Teuchos::RCP<panzer::ConnManager> conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    auto faceToElement = Teuchos::rcp(new panzer::FaceToElement<panzer::LocalOrdinal,panzer::GlobalOrdinal>());
    faceToElement->initialize(*conn_manager);

    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif

  TEUCHOS_UNIT_TEST(FaceToElement, test_comm_world)
  {
    Teuchos::RCP<const Teuchos::MpiComm<int>> comm_world(new Teuchos::MpiComm< int>(MPI_COMM_WORLD));

    RCP<Teuchos::ParameterList> mesh_pl(new Teuchos::ParameterList("Mesh"));
    int nx=4, ny=8, nz=-1, bx=1, by=2, bz=-1;
    setParamList(mesh_pl, nx, ny, nz, bx, by, bz);

    RCP<panzer_stk::STK_MeshFactory> mesh_factory = rcp(new panzer_stk::SquareTriMeshFactory);
    mesh_factory->setParameterList(mesh_pl);

    RCP<panzer_stk::STK_Interface> mesh = mesh_factory->buildUncommitedMesh(*comm_world->getRawMpiComm());
    mesh_factory->completeMeshConstruction(*mesh,*comm_world->getRawMpiComm());
    const Teuchos::RCP<panzer::ConnManager> conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    auto faceToElement = Teuchos::rcp(new panzer::FaceToElement<panzer::LocalOrdinal,panzer::GlobalOrdinal>());
    faceToElement->initialize(*conn_manager, comm_world);

    MPI_Barrier(MPI_COMM_WORLD);
  }

  TEUCHOS_UNIT_TEST(FaceToElement, test_sub_communicator)
  {
    std::vector<int> this_ranks = {0,1};
    Teuchos::ArrayView<int> array_view(this_ranks);
    Teuchos::RCP<const Teuchos::MpiComm<int>> comm_world(new Teuchos::MpiComm< int>(MPI_COMM_WORLD));
    auto sub_comm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(comm_world->createSubcommunicator(array_view),true);

    if (comm_world->getRank() < 2) {
      RCP<Teuchos::ParameterList> mesh_pl(new Teuchos::ParameterList("Mesh"));
      int nx=4, ny=8, nz=-1, bx=1, by=2, bz=-1;
      setParamList(mesh_pl, nx, ny, nz, bx, by, bz);

      RCP<panzer_stk::STK_MeshFactory> mesh_factory = rcp(new panzer_stk::SquareTriMeshFactory);
      mesh_factory->setParameterList(mesh_pl);

      RCP<panzer_stk::STK_Interface> mesh = mesh_factory->buildUncommitedMesh(*sub_comm->getRawMpiComm());
      mesh_factory->completeMeshConstruction(*mesh,*sub_comm->getRawMpiComm());
      const Teuchos::RCP<panzer::ConnManager> conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

      auto faceToElement = Teuchos::rcp(new panzer::FaceToElement<panzer::LocalOrdinal,panzer::GlobalOrdinal>());
      faceToElement->initialize(*conn_manager, sub_comm);
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }
}
}
}
