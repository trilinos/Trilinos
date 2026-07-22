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

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Tuple.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_Utilities.hpp"
#include "Panzer_STK_PeriodicBC_Matcher.hpp"
#include "Panzer_STK_PeriodicBC_Parser.hpp"
#include "Panzer_STK_PeriodicBC_MatchConditions.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_DOFManager.hpp"

#include "Epetra_MpiComm.h"

#include "Kokkos_DynRankView.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"

#include <string>

typedef Kokkos::DynRankView<double,PHX::Device> FieldContainer;

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer_stk {


  template <typename Intrepid2Type>
  RCP<const panzer::FieldPattern> buildFieldPattern()
  {
     // build a geometric pattern from a single basis
     RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2Type);
     RCP<const panzer::FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
     return pattern;
  }

  TEUCHOS_UNIT_TEST(periodic_mesh_nosubcells, serial)
  {
    using Teuchos::RCP;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    TEUCHOS_ASSERT(Comm.NumProc()==1);

    // panzer::pauseToAttach();

    panzer_stk::CubeHexMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",1);
       pl->set("Y Blocks",1);
       pl->set("Z Blocks",1);
       pl->set("X Elements",2);
       pl->set("Y Elements",2);
       pl->set("Z Elements",2);
       pl->set("X Procs",1);
       pl->set("Y Procs",1);
       pl->set("Z Procs",1);
       pl->set("X0",0.0);
       pl->set("Y0",0.0);
       pl->set("Z0",0.0);
       pl->set("Xf",6.0);
       pl->set("Yf",3.0);
       pl->set("Zf",4.0);
       mesh_factory.setParameterList(pl);
       mesh_factory.buildSubcells(false);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    panzer_stk::PlaneMatcher side_matcher(1,2);
    panzer_stk::PlaneMatcher front_matcher(0,1);
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("right","left",side_matcher));
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("front","back",front_matcher));

    mesh->writeToExodus("file.exo");

    // connection manager
    /////////////////////////////////////////////
    RCP<const panzer::FieldPattern> fp
         = buildFieldPattern<Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::exec_space,double,double> >();

    Teuchos::RCP<panzer::ConnManager> connMngr
          = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    Teuchos::RCP<panzer::DOFManager> dofManager
          = Teuchos::rcp(new panzer::DOFManager(connMngr,MPI_COMM_WORLD));
    dofManager->addField("VAR",fp);
    dofManager->buildGlobalUnknowns();

    std::vector<panzer::GlobalOrdinal> owned;
    dofManager->getOwnedIndices(owned);

    std::size_t unkCount = owned.size();
    std::size_t nodeCount = mesh->getEntityCounts(mesh->getNodeRank());
    out << "Unknown Count = " << unkCount
        << ", Node Count = "  << nodeCount << std::endl;

    TEST_ASSERT(unkCount < nodeCount);
  }
}
