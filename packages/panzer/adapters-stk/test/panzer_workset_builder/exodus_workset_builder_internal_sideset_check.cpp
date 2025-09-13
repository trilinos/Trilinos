// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_WorksetFactory.hpp"

namespace panzer {

  // This test checks that worksets for internal sidesets are created
  // correctly. If a sideset is internal, the side entities (edges in
  // 2D and face in 3D) can be owned by a different mpi process than
  // the cell associated with a particular boundary condition. The
  // search needs to consider ghosted edges. A fix was made in PR
  // #14006:
  //
  //   https://github.com/trilinos/Trilinos/pull/14006
  //
  // A specific exodus mesh was created to reproduce the issue. It
  // requires rib/rcb runtime parallel decompositon on exactly 2 mpi
  // processes to create the interface:
  //
  //   export IOSS_PROPERTIES="DECOMPOSITION_METHOD=rib"
  //

  TEUCHOS_UNIT_TEST(workset_builder, internal_sideset_issue_13717)
  {
    using Teuchos::RCP;
    using Teuchos::make_rcp;

    auto comm = Teuchos::make_rcp<const Teuchos::MpiComm<int>>(MPI_COMM_WORLD);
    TEST_EQUALITY(comm->getSize(),2);
    setenv("IOSS_PROPERTIES", "DECOMPOSITION_METHOD=rib", 1);

    RCP<panzer_stk::STK_Interface> mesh;
    {
      RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
      pl->set("File Name","energy-ss-tp-multiblock-ic-bc-issue.gen");
      panzer_stk::STK_ExodusReaderFactory factory;
      factory.setParameterList(pl);
      mesh = factory.buildMesh(MPI_COMM_WORLD);
    }

    const std::string element_block = "left";
    const std::string sideset = "center";
    const int workset_size = 20;

    auto conn_manager = make_rcp<panzer_stk::STKConnManager>(mesh);
    auto dof_manager = make_rcp<panzer::DOFManager>(conn_manager,MPI_COMM_WORLD);
    {
      auto intrepid_basis =
        panzer::createIntrepid2Basis<PHX::Device::execution_space,double,double>
        ("HGrad",1,*mesh->getCellTopology(element_block));

      auto field_pattern = make_rcp<panzer::Intrepid2FieldPattern>(intrepid_basis);

      dof_manager->addField(element_block, "T", field_pattern);
      dof_manager->buildGlobalUnknowns();
    }

    auto wkstFactory = make_rcp<panzer_stk::WorksetFactory>(mesh);
    auto wkstContainer = make_rcp<panzer::WorksetContainer>();

    {
      BasisDescriptor basis_desc(1,"HGrad");
      WorksetNeeds needs;
      needs.cellData = CellData(workset_size,mesh->getCellTopology(element_block));
      needs.addBasis(basis_desc);
      wkstContainer->setNeeds(element_block,needs);
    }

    wkstContainer->setFactory(wkstFactory);
    wkstContainer->setGlobalIndexer(dof_manager);
    wkstContainer->setWorksetSize(workset_size);

    Teuchos::RCP<std::map<unsigned,Workset> > rcp_worksets = wkstContainer->getSideWorksets(sidesetDescriptor(element_block,sideset));

    int local_side_element_count = 0;
    int global_side_element_count = 0;

    if (nonnull(rcp_worksets)) {
      for (auto side : *rcp_worksets) {
        local_side_element_count += side.second.num_cells;
      }
    }

    MPI_Allreduce(&local_side_element_count,&global_side_element_count,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    TEST_EQUALITY(global_side_element_count,7);

    unsetenv("IOSS_PROPERTIES");
  }

}
