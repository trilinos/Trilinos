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

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"

#include "Panzer_WorksetContainer.hpp"
#include "Panzer_WorksetDescriptor.hpp"
#include "Panzer_WorksetNeeds.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"
#include "Panzer_IntegrationDescriptor.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_SubcellConnectivity.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_WorksetFactory.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(check_rotation_matrices, basic)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef Intrepid2::Basis<PHX::Device::execution_space,double,double> IntrepidBasis;

    std::string element_block = "eblock-0_0_0";
    std::string sideset = "left";

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Elements",1);
    pl->set("Y Elements",1);
    pl->set("Z Elements",1);
    pl->set("X0",-0.001);
    pl->set("Xf", 0.001);
    pl->set("Y0",-0.0005);
    pl->set("Yf", 0.0005);
    pl->set("Z0",-0.0003333333333333333);
    pl->set("Zf", 0.0003333333333333333);

    panzer_stk::CubeHexMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    RCP<const Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));

    // build DOF Manager (with a single HDiv basis)
    /////////////////////////////////////////////////////////////

    // build the connection manager
    const RCP<panzer::ConnManager>
      conn_manager = rcp(new panzer_stk::STKConnManager(mesh));

    RCP<panzer::DOFManager> dof_manager
        = rcp(new panzer::DOFManager(conn_manager,MPI_COMM_WORLD));

    // build an intrepid basis and a related field pattern for seeding the DOFManager
    {
       RCP<IntrepidBasis> intrepid_basis
           = panzer::createIntrepid2Basis<PHX::Device::execution_space,double,double>("HGrad",1,
                                                                                      *mesh->getCellTopology(element_block));
      RCP<panzer::Intrepid2FieldPattern> field_pattern = rcp(new panzer::Intrepid2FieldPattern(intrepid_basis));

      dof_manager->addField(element_block, "T", field_pattern);
    }
    dof_manager->buildGlobalUnknowns();

    // build WorksetContainer
    //////////////////////////////////////////////////////////////

    panzer::IntegrationDescriptor sid(2*2, panzer::IntegrationDescriptor::SURFACE);

    RCP<panzer_stk::WorksetFactory> wkstFactory
       = rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
    RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
       = rcp(new panzer::WorksetContainer(wkstFactory));

    wkstContainer->setGlobalIndexer(dof_manager);

    panzer::WorksetDescriptor workset_descriptor(element_block, panzer::WorksetSizeType::ALL_ELEMENTS, true,false);

    auto worksets = wkstContainer->getWorksets(workset_descriptor);

    TEST_ASSERT(worksets->size()==1);

    auto & workset = (*worksets)[0];

    auto rot_matrices = workset.getIntegrationValues(sid).getSurfaceRotationMatrices();
    auto normals = workset.getIntegrationValues(sid).getSurfaceNormals();

    const int num_owned_cells   = workset.numOwnedCells();
    const int num_ghost_cells   = workset.numGhostCells();
    const int num_real_cells    = num_owned_cells + num_ghost_cells;
    const int num_virtual_cells = workset.numVirtualCells();
    const int num_cells         = num_real_cells + num_virtual_cells;

    const int faces_per_cell    = 6; // hexahedron
    auto & face_connectivity    = workset.getFaceConnectivity();

    // sanity check on cell counts: should have 6 virtual, 1 owned
    TEST_EQUALITY(num_owned_cells,   1);
    TEST_EQUALITY(num_ghost_cells,   0);
    TEST_EQUALITY(num_virtual_cells, 6);

    TEST_ASSERT(rot_matrices.rank()==4);
    TEST_ASSERT(rot_matrices.extent_int(0)==num_cells);
    TEST_ASSERT(rot_matrices.extent_int(2)==3);
    TEST_ASSERT(rot_matrices.extent_int(3)==3);
    out << "SIZES = "
        << rot_matrices.extent(0) << " "
        << rot_matrices.extent(1) << " "
        << rot_matrices.extent(2) << " "
        << rot_matrices.extent(3) << std::endl;

    out << std::endl;

    const int num_points      = rot_matrices.extent_int(1);
    const int points_per_face = num_points / faces_per_cell;
    for(int c=0;c<rot_matrices.extent_int(0);c++) {
      // rule for virtual cells is that they have rotation matrices defined only on local face ordinal 0;
      // the rotation matrices for the other faces of the virtual cell should be all 0s.
      const bool is_virtual = (c >= num_real_cells);
      int virtual_local_face_id = -1; // the virtual cell face that adjoins the real cell

      if (is_virtual)
      {
        // determine which face adjoins the real cell:
        int face_ordinal = -1;
        for (int local_face_id=0; local_face_id<faces_per_cell; local_face_id++)
        {
          face_ordinal = face_connectivity.subcellForCellHost(c, local_face_id);
          if (face_ordinal >= 0)
          {
            virtual_local_face_id = local_face_id;
            break;
          }
        }
      }

      auto rot_matrices_view = PHX::as_view(rot_matrices);
      auto rot_matrices_h = Kokkos::create_mirror_view(rot_matrices_view);
      Kokkos::deep_copy(rot_matrices_h, rot_matrices_view);

      for(int p=0;p<num_points;p++) {
        const int local_face_ordinal = p / points_per_face;
        bool expect_rotation_matrix = true; // if false, we expect all 0s
        if (is_virtual)
        {
          expect_rotation_matrix = (local_face_ordinal == virtual_local_face_id);
        }
        out << "Cell,Point = " << c << "," << p << std::endl;
        out << "   " << rot_matrices_h(c,p,0,0) << " " << rot_matrices_h(c,p,0,1) << " " << rot_matrices_h(c,p,0,2) << std::endl;
        out << "   " << rot_matrices_h(c,p,1,0) << " " << rot_matrices_h(c,p,1,1) << " " << rot_matrices_h(c,p,1,2) << std::endl;
        out << "   " << rot_matrices_h(c,p,2,0) << " " << rot_matrices_h(c,p,2,1) << " " << rot_matrices_h(c,p,2,2) << std::endl;
        out << std::endl;

        if (expect_rotation_matrix)
        {
          // mutually orthogonal
          TEST_ASSERT(std::fabs(rot_matrices_h(c,p,0,0) * rot_matrices_h(c,p,1,0) +
                                rot_matrices_h(c,p,0,1) * rot_matrices_h(c,p,1,1) +
                                rot_matrices_h(c,p,0,2) * rot_matrices_h(c,p,1,2))<=1e-14);
          TEST_ASSERT(std::fabs(rot_matrices_h(c,p,0,0) * rot_matrices_h(c,p,2,0) +
                                rot_matrices_h(c,p,0,1) * rot_matrices_h(c,p,2,1) +
                                rot_matrices_h(c,p,0,2) * rot_matrices_h(c,p,2,2))<=1e-14);
          TEST_ASSERT(std::fabs(rot_matrices_h(c,p,1,0) * rot_matrices_h(c,p,2,0) +
                                rot_matrices_h(c,p,1,1) * rot_matrices_h(c,p,2,1) +
                                rot_matrices_h(c,p,1,2) * rot_matrices_h(c,p,2,2))<=1e-14);

          // normalized
          TEST_FLOATING_EQUALITY(std::sqrt(rot_matrices_h(c,p,0,0) * rot_matrices_h(c,p,0,0) +
                                           rot_matrices_h(c,p,0,1) * rot_matrices_h(c,p,0,1) +
                                           rot_matrices_h(c,p,0,2) * rot_matrices_h(c,p,0,2)),1.0,1e-14);
          TEST_FLOATING_EQUALITY(std::sqrt(rot_matrices_h(c,p,1,0) * rot_matrices_h(c,p,1,0) +
                                           rot_matrices_h(c,p,1,1) * rot_matrices_h(c,p,1,1) +
                                           rot_matrices_h(c,p,1,2) * rot_matrices_h(c,p,1,2)),1.0,1e-14);
          TEST_FLOATING_EQUALITY(std::sqrt(rot_matrices_h(c,p,2,0) * rot_matrices_h(c,p,2,0) +
                                           rot_matrices_h(c,p,2,1) * rot_matrices_h(c,p,2,1) +
                                           rot_matrices_h(c,p,2,2) * rot_matrices_h(c,p,2,2)),1.0,1e-14);
        }
        else
        {
          // virtual cell on a face other than local face 0: expect zeros
          for (int i=0; i<3; i++)
          {
            for (int j=0; j<3; j++)
            {
              // since these should be filled in as exact machine 0s, we do non-floating comparison
              TEST_EQUALITY(rot_matrices_h(c,p,i,j), 0.0);
            }
          }
        }
      }
    }
  }
}
