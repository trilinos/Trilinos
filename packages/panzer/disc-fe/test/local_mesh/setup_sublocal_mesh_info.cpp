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

#include "PanzerCore_config.hpp"

#include "Panzer_LocalMeshInfo.hpp"
#include "Panzer_LocalPartitioningUtilities.hpp"

#include "Panzer_UnitTest_LocalMeshUtilities.hpp"

#include <vector>
#include <string>
#include <iostream>

namespace panzer {

TEUCHOS_UNIT_TEST(setupSubLocalMeshInfo, basic)
{

  // Make sure passing the function an empty parent info throws an error
  {
    Teuchos::RCP<panzer::LocalMeshInfoBase> mesh(new panzer::LocalMeshInfoBase);
    std::vector<panzer::LocalOrdinal> local_cells = {0};
    panzer::LocalMeshInfoBase sub_mesh;

    TEST_THROW(partitioning_utilities::setupSubLocalMeshInfo(*mesh,local_cells,sub_mesh),std::logic_error);
  }

  // Make sure passing the function an empty set of cells throws an error
  {
    Teuchos::RCP<panzer::LocalMeshInfoBase> mesh = generateLocalMeshInfoBase();
    std::vector<panzer::LocalOrdinal> local_cells = {};
    panzer::LocalMeshInfoBase sub_mesh;

    TEST_THROW(partitioning_utilities::setupSubLocalMeshInfo(*mesh,local_cells,sub_mesh),std::logic_error);
  }

  // Make sure that we can grab a couple of the cells as a local mesh
  {
    Teuchos::RCP<panzer::LocalMeshInfoBase> mesh = generateLocalMeshInfoBase();

    // Skip cell 1 to make it a ghost cell
    std::vector<panzer::LocalOrdinal> n_local_cells = {0,2};
    panzer::LocalMeshInfoBase sub_mesh;

    partitioning_utilities::setupSubLocalMeshInfo(*mesh,n_local_cells,sub_mesh);

    TEST_EQUALITY(sub_mesh.num_owned_cells,2);
    TEST_EQUALITY(sub_mesh.num_ghstd_cells,2);
    TEST_EQUALITY(sub_mesh.num_virtual_cells,1);

    // The following is just an analysis of the sub mesh
    // If any of these fail after changing the setupSubLocalMeshInfo call, it may not be because the changes were wrong.
    // It may just be a different ordering compared to what is currently expected (that is 100% OK, this test just needs to be updated).


    // What we had in 'mesh'
    //    g   o   o   o   v
    //  ---------------------
    //  | 3 | 0 | 1 | 2 | 4 |
    //  ---------------------
    //      0   1   2   3

    // What we expect in 'sub_mesh'
    //    g   o   g   o   v
    //  ---------------------
    //  | 3 | 0 | 2 | 1 | 4 |
    //  ---------------------
    //      3   2   1   0

    auto local_cells   = Kokkos::create_mirror_view(sub_mesh.local_cells);
    auto global_cells  = Kokkos::create_mirror_view(sub_mesh.global_cells);
    auto face_to_cells = Kokkos::create_mirror_view(sub_mesh.face_to_cells);
    auto face_to_lidx  = Kokkos::create_mirror_view(sub_mesh.face_to_lidx);
    auto cell_to_faces = Kokkos::create_mirror_view(sub_mesh.cell_to_faces);
    Kokkos::deep_copy(local_cells  , sub_mesh.local_cells);
    Kokkos::deep_copy(global_cells , sub_mesh.global_cells);
    Kokkos::deep_copy(face_to_cells, sub_mesh.face_to_cells);
    Kokkos::deep_copy(face_to_lidx , sub_mesh.face_to_lidx);
    Kokkos::deep_copy(cell_to_faces, sub_mesh.cell_to_faces);
    TEST_EQUALITY(local_cells(0),0);
    TEST_EQUALITY(local_cells(1),2);
    TEST_EQUALITY(local_cells(2),1);
    TEST_EQUALITY(local_cells(3),3);
    TEST_EQUALITY(local_cells(4),4);

    TEST_EQUALITY(global_cells(0),8);
    TEST_EQUALITY(global_cells(1),2);
    TEST_EQUALITY(global_cells(2),6);
    TEST_EQUALITY(global_cells(3),1);
    TEST_EQUALITY(global_cells(4),9);

    TEST_EQUALITY(face_to_cells(0,0),1);
    TEST_EQUALITY(face_to_cells(0,1),4);
    TEST_EQUALITY(face_to_cells(1,0),1);
    TEST_EQUALITY(face_to_cells(1,1),2);
    TEST_EQUALITY(face_to_cells(2,0),0);
    TEST_EQUALITY(face_to_cells(2,1),2);
    TEST_EQUALITY(face_to_cells(3,0),0);
    TEST_EQUALITY(face_to_cells(3,1),3);

    TEST_EQUALITY(face_to_lidx(0,0),1);
    TEST_EQUALITY(face_to_lidx(0,1),0);
    TEST_EQUALITY(face_to_lidx(1,0),0);
    TEST_EQUALITY(face_to_lidx(1,1),1);
    TEST_EQUALITY(face_to_lidx(2,0),1);
    TEST_EQUALITY(face_to_lidx(2,1),0);
    TEST_EQUALITY(face_to_lidx(3,0),0);
    TEST_EQUALITY(face_to_lidx(3,1),0);

    TEST_EQUALITY(cell_to_faces(0,0),3);
    TEST_EQUALITY(cell_to_faces(0,1),2);
    TEST_EQUALITY(cell_to_faces(1,0),1);
    TEST_EQUALITY(cell_to_faces(1,1),0);
    TEST_EQUALITY(cell_to_faces(2,0),2);
    TEST_EQUALITY(cell_to_faces(2,1),1);
    TEST_EQUALITY(cell_to_faces(3,0),3);
    TEST_EQUALITY(cell_to_faces(3,1),-1);
    TEST_EQUALITY(cell_to_faces(4,0),0);
    TEST_EQUALITY(cell_to_faces(4,1),-1);

  }
}

} // end namespace panzer
