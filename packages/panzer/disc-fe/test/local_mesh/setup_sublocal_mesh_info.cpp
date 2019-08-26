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
    std::vector<panzer::LocalOrdinal> local_cells = {0,2};
    panzer::LocalMeshInfoBase sub_mesh;

    partitioning_utilities::setupSubLocalMeshInfo(*mesh,local_cells,sub_mesh);

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

    TEST_EQUALITY(sub_mesh.local_cells(0),0);
    TEST_EQUALITY(sub_mesh.local_cells(1),2);
    TEST_EQUALITY(sub_mesh.local_cells(2),1);
    TEST_EQUALITY(sub_mesh.local_cells(3),3);
    TEST_EQUALITY(sub_mesh.local_cells(4),4);

    TEST_EQUALITY(sub_mesh.global_cells(0),8);
    TEST_EQUALITY(sub_mesh.global_cells(1),2);
    TEST_EQUALITY(sub_mesh.global_cells(2),6);
    TEST_EQUALITY(sub_mesh.global_cells(3),1);
    TEST_EQUALITY(sub_mesh.global_cells(4),9);

    TEST_EQUALITY(sub_mesh.face_to_cells(0,0),1);
    TEST_EQUALITY(sub_mesh.face_to_cells(0,1),4);
    TEST_EQUALITY(sub_mesh.face_to_cells(1,0),1);
    TEST_EQUALITY(sub_mesh.face_to_cells(1,1),2);
    TEST_EQUALITY(sub_mesh.face_to_cells(2,0),0);
    TEST_EQUALITY(sub_mesh.face_to_cells(2,1),2);
    TEST_EQUALITY(sub_mesh.face_to_cells(3,0),0);
    TEST_EQUALITY(sub_mesh.face_to_cells(3,1),3);

    TEST_EQUALITY(sub_mesh.face_to_lidx(0,0),1);
    TEST_EQUALITY(sub_mesh.face_to_lidx(0,1),0);
    TEST_EQUALITY(sub_mesh.face_to_lidx(1,0),0);
    TEST_EQUALITY(sub_mesh.face_to_lidx(1,1),1);
    TEST_EQUALITY(sub_mesh.face_to_lidx(2,0),1);
    TEST_EQUALITY(sub_mesh.face_to_lidx(2,1),0);
    TEST_EQUALITY(sub_mesh.face_to_lidx(3,0),0);
    TEST_EQUALITY(sub_mesh.face_to_lidx(3,1),0);

    TEST_EQUALITY(sub_mesh.cell_to_faces(0,0),3);
    TEST_EQUALITY(sub_mesh.cell_to_faces(0,1),2);
    TEST_EQUALITY(sub_mesh.cell_to_faces(1,0),1);
    TEST_EQUALITY(sub_mesh.cell_to_faces(1,1),0);
    TEST_EQUALITY(sub_mesh.cell_to_faces(2,0),2);
    TEST_EQUALITY(sub_mesh.cell_to_faces(2,1),1);
    TEST_EQUALITY(sub_mesh.cell_to_faces(3,0),3);
    TEST_EQUALITY(sub_mesh.cell_to_faces(3,1),-1);
    TEST_EQUALITY(sub_mesh.cell_to_faces(4,0),0);
    TEST_EQUALITY(sub_mesh.cell_to_faces(4,1),-1);

  }
}

} // end namespace panzer
