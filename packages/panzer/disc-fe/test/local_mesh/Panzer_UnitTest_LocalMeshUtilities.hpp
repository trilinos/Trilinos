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

#ifndef __PANZER_UNITTEST_LOCALMESHUTILITIES_H__
#define __PANZER_UNITTEST_LOCALMESHUTILITIES_H__

#include "PanzerCore_config.hpp"

#include "Panzer_WorksetDescriptor.hpp"

#include "Shards_CellTopology.hpp"
#include "Shards_BasicTopologies.hpp"

#include "Panzer_LocalMeshInfo.hpp"

#include <vector>
#include <string>
#include <iostream>

namespace panzer {

// Generate a mesh info object
inline
Teuchos::RCP<panzer::LocalMeshInfo>
generateLocalMeshInfo()
{
  Teuchos::RCP<panzer::LocalMeshInfo> mesh_info(new panzer::LocalMeshInfo);

  // Local (=global) mesh indexes (4 owned, 2 virtual)
  //     v    o    o    o    o    v
  //  -------------------------------
  //  |  4 |  0 |  1 |  2 |  3 |  5 |
  //  -------------------------------

  // Element blocks
  //  -------------------------------
  //  |    | b0 | b0 | b1 | b1 |    |
  //  -------------------------------

  // Sidesets
  //  -------------------------------
  //  |    |    |    |    |    |    |
  //  -------------------------------
  //      s0        s2        s1

  Teuchos::RCP<const shards::CellTopology> cell_topology(new shards::CellTopology(shards::getCellTopologyData<shards::Line<2>>()));

  // Create the mesh info
  {

    // Mesh cell and face indexes
    //     v    o    o    o    o    v
    //  -------------------------------
    //  |  4 |  0 |  1 |  2 |  3 |  5 |
    //  -------------------------------
    //       0    1    2    3    4


    auto & mesh = *mesh_info;

    // Set some numbers
    mesh.num_owned_cells = 4;
    mesh.num_ghstd_cells = 0;
    mesh.num_virtual_cells = 2;

    // Set vertices
    mesh.cell_vertices = PHX::View<double***>("vertices",6,2,1);

    // Set local cells
    mesh.local_cells = PHX::View<panzer::LocalOrdinal*>("local_cells",6);
    mesh.local_cells(0) = 0;
    mesh.local_cells(1) = 1;
    mesh.local_cells(2) = 2;
    mesh.local_cells(3) = 3;
    mesh.local_cells(4) = 4;
    mesh.local_cells(5) = 5;

    // Set global cells
    mesh.global_cells = PHX::View<panzer::GlobalOrdinal*>("global_cells",6);
    mesh.global_cells(0) = 0;
    mesh.global_cells(1) = 1;
    mesh.global_cells(2) = 2;
    mesh.global_cells(3) = 3;
    mesh.global_cells(4) = 4;
    mesh.global_cells(5) = 5;

    // Set face to cells
    mesh.face_to_cells = PHX::View<panzer::LocalOrdinal*[2]>("face_to_cells",5);
    mesh.face_to_cells(0,0) = 0; mesh.face_to_cells(0,1) = 4;
    mesh.face_to_cells(1,0) = 0; mesh.face_to_cells(1,1) = 1;
    mesh.face_to_cells(2,0) = 1; mesh.face_to_cells(2,1) = 2;
    mesh.face_to_cells(3,0) = 2; mesh.face_to_cells(3,1) = 3;
    mesh.face_to_cells(4,0) = 3; mesh.face_to_cells(4,1) = 5;

    // Set face to local face indexes
    mesh.face_to_lidx = PHX::View<panzer::LocalOrdinal*[2]>("face_to_lidx",5);
    mesh.face_to_lidx(0,0) = 0; mesh.face_to_lidx(0,1) = 0;
    mesh.face_to_lidx(1,0) = 1; mesh.face_to_lidx(1,1) = 0;
    mesh.face_to_lidx(2,0) = 1; mesh.face_to_lidx(2,1) = 0;
    mesh.face_to_lidx(3,0) = 1; mesh.face_to_lidx(3,1) = 0;
    mesh.face_to_lidx(4,0) = 1; mesh.face_to_lidx(4,1) = 0;

    // Set cell to faces
    mesh.cell_to_faces = PHX::View<panzer::LocalOrdinal**>("cell_to_faces",6,2);
    mesh.cell_to_faces(0,0) = 0; mesh.cell_to_faces(0,1) = 1;
    mesh.cell_to_faces(1,0) = 1; mesh.cell_to_faces(1,1) = 2;
    mesh.cell_to_faces(2,0) = 2; mesh.cell_to_faces(2,1) = 3;
    mesh.cell_to_faces(3,0) = 3; mesh.cell_to_faces(3,1) = 4;
    mesh.cell_to_faces(4,0) = 0; mesh.cell_to_faces(4,1) =-1;
    mesh.cell_to_faces(5,0) = 4; mesh.cell_to_faces(5,1) =-1;

  }

  // Create block 0
  {

    // Element block cell and face indexes
    //    v    o     o    g
    //  -------------------------------
    //  |  3 |  0 |  1 |  2 |    |    |
    //  -------------------------------
    //       0    1    2

    LocalMeshBlockInfo & block = mesh_info->element_blocks["block0"];
    block.element_block_name = "block0";
    block.cell_topology = cell_topology;

    // Set some numbers
    block.num_owned_cells = 2;
    block.num_ghstd_cells = 1;
    block.num_virtual_cells = 1;

    // Set vertices
    block.cell_vertices = PHX::View<double***>("vertices",4,2,1);

    // Set local cells
    block.local_cells = PHX::View<panzer::LocalOrdinal*>("local_cells",4);
    block.local_cells(0) = 0;
    block.local_cells(1) = 1;
    block.local_cells(2) = 2;
    block.local_cells(3) = 3;

    // Set global cells
    block.global_cells = PHX::View<panzer::GlobalOrdinal*>("global_cells",4);
    block.global_cells(0) = 0;
    block.global_cells(1) = 1;
    block.global_cells(2) = 2;
    block.global_cells(3) = 4;

    // Set face to cells
    block.face_to_cells = PHX::View<panzer::LocalOrdinal*[2]>("face_to_cells",3);
    block.face_to_cells(0,0) = 0; block.face_to_cells(0,1) = 3;
    block.face_to_cells(1,0) = 0; block.face_to_cells(1,1) = 1;
    block.face_to_cells(2,0) = 1; block.face_to_cells(2,1) = 2;

    // Set face to local face indexes
    block.face_to_lidx = PHX::View<panzer::LocalOrdinal*[2]>("face_to_lidx",3);
    block.face_to_lidx(0,0) = 0; block.face_to_lidx(0,1) = 0;
    block.face_to_lidx(1,0) = 1; block.face_to_lidx(1,1) = 0;
    block.face_to_lidx(2,0) = 1; block.face_to_lidx(2,1) = 0;

    // Set cell to faces
    block.cell_to_faces = PHX::View<panzer::LocalOrdinal**>("cell_to_faces",4,2);
    block.cell_to_faces(0,0) = 0; block.cell_to_faces(0,1) = 1;
    block.cell_to_faces(1,0) = 1; block.cell_to_faces(1,1) = 2;
    block.cell_to_faces(2,0) = 2; block.cell_to_faces(2,1) =-1;
    block.cell_to_faces(3,0) = 0; block.cell_to_faces(3,1) =-1;

  }


  // Create block 1
  {

    // Element block cell and face indexes
    //              g    o     o    v
    //  -------------------------------
    //  |    |    |  2 |  0 |  1 |  3 |
    //  -------------------------------
    //                 0    1    2

    LocalMeshBlockInfo & block = mesh_info->element_blocks["block1"];
    block.element_block_name = "block1";
    block.cell_topology = cell_topology;

    // Set some numbers
    block.num_owned_cells = 2;
    block.num_ghstd_cells = 1;
    block.num_virtual_cells = 1;

    // Set vertices
    block.cell_vertices = PHX::View<double***>("vertices",4,2,1);

    // Set local cells
    block.local_cells = PHX::View<panzer::LocalOrdinal*>("local_cells",4);
    block.local_cells(0) = 0;
    block.local_cells(1) = 1;
    block.local_cells(2) = 2;
    block.local_cells(3) = 3;

    // Set global cells
    block.global_cells = PHX::View<panzer::GlobalOrdinal*>("global_cells",4);
    block.global_cells(0) = 2;
    block.global_cells(1) = 3;
    block.global_cells(2) = 1;
    block.global_cells(3) = 5;

    // Set face to cells
    block.face_to_cells = PHX::View<panzer::LocalOrdinal*[2]>("face_to_cells",3);
    block.face_to_cells(0,0) = 0; block.face_to_cells(0,1) = 2;
    block.face_to_cells(1,0) = 0; block.face_to_cells(1,1) = 1;
    block.face_to_cells(2,0) = 1; block.face_to_cells(2,1) = 3;

    // Set face to local face indexes
    block.face_to_lidx = PHX::View<panzer::LocalOrdinal*[2]>("face_to_lidx",3);
    block.face_to_lidx(0,0) = 0; block.face_to_lidx(0,1) = 0;
    block.face_to_lidx(1,0) = 1; block.face_to_lidx(1,1) = 0;
    block.face_to_lidx(2,0) = 1; block.face_to_lidx(2,1) = 0;

    // Set cell to faces
    block.cell_to_faces = PHX::View<panzer::LocalOrdinal**>("cell_to_faces",4,2);
    block.cell_to_faces(0,0) = 0; block.cell_to_faces(0,1) = 1;
    block.cell_to_faces(1,0) = 1; block.cell_to_faces(1,1) = 2;
    block.cell_to_faces(2,0) = 0; block.cell_to_faces(2,1) =-1;
    block.cell_to_faces(3,0) = 2; block.cell_to_faces(3,1) =-1;
  }

  // Create sideset 0
  {

    // Sideset 0 cell and face indexes
    //    v    o
    //  -------------------------------
    //  |  1 |  0 |    |    |    |    |
    //  -------------------------------
    //       0

    LocalMeshSidesetInfo & sideset = mesh_info->sidesets["block0"]["sideset0"];
    sideset.element_block_name = "block0";
    sideset.sideset_name = "sideset0";
    sideset.cell_topology = cell_topology;

    // Set some numbers
    sideset.num_owned_cells = 1;
    sideset.num_ghstd_cells = 0;
    sideset.num_virtual_cells = 1;

    // Set vertices
    sideset.cell_vertices = PHX::View<double***>("vertices",2,2,1);

    // Set local cells
    sideset.local_cells = PHX::View<panzer::LocalOrdinal*>("local_cells",2);
    sideset.local_cells(0) = 0;
    sideset.local_cells(1) = 1;

    // Set global cells
    sideset.global_cells = PHX::View<panzer::GlobalOrdinal*>("global_cells",2);
    sideset.global_cells(0) = 0;
    sideset.global_cells(1) = 4;

    // Set face to cells
    sideset.face_to_cells = PHX::View<panzer::LocalOrdinal*[2]>("face_to_cells",1);
    sideset.face_to_cells(0,0) = 0; sideset.face_to_cells(0,1) = 1;

    // Set face to local face indexes
    sideset.face_to_lidx = PHX::View<panzer::LocalOrdinal*[2]>("face_to_lidx",1);
    sideset.face_to_lidx(0,0) = 0; sideset.face_to_lidx(0,1) = 0;

    // Set cell to faces
    sideset.cell_to_faces = PHX::View<panzer::LocalOrdinal**>("cell_to_faces",2,2);
    sideset.cell_to_faces(0,0) = 0; sideset.cell_to_faces(0,1) =-1;
    sideset.cell_to_faces(1,0) = 0; sideset.cell_to_faces(1,1) =-1;
  }

  // Create sideset 1
  {

    // Sideset 1 cell and face indexes
    //                         o    v
    //  -------------------------------
    //  |    |    |    |    |  0 |  1 |
    //  -------------------------------
    //                           0

    LocalMeshSidesetInfo & sideset = mesh_info->sidesets["block1"]["sideset1"];
    sideset.element_block_name = "block1";
    sideset.sideset_name = "sideset1";
    sideset.cell_topology = cell_topology;

    // Set some numbers
    sideset.num_owned_cells = 1;
    sideset.num_ghstd_cells = 0;
    sideset.num_virtual_cells = 1;

    // Set vertices
    sideset.cell_vertices = PHX::View<double***>("vertices",2,2,1);

    // Set local cells
    sideset.local_cells = PHX::View<panzer::LocalOrdinal*>("local_cells",2);
    sideset.local_cells(0) = 0;
    sideset.local_cells(1) = 1;

    // Set global cells
    sideset.global_cells = PHX::View<panzer::GlobalOrdinal*>("global_cells",2);
    sideset.global_cells(0) = 3;
    sideset.global_cells(1) = 5;

    // Set face to cells
    sideset.face_to_cells = PHX::View<panzer::LocalOrdinal*[2]>("face_to_cells",1);
    sideset.face_to_cells(0,0) = 0; sideset.face_to_cells(0,1) = 1;

    // Set face to local face indexes
    sideset.face_to_lidx = PHX::View<panzer::LocalOrdinal*[2]>("face_to_lidx",1);
    sideset.face_to_lidx(0,0) = 0; sideset.face_to_lidx(0,1) = 0;

    // Set cell to faces
    sideset.cell_to_faces = PHX::View<panzer::LocalOrdinal**>("cell_to_faces",2,2);
    sideset.cell_to_faces(0,0) = 0; sideset.cell_to_faces(0,1) =-1;
    sideset.cell_to_faces(1,0) = 0; sideset.cell_to_faces(1,1) =-1;
  }

  // Create sideset 2
  {

    // Sideset 2 cell and face indexes
    //               o    g
    //  -------------------------------
    //  |    |    |  0 |  1 |    |    |
    //  -------------------------------
    //                 0

    LocalMeshSidesetInfo & sideset = mesh_info->sidesets["block0"]["sideset2"];
    sideset.element_block_name = "block0";
    sideset.sideset_name = "sideset2";
    sideset.cell_topology = cell_topology;

    // Set some numbers
    sideset.num_owned_cells = 1;
    sideset.num_ghstd_cells = 1;
    sideset.num_virtual_cells = 0;

    // Set vertices
    sideset.cell_vertices = PHX::View<double***>("vertices",2,2,1);

    // Set local cells
    sideset.local_cells = PHX::View<panzer::LocalOrdinal*>("local_cells",2);
    sideset.local_cells(0) = 0;
    sideset.local_cells(1) = 1;

    // Set global cells
    sideset.global_cells = PHX::View<panzer::GlobalOrdinal*>("global_cells",2);
    sideset.global_cells(0) = 1;
    sideset.global_cells(1) = 2;

    // Set face to cells
    sideset.face_to_cells = PHX::View<panzer::LocalOrdinal*[2]>("face_to_cells",1);
    sideset.face_to_cells(0,0) = 0; sideset.face_to_cells(0,1) = 1;

    // Set face to local face indexes
    sideset.face_to_lidx = PHX::View<panzer::LocalOrdinal*[2]>("face_to_lidx",1);
    sideset.face_to_lidx(0,0) = 1; sideset.face_to_lidx(0,1) = 0;

    // Set cell to faces
    sideset.cell_to_faces = PHX::View<panzer::LocalOrdinal**>("cell_to_faces",2,2);
    sideset.cell_to_faces(0,0) =-1; sideset.cell_to_faces(0,1) = 0;
    sideset.cell_to_faces(1,0) = 0; sideset.cell_to_faces(1,1) =-1;
  }
  {

    // Sideset 2 cell and face indexes
    //               g    o
    //  -------------------------------
    //  |    |    |  1 |  0 |    |    |
    //  -------------------------------
    //                 0

    LocalMeshSidesetInfo & sideset = mesh_info->sidesets["block1"]["sideset2"];
    sideset.element_block_name = "block1";
    sideset.sideset_name = "sideset2";
    sideset.cell_topology = cell_topology;

    // Set some numbers
    sideset.num_owned_cells = 1;
    sideset.num_ghstd_cells = 1;
    sideset.num_virtual_cells = 0;

    // Set vertices
    sideset.cell_vertices = PHX::View<double***>("vertices",2,2,1);

    // Set local cells
    sideset.local_cells = PHX::View<panzer::LocalOrdinal*>("local_cells",2);
    sideset.local_cells(0) = 0;
    sideset.local_cells(1) = 1;

    // Set global cells
    sideset.global_cells = PHX::View<panzer::GlobalOrdinal*>("global_cells",2);
    sideset.global_cells(0) = 2;
    sideset.global_cells(1) = 1;

    // Set face to cells
    sideset.face_to_cells = PHX::View<panzer::LocalOrdinal*[2]>("face_to_cells",1);
    sideset.face_to_cells(0,0) = 0; sideset.face_to_cells(0,1) = 1;

    // Set face to local face indexes
    sideset.face_to_lidx = PHX::View<panzer::LocalOrdinal*[2]>("face_to_lidx",1);
    sideset.face_to_lidx(0,0) = 0; sideset.face_to_lidx(0,1) = 1;

    // Set cell to faces
    sideset.cell_to_faces = PHX::View<panzer::LocalOrdinal**>("cell_to_faces",2,2);
    sideset.cell_to_faces(0,0) = 0; sideset.cell_to_faces(0,1) =-1;
    sideset.cell_to_faces(1,0) =-1; sideset.cell_to_faces(1,1) = 0;
  }

  return mesh_info;
}


// Generate a mesh info object
inline
Teuchos::RCP<panzer::LocalMeshInfoBase>
generateLocalMeshInfoBase()
{
  Teuchos::RCP<panzer::LocalMeshInfoBase> block_rcp(new panzer::LocalMeshInfoBase);
  panzer::LocalMeshInfoBase & block = *block_rcp;

  // Local mesh indexes (3 owned, 1 ghost, 1 virtual)
  //    g   o   o   o   v
  //  ---------------------
  //  | 3 | 0 | 1 | 2 | 4 |
  //  ---------------------
  //      0   1   2   3

  // Set some numbers
  block.num_owned_cells = 3;
  block.num_ghstd_cells = 1;
  block.num_virtual_cells = 1;

  // Set vertices
  block.cell_vertices = PHX::View<double***>("vertices",5,2,1);

  // Set local cells
  block.local_cells = PHX::View<panzer::LocalOrdinal*>("local_cells",5);
  block.local_cells(0) = 0;
  block.local_cells(1) = 1;
  block.local_cells(2) = 2;
  block.local_cells(3) = 3;
  block.local_cells(4) = 4;

  // Set global cells - arbitrary
  block.global_cells = PHX::View<panzer::GlobalOrdinal*>("global_cells",5);
  block.global_cells(0) = 8;
  block.global_cells(1) = 6;
  block.global_cells(2) = 2;
  block.global_cells(3) = 1;
  block.global_cells(4) = 9;

  // Set face to cells
  block.face_to_cells = PHX::View<panzer::LocalOrdinal*[2]>("face_to_cells",4);
  block.face_to_cells(0,0) = 0; block.face_to_cells(0,1) = 3;
  block.face_to_cells(1,0) = 0; block.face_to_cells(1,1) = 1;
  block.face_to_cells(2,0) = 1; block.face_to_cells(2,1) = 2;
  block.face_to_cells(3,0) = 2; block.face_to_cells(3,1) = 4;

  // Set face to local face indexes
  block.face_to_lidx = PHX::View<panzer::LocalOrdinal*[2]>("face_to_lidx",4);
  block.face_to_lidx(0,0) = 0; block.face_to_lidx(0,1) = 0;
  block.face_to_lidx(1,0) = 1; block.face_to_lidx(1,1) = 0;
  block.face_to_lidx(2,0) = 1; block.face_to_lidx(2,1) = 0;
  block.face_to_lidx(3,0) = 1; block.face_to_lidx(3,1) = 0;

  // Set cell to faces
  block.cell_to_faces = PHX::View<panzer::LocalOrdinal**>("cell_to_faces",5,2);
  block.cell_to_faces(0,0) = 0; block.cell_to_faces(0,1) = 1;
  block.cell_to_faces(1,0) = 1; block.cell_to_faces(1,1) = 2;
  block.cell_to_faces(2,0) = 2; block.cell_to_faces(2,1) = 3;
  block.cell_to_faces(3,0) = 0; block.cell_to_faces(3,1) =-1;
  block.cell_to_faces(4,0) = 3; block.cell_to_faces(4,1) =-1;

  return block_rcp;
}


} // end namespace panzer

#endif
