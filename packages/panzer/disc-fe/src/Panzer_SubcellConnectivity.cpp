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

#include "Panzer_SubcellConnectivity.hpp"
#include "Panzer_LocalMeshInfo.hpp"

namespace panzer
{

void
FaceConnectivity::
setup(const panzer::LocalMeshPartition & partition)
{
  const int num_cells = partition.cell_to_faces.extent(0);
  const int num_faces = partition.face_to_cells.extent(0);
  const int num_faces_per_cell = partition.cell_to_faces.extent(1);
  const int num_cells_per_face = 2;

  _subcell_to_cells_adj = PHX::View<int*>("subcell_to_cells_adj", num_faces+1);
  _subcell_to_cells = PHX::View<int*>("subcell_to_cells", num_faces*num_cells_per_face);
  _subcell_to_local_subcells = PHX::View<int*>("subcell_to_local_subcells", num_faces*num_cells_per_face);
  _cell_to_subcells_adj = PHX::View<int*>("cell_to_subcells_adj", num_cells+1);
  _cell_to_subcells = PHX::View<int*>("cell_to_subcells", num_cells*num_faces_per_cell);

  // Host copies
  _subcell_to_cells_adj_host = PHX::View<int*>::HostMirror("subcell_to_cells_adj_host", num_faces+1);
  _subcell_to_cells_host = PHX::View<int*>::HostMirror("subcell_to_cells_host", num_faces*num_cells_per_face);
  _subcell_to_local_subcells_host = PHX::View<int*>::HostMirror("subcell_to_local_subcells_host", num_faces*num_cells_per_face);
  _cell_to_subcells_adj_host = PHX::View<int*>::HostMirror("cell_to_subcells_adj_host", num_cells+1);
  _cell_to_subcells_host = PHX::View<int*>::HostMirror("cell_to_subcells_host", num_cells*num_faces_per_cell);

  // This line not needed since kokkos initializes the arrays above to zero
  //_subcell_to_cells_adj(0)=0;

  {
    auto face_to_cells = partition.face_to_cells;
    auto face_to_lidx = partition.face_to_lidx;
    auto subcell_to_cells_adj = _subcell_to_cells_adj;
    auto subcell_to_cells = _subcell_to_cells;
    auto subcell_to_local_subcells = _subcell_to_local_subcells;
    Kokkos::parallel_for("subcell connectivity 0",num_faces,KOKKOS_LAMBDA (const int face) {
      subcell_to_cells_adj(face+1) = (face * num_cells_per_face) + num_cells_per_face;
      subcell_to_cells(num_cells_per_face*face + 0) = face_to_cells(face,0);
      subcell_to_cells(num_cells_per_face*face + 1) = face_to_cells(face,1);
      subcell_to_local_subcells(num_cells_per_face*face + 0) = face_to_lidx(face,0);
      subcell_to_local_subcells(num_cells_per_face*face + 1) = face_to_lidx(face,1);
    });
    PHX::Device::execution_space().fence();
  }

  // This line not needed since kokkos initializes the arrays above to zero
  //_cell_to_subcells_adj(0)=0;

  {
    auto cell_to_faces = partition.cell_to_faces;
    auto cell_to_subcells_adj = _cell_to_subcells_adj;
    auto cell_to_subcells = _cell_to_subcells;
    Kokkos::parallel_for("subcell connectivity 1",num_cells,KOKKOS_LAMBDA (const int cell) {
      cell_to_subcells_adj(cell+1) = (cell * num_faces_per_cell) + num_faces_per_cell;
      for(int local_face=0;local_face<num_faces_per_cell;++local_face){
        cell_to_subcells(num_faces_per_cell*cell+local_face) = cell_to_faces(cell,local_face);
      }
    });
    PHX::Device::execution_space().fence();
  }

  // Copy values to host
  Kokkos::deep_copy(_subcell_to_cells_adj_host, _subcell_to_cells_adj);
  Kokkos::deep_copy(_subcell_to_cells_host, _subcell_to_cells);
  Kokkos::deep_copy(_subcell_to_local_subcells_host,_subcell_to_local_subcells);
  Kokkos::deep_copy(_cell_to_subcells_adj_host,_cell_to_subcells_adj);
  Kokkos::deep_copy(_cell_to_subcells_host,_cell_to_subcells);
}

}
