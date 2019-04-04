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
setup(const panzer::LocalMeshPartition<int,panzer::Ordinal64> & partition)
{
  const int num_cells = partition.cell_to_faces.extent(0);
  const int num_faces = partition.face_to_cells.extent(0);
  const int num_faces_per_cell = partition.cell_to_faces.extent(1);
  const int num_cells_per_face = 2;

  _num_subcells = num_faces;
  _num_cells = num_cells;

  _subcell_to_cells_adj = Kokkos::View<int*>("subcell_to_cells_adj", num_faces+1);
  _subcell_to_cells = Kokkos::View<int*>("subcell_to_cells", num_faces*num_cells_per_face);
  _subcell_to_local_subcells = Kokkos::View<int*>("subcell_to_local_subcells", num_faces*num_cells_per_face);
  _cell_to_subcells_adj = Kokkos::View<int*>("cell_to_subcells_adj", num_cells+1);
  _cell_to_subcells = Kokkos::View<int*>("cell_to_subcells", num_cells*num_faces_per_cell);

  _subcell_to_cells_adj(0)=0;
  for(int face=0;face<num_faces;++face){
    _subcell_to_cells_adj(face+1) =_subcell_to_cells_adj(face)+num_cells_per_face;
    _subcell_to_cells(num_cells_per_face*face + 0) = partition.face_to_cells(face,0);
    _subcell_to_cells(num_cells_per_face*face + 1) = partition.face_to_cells(face,1);
    _subcell_to_local_subcells(num_cells_per_face*face + 0) = partition.face_to_lidx(face,0);
    _subcell_to_local_subcells(num_cells_per_face*face + 1) = partition.face_to_lidx(face,1);
  }

  _cell_to_subcells_adj(0)=0;
  for(int cell=0;cell<num_cells;++cell){
    _cell_to_subcells_adj(cell+1) =_cell_to_subcells_adj(cell)+num_faces_per_cell;
    for(int local_face=0;local_face<num_faces_per_cell;++local_face){
      _cell_to_subcells(num_faces_per_cell*cell+local_face) = partition.cell_to_faces(cell,local_face);
    }
  }

}

}
