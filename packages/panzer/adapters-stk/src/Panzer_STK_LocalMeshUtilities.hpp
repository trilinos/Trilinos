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

#ifndef PANZER_STK_LOCAL_MESH_UTILITIES_HPP
#define PANZER_STK_LOCAL_MESH_UTILITIES_HPP

#include "Kokkos_View.hpp"
#include "Kokkos_DynRankView.hpp"

#include "Phalanx_KokkosDeviceTypes.hpp"

#include "Panzer_STK_Interface.hpp"

#include "Tpetra_Import.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"
#include <vector>
#include <string>

namespace panzer
{

// TODO: Move this class to disc-fe
template <typename LO, typename GO>
struct LocalMeshInfo {

  // cell ids
  Kokkos::View<const GO*> owned_cells;
  Kokkos::View<const GO*> ghstd_cells;
  Kokkos::View<const LO*> virtual_cells;

  // vertices
  Kokkos::DynRankView<double,PHX::Device> owned_vertices;
  Kokkos::DynRankView<double,PHX::Device> ghstd_vertices;

  // Face to neighbors
  Kokkos::View<const LO*[2]> face_to_cells;    // this is local numbering
                                                // that indexes first into
                                                // the owned_cells and then
                                                // into ghstd_cells
  Kokkos::View<const LO*[2]> face_to_lidx;     // maps faces to the cell local
                                                // face index
  Kokkos::View<const LO**> cell_to_face;       // using cell local numbering,
                                                // produce face index

};

}

namespace panzer_stk
{

template <typename LO, typename GO>
panzer::LocalMeshInfo<LO,GO>
generateLocalMeshInfo(const panzer_stk::STK_Interface & mesh,
                      const std::string & element_block_name);


}
#endif
