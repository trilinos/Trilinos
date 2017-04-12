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


#ifndef PANZER_SUBCELL_CONNECTIVITY_HPP
#define PANZER_SUBCELL_CONNECTIVITY_HPP

#include "Kokkos_View.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"

namespace panzer {

template <typename LO, typename GO>
class LocalMeshPartition;

class SubcellConnectivity
{
public:

  SubcellConnectivity():_num_subcells(-1){}
  ~SubcellConnectivity() = default;

  /// Number of faces associated with workset
  int numSubcells() const {return _num_subcells;}

  /// Number of faces on a given cell
  int numSubcellsOnCell(const int cell) const {return _cell_to_subcells_adj(cell+1)-_cell_to_subcells_adj(cell);}

  int numCellsOnSubcell(const int subcell) const {return _subcell_to_cells_adj(subcell+1)-_subcell_to_cells_adj(subcell);}

  int subcellForCell(const int cell, const int local_subcell_index) const
  {
    const int index = _cell_to_subcells_adj(cell)+local_subcell_index;
    return _cell_to_subcells(index);
  };

  int cellForSubcell(const int subcell, const int local_cell_index) const
  {
    const int index = _subcell_to_cells_adj(subcell)+local_cell_index;
    return _subcell_to_cells(index);
  };

  int localSubcellForSubcell(const int subcell, const int local_cell_index) const
  {
    const int index = _subcell_to_cells_adj(subcell)+local_cell_index;
    return _subcell_to_local_subcells(index);
  };

protected:
  int _num_subcells;
  Kokkos::View<int*, PHX::Device> _subcell_to_cells_adj;
  Kokkos::View<int*, PHX::Device> _subcell_to_cells;
  Kokkos::View<int*, PHX::Device> _subcell_to_local_subcells;
  Kokkos::View<int*, PHX::Device> _cell_to_subcells_adj;
  Kokkos::View<int*, PHX::Device> _cell_to_subcells;

};


class FaceConnectivity:
    public SubcellConnectivity
{
public:

  FaceConnectivity() = default;
  ~FaceConnectivity() = default;

  void setup(const panzer::LocalMeshPartition<int,int> & partition);

protected:

};

} // namespace panzer

#endif
