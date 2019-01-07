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

#include "PanzerCore_config.hpp"
#include "Kokkos_View.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Teuchos_Assert.hpp"

namespace panzer {

template <typename LO, typename GO>
struct LocalMeshPartition;

class SubcellConnectivity
{
public:

  /// Default constructor
  SubcellConnectivity():_num_subcells(-1),_num_cells(-1){}

  /// Default destructor
  ~SubcellConnectivity() = default;

  /**
   * \brief Gives number of subcells (e.g. faces) in connectivity
   *
   * \return Number of subcells associated with the cells
   */
  KOKKOS_INLINE_FUNCTION
  int numSubcells() const {return _num_subcells;}

  /**
   * \brief Gives number of cells in connectivity
   *
   * \return Number of subcells associated with the cells
   */
  KOKKOS_INLINE_FUNCTION
  int numCells() const {return _num_cells;}

  /**
   * \brief gives number of subcells (e.g. faces) found on a given cell
   *
   * \throw If cell is out of range
   *
   * param[in] Cell index
   *
   * \return Number of subcells on a given cell
   */
  KOKKOS_INLINE_FUNCTION
  int numSubcellsOnCell(const int cell) const;

  /**
   * \brief Returns the number of cells attached to a given subcell
   *
   * For example:
   * 1) A face joins two cells.
   * 2) A node joins four cells on a 2D rectilinear mesh (eight in 3D)
   *
   * \param[in] subcell Subcell index
   *
   * \return Number of cells connected to subcell
   */
  KOKKOS_INLINE_FUNCTION
  int numCellsOnSubcell(const int subcell) const;

  /**
   * \brief Get the subcell index for a given cell and local subcell index
   *
   * A local subcell index is the local indexing scheme for the cell.
   *
   * For example:
   * 1) A quad cell has four faces indexed by 0,1,2,3 (local subcell indexes)
   * 2) A hex cell has eight nodes indexed by 0,1,2,3,4,5,6,7 (local subcell indexes)
   *
   * \param[in] cell Cell index
   * \param[in] local_subcell_index Index of subcell in cell local indexing
   *
   * \return Subcell index
   */
  KOKKOS_INLINE_FUNCTION
  int subcellForCell(const int cell, const int local_subcell_index) const;

  /**
   * \brief Get the cell for a given subcell and a local_cell_index
   *
   * A local cell index is the indexing scheme local to a subcell.
   *
   * For example:
   * 1) A 1D mesh has subcells (nodes) connects two cells (lines) with local cell indexes 0,1
   * 2) A 2D quad can have nodal subcells (on structured mesh) that connect four cells with local cell indexes 0,1,2,3
   *
   * \param[in] subcell Subcell index
   * \param[in] local_cell_index
   *
   * \return Cell index
   */
  KOKKOS_INLINE_FUNCTION
  int cellForSubcell(const int subcell, const int local_cell_index) const;

  /**
   * \brief Get the local subcell index given a subcell and a local cell index
   *
   * This is the mapping between local subcell indexes and local cell indexes
   *
   * \param[in] subcell Subcell index
   * \param[in] local_cell_index Local cell index on subcell
   *
   * \return Local subcell index for cell identified by subcell index and local_cell_index
   */
  KOKKOS_INLINE_FUNCTION
  int localSubcellForSubcell(const int subcell, const int local_cell_index) const;

protected:

  /// Number of subcells for a given number of cells
  int _num_subcells;

  /// Number of cells
  int _num_cells;

  /// Adjacency array for indexing into subcell_to_cells array
  Kokkos::View<int*, PHX::Device> _subcell_to_cells_adj;

  /// Mapping from subcells to cells
  Kokkos::View<int*, PHX::Device> _subcell_to_cells;

  /// Mapping from subcell indexes to local subcell indexes
  Kokkos::View<int*, PHX::Device> _subcell_to_local_subcells;

  /// Adjacency array for indexing into cell_to_subcells array
  Kokkos::View<int*, PHX::Device> _cell_to_subcells_adj;

  /// Mapping from cells to subcells
  Kokkos::View<int*, PHX::Device> _cell_to_subcells;

};

/**
 * \class FaceConnectivity
 *
 * \brief Generates a SubcellConnectivity associated with faces and cells given a partition of the local mesh
 */
class FaceConnectivity:
    public SubcellConnectivity
{
public:

  /// Default constructor
  FaceConnectivity() = default;

  /// Default destructor
  ~FaceConnectivity() = default;

  /**
   * \brief Setup the face connectivity from a partition of the local mesh
   *
   * \param[in] partition Partition of mesh
   */
  void setup(const panzer::LocalMeshPartition<int,panzer::Ordinal64> & partition);

};

// **********************************
// Inlined functions
// **********************************

int
SubcellConnectivity::
numSubcellsOnCell(const int cell) const
{
#ifdef PANZER_DEBUG
#ifndef KOKKOS_ENABLE_CUDA
  TEUCHOS_ASSERT(cell >= 0 and cell < _num_cells);
#endif
#endif
  return _cell_to_subcells_adj(cell+1)-_cell_to_subcells_adj(cell);
}

int
SubcellConnectivity::
numCellsOnSubcell(const int subcell) const
{
#ifdef PANZER_DEBUG
#ifndef KOKKOS_ENABLE_CUDA
  TEUCHOS_ASSERT(subcell >= 0 and subcell < _num_subcells);
#endif
#endif
  return _subcell_to_cells_adj(subcell+1)-_subcell_to_cells_adj(subcell);
}

int
SubcellConnectivity::
subcellForCell(const int cell, const int local_subcell_index) const
{
#ifdef PANZER_DEBUG
#ifndef KOKKOS_ENABLE_CUDA
  TEUCHOS_ASSERT(cell >= 0 and cell < _num_cells);
  TEUCHOS_ASSERT(local_subcell_index < numSubcellsOnCell(cell));
#endif
#endif
  const int index = _cell_to_subcells_adj(cell)+local_subcell_index;
  return _cell_to_subcells(index);
}

int
SubcellConnectivity::
cellForSubcell(const int subcell, const int local_cell_index) const
{
#ifdef PANZER_DEBUG
#ifndef KOKKOS_ENABLE_CUDA
  TEUCHOS_ASSERT(subcell >= 0 and subcell < _num_subcells);
  TEUCHOS_ASSERT(local_cell_index < numCellsOnSubcell(subcell));
#endif
#endif
  const int index = _subcell_to_cells_adj(subcell)+local_cell_index;
  return _subcell_to_cells(index);
}

int
SubcellConnectivity::
localSubcellForSubcell(const int subcell, const int local_cell_index) const
{
#ifdef PANZER_DEBUG
#ifndef KOKKOS_ENABLE_CUDA
  TEUCHOS_ASSERT(subcell >= 0 and subcell < _num_subcells);
  TEUCHOS_ASSERT(local_cell_index < numCellsOnSubcell(subcell));
#endif
#endif
  const int index = _subcell_to_cells_adj(subcell)+local_cell_index;
  return _subcell_to_local_subcells(index);
}

} // namespace panzer

#endif
