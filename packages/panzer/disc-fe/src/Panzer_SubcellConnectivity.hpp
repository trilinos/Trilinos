// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_SUBCELL_CONNECTIVITY_HPP
#define PANZER_SUBCELL_CONNECTIVITY_HPP

#include "PanzerCore_config.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Teuchos_Assert.hpp"

namespace panzer {

struct LocalMeshPartition;

class SubcellConnectivity
{
public:

  /// Default constructor
  SubcellConnectivity() = default;

  /// Default destructor
  ~SubcellConnectivity() = default;

  /**
   * \brief Gives number of subcells (e.g. faces) in connectivity
   *
   * \return Number of subcells associated with the cells
   */
  KOKKOS_INLINE_FUNCTION
  int numSubcells() const {return _subcell_to_cells_adj.extent(0)-1;}

  /**
   * \brief Gives number of cells in connectivity
   *
   * \return Number of subcells associated with the cells
   */
  KOKKOS_INLINE_FUNCTION
  int numCells() const {return _cell_to_subcells_adj.extent(0)-1;}

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
  inline
  int numSubcellsOnCellHost(const int cell) const;

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
  inline
  int numCellsOnSubcellHost(const int subcell) const;

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
  inline
  int subcellForCellHost(const int cell, const int local_subcell_index) const;

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
  inline
  int cellForSubcellHost(const int subcell, const int local_cell_index) const;

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
  inline
  int localSubcellForSubcellHost(const int subcell, const int local_cell_index) const;

protected:

  /// Adjacency array for indexing into subcell_to_cells array
  PHX::View<int*> _subcell_to_cells_adj;
  PHX::View<int*>::HostMirror _subcell_to_cells_adj_host;

  /// Mapping from subcells to cells
  PHX::View<int*> _subcell_to_cells;
  PHX::View<int*>::HostMirror _subcell_to_cells_host;

  /// Mapping from subcell indexes to local subcell indexes
  PHX::View<int*> _subcell_to_local_subcells;
  PHX::View<int*>::HostMirror _subcell_to_local_subcells_host;

  /// Adjacency array for indexing into cell_to_subcells array
  PHX::View<int*> _cell_to_subcells_adj;
  PHX::View<int*>::HostMirror _cell_to_subcells_adj_host;

  /// Mapping from cells to subcells
  PHX::View<int*> _cell_to_subcells;
  PHX::View<int*>::HostMirror _cell_to_subcells_host;

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
  void setup(const panzer::LocalMeshPartition & partition);

};

// **********************************
// Inlined functions
// **********************************

PHALANX_HIP_HACK_KOKKOS_FUNCTION
int
SubcellConnectivity::
numSubcellsOnCell(const int cell) const
{
#ifdef PANZER_DEBUG
  KOKKOS_ASSERT(cell >= 0 and cell < numCells());
#endif
  return _cell_to_subcells_adj(cell+1)-_cell_to_subcells_adj(cell);
}

int
SubcellConnectivity::
numSubcellsOnCellHost(const int cell) const
{
#ifdef PANZER_DEBUG
  KOKKOS_ASSERT(cell >= 0 and cell < numCells());
#endif
  return _cell_to_subcells_adj_host(cell+1)-_cell_to_subcells_adj_host(cell);
}

PHALANX_HIP_HACK_KOKKOS_FUNCTION
int
SubcellConnectivity::
numCellsOnSubcell(const int subcell) const
{
#ifdef PANZER_DEBUG
  KOKKOS_ASSERT(subcell >= 0 and subcell < numSubcells());
#endif
  return _subcell_to_cells_adj(subcell+1)-_subcell_to_cells_adj(subcell);
}

int
SubcellConnectivity::
numCellsOnSubcellHost(const int subcell) const
{
#ifdef PANZER_DEBUG
  KOKKOS_ASSERT(subcell >= 0 and subcell < numSubcells());
#endif
  return _subcell_to_cells_adj_host(subcell+1)-_subcell_to_cells_adj_host(subcell);
}

PHALANX_HIP_HACK_KOKKOS_FUNCTION
int
SubcellConnectivity::
subcellForCell(const int cell, const int local_subcell_index) const
{
#ifdef PANZER_DEBUG
  KOKKOS_ASSERT(cell >= 0 and cell < numCells());
  KOKKOS_ASSERT(local_subcell_index < numSubcellsOnCell(cell));
#endif
  const int index = _cell_to_subcells_adj(cell)+local_subcell_index;
  return _cell_to_subcells(index);
}

int
SubcellConnectivity::
subcellForCellHost(const int cell, const int local_subcell_index) const
{
#ifdef PANZER_DEBUG
  KOKKOS_ASSERT(cell >= 0 and cell < numCells());
  KOKKOS_ASSERT(local_subcell_index < numSubcellsOnCellHost(cell));
#endif
  const int index = _cell_to_subcells_adj_host(cell)+local_subcell_index;
  return _cell_to_subcells_host(index);
}

PHALANX_HIP_HACK_KOKKOS_FUNCTION
int
SubcellConnectivity::
cellForSubcell(const int subcell, const int local_cell_index) const
{
#ifdef PANZER_DEBUG
  KOKKOS_ASSERT(subcell >= 0 and subcell < numSubcells());
  KOKKOS_ASSERT(local_cell_index < numCellsOnSubcell(subcell));
#endif
  const int index = _subcell_to_cells_adj(subcell)+local_cell_index;
  return _subcell_to_cells(index);
}

int
SubcellConnectivity::
cellForSubcellHost(const int subcell, const int local_cell_index) const
{
#ifdef PANZER_DEBUG
  KOKKOS_ASSERT(subcell >= 0 and subcell < numSubcells());
  KOKKOS_ASSERT(local_cell_index < numCellsOnSubcellHost(subcell));
#endif
  const int index = _subcell_to_cells_adj_host(subcell)+local_cell_index;
  return _subcell_to_cells_host(index);
}

PHALANX_HIP_HACK_KOKKOS_FUNCTION
int
SubcellConnectivity::
localSubcellForSubcell(const int subcell, const int local_cell_index) const
{
#ifdef PANZER_DEBUG
  KOKKOS_ASSERT(subcell >= 0 and subcell < numSubcells());
  KOKKOS_ASSERT(local_cell_index < numCellsOnSubcell(subcell));
#endif
  const int index = _subcell_to_cells_adj(subcell)+local_cell_index;
  return _subcell_to_local_subcells(index);
}

int
SubcellConnectivity::
localSubcellForSubcellHost(const int subcell, const int local_cell_index) const
{
#ifdef PANZER_DEBUG
  KOKKOS_ASSERT(subcell >= 0 and subcell < numSubcells());
  KOKKOS_ASSERT(local_cell_index < numCellsOnSubcellHost(subcell));
#endif
  const int index = _subcell_to_cells_adj_host(subcell)+local_cell_index;
  return _subcell_to_local_subcells_host(index);
}

} // namespace panzer

#endif
