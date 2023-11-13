// Copyright(C) 2021, 2022, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once

#include <vector>

#include "Cell.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_Region.h"
#include "UnitCell.h"
#include "ZE_SystemInterface.h"

//! \file

//! `grid` stores the data for the tessellation of size `IxJ`
//!  * gridI -- extent in I direction
//!  * gridJ -- extend in J direction
//!  * grid(i,j) -- return database information at i,j location
//!  * grid(index) -- return database information at i,j location corresponding to `index`
//!  * max number of element blocks in a unit_cell mesh
//!  * range of x, y size of unit_cell mesh. (assumes unit_cell minx, miny == 0.0)
//!  * Ti, Tj, Tk -- size of regular mesh on boundary of unit_cell (? may need the 3 `Tk` values for
//!  padding)
//!  * std::vector<Cell> m_grid -- contains database information...
class Grid
{
public:
  //! Create an empty grid of size `extent_i` x `extent_j`.  The output mesh will
  //! be written to the exodus database in the Ioss::Region `region`
  explicit Grid(SystemInterface &interFace);

  void set_extent(size_t extent_i, size_t extent_j, size_t /* unused */);

  UnitCellMap &unit_cells() { return m_unitCells; }

  //! Return a reference to the Cell cell at location `(i,j)`.
  //! Does not check that `i` and `j` are in bounds.
  Cell &get_cell(size_t i, size_t j)
  {
    size_t idx = i * m_gridJ + j;
    return m_grid[idx];
  }

  //! Return `I` extent of the grid / lattice
  size_t II() const { return m_gridI; }
  //! Return `J` extent of the grid / lattice
  size_t JJ() const { return m_gridJ; }
  //! Return total number of cells in the grid / lattice
  size_t size() const { return m_gridI * m_gridJ; }
  int    parallel_size() const { return m_parallelSize; }

  //! Set the rank that the should start processing and outputting (used for subcycling)
  void set_start_rank(int start_rank) { m_startRank = start_rank; }

  //! Determine if can keep all files open at all times,
  //! or if we need to close some/all after access...
  void handle_file_count();

  //! Are nodes at the boundaries of the unit cells equivalenced.
  bool equivalence_nodes() const { return m_equivalenceNodes; }

  //! Create a Cell object referencing the UnitCell specified by `key` at location `(i,j)`
  bool initialize(size_t i, size_t j, const std::string &key);

  void add_unit_cell(const std::string &key, const std::string &unit_filename, bool ints32bit);

  //! Specify the X and Y location of each grid cell in the overall grid space.
  void set_coordinate_offsets();

  //!
  void generate_sidesets();

  //! Once all Cell objects have been initialized, Determine the coordinate extents and
  //! offsets of each cell, the size of the output mesh, the node and element id offsets
  //! for each cell, the number of nodes and elements in the output mesh and initialize
  //! the output mesh.
  template <typename INT> void process(SystemInterface &interFace, INT /* dummy */);

  void decompose(const std::string &method);

  const Ioss::ParallelUtils &util() const { return m_pu; }

  Ioss::Region *output_region(int rank = 0) { return m_outputRegions[rank].get(); }

  bool         minimize_open_files(Minimize type) { return (int)m_minimizeOpenFiles & int(type); }
  unsigned int get_generated_sidesets() { return m_generatedSideSets; }

  void set_sideset_names(const std::string &names);

  std::array<std::string, 6> generated_surface_names{
      {"min_i", "max_i", "min_j", "max_j", "min_k", "max_k"}};

private:
  //! Output node coordinates and element block connectivities for the output mesh.
  template <typename INT> void output_model(INT /*dummy*/);
  void                         internal_process();

  void create_output_regions(SystemInterface &interFace);
  void categorize_processor_boundaries();

  void output_nodal_coordinates(const Cell &cell);
  template <typename INT>
  void output_block_connectivity(Cell &cell, const std::vector<INT> &node_map);
  template <typename INT>
  void output_nodal_communication_map(Cell &cell, const std::vector<INT> &node_map);
  template <typename INT> void output_element_map(Cell &cell, INT /*dummy*/);
  template <typename INT> void output_node_map(const Cell &cell, INT /*dummy*/);

  template <typename INT> void output_surfaces(Cell &cell, INT /*dummy*/);
  template <typename INT> void output_generated_surfaces(Cell &cell, INT /*dummy*/);

  UnitCellMap                                m_unitCells;
  std::vector<std::unique_ptr<Ioss::Region>> m_outputRegions;
  std::vector<Cell>                          m_grid{};
  Ioss::ParallelUtils                        m_pu{};
  size_t                                     m_gridI{0};
  size_t                                     m_gridJ{0};
  vector3d                                   m_offset{0.0, 0.0, 0.0};
  double                                     m_scaleFactor{1.0};
  int                                        m_parallelSize{1}; //! Number of ranks to decompose for
  int          m_rankCount{1}; //! Number of ranks to process at a time.
  int          m_startRank{0}; //! Which rank to start outputting...
  bool         m_equivalenceNodes{true};
  bool         m_useInternalSidesets{true};
  bool         m_subCycle{false};
  Minimize     m_minimizeOpenFiles{Minimize::NONE}; // 1: Unit, 2: output, 3: all
  unsigned int m_generatedSideSets{0};
};
