// Copyright(C) 2021, 2022, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <cstdlib>
#include <numeric>

#include "Cell.h"
#include "Decompose.h"
#include "Grid.h"
#include "ZE_SystemInterface.h"
#include "ZE_Version.h"
#include <Ioss_CommSet.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_IOFactory.h>
#include <Ioss_NodeBlock.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_Region.h>
#include <Ioss_SideBlock.h>
#include <Ioss_SideSet.h>
#include <Ioss_SmartAssert.h>

#include <exodusII.h>
#include <fmt/chrono.h>
#include <fmt/color.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <open_file_limit.h>
#include <time_stamp.h>
#include <tokenize.h>

//! \file

extern unsigned int debug_level;

#if __cplusplus == 201103L
// From: https://stackoverflow.com/questions/17902405/how-to-implement-make-unique-function-in-c11
namespace std {
  template <class T> struct _Unique_if
  {
    typedef unique_ptr<T> _Single_object;
  };

  template <class T> struct _Unique_if<T[]>
  {
    typedef unique_ptr<T[]> _Unknown_bound;
  };

  template <class T, size_t N> struct _Unique_if<T[N]>
  {
    typedef void _Known_bound;
  };

  template <class T, class... Args>
  typename _Unique_if<T>::_Single_object make_unique(Args &&...args)
  {
    return unique_ptr<T>(new T(std::forward<Args>(args)...));
  }

  template <class T> typename _Unique_if<T>::_Unknown_bound make_unique(size_t n)
  {
    typedef typename remove_extent<T>::type U;
    return unique_ptr<T>(new U[n]());
  }

  template <class T, class... Args>
  typename _Unique_if<T>::_Known_bound make_unique(Args &&...) = delete;
} // namespace std
#endif

namespace {
  std::string tsFormat = "[{:%H:%M:%S}]";
  int         axis_index(const std::string &axis_str)
  {
    char axis = axis_str[0];
    if (axis == 'x' || axis == 'i') {
      return 0;
    }
    if (axis == 'X' || axis == 'I') {
      return 1;
    }
    if (axis == 'y' || axis == 'j') {
      return 2;
    }
    if (axis == 'Y' || axis == 'J') {
      return 3;
    }
    if (axis == 'z' || axis == 'k') {
      return 4;
    }
    if (axis == 'Z' || axis == 'K') {
      return 5;
    }
    return -1;
  }

  unsigned int which_sidesets(const std::string &sideset_string)
  {
    unsigned generate = 0;
    for (const auto &c : sideset_string) {
      if (c == 'x' || c == 'i') {
        generate |= (unsigned)Flg::MIN_I;
      }
      if (c == 'y' || c == 'j') {
        generate |= (unsigned)Flg::MIN_J;
      }
      if (c == 'z' || c == 'k') {
        generate |= (unsigned)Flg::MIN_K;
      }
      if (c == 'X' || c == 'I') {
        generate |= (unsigned)Flg::MAX_I;
      }
      if (c == 'Y' || c == 'J') {
        generate |= (unsigned)Flg::MAX_J;
      }
      if (c == 'Z' || c == 'K') {
        generate |= (unsigned)Flg::MAX_K;
      }
    }
    return generate;
  }

  Ioss::PropertyManager parse_properties(SystemInterface &interFace, int int_size);

  std::shared_ptr<Ioss::Region> create_input_region(const std::string &key, std::string filename,
                                                    bool ints_32_bits);

  template <typename INT>
  std::vector<INT> generate_node_map(Grid &grid, const Cell &cell, Mode mode, INT /*dummy*/);

  void   output_summary(Ioss::Region *region, std::ostream &strm);
  size_t handle_elements(Grid &grid, int start_rank, int rank_count);
  size_t handle_nodes(Grid &grid, int start_rank, int rank_count);
  void   handle_communications(Grid &grid, int start_rank, int rank_count);
  void   handle_surfaces(Grid &grid, int start_rank, int rank_count);
  void   generate_surfaces(Grid &grid, int start_rank, int rank_count);
} // namespace

Grid::Grid(SystemInterface &interFace)
    : m_offset(interFace.offset()), m_scaleFactor(interFace.scale_factor()),
      m_parallelSize(interFace.ranks()), m_rankCount(interFace.rank_count()),
      m_startRank(interFace.start_rank()), m_equivalenceNodes(interFace.equivalence_nodes()),
      m_useInternalSidesets(!interFace.ignore_internal_sidesets()),
      m_subCycle(interFace.subcycle()), m_minimizeOpenFiles(interFace.minimize_open_files()),
      m_generatedSideSets(which_sidesets(interFace.sideset_surfaces()))
{
  set_sideset_names(interFace.sideset_names());
}

void Grid::set_extent(size_t extent_i, size_t extent_j, size_t /* unused */)
{
  m_gridI = extent_i;
  m_gridJ = extent_j;
  m_grid.resize(m_gridI * m_gridJ);
}

bool Grid::initialize(size_t i, size_t j, const std::string &key)
{
  if (unit_cells().find(key) == unit_cells().end()) {
    return false;
  }
  const auto &unit_cell = unit_cells()[key];
  SMART_ASSERT(unit_cell->m_region != nullptr)(i)(j)(key);

  auto &cell = get_cell(i, j);
  cell.initialize(i, j, unit_cell);
  return true;
}

void Grid::add_unit_cell(const std::string &key, const std::string &unit_filename, bool ints32bit)
{
  static size_t open_files = open_file_limit();
  if (!minimize_open_files(Minimize::UNIT) && unit_cells().size() >= open_files) {
    // Just hit the limit...  Close all previous unit_cell files and set the minimize_open_files
    // behavior to UNIT.
    for (const auto &unit_cell : m_unitCells) {
      unit_cell.second->m_region->get_database()->closeDatabase();
    }
    fmt::print(stderr, fmt::fg(fmt::color::yellow),
               "\nWARNING: Number of unit cells exceeds open file limit. Setting "
               "'mimimize_open_file' mode to UNIT.\n\n");
    m_minimizeOpenFiles = Minimize((unsigned)m_minimizeOpenFiles | (unsigned)Minimize::UNIT);
  }

  if (unit_cells().find(key) != unit_cells().end()) {
    fmt::print(stderr, fmt::fg(fmt::color::red),
               "\nERROR: There is a duplicate `unit cell` ({}) in the lattice dictionary.\n\n",
               key);
    exit(EXIT_FAILURE);
  }

  auto region = create_input_region(key, unit_filename, ints32bit);

  if (region == nullptr) {
    fmt::print(stderr, fmt::fg(fmt::color::red),
               "\nERROR: Unable to open the database '{}' associated with the unit cell '{}'.\n\n",
               unit_filename, key);
    exit(EXIT_FAILURE);
  }

  unit_cells().emplace(key, std::make_shared<UnitCell>(region));
  if (minimize_open_files(Minimize::UNIT)) {
    unit_cells()[key]->m_region->get_database()->closeDatabase();
  }

  if (debug_level & 2) {
    m_pu.progress(fmt::format("\tCreated Unit Cell {}", key));
  }
}

void Grid::set_coordinate_offsets()
{
  // All unit_cells have same X, Y (and Z) extent.  Only need X and Y
  auto x_range = get_cell(0, 0).get_coordinate_range(Axis::X);
  auto y_range = get_cell(0, 0).get_coordinate_range(Axis::Y);
  auto delta_x = x_range.second - x_range.first;
  auto delta_y = y_range.second - y_range.first;

  for (size_t j = 0; j < JJ(); j++) {
    for (size_t i = 0; i < II(); i++) {
      auto &cell  = get_cell(i, j);
      cell.m_offX = delta_x * (double)i;
      cell.m_offY = delta_y * (double)j;

      if (debug_level & 2) {
        util().progress(
            fmt::format("\tCell({}, {}) X = {}, Y = {}", i, j, cell.m_offX, cell.m_offY));
      }
    }
  }
}

void Grid::create_output_regions(SystemInterface &interFace)
{
  if (debug_level & 2) {
    util().progress(__func__);
  }
  m_outputRegions.resize(interFace.ranks());

  int                   int_size   = interFace.ints32bit() ? 4 : 8;
  Ioss::PropertyManager properties = parse_properties(interFace, int_size);
  if (parallel_size() == 1) {
    properties.add(Ioss::Property("OMIT_EXODUS_NUM_MAPS", 1));
  }
  // Disable this for now.  Readers need to be modified and propagated to allow this.
  //   properties.add(Ioss::Property("MINIMAL_NEMESIS_DATA", 1));

  if (debug_level & 2) {
    properties.add(Ioss::Property("ENABLE_TRACING", 1));
  }

  for (int i = m_startRank; i < m_startRank + m_rankCount; i++) {
    // Define the output database(s)...
    if (m_parallelSize > 1) {
      properties.add(Ioss::Property("processor_count", parallel_size()));
      properties.add(Ioss::Property("my_processor", i));
    }
    Ioss::DatabaseIO *dbo =
        Ioss::IOFactory::create("exodus", interFace.outputName_, Ioss::WRITE_RESTART,
                                Ioss::ParallelUtils::comm_self(), properties);
    if (dbo == nullptr || !dbo->ok(true)) {
      std::exit(EXIT_FAILURE);
    }
    m_outputRegions[i] = std::make_unique<Ioss::Region>(dbo, "zellij_output_region");
    output_region(i)->begin_mode(Ioss::STATE_DEFINE_MODEL);
    output_region(i)->property_add(Ioss::Property("code_name", qainfo[0]));
    output_region(i)->property_add(Ioss::Property("code_version", qainfo[2]));

    if (interFace.minimize_open_files() == Minimize::OUTPUT ||
        interFace.minimize_open_files() == Minimize::ALL) {
      output_region(i)->get_database()->closeDatabase();
    }
  }
}

void Grid::set_sideset_names(const std::string &names)
{
  // 'names' is a string of the form "Axis:name,Axis:name, ..."
  // For example "x:left,y:top,X:right,Y:bottom"
  // Parse the list and only update the names specified in the list...

  if (names.empty()) {
    return;
  }

  auto tokens = Ioss::tokenize(names, ",");
  for (auto &token : tokens) {
    auto axis  = token.substr(0, 1);
    bool valid = axis.find_first_not_of("ijkIJKxyzXYZ") == std::string::npos;
    if (!valid) {
      fmt::print(stderr, fmt::fg(fmt::color::red),
                 "\nERROR: Invalid axis '{}' specified for sideset name.  Valid is one of "
                 "'ijkIJKxyzXYZ'.\n\n",
                 axis);
      exit(EXIT_FAILURE);
    }

    // Now get the name...
    auto ss_name = token.substr(2);

    // Update the name in the list of generated sideset names...
    auto index = axis_index(axis);
    SMART_ASSERT(index >= 0)(axis)(index);
    generated_surface_names[index] = std::move(ss_name);
  }
}

void Grid::decompose(const std::string &method)
{
  if (debug_level & 2) {
    util().progress(__func__);
  }
  decompose_grid(*this, m_parallelSize, method);

  categorize_processor_boundaries();
}

void Grid::categorize_processor_boundaries()
{
  if (debug_level & 2) {
    util().progress(__func__);
  }

  // Now iterate the cells and tell each cell the rank of all neighboring cells...
  // boundary with its "left" or "lower" neighboring cell.

  for (size_t j = 0; j < JJ(); j++) {
    for (size_t i = 0; i < II(); i++) {
      auto &cell = get_cell(i, j);
      if (i > 0) {
        const auto &left = get_cell(i - 1, j);
        cell.set_rank(Loc::L, left.rank(Loc::C));
        if (j > 0) {
          const auto &BL = get_cell(i - 1, j - 1);
          cell.set_rank(Loc::BL, BL.rank(Loc::C));
        }
        if (j < JJ() - 1) {
          const auto &TL = get_cell(i - 1, j + 1);
          cell.set_rank(Loc::TL, TL.rank(Loc::C));
        }
      }
      if (i < II() - 1) {
        const auto &right = get_cell(i + 1, j);
        cell.set_rank(Loc::R, right.rank(Loc::C));
        if (j > 0) {
          const auto &BR = get_cell(i + 1, j - 1);
          cell.set_rank(Loc::BR, BR.rank(Loc::C));
        }
        if (j < JJ() - 1) {
          const auto &TR = get_cell(i + 1, j + 1);
          cell.set_rank(Loc::TR, TR.rank(Loc::C));
        }
      }
      if (j > 0) {
        const auto &B = get_cell(i, j - 1);
        cell.set_rank(Loc::B, B.rank(Loc::C));
      }
      if (j < JJ() - 1) {
        const auto &T = get_cell(i, j + 1);
        cell.set_rank(Loc::T, T.rank(Loc::C));
      }
    }
  }

  if (debug_level & 32) {
    auto width = Ioss::Utils::number_width(parallel_size());
    for (size_t j = 0; j < JJ(); j++) {
      for (size_t i = 0; i < II(); i++) {
        const auto &cell  = get_cell(i, j);
        auto        left  = cell.processor_boundary(Loc::L) ? '<' : ' ';
        auto        below = cell.processor_boundary(Loc::B) ? '^' : ' ';
        fmt::print(" {0}{1:{3}}{2}", left, cell.rank(Loc::C), below, width);
      }
      fmt::print("\n");
    }
  }
}

void Grid::generate_sidesets()
{
  if (m_generatedSideSets != 0) {
    for (const auto &unit_cell : m_unitCells) {
      unit_cell.second->generate_boundary_faces(m_generatedSideSets);
    }
  }
}

void Grid::internal_process()
{
  if (debug_level & 2) {
    util().progress(__func__);
  }
  auto node_count    = handle_nodes(*this, m_startRank, m_rankCount);
  auto element_count = handle_elements(*this, m_startRank, m_rankCount);
  handle_communications(*this, m_startRank, m_rankCount);
  if (m_useInternalSidesets) {
    handle_surfaces(*this, m_startRank, m_rankCount);
  }
  generate_surfaces(*this, m_startRank, m_rankCount);

  for (int i = m_startRank; i < m_startRank + m_rankCount; i++) {
    output_region(i)->end_mode(Ioss::STATE_DEFINE_MODEL);
    if (debug_level & 64) {
      output_summary(output_region(i), std::cerr);
    }
  }
  if (util().parallel_rank() == 0) {
    fmt::print("                {} Nodes; {} Elements.\n", fmt::group_digits(node_count),
               fmt::group_digits(element_count));
  }
}

template void Grid::process(SystemInterface &, int64_t);
template void Grid::process(SystemInterface &, int);

template <typename INT> void Grid::process(SystemInterface &interFace, INT /* dummy */)
{
  bool subcycle   = m_subCycle;
  int  start_rank = m_startRank;
  int  rank_count = m_rankCount;
  int  end_rank   = subcycle ? m_parallelSize : m_startRank + m_rankCount;

  if (end_rank > m_parallelSize) {
    end_rank    = m_parallelSize;
    m_rankCount = m_parallelSize - m_startRank;
  }

  for (int begin = start_rank; begin < end_rank; begin += rank_count) {
    m_startRank = begin;
    if (m_startRank + m_rankCount > m_parallelSize) {
      m_rankCount = m_parallelSize - m_startRank;
    }
    if (debug_level & 2) {
      fmt::print(stderr, "{} Processing Ranks {} to {}\n", time_stamp(tsFormat), begin,
                 begin + rank_count - 1);
    }

    create_output_regions(interFace);

    internal_process();
    if (debug_level & 2) {
      fmt::print(stderr, "{} Lattice Processing Finalized\n", time_stamp(tsFormat));
    }

    output_model(INT(0));
    if (debug_level & 2) {
      fmt::print(stderr, "{} Model Output\n", time_stamp(tsFormat));
    }
  }
}

template void Grid::output_model(int64_t);
template void Grid::output_model(int);

// Write the output database(s) for all ranks...
template <typename INT> void Grid::output_model(INT /*dummy*/)
{
  if (debug_level & 2) {
    util().progress(__func__);
  }
  // Coordinates do not depend on order of operations, so can do these
  // a rank at a time which can be more efficient once we exceed
  // maximum number of open files...  Also easier to parallelize...
  for (int r = m_startRank; r < m_startRank + m_rankCount; r++) {
    for (size_t j = 0; j < JJ(); j++) {
      for (size_t i = 0; i < II(); i++) {
        const auto &cell = get_cell(i, j);
        if (cell.rank(Loc::C) == r) {
          output_nodal_coordinates(cell);
        }
      }
    }
    if (minimize_open_files(Minimize::OUTPUT)) {
      output_region(r)->get_database()->closeDatabase();
    }
  }
  if (debug_level & 2) {
    util().progress("\tEnd Nodal Coordinate Output");
  }

  // Surfaces do not depend on order of operations, so can do these
  // a rank at a time which can be more efficient once we exceed
  // maximum number of open files...  Also easier to parallelize...
  for (int r = m_startRank; r < m_startRank + m_rankCount; r++) {
    for (size_t j = 0; j < JJ(); j++) {
      for (size_t i = 0; i < II(); i++) {
        auto &cell = get_cell(i, j);
        if (cell.rank(Loc::C) == r) {
          if (m_useInternalSidesets) {
            output_surfaces(cell, INT(0));
          }
          output_generated_surfaces(cell, INT(0));
        }
      }
    }
    if (minimize_open_files(Minimize::OUTPUT)) {
      output_region(r)->get_database()->closeDatabase();
    }
  }
  if (debug_level & 2) {
    util().progress("\tEnd Surface Output");
  }

  // All the rest of these depend on progressing through the cells in
  // correct order so that can pass information correctly from cell to
  // cell.  Need to figure out how to eliminate this ordering so can
  // parallelize and not worry about open file count.

  for (size_t j = 0; j < JJ(); j++) {
    for (size_t i = 0; i < II(); i++) {
      auto &cell     = get_cell(i, j);
      auto  node_map = generate_node_map(*this, cell, Mode::PROCESSOR, INT(0));
      output_block_connectivity(cell, node_map);
      if (parallel_size() > 1) {
        output_nodal_communication_map(cell, node_map);
      }
    }
  }
  if (debug_level & 2) {
    util().progress("\tEnd Nodal Communication Map Output");
  }

  if (parallel_size() > 1) {
    for (size_t j = 0; j < JJ(); j++) {
      for (size_t i = 0; i < II(); i++) {
        auto &cell = get_cell(i, j);
        output_node_map(cell, INT(0));
        output_element_map(cell, INT(0));
      }
    }
    if (debug_level & 2) {
      util().progress("\tEnd Node/Element Map Output");
    }
  }
}

void Grid::output_nodal_coordinates(const Cell &cell)
{
  int rank = cell.rank(Loc::C);

  auto               *nb = cell.region()->get_node_blocks()[0];
  std::vector<double> coord_x;
  std::vector<double> coord_y;
  std::vector<double> coord_z;

  // Are we modifying the coordinates ... scale and/or offset and/or offset_unit_cell...
  bool mod_x = cell.m_offX != 0.0 || m_scaleFactor != 1.0 || m_offset[0] != 0.0;
  bool mod_y = cell.m_offY != 0.0 || m_scaleFactor != 1.0 || m_offset[1] != 0.0;
  bool mod_z = m_scaleFactor != 1.0 || m_offset[2] != 0.0;

  double scale = m_scaleFactor;
  nb->get_field_data("mesh_model_coordinates_x", coord_x);
  if (mod_x) {
    double offset = m_offset[0];
    std::for_each(coord_x.begin(), coord_x.end(),
                  [&cell, scale, offset](double &d) { d = (d + cell.m_offX) * scale + offset; });
  }

  nb->get_field_data("mesh_model_coordinates_y", coord_y);
  if (mod_y) {
    double offset = m_offset[1];
    std::for_each(coord_y.begin(), coord_y.end(),
                  [&cell, scale, offset](double &d) { d = (d + cell.m_offY) * scale + offset; });
  }

  nb->get_field_data("mesh_model_coordinates_z", coord_z);
  if (mod_z) {
    double offset = m_offset[2];
    std::for_each(coord_z.begin(), coord_z.end(),
                  [scale, offset](double &d) { d = d * scale + offset; });
  }

  // Filter coordinates down to only "new nodes"...
  if (m_equivalenceNodes && (cell.has_neighbor_i() || cell.has_neighbor_j())) {
    auto   mode              = parallel_size() > 1 ? Mode::PROCESSOR : Mode::GLOBAL;
    auto   categorized_nodes = cell.categorize_nodes(mode);
    size_t nn                = 0;
    for (size_t n = 0; n < categorized_nodes.size(); n++) {
      if (categorized_nodes[n] == 0) {
        coord_x[nn] = coord_x[n];
        coord_y[nn] = coord_y[n];
        coord_z[nn] = coord_z[n];
        nn++;
      }
    }
  }

  int  exoid = output_region(rank)->get_database()->get_file_pointer();
  auto start = cell.m_localNodeIdOffset + 1;
  auto count = cell.added_node_count(Mode::PROCESSOR, m_equivalenceNodes);
  ex_put_partial_coord(exoid, start, count, Data(coord_x), Data(coord_y), Data(coord_z));

  if (minimize_open_files(Minimize::UNIT)) {
    cell.region()->get_database()->closeDatabase();
  }
}

template <typename INT> void Grid::output_generated_surfaces(Cell &cell, INT /*dummy*/)
{
  auto generated_ijk = get_generated_sidesets();
  if (generated_ijk == 0) {
    return;
  }

  int rank = cell.rank(Loc::C);

  const std::array<int, 6> boundary_rank{
      cell.rank(Loc::L), cell.rank(Loc::R), cell.rank(Loc::B), cell.rank(Loc::T), -1, -1};
  std::array<enum Flg, 6> boundary_flag{Flg::MIN_I, Flg::MAX_I, Flg::MIN_J,
                                        Flg::MAX_J, Flg::MIN_K, Flg::MAX_K};

  int exoid = output_region(rank)->get_database()->get_file_pointer();

  for (int face = 0; face < 6; face++) {
    if ((generated_ijk & (unsigned)boundary_flag[face]) && boundary_rank[face] == -1) {
      // Find surface on output mesh...
      auto *osurf = output_region(rank)->get_sideset(generated_surface_names[face]);
      SMART_ASSERT(osurf != nullptr);
      const auto &oblocks = osurf->get_side_blocks();
      SMART_ASSERT(oblocks.size() == 1)(oblocks.size());

      auto            &boundary = cell.unit()->boundary_blocks[face];
      auto             count    = boundary.size();
      std::vector<INT> elements;
      std::vector<INT> faces;
      elements.reserve(count);
      faces.reserve(count);

      for (const auto &block_faces : boundary.m_faces) {
        const auto &block_name = block_faces.first;
        const auto &bnd_faces  = block_faces.second;

        // This is the offset within this element block -- i.e., the 'element_offsetth' element in
        // this block.
        size_t element_offset = cell.m_localElementIdOffset[block_name];
        // This is the offset of the elements in this element block within all of the elements in
        // the output file.
        size_t global_offset = output_region(rank)->get_element_block(block_name)->get_offset();

        for (const auto &bface : bnd_faces) {
          elements.push_back(bface / 10 + element_offset + global_offset);
          faces.push_back(bface % 10 + 1);
        }
      }

      auto id    = osurf->get_property("id").get_int();
      auto start = cell.m_localSurfaceOffset[generated_surface_names[face]];
      ex_put_partial_set(exoid, EX_SIDE_SET, id, start + 1, count, Data(elements), Data(faces));
    }
  }
}

template <typename INT> void Grid::output_surfaces(Cell &cell, INT /*dummy*/)
{
  int rank  = cell.rank(Loc::C);
  int exoid = output_region(rank)->get_database()->get_file_pointer();

  // Get the surfaces on this cell...
  const auto &surfaces = cell.region()->get_sidesets();
  for (const auto *surface : surfaces) {

    // Find corresponding surface on output mesh...
    auto *osurf = output_region(rank)->get_sideset(surface->name());
    SMART_ASSERT(osurf != nullptr);
    const auto &oblocks = osurf->get_side_blocks();

    std::vector<INT> elements;
    std::vector<INT> faces;
    elements.reserve(oblocks[0]->entity_count());
    faces.reserve(oblocks[0]->entity_count());

    const auto &blocks = surface->get_side_blocks();
    for (const auto *block : blocks) {
      // Get the element/face pairs for the SideBlock in this surface...
      std::vector<INT> element_side;
      block->get_field_data("element_side_raw", element_side);

      // Now separate the element and the face into different vectors and
      // update the element number to match the elements in this cell...
      const auto *parent = block->parent_element_block();
      SMART_ASSERT(parent != nullptr);

      // This is the offset of the elements in the unit cell parent element block
      size_t unit_block_offset = parent->get_offset();

      // This is the offset within this element block -- i.e., the 'element_offsetth' element in
      // this block.
      size_t element_offset = cell.m_localElementIdOffset[parent->name()];
      // This is the offset of the elements in this element block within all of the elements in the
      // output file.
      size_t global_offset = output_region(rank)->get_element_block(parent->name())->get_offset();

      size_t entity_count = block->entity_count();
      for (size_t i = 0; i < entity_count; i++) {
        auto element        = element_side[2 * i + 0];
        auto output_element = element - unit_block_offset + element_offset + global_offset;
        elements.push_back(output_element);
        faces.push_back(element_side[2 * i + 1]);
      }
    }

    auto id    = osurf->get_property("id").get_int();
    auto start = cell.m_localSurfaceOffset[osurf->name()];
    auto count = elements.size();
    ex_put_partial_set(exoid, EX_SIDE_SET, id, start + 1, count, Data(elements), Data(faces));
  }

  if (minimize_open_files(Minimize::UNIT)) {
    cell.region()->get_database()->closeDatabase();
  }
}

template <typename INT>
void Grid::output_block_connectivity(Cell &cell, const std::vector<INT> &node_map)
{
  int rank = cell.rank(Loc::C);
  if (rank >= m_startRank && rank < m_startRank + m_rankCount) {
    int exoid = output_region(rank)->get_database()->get_file_pointer();

    const auto      &blocks = cell.region()->get_element_blocks();
    std::vector<INT> connect;
    for (const auto *block : blocks) {
      block->get_field_data("connectivity_raw", connect);
      for (size_t k = 0; k < connect.size(); k++) {
        connect[k] = node_map[connect[k]];
      }
      auto start = cell.m_localElementIdOffset[block->name()] + 1;
      auto count = block->entity_count();
      auto id    = block->get_property("id").get_int();
      if (debug_level & 8) {
        fmt::print(stderr, "Rank: {}, Cell({}, {}), Block {}, id {}, start {}, count {}\n", rank,
                   cell.m_i, cell.m_j, block->name(), id, start, count);
      }
      ex_put_partial_conn(exoid, EX_ELEM_BLOCK, id, start, count, Data(connect), nullptr, nullptr);
    }

    if (debug_level & 2) {
      util().progress(fmt::format("Generated Node Map / Output Connectivity for Cell({}, {})",
                                  cell.m_i, cell.m_j));
    }
    if (minimize_open_files(Minimize::UNIT)) {
      cell.region()->get_database()->closeDatabase();
    }
    if (minimize_open_files(Minimize::OUTPUT)) {
      output_region(rank)->get_database()->closeDatabase();
    }
  }
}

template <typename INT>
void Grid::output_nodal_communication_map(Cell &cell, const std::vector<INT> &node_map)
{
  // The `node_map` has the processor-local node ids for all nodes on this cell.
  // Note that the `node_map` starts at index 1 in the vector...
  // Need to check the boundaries of the cell and determine which boundaries are on a different rank
  // and output the nodes and processors to the communication map.
  int rank = cell.rank(Loc::C);
  if (rank >= m_startRank && rank < m_startRank + m_rankCount) {
    std::vector<INT> nodes;
    std::vector<INT> procs;
    cell.populate_node_communication_map(node_map, nodes, procs);

    int exoid = output_region(rank)->get_database()->get_file_pointer();

    auto start = cell.m_communicationNodeOffset + 1;
    auto count = cell.m_communicationNodeCount;

    ex_put_partial_node_cmap(exoid, 1, start, count, Data(nodes), Data(procs), rank);

    if (minimize_open_files(Minimize::OUTPUT)) {
      output_region(rank)->get_database()->closeDatabase();
    }

    if (debug_level & 32) {
      fmt::print(stderr, "Rank: {}, Cell({}, {}), Node Comm Map: start {}, count {}\n", rank,
                 cell.m_i, cell.m_j, start, count);
    }

    if (debug_level & 2) {
      util().progress(
          fmt::format("Output Nodal Communication Map for Cell({}, {})", cell.m_i, cell.m_j));
    }
  }
}

template <typename INT> void Grid::output_node_map(const Cell &cell, INT /*dummy*/)
{
  int rank = cell.rank(Loc::C);

  auto start = cell.m_localNodeIdOffset + 1;
  auto count = cell.added_node_count(Mode::PROCESSOR, m_equivalenceNodes);

  if (parallel_size() == 1) {
    auto             gid = cell.m_globalNodeIdOffset + 1;
    std::vector<INT> map(count);
    std::iota(map.begin(), map.end(), gid);
    int exoid = output_region(rank)->get_database()->get_file_pointer();
    ex_put_partial_id_map(exoid, EX_NODE_MAP, start, count, Data(map));
  }
  else {
    auto map = generate_node_map(*this, cell, Mode::GLOBAL, INT(0));

    if (rank >= m_startRank && rank < m_startRank + m_rankCount) {

      // Filter nodes down to only "new nodes"...
      if (m_equivalenceNodes && (cell.has_neighbor_i() || cell.has_neighbor_j())) {
        auto   mode              = Mode::PROCESSOR;
        auto   categorized_nodes = cell.categorize_nodes(mode);
        size_t nn                = 0;
        for (size_t n = 0; n < categorized_nodes.size(); n++) {
          if (categorized_nodes[n] == 0) {
            map[nn + 1] = map[n + 1];
            nn++;
          }
        }
      }
      if (debug_level & 8) {
        fmt::print("Cell({}, {}), start {}, count {}\n", cell.m_i, cell.m_j, start, count);
      }
      int exoid = output_region(rank)->get_database()->get_file_pointer();
      ex_put_partial_id_map(exoid, EX_NODE_MAP, start, count, &map[1]);
      if (minimize_open_files(Minimize::OUTPUT)) {
        output_region(rank)->get_database()->closeDatabase();
      }
    }
  }

  if (debug_level & 2) {
    util().progress(fmt::format("Generated Node Map for Rank {}, Cell({}, {}): start {}, count "
                                "{}\n",
                                rank, cell.m_i, cell.m_j, start, count));
  }
}

template <typename INT> void Grid::output_element_map(Cell &cell, INT /*dummy*/)
{
  int rank = cell.rank(Loc::C);
  if (rank >= m_startRank && rank < m_startRank + m_rankCount) {
    int exoid = output_region(rank)->get_database()->get_file_pointer();

    const auto &output_blocks = output_region(rank)->get_element_blocks();

    // This is the element block offset for the "single output file"
    // for the block being output For example, if the total mesh has
    // 3 blocks, with 100, 200, 100 elements, then the element block
    // offset would be 0, 100, 300 (Ioss::ElementBlock::get_offset())
    size_t global_id_offset = 0;

    for (const auto *output_element_block : output_blocks) {
      auto *block = cell.region()->get_element_block(output_element_block->name());
      if (block != nullptr) {

        auto             gid = cell.m_globalElementIdOffset[block->name()] + 1 + global_id_offset;
        auto             element_count = block->entity_count();
        std::vector<INT> map(element_count);

        std::iota(map.begin(), map.end(), gid);

        auto output_block_offset = output_element_block->get_offset();

        // This cells element block ids start this far into the portion of the map for this
        // element block
        auto local_offset = cell.m_localElementIdOffset[block->name()];

        auto start = output_block_offset + local_offset + 1;
        ex_put_partial_id_map(exoid, EX_ELEM_MAP, start, element_count, Data(map));

        if (debug_level & 8) {
          fmt::print("Rank {}: Cell({}, {}), Block {}, start {}, element_count {}, gid {}\n", rank,
                     cell.m_i, cell.m_j, block->name(), start, element_count, gid);
        }
      }
      // If we were outputting a single file, then this element
      // block in that file would have this many elements.
      auto global_block_element_count =
          output_element_block->get_property("global_entity_count").get_int();
      global_id_offset += global_block_element_count;
    }
    if (minimize_open_files(Minimize::OUTPUT)) {
      output_region(rank)->get_database()->closeDatabase();
    }
  }
}

void Grid::handle_file_count()
{
  // If the user has specified a 'minimize_open_files' behavior on the command line,
  // honor that mode here.  Note that we have already processed the unit cells and
  // if the unit cell count exceeds the open file limit, then that code has already
  // automatically set the minimize_open_files to include UNIT...
  if (m_minimizeOpenFiles == Minimize::ALL) {
    return;
  }

  size_t open_files = open_file_limit();
  if (util().parallel_rank() == 0) {
    fmt::print("\n Maximum Open File Count = {}\n", open_file_limit());
  }

  auto unit_cell_size = unit_cells().size();
  if (minimize_open_files(Minimize::UNIT)) {
    unit_cell_size = 1; // There will only be a single unit cell file open at a time.
  }

  if (unit_cell_size + m_rankCount > open_files) {
    // Too many files, need to close some of them after accessing...
    // We have already checked the unit_cell() size earlier.
    // If user has specified 'minimize_open_files == OUTPUT', then honor that...
    if (minimize_open_files(Minimize::OUTPUT)) {
      return;
    }

    // Have decision on whether to keep all unit_files open
    // and reduce open output files further, or to close
    // unit files and have more output files open...
    int output_open = open_files - (int)unit_cell_size;
    if (output_open < (int)(.2 * m_rankCount)) {
      // Close the unit files.  Prefer to have output files open instead.
      m_minimizeOpenFiles = Minimize((unsigned)m_minimizeOpenFiles | (unsigned)Minimize::UNIT);
      unit_cell_size      = 1;
    }

    // Set the rank count to the number of files left after subtracting off 'unit_cell_size'
    // and set the 'subcycle' mode.
    auto max_rank_count = open_files - unit_cell_size;
    if (max_rank_count < (size_t)m_rankCount) {
      m_rankCount = max_rank_count;
    }
    m_subCycle = true;
  }

  if (util().parallel_rank() == 0 && m_minimizeOpenFiles != Minimize::NONE) {
    std::array<std::string, 4> smode{"NONE", "UNIT", "OUTPUT", "ALL"};
    fmt::print(" Setting `minimize_open_files` mode to {}.\n", smode[(int)m_minimizeOpenFiles]);
  }
}

namespace {
  size_t handle_elements(Grid &grid, int start_rank, int rank_count)
  {
    // Not all unit cells have the same element blocks and the output
    // grid will contain the union of the element blocks on each unit
    // cell...
    //
    // While finalizing the cells, we will create this union for use in
    // the output region.  Would be quicker to do iteration on unit_cell
    // map which is smaller, but for now do it here and see if becomes
    // bottleneck.

    std::map<std::string, std::unique_ptr<Ioss::ElementBlock>> output_element_blocks;
    std::vector<std::map<std::string, size_t>> element_block_elem_count(grid.parallel_size());
    std::map<std::string, size_t>              global_element_block_elem_count;

    for (size_t j = 0; j < grid.JJ(); j++) {
      for (size_t i = 0; i < grid.II(); i++) {
        auto &cell = grid.get_cell(i, j);
        auto  rank = cell.rank(Loc::C);

        const auto &element_blocks = cell.region()->get_element_blocks();
        for (const auto *block : element_blocks) {
          const auto &blk                   = block->name();
          cell.m_globalElementIdOffset[blk] = global_element_block_elem_count[blk];
          cell.m_localElementIdOffset[blk]  = element_block_elem_count[rank][blk];

          element_block_elem_count[rank][blk] += block->entity_count();
          global_element_block_elem_count[blk] += block->entity_count();

          if (debug_level & 8) {
            fmt::print("rank, i, j, blk, loffset, goffset, block_count: {}: {} {} {} {} {} {}\n",
                       rank, i, j, blk, cell.m_localElementIdOffset[blk],
                       cell.m_globalElementIdOffset[blk], element_block_elem_count[rank][blk]);
          }

          // Create output element block if does not exist yet...
          if (output_element_blocks.find(blk) == output_element_blocks.end()) {
            output_element_blocks.emplace(
                blk, std::make_unique<Ioss::ElementBlock>(Ioss::ElementBlock(*block)));
          }
        }
      }
    }

    // Calculate values needed to set the "global_entity_count" property on the output element
    // blocks.
    size_t                               global_element_count = 0;
    std::map<const std::string, int64_t> global_block_element_count;
    for (auto &blk : output_element_blocks) {
      for (int rank = 0; rank < grid.parallel_size(); rank++) {
        auto *block = blk.second.get();
        global_block_element_count[block->name()] += element_block_elem_count[rank][block->name()];
        global_element_count += element_block_elem_count[rank][block->name()];
      }
    }

    // Define the element blocks in the output database...
    for (int rank = start_rank; rank < start_rank + rank_count; rank++) {
      for (auto &blk : output_element_blocks) {
        auto *block = new Ioss::ElementBlock(*blk.second);
        block->property_update("entity_count", element_block_elem_count[rank][block->name()]);
        block->property_update("global_entity_count", global_block_element_count[block->name()]);
        grid.output_region(rank)->property_add(
            Ioss::Property("global_element_count", (int64_t)global_element_count));
        grid.output_region(rank)->add(block);
        if (debug_level & 8) {
          fmt::print("rank, blk, element_count: {}: {} {}\n", rank, block->name(),
                     block->entity_count());
        }
      }
    }
    return global_element_count;
  }

  size_t handle_nodes(Grid &grid, int start_rank, int rank_count)
  {
    size_t              global_node_count = 0;
    std::vector<size_t> local_node_count(grid.parallel_size());

    for (size_t j = 0; j < grid.JJ(); j++) {
      for (size_t i = 0; i < grid.II(); i++) {
        auto &cell                = grid.get_cell(i, j);
        auto  rank                = cell.rank(Loc::C);
        cell.m_globalNodeIdOffset = global_node_count;
        cell.m_localNodeIdOffset  = local_node_count[rank];
        SMART_ASSERT(cell.region() != nullptr)(i)(j);

        auto new_global_nodes    = cell.added_node_count(Mode::GLOBAL, grid.equivalence_nodes());
        auto new_processor_nodes = cell.added_node_count(Mode::PROCESSOR, grid.equivalence_nodes());
        if (debug_level & 8) {
          fmt::print("rank: i, j, node_offset, added_nodes: {}: {} {} {} {} {}\n", rank, i, j,
                     local_node_count[rank], new_global_nodes, new_processor_nodes);
        }
        local_node_count[rank] += new_processor_nodes;
        global_node_count += new_global_nodes;
      }
    }
    // Define the output database node block...
    for (int i = start_rank; i < start_rank + rank_count; i++) {
      std::string block_name        = "nodeblock_1";
      int         spatial_dimension = 3;
      auto       *block = new Ioss::NodeBlock(grid.output_region(i)->get_database(), block_name,
                                              local_node_count[i], spatial_dimension);
      block->property_add(Ioss::Property("id", 1));
      grid.output_region(i)->add(block);
      grid.output_region(i)->property_add(
          Ioss::Property("global_node_count", (int64_t)global_node_count));
    }
    return global_node_count;
  }

  void handle_communications(Grid &grid, int start_rank, int rank_count)
  {
    // Need to determine the number of nodes on the processor
    // boundaries for each output database and create a CommSet...
    std::vector<size_t> bnode(grid.parallel_size());
    for (size_t j = 0; j < grid.JJ(); j++) {
      for (size_t i = 0; i < grid.II(); i++) {
        auto &cell                     = grid.get_cell(i, j);
        auto  rank                     = cell.rank(Loc::C);
        cell.m_communicationNodeOffset = bnode[rank];
        auto bnd                       = cell.processor_boundary_node_count();
        bnode[rank] += bnd;
        if (debug_level & 32) {
          fmt::print("rank: {}, Cell({}, {}): Boundary Count = {}, Total = {}\n", rank, i, j, bnd,
                     bnode[rank]);
        }
      }
    }

    // Now add communication sets to all output databases...
    for (int rank = start_rank; rank < start_rank + rank_count; rank++) {
      auto *cs = new Ioss::CommSet(grid.output_region(rank)->get_database(), "commset_node", "node",
                                   bnode[rank]);
      grid.output_region(rank)->add(cs);
    }
  }

  void handle_surfaces(Grid &grid, int start_rank, int rank_count)
  {
    // Not all unit cells have the same surfaces and the output
    // grid will contain the union of the surfaces on each unit
    // cell that has a surface on the boundary of the output mesh.
    //
    // While finalizing the cells, we will create this union for use in
    // the output region.  Would be quicker to do iteration on unit_cell
    // map which is smaller, but for now do it here and see if becomes
    // bottleneck.
    std::map<std::string, std::unique_ptr<Ioss::SideSet>> output_surfaces;
    std::map<std::string, size_t>                         global_surface_face_count;

    // Per rank...
    std::vector<std::map<std::string, size_t>> surface_face_count(grid.parallel_size());
    std::vector<std::map<std::string, size_t>> local_surface_offset(grid.parallel_size());

    for (size_t j = 0; j < grid.JJ(); j++) {
      for (size_t i = 0; i < grid.II(); i++) {
        auto &cell = grid.get_cell(i, j);
        auto  rank = cell.rank(Loc::C);

        const auto &surfaces = cell.region()->get_sidesets();
        for (const auto *surface : surfaces) {
          const auto &surf                = surface->name();
          cell.m_localSurfaceOffset[surf] = local_surface_offset[rank][surf];

          const auto &blocks = surface->get_side_blocks();
          for (const auto *blk : blocks) {
            surface_face_count[rank][surf] += blk->entity_count();
            global_surface_face_count[surf] += blk->entity_count();
            local_surface_offset[rank][surf] += blk->entity_count();
          }

          // Create output surface if does not exist yet...
          if (output_surfaces.find(surf) == output_surfaces.end()) {
            output_surfaces.emplace(surf, std::make_unique<Ioss::SideSet>(Ioss::SideSet(*surface)));
          }
        }
      }
    }

    // Define the surfaces in the output database...
    for (int rank = start_rank; rank < start_rank + rank_count; rank++) {
      for (auto &surf : output_surfaces) {
        auto *surface =
            new Ioss::SideSet(grid.output_region(rank)->get_database(), surf.second->name());
        auto *block =
            new Ioss::SideBlock(grid.output_region(rank)->get_database(), surf.second->name(),
                                "quad4", "hex8", surface_face_count[rank][surface->name()]);
        surface->add(block);
        block->property_update("global_entity_count", global_surface_face_count[surface->name()]);
        surface->property_update("global_entity_count", global_surface_face_count[surface->name()]);
        block->property_update("distribution_factor_count", 0);
        grid.output_region(rank)->add(surface);
      }
    }
  }

  void generate_surfaces(Grid &grid, int start_rank, int rank_count)
  {
    // The user may have requested that sidesets on one or more
    // of the models boundary surfaces (min i,j,k and/or max i,j,k)
    // be automatically generated (they don't exist on the input
    // unit cells).  Handle those sidesets here.

    auto generated_ijk = grid.get_generated_sidesets();
    if (generated_ijk == 0) {
      return;
    }

    std::array<size_t, 6> global_surface_face_count{};

    // Per rank...
    std::vector<std::array<size_t, 6>> surface_face_count(grid.parallel_size());
    std::vector<std::array<size_t, 6>> local_surface_offset(grid.parallel_size());

    std::array<enum Flg, 6> boundary_flag{Flg::MIN_I, Flg::MAX_I, Flg::MIN_J,
                                          Flg::MAX_J, Flg::MIN_K, Flg::MAX_K};

    for (size_t j = 0; j < grid.JJ(); j++) {
      for (size_t i = 0; i < grid.II(); i++) {
        auto                    &cell = grid.get_cell(i, j);
        auto                     rank = cell.rank(Loc::C);
        const std::array<int, 6> boundary_rank{
            cell.rank(Loc::L), cell.rank(Loc::R), cell.rank(Loc::B), cell.rank(Loc::T), -1, -1};

        for (int face = 0; face < 6; face++) {
          if ((generated_ijk & (unsigned)boundary_flag[face]) && boundary_rank[face] == -1) {
            cell.m_localSurfaceOffset[grid.generated_surface_names[face]] =
                local_surface_offset[rank][face];

            auto count = cell.unit()->boundary_blocks[face].size();
            global_surface_face_count[face] += count;
            surface_face_count[rank][face] += count;
            local_surface_offset[rank][face] += count;
          }
        }
      }
    }

    // Define the surfaces in the output database...
    for (int rank = start_rank; rank < start_rank + rank_count; rank++) {
      for (size_t i = 0; i < 6; i++) {
        if (global_surface_face_count[i] > 0) {
          auto *surface = new Ioss::SideSet(grid.output_region(rank)->get_database(),
                                            grid.generated_surface_names[i]);
          auto *block   = new Ioss::SideBlock(grid.output_region(rank)->get_database(),
                                              grid.generated_surface_names[i], "quad4", "hex8",
                                              surface_face_count[rank][i]);
          surface->property_update("global_entity_count", global_surface_face_count[i]);
          block->property_update("global_entity_count", global_surface_face_count[i]);
          block->property_update("distribution_factor_count", 0);
          surface->add(block);
          grid.output_region(rank)->add(surface);
        }
      }
    }
  }

  template <typename INT>
  std::vector<INT> generate_node_map(Grid &grid, const Cell &cell, Mode mode, INT /*dummy*/)
  {
    // Generate a "map" from nodes in the input connectivity to the
    // output connectivity in global nodes.  If no neighbors, then
    // this would just be adding `cell.m_globalNodeIdOffset` to each
    // connectivity entry.

    std::vector<INT> map = cell.generate_node_map(mode, grid.equivalence_nodes(), INT(0));
    if (debug_level & 8) {
      fmt::print("Cell({},{}) PROCESSOR MAP: {}\n", cell.m_i, cell.m_j, fmt::join(map, " "));
    }

    // Now that we have the node map for this cell, we need to save
    // the mappings for the max_I and max_J faces and max_I-max_J edge
    // for use by later neighbors...  Check whether cell has neighbors
    // on max_I or max_J faces...
    if (grid.equivalence_nodes()) {
      if ((mode == Mode::GLOBAL ||
           (mode == Mode::PROCESSOR && cell.rank(Loc::C) == cell.rank(Loc::R))) &&
          (cell.m_i + 1 < grid.II())) {
        const auto &neighbor = grid.get_cell(cell.m_i + 1, cell.m_j);
        cell.populate_neighbor(Loc::L, map, neighbor);
      }

      if ((mode == Mode::GLOBAL ||
           (mode == Mode::PROCESSOR && cell.rank(Loc::C) == cell.rank(Loc::T))) &&
          (cell.m_j + 1 < grid.JJ())) {
        const auto &neighbor = grid.get_cell(cell.m_i, cell.m_j + 1);
        cell.populate_neighbor(Loc::B, map, neighbor);
      }

      if (mode == Mode::PROCESSOR) {
        if (cell.processor_boundary(Loc::L) && (cell.rank(Loc::TL) == cell.rank(Loc::C))) {
          const auto &tl_corner = grid.get_cell(cell.m_i - 1, cell.m_j + 1);
          cell.populate_neighbor(Loc::BR, map, tl_corner);
        }
        // Now the other "corner case"
        if (cell.processor_boundary(Loc::R) && (cell.rank(Loc::TR) == cell.rank(Loc::C))) {
          const auto &tr_corner = grid.get_cell(cell.m_i + 1, cell.m_j + 1);
          cell.populate_neighbor(Loc::BL, map, tr_corner);
        }
      }
    }
    return map;
  }

  Ioss::PropertyManager parse_properties(SystemInterface &interFace, int int_size)
  {
    Ioss::PropertyManager properties;
    if (int_size == 8) {
      properties.add(Ioss::Property("INTEGER_SIZE_DB", 8));
      properties.add(Ioss::Property("INTEGER_SIZE_API", 8));
    }

    if (interFace.use_netcdf4()) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf4"));
    }

    if (interFace.use_netcdf5()) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf5"));
    }

    if (interFace.compression_level() > 0 || interFace.szip()) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf4"));
      properties.add(Ioss::Property("COMPRESSION_LEVEL", interFace.compression_level()));
      properties.add(Ioss::Property("COMPRESSION_SHUFFLE", true));
      if (interFace.szip()) {
        properties.add(Ioss::Property("COMPRESSION_METHOD", "szip"));
      }
      else if (interFace.zlib()) {
        properties.add(Ioss::Property("COMPRESSION_METHOD", "zlib"));
      }
    }
    return properties;
  }

  std::shared_ptr<Ioss::Region> create_input_region(const std::string &key, std::string filename,
                                                    bool ints_32_bits)
  {
    // Check that 'filename' does not contain a starting/ending double quote...
    filename.erase(remove(filename.begin(), filename.end(), '\"'), filename.end());
    Ioss::DatabaseIO *dbi = Ioss::IOFactory::create("exodus", filename, Ioss::READ_RESTART,
                                                    Ioss::ParallelUtils::comm_self());
    if (dbi == nullptr || !dbi->ok(true)) {
      std::exit(EXIT_FAILURE);
    }

    dbi->set_surface_split_type(Ioss::SPLIT_BY_DONT_SPLIT);
    if (ints_32_bits) {
      dbi->set_int_byte_size_api(Ioss::USE_INT32_API);
    }
    else {
      dbi->set_int_byte_size_api(Ioss::USE_INT64_API);
    }

    // Splitting surfaces by element block makes it easier to transform the input
    // element id into the output element ids.
    dbi->set_surface_split_type(Ioss::SPLIT_BY_ELEMENT_BLOCK);

    // Generate a name for the region based on the key...
    std::string name = "Region_" + key;
    // NOTE: region owns database pointer at this time...
    return std::make_shared<Ioss::Region>(dbi, name);
  }

  void output_summary(Ioss::Region *region, std::ostream &strm)
  {
    int64_t nodes    = region->get_property("node_count").get_int();
    int64_t elements = region->get_property("element_count").get_int();

    fmt::print(strm, " Database: {}\tNodes = {} \tElements = {}\n",
               region->get_database()->get_filename(), fmt::group_digits(nodes),
               fmt::group_digits(elements));
  }

} // namespace
