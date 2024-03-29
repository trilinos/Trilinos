// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#include <algorithm>
#include <exception>
#include <fmt/chrono.h>
#include <fmt/ostream.h>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <unistd.h>

#include <cctype>
#include <cstring>
#include <ctime>

#include "add_to_log.h"
#include "adler.h"
#include "copy_string_cpp.h"
#include "format_time.h"
#include "sys_info.h"
#include "vector_data.h"
#if !USE_STD_SORT
#include "pdqsort.h"
#endif
#include "smart_assert.h"
#include "time_stamp.h"
#include <exodusII.h>

#if EX_API_VERS_NODOT <= 467
#error "Requires exodusII version 4.68 or later"
#endif

#include "CJ_ExodusEntity.h"
#include "CJ_ExodusFile.h"
#include "CJ_Internals.h"
#include "CJ_ObjectType.h"
#include "CJ_SystemInterface.h"
#include "CJ_Variables.h"
#include "CJ_Version.h"

namespace {
  bool                       check_variable_params(size_t p, const Excn::Variables &vars);
  template <typename T> void clear(std::vector<T> &vec)
  {
    vec.clear();
    vec.shrink_to_fit();
    SMART_ASSERT(vec.capacity() == 0);
  }

  template <typename T> bool approx_equal(T v1, T v2)
  {
#if 0
  static const T tolerance = 100.0 * std::numeric_limits<T>::epsilon();
  return std::fabs(v1 - v2) <= std::fabs(v1+v2)*tolerance;
#else
    return (float)v1 == (float)v2;
#endif
  }
} // namespace

struct NodeInfo
{
  NodeInfo() = default;
  NodeInfo(size_t id_, double x_, double y_, double z_) : id(id_), x(x_), y(y_), z(z_) {}
  size_t id{0};
  double x{0.0};
  double y{0.0};
  double z{0.0};

  bool operator==(const NodeInfo &other) const
  {
    return id == other.id && approx_equal(x, other.x) && approx_equal(y, other.y) &&
           approx_equal(z, other.z);
  }

  bool operator!=(const NodeInfo &other) const { return !(*this == other); }

  bool operator<(const NodeInfo &other) const
  {
    if (id < other.id) {
      return true;
    }
    if (id > other.id) {
      return false;
    }
    SMART_ASSERT(id == other.id);
    if (!approx_equal(x, other.x) && x < other.x) {
      return true;
    }
    if (!approx_equal(x, other.x) && x > other.x) {
      return false;
    }
    SMART_ASSERT(approx_equal(x, other.x));
    if (!approx_equal(y, other.y) && y < other.y) {
      return true;
    }
    if (!approx_equal(y, other.y) && y > other.y) {
      return false;
    }
    SMART_ASSERT(approx_equal(y, other.y));
    if (!approx_equal(z, other.z) && z < other.z) {
      return true;
    }
    return false;
  }
};

using GlobalMap = std::vector<NodeInfo>;
using GMapIter  = GlobalMap::iterator;

extern double seacas_timer();

namespace {
  template <class T> struct TimeStepMap
  {
    TimeStepMap(size_t part, int step, T time)
        : partNumber(part), localStepNumber(step), timeValue(time)
    {
    }
    size_t partNumber;
    int    localStepNumber;
    T      timeValue;
  };

  int get_width(size_t max_value);

  ex_entity_type exodus_object_type(const Excn::ObjectType &conjoin_type)
  {
    switch (conjoin_type) {
    case Excn::ObjectType::EBLK: return EX_ELEM_BLOCK;
    case Excn::ObjectType::SSET: return EX_SIDE_SET;
    case Excn::ObjectType::NSET: return EX_NODE_SET;
    default:
      throw std::runtime_error("Invalid Object Type in exodus_object_type: " +
                               std::to_string(static_cast<int>(conjoin_type)));
    }
  }

  char **get_name_array(int size, size_t length)
  {
    char **names = nullptr;
    if (size > 0) {
      names = new char *[size];
      for (int i = 0; i < size; i++) {
        names[i] = new char[length + 1];
        std::memset(names[i], '\0', length + 1);
      }
    }
    return names;
  }

  void free_name_array(char **names, int size)
  {
    for (int i = 0; i < size; i++) {
      delete[] names[i];
    }
    delete[] names;
    names = nullptr;
  }
} // namespace

std::string tsFormat = "[{:%H:%M:%S}] ";

// prototypes

// The main program templated to permit float/double transfer.
template <typename T, typename INT>
int conjoin(Excn::SystemInterface &interFace, T /* dummy */, INT /* dummy int */);

namespace {
  void sort_file_times(StringVector &input_files);

  void                         compress_white_space(char *str);
  void                         add_info_record(char *info_record, int size);
  template <typename INT> void put_mesh_summary(const Excn::Mesh<INT> &mesh);

  template <typename T>
  void get_put_qa(int id, int id_out, const std::vector<TimeStepMap<T>> &global_times,
                  Excn::SystemInterface &interFace);
  int  get_put_coordinate_names(int in, int out, int dimensionality);

  template <typename T, typename INT>
  int get_put_coordinates(Excn::Mesh<INT> &global, size_t part_count,
                          std::vector<Excn::Mesh<INT>> &local_mesh, T dummy);

  template <typename T, typename INT>
  int get_coordinates(int id, int dimensionality, size_t num_nodes,
                      std::vector<INT> local_node_to_global, size_t part, std::vector<T> &x,
                      std::vector<T> &y, std::vector<T> &z);

  StringVector get_exodus_variable_names(int id, ex_entity_type elType, size_t var_count);

  template <typename T, typename INT>
  void filter_truth_table(int id, Excn::Mesh<INT> &global, std::vector<T> &glob_blocks,
                          Excn::Variables &vars, const StringIdVector &variable_names);

  template <typename T, typename INT>
  void get_truth_table(Excn::Mesh<INT> &global, std::vector<std::vector<T>> &blocks,
                       std::vector<T> &glob_blocks, Excn::Variables &vars, int debug);

  template <typename U>
  void create_output_truth_table(std::vector<U> &global_sets, Excn::Variables &vars,
                                 std::vector<int> &truth_table);

  template <typename T, typename U, typename INT>
  int read_write_master_values(Excn::Variables &vars, const Excn::Mesh<INT> &global,
                               std::vector<U>               &global_sets,
                               std::vector<Excn::Mesh<INT>> &local_mesh,
                               std::vector<std::vector<U>> &local_sets, std::vector<T> &values,
                               size_t p, Excn::ExodusFile &id, int time_step, int time_step_out);

  void get_variable_params(int id, Excn::Variables &vars, const StringIdVector &variable_list);

  void get_put_variable_names(int id, int out, Excn::Variables &vars, Excn::SystemInterface &si,
                              int *combined_status_variable_index = nullptr);

  template <typename INT>
  void build_reverse_node_map(std::vector<Excn::Mesh<INT>> &local_mesh, Excn::Mesh<INT> *global,
                              size_t part_count, GlobalMap &global_node_map);

  template <typename INT>
  void build_reverse_node_map(std::vector<Excn::Mesh<INT>> &local_mesh, Excn::Mesh<INT> *global,
                              size_t part_count, std::vector<INT> &global_node_map);

  template <typename INT>
  void build_reverse_element_map(std::vector<Excn::Mesh<INT>>          &local_mesh,
                                 std::vector<std::vector<Excn::Block>> &blocks,
                                 std::vector<Excn::Block> &glob_blocks, Excn::Mesh<INT> *global,
                                 size_t                               part_count,
                                 std::vector<std::pair<INT, size_t>> &global_element_map);

  template <typename INT>
  void get_nodesets(size_t total_node_count, std::vector<Excn::Mesh<INT>> &local_mesh,
                    std::vector<std::vector<Excn::NodeSet<INT>>> &nodesets,
                    std::vector<Excn::NodeSet<INT>>              &glob_sets);

  template <typename INT>
  void get_element_blocks(const std::vector<Excn::Mesh<INT>>    &local_mesh,
                          const Excn::Mesh<INT>                 &global,
                          std::vector<std::vector<Excn::Block>> &blocks,
                          std::vector<Excn::Block>              &glob_blocks);
  template <typename T, typename INT>
  void put_element_blocks(std::vector<Excn::Mesh<INT>>          &local_mesh,
                          std::vector<std::vector<Excn::Block>> &blocks,
                          std::vector<Excn::Block> &glob_blocks, T single_or_double);

  template <typename INT> void put_nodesets(std::vector<Excn::NodeSet<INT>> &glob_sets);

  template <typename INT>
  void                         get_sideset_metadata(std::vector<Excn::Mesh<INT>>                 &local_mesh,
                                                    std::vector<std::vector<Excn::SideSet<INT>>> &sets,
                                                    std::vector<Excn::SideSet<INT>>              &glob_ssets);
  template <typename INT> void get_put_sidesets(std::vector<Excn::SideSet<INT>> &glob_ssets);

  template <typename T, typename INT>
  void add_status_variable(int id_out, const Excn::Mesh<INT> &global,
                           const std::vector<Excn::Block> &blocks,
                           const std::vector<Excn::Block> &glob_blocks,
                           const std::vector<INT> &local_element_to_global, int step, int variable,
                           T alive, int combined_variable_index);

  template <typename INT>
  size_t find_max_entity_count(size_t part_count, std::vector<Excn::Mesh<INT>> &local_mesh,
                               const Excn::Mesh<INT>                        &global,
                               std::vector<std::vector<Excn::Block>>        &blocks,
                               std::vector<std::vector<Excn::NodeSet<INT>>> &nodesets,
                               std::vector<std::vector<Excn::SideSet<INT>>> &sidesets);

  bool case_compare(const std::string &s1, const std::string &s2);

  template <typename T>
  void verify_set_position_mapping(const std::string &type, size_t part_count,
                                   const std::vector<T>              &global_sets,
                                   const std::vector<std::vector<T>> &sets)
  {
    bool problem = false;
    for (size_t p = 0; p < part_count; p++) {
      for (size_t i = 0; i < sets[p].size(); i++) {
        if (sets[p][i].id == 0) {
          continue;
        }
        auto glob_pos = sets[p][i].position_;
        if (global_sets[glob_pos].id != sets[p][i].id ||
            !case_compare(global_sets[glob_pos].name_, sets[p][i].name_)) {
          problem = true;
          fmt::print(stderr,
                     "\nERROR: {0} Mismatch on part {1}:\n"
                     "\tpart {0} at position {2} has id {3} and name {4}\n"
                     "\tglobal {0} at position {5} has id {6} and name {7}\n",
                     type, p + 1, i, sets[p][i].id, sets[p][i].name_, glob_pos,
                     global_sets[glob_pos].id, global_sets[glob_pos].name_);
          global_sets[glob_pos].dump();
          sets[p][i].dump();
        }
      }
    }

    if (problem) {
      throw std::runtime_error(type + " mismatch");
    }
  }

  // SEE: http://lemire.me/blog/2017/04/10/removing-duplicates-from-lists-quickly
  template <typename T> size_t unique(std::vector<T> &out)
  {
    if (out.empty()) {
      return 0;
    }
    size_t pos  = 1;
    T      oldv = out[0];
    for (size_t i = 1; i < out.size(); ++i) {
      T newv   = out[i];
      out[pos] = newv;
      pos += (newv != oldv);
      oldv = newv;
    }
    return pos;
  }

  template <typename T> void uniquify(std::vector<T> &vec)
  {
#if USE_STD_SORT
    std::sort(vec.begin(), vec.end());
#else
    pdqsort(vec.begin(), vec.end());
#endif
    vec.resize(unique(vec));
    vec.shrink_to_fit();
  }
} // namespace

unsigned int debug_level = 0;
int          main(int argc, char *argv[])
{
  try {
    time_t begin_time = time(nullptr);
    Excn::SystemInterface::show_version();

    Excn::SystemInterface interFace;
    bool                  ok = interFace.parse_options(argc, argv);

    if (!ok) {
      fmt::print(stderr, "\nERROR: Problems parsing command line arguments.\n\n");
      exit(EXIT_FAILURE);
    }

    debug_level = interFace.debug();

    if ((debug_level & 64) != 0U) {
      ex_opts(EX_VERBOSE | EX_DEBUG);
    }
    else {
      ex_opts(0);
    }

    int error = 0;

    if (interFace.sort_times()) {
      sort_file_times(interFace.inputFiles_);
    }

    if (!Excn::ExodusFile::initialize(interFace)) {
      fmt::print(stderr, "ERROR: Problem initializing input and/or output files.\n");
      exit(EXIT_FAILURE);
    }

    int int_byte_size = 4;
    if (interFace.ints_64_bit()) {
      int_byte_size = 8;
    }

    if (Excn::ExodusFile::io_word_size() == 4) {
      if (int_byte_size == 4) {
        error = conjoin(interFace, static_cast<float>(0.0), 0);
      }
      else {
        error = conjoin(interFace, static_cast<float>(0.0), static_cast<int64_t>(0));
      }
    }
    else {
      if (int_byte_size == 4) {
        error = conjoin(interFace, 0.0, 0);
      }
      else {
        error = conjoin(interFace, 0.0, static_cast<int64_t>(0));
      }
    }

    Excn::ExodusFile::close_all();

    time_t end_time = time(nullptr);
    add_to_log(argv[0], end_time - begin_time);
    return error;
  }
  catch (std::exception &e) {
    fmt::print(stderr, "ERROR: Standard exception: {}\n", e.what());
  }
}

template <typename T, typename INT>
int conjoin(Excn::SystemInterface &interFace, T /* dummy */, INT /* dummy int */)
{
  SMART_ASSERT(sizeof(T) == Excn::ExodusFile::io_word_size());

  const T alive      = interFace.alive_value();
  size_t  part_count = interFace.inputFiles_.size();

  std::array<char, MAX_LINE_LENGTH + 1> mytitle{};

  Excn::Mesh<INT> global;

  std::vector<Excn::Mesh<INT>> local_mesh(part_count);

  // ******************************************************************
  // 1. Read global info

  int error = 0;

  if (debug_level & 1) {
    fmt::print("{}", time_stamp(tsFormat));
  }

  std::string title0;

  for (size_t p = 0; p < part_count; p++) {
    ex_init_params info{};
    error += ex_get_init_ext(Excn::ExodusFile(p), &info);

    local_mesh[p].title          = info.title;
    local_mesh[p].dimensionality = info.num_dim;
    local_mesh[p].nodeCount      = info.num_nodes;
    local_mesh[p].elementCount   = info.num_elem;
    local_mesh[p].blockCount     = info.num_elem_blk;
    local_mesh[p].nodesetCount   = info.num_node_sets;
    local_mesh[p].sidesetCount   = info.num_side_sets;

    if (p == 0) {
      global.title          = mytitle.data();
      global.dimensionality = local_mesh[p].count(Excn::ObjectType::DIM);
      global.blockCount     = local_mesh[p].count(Excn::ObjectType::EBLK);
    }
    else {
      SMART_ASSERT(global.count(Excn::ObjectType::DIM) ==
                   local_mesh[p].count(Excn::ObjectType::DIM))
      (global.count(Excn::ObjectType::DIM))(local_mesh[p].count(Excn::ObjectType::DIM))(p);
      SMART_ASSERT(global.count(Excn::ObjectType::EBLK) ==
                   local_mesh[p].count(Excn::ObjectType::EBLK))
      (global.count(Excn::ObjectType::EBLK))(local_mesh[p].count(Excn::ObjectType::EBLK))(p);
    }

    local_mesh[p].localNodeToGlobal.resize(local_mesh[p].count(Excn::ObjectType::NODE));
    local_mesh[p].localElementToGlobal.resize(local_mesh[p].count(Excn::ObjectType::ELEM));

  } // end for (p=0..part_count)

  if (interFace.omit_nodesets()) {
    global.nodesetCount = 0;
  }

  if (interFace.omit_sidesets()) {
    global.sidesetCount = 0;
  }

  // Get database times...  Save mapping from time step to part providing that time step.
  std::vector<TimeStepMap<T>> global_times;

  double t_min = FLT_MAX;
  for (size_t p = part_count; p > 0; p--) {

    Excn::ExodusFile id(p - 1);
    bool             used = false;
    int              nts  = ex_inquire_int(id, EX_INQ_TIME);
    if (nts == 0) {
      std::string part = "Part " + std::to_string(p) + ": ";
      part += interFace.inputFiles_[p - 1];
      fmt::print(stderr,
                 "\nWARNING: '{}'\n\tdoes not contain any time steps so it will not be in the "
                 "conjoined file.\n",
                 part);
      local_mesh[p - 1].isActive = false;
    }
    else {
      std::vector<T> times(nts);
      ex_get_all_times(id, Data(times));

      // A database will include all times from step 0 up to the
      // last time that is less than t_min on the following database.
      // Note that we are iterating through the databases from last to first...

      // Check delta between the last used timestep on this database
      // and the first time on the previous (later in time) database
      // (which is t_min). It must be > user-specified minimum delta.
      // Set by -interpart_minimum_time_delta {delta} command line option.

      // NOTE that this is not the delta between individual timesteps
      // within a database, it is the delta between the last timestep
      // on one database and the first on the next database.

      int i = 0;
      for (i = nts; i > 0; i--) {
        if (times[i - 1] < t_min) {
          if (used || t_min - times[i - 1] >= interFace.interpart_minimum_time_delta()) {
            used = true;
            global_times.push_back(TimeStepMap<T>(p - 1, i - 1, times[i - 1]));
          }
        }
      }
      local_mesh[p - 1].timestepCount = i;
      t_min                           = t_min < times[0] ? t_min : times[0];

      if (!used) {
        std::string part = "Part " + std::to_string(p) + ": ";
        part += interFace.inputFiles_[p - 1];
        fmt::print(stderr,
                   "\nWARNING: '{}'\n\tdoes not contain any time steps which will be used in "
                   "conjoined file.\n"
                   "\tCurrent minimum time = {}, timestep range on this part is {} to {}\n",
                   part, t_min, times[0], times[nts - 1]);
        local_mesh[p - 1].isActive = false;
      }
    }
  }
  global.timestepCount = global_times.size();
  std::reverse(global_times.begin(), global_times.end());

  // Need these throughout run, so declare outside of this block...
  // TODO(gdsjaar): Add these to the "Mesh" class.
  std::vector<Excn::Block>              glob_blocks(global.count(Excn::ObjectType::EBLK));
  std::vector<std::vector<Excn::Block>> blocks(part_count);

  std::vector<Excn::SideSet<INT>>              glob_ssets;
  std::vector<std::vector<Excn::SideSet<INT>>> sidesets(part_count);

  std::vector<Excn::NodeSet<INT>>              glob_nsets;
  std::vector<std::vector<Excn::NodeSet<INT>>> nodesets(part_count);
  {
    // Now, build the reverse global node map which permits access of the
    // local id given the global id.
    std::vector<INT> global_node_map;
    if (interFace.ignore_coordinates()) {
      build_reverse_node_map(local_mesh, &global, part_count, global_node_map);
    }
    else {
      GlobalMap global_node_coord_map;
      build_reverse_node_map(local_mesh, &global, part_count, global_node_coord_map);
      global_node_map.reserve(global_node_coord_map.size());
      for (const auto &ni : global_node_coord_map) {
        global_node_map.push_back(ni.id);
      }
    }

    if (debug_level & 1) {
      fmt::print("{}", time_stamp(tsFormat));
    }

    // ****************************************************************************
    // 5. Get Block information including element attributes
    // must check for zero length blocks
    get_element_blocks(local_mesh, global, blocks, glob_blocks);

    std::vector<std::pair<INT, size_t>> global_element_map(global.count(Excn::ObjectType::ELEM));
    build_reverse_element_map(local_mesh, blocks, glob_blocks, &global, part_count,
                              global_element_map);

    //
    //    NOTE:  Node set/side set information can be different for each part
    /************************************************************************/
    // 7. Get Side sets
    if (!interFace.omit_sidesets()) {
      if (debug_level & 1) {
        fmt::print("{}", time_stamp(tsFormat));
      }
      get_sideset_metadata(local_mesh, sidesets, glob_ssets);
      if (global.count(Excn::ObjectType::SSET) != glob_ssets.size()) {
        global.sidesetCount = glob_ssets.size();
      }
    }

    /************************************************************************/
    // 6. Get Node sets
    if (!interFace.omit_nodesets()) {
      if (debug_level & 1) {
        fmt::print("{}", time_stamp(tsFormat));
      }
      get_nodesets(global.count(Excn::ObjectType::NODE), local_mesh, nodesets, glob_nsets);
      if (global.count(Excn::ObjectType::NSET) != glob_nsets.size()) {
        global.nodesetCount = glob_nsets.size();
      }
    }

    /************************************************************************/
    // Start writing the output file...

    if (debug_level & 1) {
      fmt::print("{}", time_stamp(tsFormat));
    }
    Excn::CommunicationMetaData comm_data;

    // Create the output file...
    Excn::ExodusFile::create_output(interFace);

    put_mesh_summary(global);

    get_put_qa(Excn::ExodusFile(0), Excn::ExodusFile::output(), global_times, interFace);

    Excn::Internals exodus(Excn::ExodusFile::output(), Excn::ExodusFile::max_name_length());

    exodus.write_meta_data(global, glob_blocks, glob_nsets, glob_ssets, comm_data);

    // Output bulk mesh data....
    put_nodesets(glob_nsets);

    // c.2.  Write Global Node Number Map
    if (debug_level & 1) {
      fmt::print("{}", time_stamp(tsFormat));
    }

    error = ex_put_id_map(Excn::ExodusFile::output(), EX_NODE_MAP, Data(global_node_map));

    if (debug_level & 1) {
      fmt::print("{}", time_stamp(tsFormat));
    }
    if (!global_element_map.empty()) {
      std::vector<INT> global_map(global.count(Excn::ObjectType::ELEM));
      for (size_t i = 0; i < global.count(Excn::ObjectType::ELEM); i++) {
        global_map[i] = global_element_map[i].first;
      }
      ex_put_id_map(Excn::ExodusFile::output(), EX_ELEM_MAP, Data(global_map));
    }

    T dummy = 0;
    put_element_blocks(local_mesh, blocks, glob_blocks, dummy);
    get_put_sidesets(glob_ssets);
  }
  // ************************************************************************
  // 2. Get Coordinate Info.
  {
    if (debug_level & 1) {
      fmt::print("{}", time_stamp(tsFormat));
    }

    error += get_put_coordinates(global, part_count, local_mesh, (T)0);

    if (debug_level & 1) {
      fmt::print("{}", time_stamp(tsFormat));
    }
    fmt::print("Wrote coordinate information...\n");
  }
  // ####################TRANSIENT DATA SECTION###########################
  // ***********************************************************************
  // 9. Get Variable Information and names

  if (debug_level & 1) {
    fmt::print("{}", time_stamp(tsFormat));
  }

  //  I. read number of variables for each type. Note that exodusII does not
  //     provide for history variables
  //  NOTE: it is assumed that every part has the same global, nodal,
  //        and element lists

  bool            add_n_status = interFace.nodal_status_variable() != "NONE";
  bool            add_e_status = interFace.element_status_variable() != "NONE";
  Excn::Variables global_vars(Excn::ObjectType::GLOBAL);
  Excn::Variables nodal_vars(Excn::ObjectType::NODE, add_n_status);
  Excn::Variables element_vars(Excn::ObjectType::EBLK, add_e_status);
  Excn::Variables nodeset_vars(Excn::ObjectType::NSET);
  Excn::Variables sideset_vars(Excn::ObjectType::SSET);

  {
    Excn::ExodusFile id(0);

    get_variable_params(id, global_vars, interFace.global_var_names());
    get_variable_params(id, nodal_vars, interFace.node_var_names());
    get_variable_params(id, element_vars, interFace.elem_var_names());
    get_variable_params(id, nodeset_vars, interFace.nset_var_names());
    get_variable_params(id, sideset_vars, interFace.sset_var_names());

    get_truth_table(global, blocks, glob_blocks, element_vars, 4);
    filter_truth_table(id, global, glob_blocks, element_vars, interFace.elem_var_names());

    get_truth_table(global, nodesets, glob_nsets, nodeset_vars, 32);
    filter_truth_table(id, global, glob_nsets, nodeset_vars, interFace.nset_var_names());

    get_truth_table(global, sidesets, glob_ssets, sideset_vars, 16);
    filter_truth_table(id, global, glob_ssets, sideset_vars, interFace.sset_var_names());
  }

  // Check that the variable counts are the same on the subsequent files...
  // Error out if there is a difference...
  bool found_error = false;
  for (size_t p = 1; p < part_count; p++) {
    found_error |= check_variable_params(p, global_vars);
    found_error |= check_variable_params(p, nodal_vars);
    found_error |= check_variable_params(p, element_vars);
    found_error |= check_variable_params(p, nodeset_vars);
    found_error |= check_variable_params(p, sideset_vars);
  }
  if (found_error) {
    return 1;
  }

  // There is a slightly tricky situation here. The truthTable block order
  // is based on the ordering of the blocks on the input databases.
  // These blocks may have been reordered on output to make the 'offset'
  // variables line up correctly when the element ids are mapped back to
  // the "overall global" order.  There is not much problem since most
  // calls  outputting block-related items pass the id and don't rely
  // on ordering. However, the truth table is one of the exceptions
  // and we need to reorder the truth table to match the output block
  // order. After this call, we can use the original ordering, so just
  // need a temporary vector here...
  if (global_vars.count(Excn::InOut::OUT_) + nodal_vars.count(Excn::InOut::OUT_) +
          element_vars.count(Excn::InOut::OUT_) + nodeset_vars.count(Excn::InOut::OUT_) +
          sideset_vars.count(Excn::InOut::OUT_) >
      0) {

    std::vector<int> elem_truth_table;
    std::vector<int> nset_truth_table;
    std::vector<int> sset_truth_table;
    create_output_truth_table(glob_blocks, element_vars, elem_truth_table);
    create_output_truth_table(glob_nsets, nodeset_vars, nset_truth_table);
    create_output_truth_table(glob_ssets, sideset_vars, sset_truth_table);

    ex_put_all_var_param(Excn::ExodusFile::output(), global_vars.count(Excn::InOut::OUT_),
                         nodal_vars.count(Excn::InOut::OUT_), element_vars.count(Excn::InOut::OUT_),
                         Data(elem_truth_table), nodeset_vars.count(Excn::InOut::OUT_),
                         Data(nset_truth_table), sideset_vars.count(Excn::InOut::OUT_),
                         Data(sset_truth_table));
  }

  // II. read/write the variable names
  int combined_status_variable_index = 0;
  {
    Excn::ExodusFile id(0);
    get_put_variable_names(id, Excn::ExodusFile::output(), global_vars, interFace);
    get_put_variable_names(id, Excn::ExodusFile::output(), nodal_vars, interFace);
    get_put_variable_names(id, Excn::ExodusFile::output(), element_vars, interFace,
                           &combined_status_variable_index);
    get_put_variable_names(id, Excn::ExodusFile::output(), nodeset_vars, interFace);
    get_put_variable_names(id, Excn::ExodusFile::output(), sideset_vars, interFace);
  }
  ex_update(Excn::ExodusFile::output());

  /**********************************************************************/
  // 10. Get Transient Data
  //     This routine reads in a time dump from an EXODUSII file

  size_t num_time_steps = global.count(Excn::ObjectType::TIME);

  if (debug_level & 1) {
    fmt::print("{}", time_stamp(tsFormat));
  }

  std::vector<T> global_values(global_vars.count(Excn::InOut::IN_));
  std::vector<T> output_global_values(global_vars.count(Excn::InOut::OUT_));

  // Determine maximum number of entities on any part
  auto max_ent = find_max_entity_count(part_count, local_mesh, global, blocks, nodesets, sidesets);
  std::vector<T> values(max_ent);

  // Stage II.  Extracting transient variable data.
  //            loop over time steps

  // Determine if user wants a subset of timesteps transferred to the output file.
  // Time steps for output file
  double start_time = seacas_timer();

  // Used for output formatting at end of loop.
  int percent_width = 0;
  int field_width   = 3;
  if (num_time_steps > 100) {
    percent_width = 1;
    field_width   = 5;
  }

  int element_width = get_width(global.count(Excn::ObjectType::ELEM));
  int step_width    = get_width(num_time_steps);
  step_width += step_width / 3; // For commas -- 1,234,456

  int part_width     = get_width(part_count + 1);
  int loc_step_width = 0;
  for (auto &global_time : global_times) {
    if (loc_step_width < global_time.localStepNumber + 1) {
      loc_step_width = global_time.localStepNumber + 1;
    }
  }
  loc_step_width = get_width(loc_step_width);
  loc_step_width += loc_step_width / 3;

  size_t time_step_out = 0;
  for (size_t time_step = 0; time_step < num_time_steps; time_step++) {
    time_step_out++;

    T                time_val = 0;
    size_t           p        = global_times[time_step].partNumber;
    Excn::ExodusFile id(p);

    // read in and write out the time step information
    error += ex_get_time(id, global_times[time_step].localStepNumber + 1, &time_val);
    SMART_ASSERT(time_val == global_times[time_step].timeValue)
    (time_step)(time_val)(global_times[time_step].timeValue);
    error += ex_put_time(Excn::ExodusFile::output(), time_step_out, &time_val);

    if (global_vars.count(Excn::InOut::OUT_) > 0) {
      if (debug_level & 1) {
        fmt::print("{}Global Variables...\n", time_stamp(tsFormat));
      }
      error += ex_get_var(id, global_times[time_step].localStepNumber + 1, EX_GLOBAL, 0, 0,
                          global_vars.count(), (void *)Data(global_values));
      // Map ...
      for (int ig = 0; ig < global_vars.count(Excn::InOut::IN_); ig++) {
        if (global_vars.index_[ig] > 0) {
          SMART_ASSERT(ig < (int)global_values.size());
          output_global_values[global_vars.index_[ig] - 1] = global_values[ig];
        }
      }
      error += ex_put_var(Excn::ExodusFile::output(), time_step_out, EX_GLOBAL, 1, 0,
                          global_vars.count(Excn::InOut::OUT_), Data(output_global_values));
    }

    // ========================================================================
    // Nodal Values...
    if (debug_level & 1) {
      fmt::print("{}Nodal Variables...\n", time_stamp(tsFormat));
    }

    if (nodal_vars.count(Excn::InOut::OUT_) > 0) {
      size_t         node_count = local_mesh[p].count(Excn::ObjectType::NODE);
      std::vector<T> master_nodal_values(global.count(Excn::ObjectType::NODE));

      int offset = nodal_vars.addStatus ? 1 : 0;
      for (int i = 0; i < nodal_vars.count(Excn::InOut::OUT_) - offset;
           i++) { // Last output variable may be status
        for (int j = 0; j < nodal_vars.count(Excn::InOut::IN_); j++) {
          if (nodal_vars.index_[j] - 1 == i) {
            std::fill(master_nodal_values.begin(), master_nodal_values.end(), T(0.0));

            error += ex_get_var(id, global_times[time_step].localStepNumber + 1, EX_NODAL, j + 1, 0,
                                node_count, Data(values));

            // copy values to master nodal value information
            for (size_t jj = 0; jj < node_count; jj++) {
              // Map local nodal value to global location...
              size_t nodal_value               = local_mesh[p].localNodeToGlobal[jj];
              master_nodal_values[nodal_value] = values[jj];
            }
            error += ex_put_var(Excn::ExodusFile::output(), time_step_out, EX_NODAL, i + 1, 0,
                                global.count(Excn::ObjectType::NODE), Data(master_nodal_values));
            break;
          }
        }

        // Fill "node_status" variable -- 'alive' for alive; 1-alive for dead.
        // It is the last output variable...
        if (nodal_vars.addStatus) {
          SMART_ASSERT(alive == 0.0 || alive == 1.0)(alive);
          std::fill(master_nodal_values.begin(), master_nodal_values.end(), T((1.0 - alive)));
          for (size_t j = 0; j < node_count; j++) {
            // Map local nodal value to global location...
            size_t nodal_value               = local_mesh[p].localNodeToGlobal[j];
            master_nodal_values[nodal_value] = alive;
          }

          error += ex_put_var(Excn::ExodusFile::output(), time_step_out, EX_NODAL,
                              nodal_vars.count(Excn::InOut::OUT_), 0,
                              global.count(Excn::ObjectType::NODE), Data(master_nodal_values));
        }
      }
    }

    // ========================================================================
    // Extracting element transient variable data
    if (debug_level & 1) {
      fmt::print("{}Element Variables...\n", time_stamp(tsFormat));
    }
    if (element_vars.count(Excn::InOut::IN_) > 0) {
      read_write_master_values(element_vars, global, glob_blocks, local_mesh, blocks, values, p, id,
                               global_times[time_step].localStepNumber, time_step_out);
    }

    // Add element status variable...
    // Use the output time step for writing data
    if (interFace.element_status_variable() != "NONE") {
      add_status_variable(Excn::ExodusFile::output(), global, blocks[p], glob_blocks,
                          local_mesh[p].localElementToGlobal, time_step_out,
                          element_vars.index_[element_vars.count(Excn::InOut::IN_)], alive,
                          combined_status_variable_index);
    }

    // ========================================================================
    // Extracting sideset transient variable data
    if (debug_level & 1) {
      fmt::print("{}Sideset Variables...\n", time_stamp(tsFormat));
    }
    if (sideset_vars.count(Excn::InOut::IN_) > 0) {
      read_write_master_values(sideset_vars, global, glob_ssets, local_mesh, sidesets, values, p,
                               id, global_times[time_step].localStepNumber, time_step_out);
    }

    // ========================================================================
    // Extracting nodeset transient variable data
    if (debug_level & 1) {
      fmt::print("{}Nodeset Variables...\n", time_stamp(tsFormat));
    }
    if (nodeset_vars.count(Excn::InOut::IN_) > 0) {
      read_write_master_values(nodeset_vars, global, glob_nsets, local_mesh, nodesets, values, p,
                               id, global_times[time_step].localStepNumber, time_step_out);
    }

    // ========================================================================
    ex_update(Excn::ExodusFile::output());

    if (debug_level & 1) {
      fmt::print("{}", time_stamp(tsFormat));
    }

    fmt::print("Step {:{}}/{},  time {:.4e}  (Part {:{}}/{},  step {:{}})   Active Elem: {:{}}",
               fmt::group_digits(time_step + 1), step_width, fmt::group_digits(num_time_steps),
               time_val, p + 1, part_width, part_count,
               fmt::group_digits(global_times[time_step].localStepNumber + 1), loc_step_width,
               fmt::group_digits(local_mesh[p].count(Excn::ObjectType::ELEM)), element_width);

    double cur_time        = seacas_timer();
    double elapsed         = cur_time - start_time;
    double time_per_step   = elapsed / time_step_out;
    double percentage_done = (time_step_out * 100.0) / global.count(Excn::ObjectType::TIME);
    double estimated_remaining =
        time_per_step * (global.count(Excn::ObjectType::TIME) - time_step_out);

    fmt::print("  [{:{}.{}f}%, Elapsed={}, ETA={}]\n", percentage_done, field_width, percent_width,
               format_time(elapsed), format_time(estimated_remaining));
    if (debug_level & 1) {
      fmt::print("\n");
    }
  }

  /*************************************************************************/
  // EXIT program
  if (debug_level & 1) {
    fmt::print("{}", time_stamp(tsFormat));
  }
  fmt::print("******* END *******\n");
  return error;
}

namespace {
  template <typename T>
  void get_put_qa(int id, int id_out, const std::vector<TimeStepMap<T>> & /*global_times*/,
                  Excn::SystemInterface & /*interFace*/)
  {
    // NOTE: Assuming info and QA records for all parts
    int error = 0;

    int info_string_len = MAX_LINE_LENGTH;

    size_t num_info_records = ex_inquire_int(id, EX_INQ_INFO);
    //    size_t extra_info = global_times.size() + 2 + 1;
    size_t extra_info = 2 + 1;

    char **info_records = get_name_array(num_info_records + extra_info, info_string_len);

    if (num_info_records > 0) {
      error += ex_get_info(id, info_records);
    }

    // Add an info record for CONJOIN
    add_info_record(info_records[num_info_records], MAX_LINE_LENGTH);

#if 0
    // Add time/part mapping...
    for (size_t i=0; i < global_times.size(); i++) {
      std::ostringstream os;
      fmt::print(os, "Step {:2}, time {:.4e} (Part {}, step {})  File: {}",
                 i+1, global_times[i].timeValue,
                 global_times[i].partNumber+1,
                 global_times[i].localStepNumber+1,
                 interFace.inputFiles_[global_times[i].partNumber]);

      copy_string(info_records[num_info_records+1+i], os.str(), MAX_LINE_LENGTH + 1);
    }
#endif

    error += ex_put_info(id_out, num_info_records + extra_info, info_records);

    free_name_array(info_records, num_info_records + extra_info);

    // II. Get and store QA records, if they exist
    struct qa_element
    {
      char *qa_record[1][4];
    };

    int                     num_qa_records = ex_inquire_int(id, EX_INQ_QA);
    std::vector<qa_element> qaRecord(num_qa_records + 1);
    for (int i = 0; i < num_qa_records + 1; i++) {
      for (int j = 0; j < 4; j++) {
        qaRecord[i].qa_record[0][j]    = new char[MAX_STR_LENGTH + 1];
        qaRecord[i].qa_record[0][j][0] = '\0';
      }
    }
    if (num_qa_records) {
      error += ex_get_qa(id, qaRecord[0].qa_record);
    }

    copy_string(qaRecord[num_qa_records].qa_record[0][0], qainfo[0], MAX_STR_LENGTH + 1); // Code
    copy_string(qaRecord[num_qa_records].qa_record[0][1], qainfo[2], MAX_STR_LENGTH + 1); // Version

    std::time_t date_time = std::time(nullptr);

    auto date = fmt::format("{:%Y/%m/%d}", fmt::localtime(date_time));
    copy_string(qaRecord[num_qa_records].qa_record[0][2], date.c_str(), MAX_STR_LENGTH + 1);

    auto time = fmt::format("{:%T}", fmt::localtime(date_time));
    copy_string(qaRecord[num_qa_records].qa_record[0][3], time.c_str(), MAX_STR_LENGTH + 1);

    error += ex_put_qa(id_out, num_qa_records + 1, qaRecord[0].qa_record);

    for (int i = 0; i < num_qa_records + 1; i++) {
      for (int j = 0; j < 4; j++) {
        delete[] qaRecord[i].qa_record[0][j];
      }
    }
  }

  template <typename T, typename INT>
  int get_put_coordinates(Excn::Mesh<INT> &global, size_t part_count,
                          std::vector<Excn::Mesh<INT>> &local_mesh, T /* dummy */)
  {
    SMART_ASSERT(sizeof(T) == Excn::ExodusFile::io_word_size());
    std::vector<T> x(global.count(Excn::ObjectType::NODE));
    std::vector<T> y(global.count(Excn::ObjectType::NODE));
    std::vector<T> z(global.count(Excn::ObjectType::NODE));

    if (debug_level & 8) {
      const T FILL_VALUE = FLT_MAX;
      std::fill(x.begin(), x.end(), FILL_VALUE);
      std::fill(y.begin(), y.end(), FILL_VALUE);
      std::fill(z.begin(), z.end(), FILL_VALUE);
    }

    int error = 0;
    for (size_t p = 0; p < part_count; p++) {
      error += get_coordinates(Excn::ExodusFile(p), global.count(Excn::ObjectType::DIM),
                               local_mesh[p].count(Excn::ObjectType::NODE),
                               local_mesh[p].localNodeToGlobal, p, x, y, z);
    } // end for p=0..part_count

    // Get Coordinate Names
    // NOTE: Assuming coordinate names should be
    // the same for all files/parts therefore, only one
    // file/part needs to be loaded
    error += get_put_coordinate_names(Excn::ExodusFile(0), Excn::ExodusFile::output(),
                                      global.count(Excn::ObjectType::DIM));
    // Write out coordinate information
    error += ex_put_coord(Excn::ExodusFile::output(), Data(x), Data(y), Data(z));
    return error;
  }

  int get_put_coordinate_names(int in, int out, int dimensionality)
  {
    int    error            = 0;
    char **coordinate_names = get_name_array(dimensionality, Excn::ExodusFile::max_name_length());

    error += ex_get_coord_names(in, coordinate_names);
    error += ex_put_coord_names(out, coordinate_names);
    fmt::print("Wrote coordinate names...\n");

    free_name_array(coordinate_names, dimensionality);
    return error;
  }

  template <typename T, typename INT>
  int get_coordinates(int id, int dimensionality, size_t num_nodes,
                      std::vector<INT> local_node_to_global, size_t part, std::vector<T> &x,
                      std::vector<T> &y, std::vector<T> &z)
  {
    SMART_ASSERT(sizeof(T) == Excn::ExodusFile::io_word_size());
    SMART_ASSERT(local_node_to_global.size() == num_nodes);
    int            error = 0;
    std::vector<T> local_x(num_nodes);
    std::vector<T> local_y(num_nodes);
    std::vector<T> local_z(num_nodes);

    error += ex_get_coord(id, Data(local_x), Data(local_y), Data(local_z));

    // Check for 2D or 3D coordinates
    if (dimensionality == 3) {
      if (debug_level & 8) {
        for (size_t i = 0; i < num_nodes; i++) {
          INT     node       = local_node_to_global[i];
          const T FILL_VALUE = FLT_MAX;
          if (x[node] != FILL_VALUE && y[node] != FILL_VALUE && z[node] != FILL_VALUE) {
            if (!approx_equal(x[node], local_x[i]) || !approx_equal(y[node], local_y[i]) ||
                !approx_equal(z[node], local_z[i])) {
              fmt::print(
                  stderr,
                  "\nWARNING: Node {} has different coordinates in at least two parts.\n"
                  "         this may indicate that this id has been reused in the current part.\n"
                  "         cur value = {:14.6e} {:14.6e} {:14.6e}\n"
                  "         new value = {:14.6e} {:14.6e} {:14.6e} from part {}\n",
                  node + 1, x[node], y[node], z[node], local_x[i], local_y[i], local_z[i], part);
            }
          }
        }
      }

      for (size_t i = 0; i < num_nodes; i++) {
        // The following WILL overwrite x[node],y[node],z[node] if node is the
        // same for different parts
        INT node = local_node_to_global[i];
        x[node]  = local_x[i];
        y[node]  = local_y[i];
        z[node]  = local_z[i];
      }
    }
    else {
      if (debug_level & 8) {
        for (size_t i = 0; i < num_nodes; i++) {
          INT     node       = local_node_to_global[i];
          const T FILL_VALUE = FLT_MAX;
          if (x[node] != FILL_VALUE && y[node] != FILL_VALUE) {
            if (!approx_equal(x[node], local_x[i]) || !approx_equal(y[node], local_y[i])) {
              fmt::print(
                  stderr,
                  "\nWARNING: Node {} has different coordinates in at least two parts.\n"
                  "         this may indicate that this id has been reused in the current part.\n"
                  "         cur value = {:14.6e} {:14.6e}\n"
                  "         new value = {:14.6e} {:14.6e} from part {}\n",
                  node + 1, x[node], y[node], local_x[i], local_y[i], part);
            }
          }
        }
      }

      for (size_t i = 0; i < num_nodes; i++) {
        // The following WILL overwrite x[node],y[node] if node is the same for
        // different parts
        INT node = local_node_to_global[i];
        x[node]  = local_x[i];
        y[node]  = local_y[i];
      }
    }
    return error;
  }

  template <typename INT>
  void get_element_blocks(const std::vector<Excn::Mesh<INT>>    &local_mesh,
                          const Excn::Mesh<INT>                 &global,
                          std::vector<std::vector<Excn::Block>> &blocks,
                          std::vector<Excn::Block>              &glob_blocks)
  {
    size_t part_count = local_mesh.size();
    for (size_t ip = 0; ip < part_count; ip++) {
      blocks[ip].resize(local_mesh[ip].count(Excn::ObjectType::EBLK));
    }

    std::vector<INT> block_id(global.count(Excn::ObjectType::EBLK));

    int error = 0;
    for (size_t p = 0; p < part_count; p++) {
      Excn::ExodusFile id(p);

      error += ex_get_ids(id, EX_ELEM_BLOCK, Data(block_id));

      if (p == 0) {
        for (size_t b = 0; b < global.count(Excn::ObjectType::EBLK); b++) {
          glob_blocks[b].id = block_id[b];
        }
      }

      for (size_t b = 0; b < global.count(Excn::ObjectType::EBLK); b++) {
        ex_block block_param{};
        block_param.id   = glob_blocks[b].id;
        block_param.type = EX_ELEM_BLOCK;

        error += ex_get_block_param(id, &block_param);

        std::vector<char> name(Excn::ExodusFile::max_name_length() + 1);
        ex_get_name(id, EX_ELEM_BLOCK, glob_blocks[b].id, Data(name));

        blocks[p][b].id = glob_blocks[b].id;
        if (name[0] != '\0') {
          blocks[p][b].name_ = Data(name);
          if (p == 0) {
            glob_blocks[b].name_ = Data(name);
          }
        }

        // Find position of this block id in the local mesh
        for (size_t lb = 0; lb < global.count(Excn::ObjectType::EBLK); lb++) {
          if (block_id[lb] == glob_blocks[b].id) {
            blocks[p][b].position_ = lb;
            break;
          }
        }

        if (block_param.num_entry > 0) {
          blocks[p][b].elementCount    = block_param.num_entry;
          blocks[p][b].nodesPerElement = block_param.num_nodes_per_entry;
          blocks[p][b].attributeCount  = block_param.num_attribute;
          blocks[p][b].offset_         = block_param.num_entry;
          blocks[p][b].elType          = block_param.topology;

          // NOTE: This is not correct, but fixed below.
          glob_blocks[b].elementCount += block_param.num_entry;

          if (glob_blocks[b].nodesPerElement == 0) {
            glob_blocks[b].nodesPerElement = block_param.num_nodes_per_entry;
          }

          if (glob_blocks[b].attributeCount == 0) {
            glob_blocks[b].attributeCount = block_param.num_attribute;
          }

          glob_blocks[b].position_ = b;
          glob_blocks[b].elType    = block_param.topology;
        }

        if (block_param.num_attribute > 0 && glob_blocks[b].attributeNames.empty()) {
          // Get attribute names.  Assume the same on all parts
          // on which the block exists.
          char **names =
              get_name_array(block_param.num_attribute, Excn::ExodusFile::max_name_length());

          ex_get_attr_names(id, EX_ELEM_BLOCK, block_id[b], names);
          for (int i = 0; i < block_param.num_attribute; i++) {
            glob_blocks[b].attributeNames.emplace_back(names[i]);
          }
          free_name_array(names, block_param.num_attribute);
        }
      }
    } // end for p=0..part_count

    // Convert block_offset from elements/block/part to true offset
    for (size_t p = 0; p < part_count; p++) {
      if (debug_level & 4) {
        fmt::print("\nElement block info for part {}...\n", p);
      }
      // Number of elements per block in local block order...
      std::vector<size_t> local_order_entity_count(global.count(Excn::ObjectType::EBLK));
      for (size_t b = 0; b < global.count(Excn::ObjectType::EBLK); b++) {
        int local_order                       = blocks[p][b].position_;
        local_order_entity_count[local_order] = blocks[p][b].entity_count();
      }
      size_t sum = 0;
      for (size_t b = 0; b < global.count(Excn::ObjectType::EBLK); b++) {
        size_t save                 = local_order_entity_count[b];
        local_order_entity_count[b] = sum;
        sum += save;
      }
      for (size_t b = 0; b < global.count(Excn::ObjectType::EBLK); b++) {
        int local_order      = blocks[p][b].position_;
        blocks[p][b].offset_ = local_order_entity_count[local_order];

        if (debug_level & 4) {
          fmt::print("\tBlock {}, Id = {}, Name = '{}', Elements = {:12}, Nodes/element = {}, "
                     "Attributes = {}, Position = {}, Offset = {}\n",
                     b, glob_blocks[b].id, blocks[p][b].name_,
                     fmt::group_digits(blocks[p][b].entity_count()), blocks[p][b].nodesPerElement,
                     blocks[p][b].attributeCount, blocks[p][b].position_, blocks[p][b].offset_);
        }
      }
    }
  }

  template <typename T, typename INT>
  void put_element_blocks(std::vector<Excn::Mesh<INT>>          &local_mesh,
                          std::vector<std::vector<Excn::Block>> &blocks,
                          std::vector<Excn::Block>              &glob_blocks, T /* dummy */)
  {
    SMART_ASSERT(sizeof(T) == Excn::ExodusFile::io_word_size());
    int global_num_blocks = glob_blocks.size();

    if (debug_level & 1) {
      fmt::print("{}", time_stamp(tsFormat));
    }
    fmt::print("\nReading and Writing element connectivity & attributes\n");

    for (int b = 0; b < global_num_blocks; b++) {

      if (debug_level & 4) {
        fmt::print("\nOutput element block info for...\n"
                   "Block {}, Id = {}, Name = '{}', Elements = {:12}, Nodes/element = {}, "
                   "Attributes = {}\n",
                   b, glob_blocks[b].id, glob_blocks[b].name_,
                   fmt::group_digits(glob_blocks[b].entity_count()), glob_blocks[b].nodesPerElement,
                   glob_blocks[b].attributeCount);
      }

      if (debug_level & 4) {
        fmt::print("B{}:\t", b);
      }

      size_t max_nodes = glob_blocks[b].entity_count();
      max_nodes *= glob_blocks[b].nodesPerElement;
      std::vector<INT> block_linkage(max_nodes);

      // Initialize attributes list, if it exists
      std::vector<T> attributes(glob_blocks[b].attributeCount * glob_blocks[b].entity_count());

      int    error      = 0;
      size_t part_count = local_mesh.size();
      for (size_t p = 0; p < part_count; p++) {
        Excn::ExodusFile id(p);

        size_t global_pos;
        size_t global_block_pos;

        if (blocks[p][b].entity_count() > 0) { // non-zero length block

          if (debug_level & 4) {
            fmt::print("#");
          }
          size_t maximum_nodes = blocks[p][b].entity_count();
          maximum_nodes *= blocks[p][b].nodesPerElement;
          std::vector<INT> local_linkage(maximum_nodes);

          INT bid = blocks[p][b].id;
          error   = ex_get_conn(id, EX_ELEM_BLOCK, bid, Data(local_linkage), nullptr, nullptr);
          if (error < 0) {
            fmt::print(stderr,
                       "ERROR: Cannot get element block connectivity for block {} on part {}.\n",
                       bid, p);
          }
          size_t pos                     = 0;
          size_t goffset                 = glob_blocks[b].offset_;
          size_t element_count           = blocks[p][b].entity_count();
          size_t boffset                 = blocks[p][b].offset_;
          size_t npe                     = blocks[p][b].nodesPerElement;
          INT   *part_loc_elem_to_global = Data(local_mesh[p].localElementToGlobal);
          INT   *part_loc_node_to_global = Data(local_mesh[p].localNodeToGlobal);

          for (size_t e = 0; e < element_count; e++) {
            global_block_pos = part_loc_elem_to_global[(e + boffset)] - goffset;
            global_pos       = global_block_pos * npe;

            for (size_t n = 0; n < npe; n++) {
              size_t node = part_loc_node_to_global[local_linkage[pos++] - 1];
              if (debug_level & 4) {
                SMART_ASSERT(block_linkage[global_pos] == (int)node + 1 ||
                             block_linkage[global_pos] == 0);
              }
              block_linkage[global_pos++] = node + 1;
            }
          }

          // Get attributes list,  if it exists
          if (blocks[p][b].attributeCount > 0) {
            SMART_ASSERT(blocks[p][b].attributeCount == glob_blocks[b].attributeCount)
            (p)(b)(blocks[p][b].attributeCount)(glob_blocks[b].attributeCount);

            size_t         max_attr = blocks[p][b].entity_count() * blocks[p][b].attributeCount;
            std::vector<T> local_attr(max_attr);

            error += ex_get_attr(id, EX_ELEM_BLOCK, blocks[p][b].id, Data(local_attr));

            pos = 0;

            size_t att_count = blocks[p][b].attributeCount;
            for (size_t e = 0; e < element_count; e++) {
              // global_pos is global position within this element block...
              global_block_pos = local_mesh[p].localElementToGlobal[(e + boffset)] - goffset;
              global_pos       = global_block_pos * att_count;
              for (size_t n = 0; n < att_count; n++) {
                attributes[global_pos++] = local_attr[pos++];
              }
            }
          }
        } // end if blocks[p][b].entity_count() (non-zero length block)
        else if (debug_level & 4) {
          fmt::print(".");
        }
      } // end for p=0..part_count-1

      // Verify that connectivity has been set for all elements...
      for (size_t i = 0; i < block_linkage.size(); i++) {
        SMART_ASSERT(block_linkage[i] > 0)(block_linkage[i])(i)(b)(glob_blocks[b].id);
      }

      // Write out block info
      int id_out = Excn::ExodusFile::output(); // output file identifier

      if (!block_linkage.empty()) {
        error += ex_put_conn(id_out, EX_ELEM_BLOCK, glob_blocks[b].id, Data(block_linkage), nullptr,
                             nullptr);
      }

      // Write out attributes list if it exists
      if (glob_blocks[b].attributeCount > 0) {
        error += ex_put_attr(id_out, EX_ELEM_BLOCK, glob_blocks[b].id, Data(attributes));
      } // end for b=0..global_num_blocks-1
      if (debug_level & 4) {
        fmt::print("\n");
      }
    }
    fmt::print("\n");
  }

  template <typename INT>
  void build_reverse_element_map(std::vector<Excn::Mesh<INT>>          &local_mesh,
                                 std::vector<std::vector<Excn::Block>> &blocks,
                                 std::vector<Excn::Block> &glob_blocks, Excn::Mesh<INT> *global,
                                 size_t                               part_count,
                                 std::vector<std::pair<INT, size_t>> &global_element_map)
  {
    int error = 0;
    // Create the map that maps from a local part element to the
    // global map. This combines the mapping local part element to
    // 'global id' and then 'global id' to global position. The
    // mapping is now a direct lookup instead of a lookup followed by
    // a reverse map.
    //

    // We iterate through each element block a part at a time to
    // determine which elements are in each block. Note that an
    // element will possibly exist in multiple parts, so we need to do
    // the sort/uniquify/shrink on each element block.  The operations
    // for each element block will be:
    // 1. get elements for the block for each part into a vector.
    // 2. sort/uniquify/shrink that vector and copy into the
    // 'global_element_map' in the positions following the previous block.
    //
    // This is similar to what was done above in a global sense,
    // except it is done an element block at a time. The only need for
    // doing them both is that it lets us detect whether an element id
    // was reused and appears in both blocks.  Currently, just error
    // out if that happens...  In the future, change the id to an
    // unused value and continue....

    size_t                                           tot_size = 0;
    std::vector<std::vector<std::pair<INT, size_t>>> global_element_numbers(part_count);
    for (size_t p = 0; p < part_count; p++) {
      Excn::ExodusFile id(p);
      global_element_numbers[p].resize(local_mesh[p].count(Excn::ObjectType::ELEM));
      std::vector<INT> ids(local_mesh[p].count(Excn::ObjectType::ELEM));
      error += ex_get_id_map(id, EX_ELEM_MAP, Data(ids));
      for (size_t i = 0; i < local_mesh[p].count(Excn::ObjectType::ELEM); i++) {
        global_element_numbers[p][i] = std::make_pair(ids[i], size_t(0));
      }
      tot_size += local_mesh[p].count(Excn::ObjectType::ELEM);
    }
    global_element_map.resize(tot_size);

    size_t goffset = 0;
    for (auto &glob_block : glob_blocks) {

      size_t block_size = 0;
      for (size_t p = 0; p < part_count; p++) {
        size_t lb = 0;
        for (; lb < blocks[0].size(); lb++) {
          if (glob_block.id == blocks[p][lb].id) {
            break;
          }
        }
        SMART_ASSERT(glob_block.id == blocks[p][lb].id);
        block_size += blocks[p][lb].entity_count();
      }
      std::vector<std::pair<INT, size_t>> block_element_map(block_size);

      size_t poffset = 0;
      for (size_t p = 0; p < part_count; p++) {
        Excn::ExodusFile id(p);

        size_t lb = 0;
        for (; lb < blocks[0].size(); lb++) {
          if (glob_block.id == blocks[p][lb].id) {
            break;
          }
        }
        SMART_ASSERT(glob_block.id == blocks[p][lb].id);

        // Get connectivity array for this element block...
        size_t element_count  = blocks[p][lb].entity_count();
        size_t nodes_per_elem = blocks[p][lb].nodesPerElement;
        size_t maximum_nodes  = element_count * nodes_per_elem;

        std::vector<INT> local_linkage(maximum_nodes);

        size_t bid = blocks[p][lb].id;
        error      = ex_get_conn(id, EX_ELEM_BLOCK, bid, Data(local_linkage), nullptr, nullptr);
        if (error < 0) {
          fmt::print(stderr,
                     "ERROR: Cannot get element block connectivity for block {} on part {}.\n", bid,
                     p);
        }

        // Convert connectivity to global node numbers.
        for (size_t i = 0; i < element_count * nodes_per_elem; i++) {
          INT local_node   = local_linkage[i];
          INT global_node  = local_mesh[p].localNodeToGlobal[local_node - 1] + 1;
          local_linkage[i] = global_node;
        }

        // Have element global ids for all elements in this block,
        // and connectivity.  Can now create out "eleminfo" for these elements.
        size_t con_offset = 0;
        size_t boffset    = blocks[p][lb].offset_;
        for (size_t i = 0; i < element_count; i++) {
          size_t adler_crc = adler(0, &local_linkage[con_offset], nodes_per_elem * sizeof(int));
          global_element_numbers[p][boffset + i].second = adler_crc;
          block_element_map[poffset + i]                = global_element_numbers[p][boffset + i];
          con_offset += nodes_per_elem;
        }
        poffset += element_count;
      }

      // Sort, uniquify, shrink 'block_element_map' and the result
      // then contains the list of elements in this block...
      uniquify(block_element_map);

      size_t block_total_num_elements = block_element_map.size();
      glob_block.elementCount         = block_total_num_elements;
      glob_block.offset_              = goffset;

      // Copy into the global_element_map...
      std::copy(block_element_map.begin(), block_element_map.end(), &global_element_map[goffset]);
      goffset += block_total_num_elements;
    }

    global->elementCount = goffset;
    global_element_map.resize(goffset);

    size_t max_id        = global_element_map[global->elementCount - 1].first;
    bool   is_contiguous = max_id == global_element_map.size();
    // fmt::print("Element id map {}.\n", (is_contiguous ? "is" : "is not"));

    // The global_element_map may or may not be globally sorted; however, each
    // block is sorted, so if we do the iteration by blocks, we can
    // use lower_bound instead of doing global searches...
    for (auto &glob_block : glob_blocks) {

      auto gm_begin = global_element_map.begin() + glob_block.offset_;
      auto gm_end   = gm_begin + glob_block.elementCount;
      auto cur_pos  = gm_begin;
      for (size_t p = 0; p < part_count; p++) {
        size_t lb = 0;
        for (; lb < blocks[0].size(); lb++) {
          if (glob_block.id == blocks[p][lb].id) {
            break;
          }
        }
        SMART_ASSERT(glob_block.id == blocks[p][lb].id);

        size_t element_count = blocks[p][lb].entity_count();
        size_t boffset       = blocks[p][lb].offset_;
        for (size_t i = 0; i < element_count; i++) {
          std::pair<INT, size_t> global_element = global_element_numbers[p][boffset + i];

          if (cur_pos == gm_end || *cur_pos != global_element) {
            auto iter = std::lower_bound(gm_begin, gm_end, global_element);
            SMART_ASSERT(iter != gm_end);
            cur_pos = iter;
          }
          size_t element_value                            = cur_pos - gm_begin;
          local_mesh[p].localElementToGlobal[i + boffset] = element_value + glob_block.offset_;
          ++cur_pos;
        }
      }
    }

    // Update the element ids to give a unique, non-repeating set.  If
    // contiguous, then there is nothing to do.  If not contiguous,
    // then need to determine if there are any repeats (id reuse) and
    // if so, generate a new id for the repeated uses.  Note that
    // there is a possibility that elements in two or more element
    // blocks will have the same element id, so we generate a vector
    // containing a pair<id, position_in_global_element_map>, sort it,
    // and then use it to detect duplicates and map them to a new id.

    if (!is_contiguous) {
      std::vector<std::pair<size_t, size_t>> id_pos(global_element_map.size());
      for (size_t i = 0; i < global->elementCount; i++) {
        id_pos[i].first  = global_element_map[i].first;
        id_pos[i].second = i;
      }
#if USE_STD_SORT
      std::sort(id_pos.begin(), id_pos.end());
#else
      pdqsort(id_pos.begin(), id_pos.end());
#endif
      max_id = id_pos.back().first;
      // Check again for contiguous ids since we now have a sorted list...
      is_contiguous = max_id == global_element_map.size();

      if (!is_contiguous) {
        int    repeat_found = 0;
        size_t id_last      = id_pos[0].first;

        for (size_t i = 1; i < global->elementCount; i++) {
          if (id_pos[i].first == id_last) {
            global_element_map[id_pos[i].second].first = ++max_id;
            repeat_found++;
          }
          else {
            id_last = id_pos[i].first;
          }
        }
        if (repeat_found > 0) {
          fmt::print(stderr,
                     "WARNING: {} duplicate element ids were found. Their ids have been "
                     "renumbered to remove duplicates.\n",
                     repeat_found);
        }
      }
    }
  }

  template <typename INT>
  void build_reverse_node_map(std::vector<Excn::Mesh<INT>> &local_mesh, Excn::Mesh<INT> *global,
                              size_t part_count, GlobalMap &global_node_map)
  {
    // Append all local node maps to the global node map.
    // Sort the global node map
    // Remove duplicates.
    // Position within map is now the map...
    // When building the local-part node to global id, use binary_search...

    // Global node map and count.
    std::vector<std::vector<NodeInfo>> global_nodes(part_count);

    size_t tot_size = 0;
    for (size_t p = 0; p < part_count; p++) {
      tot_size += local_mesh[p].count(Excn::ObjectType::NODE);
      global_nodes[p].resize(local_mesh[p].count(Excn::ObjectType::NODE));
    }
    global_node_map.resize(tot_size);

    size_t offset = 0;
    for (size_t p = 0; p < part_count; p++) {
      std::vector<double> x(local_mesh[p].count(Excn::ObjectType::NODE));
      std::vector<double> y(local_mesh[p].count(Excn::ObjectType::NODE));
      std::vector<double> z(local_mesh[p].count(Excn::ObjectType::NODE));
      std::vector<INT>    nid(local_mesh[p].count(Excn::ObjectType::NODE));

      Excn::ExodusFile id(p);
      ex_get_id_map(id, EX_NODE_MAP, Data(nid));
      ex_get_coord(id, Data(x), Data(y), Data(z));
      for (size_t i = 0; i < local_mesh[p].count(Excn::ObjectType::NODE); i++) {
        global_nodes[p][i] = NodeInfo(nid[i], x[i], y[i], z[i]);
      }
      std::copy(global_nodes[p].begin(), global_nodes[p].end(), &global_node_map[offset]);
      offset += local_mesh[p].count(Excn::ObjectType::NODE);
    }
    // Now, sort the global_node_map array and remove duplicates...
    uniquify(global_node_map);

    global->nodeCount = global_node_map.size();

    // See whether the node numbers are contiguous.  If so, we can map
    // the nodes back to their original location. Since the nodes are
    // sorted and there are no duplicates, we just need to see if the id
    // at global_node_map.size() == global_node_map.size();
    INT  max_id        = global_node_map[global->nodeCount - 1].id;
    bool is_contiguous = (int64_t)max_id == static_cast<int64_t>(global_node_map.size());
    fmt::print("Node map {} contiguous.\n", (is_contiguous ? "is" : "is not"));

    // Create the map that maps from a local part node to the
    // global map. This combines the mapping local part node to
    // 'global id' and then 'global id' to global position. The
    // mapping is now a direct lookup instead of a lookup followed by
    // a reverse map.
    if (is_contiguous) {
      for (size_t p = 0; p < part_count; p++) {
        size_t node_count = local_mesh[p].count(Excn::ObjectType::NODE);
        for (size_t i = 0; i < node_count; i++) {
          const NodeInfo &global_node        = global_nodes[p][i];
          local_mesh[p].localNodeToGlobal[i] = global_node.id - 1;
        }
      }
    }
    else {
      auto cur_pos = global_node_map.begin();
      for (size_t p = 0; p < part_count; p++) {
        size_t node_count = local_mesh[p].count(Excn::ObjectType::NODE);
        for (size_t i = 0; i < node_count; i++) {
          NodeInfo global_node = global_nodes[p][i];

          if (cur_pos == global_node_map.end() || *cur_pos != global_node) {
            auto iter =
                std::lower_bound(global_node_map.begin(), global_node_map.end(), global_node);
            if (iter == global_node_map.end()) {
              NodeInfo n = global_node;
              fmt::print(stderr,
                         "ERROR: Bad Node in build_reverse_node_map: {}\tat location: {}\t{}\t{}\n",
                         n.id, n.x, n.y, n.z);
              exit(EXIT_FAILURE);
            }
            cur_pos = iter;
          }
          size_t nodal_value                 = cur_pos - global_node_map.begin();
          local_mesh[p].localNodeToGlobal[i] = nodal_value;
          ++cur_pos;
        }
      }
    }

    // Update the nodal ids to give a unique, non-repeating set.  If contiguous, then
    // there is nothing to do.  If not contiguous, then need to determine if there are any
    // repeats (id reuse) and if so, generate a new id for the repeated uses.
    // A duplicate id would have the same id, but different x y z position.
    if (!is_contiguous) {
      bool   repeat_found = false;
      size_t id_last      = global_node_map[0].id;
      for (size_t i = 1; i < global->nodeCount; i++) {
        if (global_node_map[i].id == id_last) {
          global_node_map[i].id = ++max_id;
          repeat_found          = true;
        }
        else {
          id_last = global_node_map[i].id;
        }
      }
      if (repeat_found) {
        fmt::print(
            stderr,
            "WARNING: Duplicate node ids were found. Their ids have been renumbered to remove "
            "duplicates. If the part meshes should be identical, maybe use the "
            "--ignore_coordinate option.\n");
      }
    }
  }

  template <typename INT>
  void build_reverse_node_map(std::vector<Excn::Mesh<INT>> &local_mesh, Excn::Mesh<INT> *global,
                              size_t part_count, std::vector<INT> &global_node_map)
  {
    // Append all local node maps to the global node map.
    // Sort the global node map
    // Remove duplicates.
    // Position within map is now the map...
    // When building the local-part node to global id, use binary_search...

    // Global node map and count.
    std::vector<std::vector<INT>> global_nodes(part_count);

    size_t tot_size = 0;
    for (size_t p = 0; p < part_count; p++) {
      tot_size += local_mesh[p].count(Excn::ObjectType::NODE);
      global_nodes[p].resize(local_mesh[p].count(Excn::ObjectType::NODE));
    }
    global_node_map.resize(tot_size);

    size_t offset = 0;
    for (size_t p = 0; p < part_count; p++) {
      Excn::ExodusFile id(p);
      ex_get_id_map(id, EX_NODE_MAP, Data(global_nodes[p]));
      std::copy(global_nodes[p].begin(), global_nodes[p].end(), &global_node_map[offset]);
      offset += local_mesh[p].count(Excn::ObjectType::NODE);
    }

    // Now, sort the global_node_map array and remove duplicates...
    uniquify(global_node_map);

    global->nodeCount = global_node_map.size();

    // See whether the node numbers are contiguous.  If so, we can map
    // the nodes back to their original location. Since the nodes are
    // sorted and there are no duplicates, we just need to see if the id
    // at global_node_map.size() == global_node_map.size();
    INT  max_id        = global_node_map[global->nodeCount - 1];
    bool is_contiguous = (int64_t)max_id == static_cast<int64_t>(global_node_map.size());
    fmt::print("Node map {} contiguous.\n", (is_contiguous ? "is" : "is not"));

    // Create the map that maps from a local part node to the
    // global map. This combines the mapping local part node to
    // 'global id' and then 'global id' to global position. The
    // mapping is now a direct lookup instead of a lookup followed by
    // a reverse map.
    if (is_contiguous) {
      for (size_t p = 0; p < part_count; p++) {
        size_t node_count = local_mesh[p].count(Excn::ObjectType::NODE);
        for (size_t i = 0; i < node_count; i++) {
          INT global_node                    = global_nodes[p][i];
          local_mesh[p].localNodeToGlobal[i] = global_node - 1;
        }
      }
    }
    else {
      auto cur_pos = global_node_map.begin();
      for (size_t p = 0; p < part_count; p++) {
        size_t node_count = local_mesh[p].count(Excn::ObjectType::NODE);
        for (size_t i = 0; i < node_count; i++) {
          INT global_node = global_nodes[p][i];

          if (cur_pos == global_node_map.end() || *cur_pos != global_node) {
            auto iter =
                std::lower_bound(global_node_map.begin(), global_node_map.end(), global_node);
            if (iter == global_node_map.end()) {
              INT n = global_node;
              fmt::print(stderr, "ERROR: Bad Node in build_reverse_node_map: {}\n", n);
              exit(EXIT_FAILURE);
            }
            cur_pos = iter;
          }
          size_t nodal_value                 = cur_pos - global_node_map.begin();
          local_mesh[p].localNodeToGlobal[i] = nodal_value;
          ++cur_pos;
        }
      }
    }
  }

  void get_put_variable_names(int id, int out, Excn::Variables &vars, Excn::SystemInterface &si,
                              int *combined_status_variable_index)
  {
    if (vars.count(Excn::InOut::OUT_) > 0) {

      char **output_name_list =
          get_name_array(vars.count(Excn::InOut::OUT_), Excn::ExodusFile::max_name_length());

      int num_vars       = vars.index_.size();
      int extra          = vars.addStatus ? 1 : 0;
      int num_input_vars = num_vars - extra;

      char **input_name_list = get_name_array(num_vars, Excn::ExodusFile::max_name_length());
      if (num_input_vars > 0) {
        int error = ex_get_variable_names(id, vars.type(), num_input_vars, input_name_list);
        if (error != EX_NOERR) {
          fmt::print(stderr, "ERROR: Cannot get {} variable names\n", vars.label());
          exit(EXIT_FAILURE);
        }
      }

      std::string status;
      if (vars.type() == EX_ELEM_BLOCK || vars.type() == EX_NODAL) {
        if (vars.type() == EX_ELEM_BLOCK) {
          status = si.element_status_variable();
          if (status != "NONE") {
            copy_string(input_name_list[num_vars - 1], status,
                        Excn::ExodusFile::max_name_length() + 1);
          }
        }
        else if (vars.type() == EX_NODAL) {
          status = si.nodal_status_variable();
          if (status != "NONE") {
            copy_string(input_name_list[num_vars - 1], status,
                        Excn::ExodusFile::max_name_length() + 1);
          }
        }
      }

      // Iterate through the 'var_index' and transfer
      // Assume that the number of pointers is limited to
      // the number of results variables
      size_t maxlen = 0;
      for (int i = 0; i < num_vars; i++) {
        if (vars.index_[i] > 0) {
          copy_string(output_name_list[vars.index_[i] - 1], input_name_list[i],
                      Excn::ExodusFile::max_name_length() + 1);
          if (strlen(input_name_list[i]) > maxlen) {
            maxlen = strlen(input_name_list[i]);
          }
        }
      }

      // See if any of the variable names conflict with the status variable name...
      if (status != "NONE") {
        for (size_t i = 0; i < static_cast<size_t>(vars.count(Excn::InOut::OUT_)) - 1; i++) {
          if (case_compare(output_name_list[i], status)) {
            // Error -- duplicate element variable names on output database.
            fmt::print(stderr,
                       "\nERROR: A {} variable already exists on the input database with the "
                       "same name as the status variable '{}'. This is not allowed.\n\n",
                       vars.label(), status);
            exit(EXIT_FAILURE);
          }
        }
      }

      maxlen += 2;
      // Assume 8 characters for initial tab...
      int width  = si.screen_width();
      int nfield = (width - 8) / maxlen;
      if (nfield < 1) {
        nfield = 1;
      }

      fmt::print("Found {} {} variables.\n\t", vars.count(Excn::InOut::OUT_), vars.label());
      {
        int i    = 0;
        int ifld = 1;
        while (i < vars.count(Excn::InOut::OUT_)) {
          fmt::print("{:<{}}", output_name_list[i++], maxlen);
          if (++ifld > nfield && i < vars.count(Excn::InOut::OUT_)) {
            fmt::print("\n\t");
            ifld = 1;
          }
        }
        fmt::print("\n\n");
      }

      ex_put_variable_names(out, vars.type(), vars.count(Excn::InOut::OUT_), output_name_list);

      // KLUGE: Handle finding combined status index variable here since it is
      // the only place we have the list of output variable names...
      if (vars.type() == EX_ELEM_BLOCK && (combined_status_variable_index != nullptr)) {
        *combined_status_variable_index = 0;
        const std::string &comb_stat    = si.combined_mesh_status_variable();
        if (!comb_stat.empty()) {
          for (int i = 0; i < vars.count(Excn::InOut::OUT_); i++) {
            if (case_compare(comb_stat, output_name_list[i])) {
              *combined_status_variable_index = i + 1;
              break;
            }
          }
        }
      }
      free_name_array(output_name_list, vars.count(Excn::InOut::OUT_));
      free_name_array(input_name_list, num_vars);
    }
  }

  void get_variable_params(int id, Excn::Variables &vars, const StringIdVector &variable_list)
  {
    // Determines the number of variables of type 'type()' that will
    // be written to the output database. The 'variable_list' vector
    // specifies a possibly empty list of variable names that the user
    // wants transferred to the output database. If 'variable_list' is
    // empty, then all variables of that type will be transferred; if
    // the 'variable_list' size is 1 and it contains the string 'NONE',
    // then no variables of that type will be transferred; if size is 1
    // and it contains the string 'ALL', then all variables of that type
    // will be transferred.
    //
    // Returns the number of variables which will be output Also creates
    // a 'var_index'.  The var_index is zero-based and of size
    // 'input_variable_count'. If:
    // var_index[i] ==0, variable not written to output database
    // var_index[i] > 0, variable written; is variable 'var_index[i]'

    // If 'type' is ELEMENT or NODE, then reserve space for the 'status' variable.
    int extra = vars.addStatus ? 1 : 0;
    int num_vars;
    ex_get_variable_param(id, vars.type(), &num_vars);

    vars.index_.resize(num_vars + extra);

    // Create initial index which defaults to no output...
    std::fill(vars.index_.begin(), vars.index_.end(), 0);
    // ...Except for the status variable (if any)
    if (extra == 1) {
      vars.index_[num_vars] = 1;
    }

    // If 'variable_list' is empty or specified 'ALL', then all
    // variables are to be output
    if (variable_list.empty() ||
        (variable_list.size() == 1 && case_compare(variable_list[0].first, "all"))) {
      std::iota(vars.index_.begin(), vars.index_.end(), 1);
      vars.outputCount = num_vars + extra;
      return;
    }

    // Another possibility is user specifies "NONE" for the variable
    // list so no variables will be written.  Just return 0.
    if (variable_list.size() == 1 && case_compare(variable_list[0].first, "none")) {
      vars.outputCount = extra;
      return;
    }

    // At this point, the variable_list specifies at least one
    // variable to be output to the database.
    // Get the list that the user entered and the list of files
    // from the input database...
    {
      StringVector exo_names = get_exodus_variable_names(id, vars.type(), num_vars);

      // Iterate 'variable_list' and find position in 'exo_names'.  If
      // not found, there is an error -- set index to -1; otherwise set
      // index to position (1-based) in variable_list.  Others are
      // already set to '0' by default initialization.

      // The variable_list may contain multiple entries for each
      // variable if the user is specifying output only on certain
      // element blocks...
      std::string var_name;
      int         var_count = 0;
      for (const auto &elem : variable_list) {
        if (var_name == elem.first) {
          continue;
        }
        var_name   = elem.first;
        bool found = false;
        for (size_t j = 0; j < exo_names.size() && !found; j++) {
          if (case_compare(exo_names[j], var_name)) {
            found          = true;
            vars.index_[j] = ++var_count;
          }
        }
        if (!found) {
          fmt::print(stderr, "ERROR: Variable '{}' is not valid.\n", elem.first);
          exit(EXIT_FAILURE);
        }
      }
      // Count non-zero entries in var_index;
      int nz_count = 0;
      for (auto &elem : vars.index_) {
        if (elem > 0) {
          nz_count++;
        }
      }
      SMART_ASSERT(nz_count == var_count + extra)(nz_count)(var_count);

      if (vars.addStatus) {
        vars.index_[num_vars] = nz_count; // Already counted above...
      }
      vars.outputCount = nz_count;
      return;
    }
  }

  bool check_variable_params(size_t p, const Excn::Variables &vars)
  {
    // Determines the number of variables of type 'type()' that will
    // be written to the output database. The 'variable_list' vector
    // specifies a possibly empty list of variable names that the user
    // wants transferred to the output database. If 'variable_list' is
    // empty, then all variables of that type will be transferred; if
    // the 'variable_list' size is 1 and it contains the string 'NONE',
    // then no variables of that type will be transferred; if size is 1
    // and it contains the string 'ALL', then all variables of that type
    // will be transferred.
    //
    // Returns the number of variables which will be output Also creates
    // a 'var_index'.  The var_index is zero-based and of size
    // 'input_variable_count'. If:
    // var_index[i] ==0, variable not written to output database
    // var_index[i] > 0, variable written; is variable 'var_index[i]'

    // If 'type' is ELEMENT or NODE, then reserve space for the 'status' variable.
    int  extra = vars.addStatus ? 1 : 0;
    int  num_vars;
    auto id = Excn::ExodusFile(p);
    ex_get_variable_param(id, vars.type(), &num_vars);
    if ((size_t)num_vars != vars.index_.size() - extra) {
      fmt::print("ERROR: Part mesh {} has a different number of {} variables ({}) than the root "
                 "part mesh ({}) which is not allowed.\n",
                 p, vars.label(), num_vars, vars.index_.size() - extra);
      return true;
    }
    return false;
  }

  template <typename INT> void put_mesh_summary(const Excn::Mesh<INT> &mesh)
  {
    // Write out Mesh info
    fmt::print(" Title: {}\n\n", mesh.title);
    fmt::print(" Number of coordinates per node ={:14}\n",
               fmt::group_digits(mesh.count(Excn::ObjectType::DIM)));
    fmt::print(" Number of nodes                ={:14}\n",
               fmt::group_digits(mesh.count(Excn::ObjectType::NODE)));
    fmt::print(" Number of elements             ={:14}\n",
               fmt::group_digits(mesh.count(Excn::ObjectType::ELEM)));
    fmt::print(" Number of element blocks       ={:14}\n",
               fmt::group_digits(mesh.count(Excn::ObjectType::EBLK)));
    fmt::print(" Number of nodal point sets     ={:14}\n",
               fmt::group_digits(mesh.count(Excn::ObjectType::NSET)));
    fmt::print(" Number of element side sets    ={:14}\n",
               fmt::group_digits(mesh.count(Excn::ObjectType::SSET)));
  }

  template <typename INT>
  void get_nodesets(size_t total_node_count, std::vector<Excn::Mesh<INT>> &local_mesh,
                    std::vector<std::vector<Excn::NodeSet<INT>>> &nodesets,
                    std::vector<Excn::NodeSet<INT>>              &glob_sets)
  {
    // Find number of nodesets in the global model...
    std::set<INT>    set_ids;
    std::vector<INT> ids;

    size_t part_count = local_mesh.size();
    int    bad_ns     = 0;
    {
      int ns_count;
      for (size_t p = 0; p < part_count; p++) {
        Excn::ExodusFile id(p);
        ns_count = ex_inquire_int(id, EX_INQ_NODE_SETS);

        // Get the ids for these
        ids.resize(ns_count);
        ex_get_ids(id, EX_NODE_SET, Data(ids));

        for (int iset = 0; iset < ns_count; iset++) {
          if (ids[iset] != 0) {
            set_ids.insert(ids[iset]);
          }
          else {
            bad_ns++;
          }
        }
      }
    }

    if (bad_ns != 0) {
      fmt::print(stderr,
                 "ERROR: There were {} nodesets (counting all files) which had an id equal to "
                 "0 which is not allowed.\n",
                 bad_ns);
    }

    if (set_ids.empty()) {
      return;
    }

    // set_ids now contains the union of all nodeset ids...
    glob_sets.resize(set_ids.size());
    size_t gnset_size = set_ids.size();
    ids.resize(gnset_size);

    {
      size_t i = 0;
      for (auto &set_id : set_ids) {
        glob_sets[i].id        = set_id;
        glob_sets[i].position_ = i;
        i++;
      }
    }

    {
      for (size_t p = 0; p < part_count; p++) {
        Excn::ExodusFile id(p);

        nodesets[p].resize(set_ids.size());

        // Get the ids again so we can map current order back to file order...
        ex_get_ids(id, EX_NODE_SET, Data(ids));
        int ns_count = ex_inquire_int(id, EX_INQ_NODE_SETS);

        for (int i = 0; i < ns_count; i++) {
          nodesets[p][i].id = ids[i];

          // Find which global nodeset this id corresponds to...
          size_t gi = gnset_size;
          for (size_t j = 0; j < gnset_size; j++) {
            if (ids[i] == glob_sets[j].id) {
              gi = j;
              break;
            }
          }
          SMART_ASSERT(gi != gnset_size);
          glob_sets[gi].position_  = i;
          nodesets[p][i].position_ = gi;

          // Get the parameters for this nodeset...
          if (ex_int64_status(id) & EX_BULK_INT64_API) {
            int64_t node_count;
            int64_t df_count;
            ex_get_set_param(id, EX_NODE_SET, nodesets[p][i].id, &node_count, &df_count);
            nodesets[p][i].nodeCount = node_count;
            nodesets[p][i].dfCount   = df_count;
          }
          else {
            int node_count;
            int df_count;
            ex_get_set_param(id, EX_NODE_SET, nodesets[p][i].id, &node_count, &df_count);
            nodesets[p][i].nodeCount = node_count;
            nodesets[p][i].dfCount   = df_count;
          }

          std::vector<char> name(Excn::ExodusFile::max_name_length() + 1);
          ex_get_name(id, EX_NODE_SET, nodesets[p][i].id, Data(name));
          if (name[0] != '\0') {
            nodesets[p][i].name_ = Data(name);
            if (glob_sets[gi].name_.empty()) {
              glob_sets[gi].name_ = Data(name);
            }
          }

          if (debug_level & 32) {
            fmt::print("Part {} ", p + 1);
            nodesets[p][i].dump();
          }
        }
      }
    }

    verify_set_position_mapping("nodeset", part_count, glob_sets, nodesets);

    {
      // Now get the nodeset nodes and df.
      // Currently ignore the DF.  Could add if people need it...

      // This is inefficient since the part loop is on
      // the inside...  The other ordering would use more memory...

      std::vector<INT> ns_nodes;
      for (size_t ns = 0; ns < set_ids.size(); ns++) {

        std::vector<INT> glob_ns_nodes(total_node_count + 1);
        std::fill(glob_ns_nodes.begin(), glob_ns_nodes.end(), 0);

        size_t lns = glob_sets[ns].position_;
        for (size_t p = 0; p < part_count; p++) {
          Excn::ExodusFile id(p);

          glob_sets[ns].name_ = nodesets[p][lns].name_;

          size_t size = nodesets[p][lns].entity_count();
          if (size > 0) {
            ns_nodes.resize(size);
            ex_get_set(id, EX_NODE_SET, nodesets[p][lns].id, Data(ns_nodes), nullptr);

            // The node ids are in local space -- map to global
            for (size_t iset = 0; iset < size; iset++) {
              size_t global_node         = local_mesh[p].localNodeToGlobal[ns_nodes[iset] - 1] + 1;
              glob_ns_nodes[global_node] = 1;
            }
          }
        }
        // Count number of nonzero entries and transfer to the
        // output nodeset
        // NOTE: global_node above is 1-based.
        glob_sets[ns].nodeCount =
            std::accumulate(glob_ns_nodes.begin(), glob_ns_nodes.end(), INT(0));
        glob_sets[ns].nodeSetNodes.resize(glob_sets[ns].entity_count());
        glob_sets[ns].dfCount = 0;

        size_t j = 0;
        for (size_t i = 1; i <= total_node_count; i++) {
          if (glob_ns_nodes[i] == 1) {
            glob_sets[ns].nodeSetNodes[j++] = i;
            glob_ns_nodes[i]                = j;
          }
        }
        SMART_ASSERT(j == glob_sets[ns].entity_count())(j)(glob_sets[ns].entity_count())(ns);

        // See if we need a map from local nodeset position to global nodeset position
        // Only needed if there are nodeset variables (or dist factors).
        // Assume all files have same number of variables...
        int num_vars;
        {
          Excn::ExodusFile id(0);
          ex_get_variable_param(id, EX_NODE_SET, &num_vars);
        }
        if (num_vars > 0) {
          for (size_t p = 0; p < part_count; p++) {
            Excn::ExodusFile id(p);
            // Get the nodelist, but store it in nodeOrderMap.
            // global_pos = nodeOrderMap[i]
            Excn::NodeSet<INT> &nset   = nodesets[p][ns];
            size_t              nnodes = nset.entity_count();
            nset.nodeOrderMap.resize(nnodes);
            ex_get_set(id, EX_NODE_SET, nset.id, Data(nset.nodeOrderMap), nullptr);

            for (size_t i = 0; i < nnodes; i++) {
              size_t local_node    = nset.nodeOrderMap[i];                                // 1-based
              size_t global_node   = local_mesh[p].localNodeToGlobal[local_node - 1] + 1; // 1-based
              size_t global_pos    = glob_ns_nodes[global_node];                          // 1-based
              nset.nodeOrderMap[i] = global_pos - 1;
            }
#if 0
            if (debug_level & 32)
              nset.dump_order();
#endif
          }
        }
      }
    }
  }

  template <typename INT> void put_nodesets(std::vector<Excn::NodeSet<INT>> &glob_sets)
  {
    int exoid = Excn::ExodusFile::output();

    if (debug_level & 32) {
      fmt::print("\nOutput NodeSets:\n");
    }
    for (auto &glob_set : glob_sets) {
      ex_put_set(exoid, EX_NODE_SET, glob_set.id, Data(glob_set.nodeSetNodes), nullptr);
      //    ex_put_node_set_dist_fact(exoid, glob_sets[ns].id, &glob_sets[ns].distFactors[0]);
      // Done with the memory; clear out the vector containing the bulk data nodes and distFactors.
      clear(glob_set.nodeSetNodes);
      clear(glob_set.distFactors);

      if (debug_level & 32) {
        glob_set.dump();
      }
    }
  }

  template <typename INT>
  void get_sideset_metadata(std::vector<Excn::Mesh<INT>>                 &local_mesh,
                            std::vector<std::vector<Excn::SideSet<INT>>> &sets,
                            std::vector<Excn::SideSet<INT>>              &glob_ssets)
  {
    // Find number of sidesets in the global model...
    std::set<int>    set_ids;
    std::vector<INT> ids;

    size_t part_count = local_mesh.size();
    int    bad_ss     = 0;
    {
      for (size_t p = 0; p < part_count; p++) {
        Excn::ExodusFile id(p);
        int              ss_count = ex_inquire_int(id, EX_INQ_SIDE_SETS);

        // Get the ids for these
        ids.resize(ss_count);
        ex_get_ids(id, EX_SIDE_SET, Data(ids));

        for (int i = 0; i < ss_count; i++) {
          if (ids[i] != 0) {
            set_ids.insert(ids[i]);
          }
          else {
            bad_ss++;
          }
        }
      }
    }

    if (bad_ss != 0) {
      fmt::print(stderr,
                 "ERROR: There were {} sidesets (counting all files) which had an id equal to 0 "
                 "which is not allowed.\n",
                 bad_ss);
    }

    if (set_ids.empty()) {
      return;
    }

    // set_ids now contains the union of all sideset ids...
    glob_ssets.resize(set_ids.size());
    size_t gsset_size = set_ids.size();
    ids.resize(gsset_size);

    {
      size_t i = 0;
      for (auto set_id : set_ids) {
        glob_ssets[i].id        = set_id;
        glob_ssets[i].position_ = i; // Not used
        i++;
      }
    }

    {
      std::vector<char> name(Excn::ExodusFile::max_name_length() + 1);
      for (size_t p = 0; p < part_count; p++) {
        Excn::ExodusFile id(p);

        sets[p].resize(set_ids.size());

        // Get the ids again so we can map current order back to file order...
        ex_get_ids(id, EX_SIDE_SET, Data(ids));

        int ss_count = ex_inquire_int(id, EX_INQ_SIDE_SETS);
        for (int i = 0; i < ss_count; i++) {
          sets[p][i].id = ids[i];

          // Find which global sideset this id corresponds to...
          size_t gi = gsset_size;
          for (size_t j = 0; j < gsset_size; j++) {
            if (ids[i] == glob_ssets[j].id) {
              gi = j;
              break;
            }
          }
          SMART_ASSERT(gi != gsset_size);
          sets[p][i].position_ = gi;

          // Get the parameters for this sideset...
          ex_set set{};
          set.type                     = EX_SIDE_SET;
          set.id                       = sets[p][i].id;
          set.entry_list               = nullptr;
          set.extra_list               = nullptr;
          set.distribution_factor_list = nullptr;
          int error                    = ex_get_sets(id, 1, &set);
          if (error != EX_NOERR) {
            fmt::print(stderr, "ERROR: Cannot get side set with id {}\n", set.id);
            exit(EXIT_FAILURE);
          }

          sets[p][i].sideCount = set.num_entry;
          sets[p][i].dfCount   = set.num_distribution_factor;
          glob_ssets[gi].sideCount += sets[p][i].entity_count();
          glob_ssets[gi].dfCount += sets[p][i].dfCount;

          ex_get_name(id, EX_SIDE_SET, sets[p][i].id, Data(name));
          if (name[0] != '\0') {
            sets[p][i].name_     = Data(name);
            glob_ssets[gi].name_ = Data(name);
          }
        }
      }

      verify_set_position_mapping("sideset", part_count, glob_ssets, sets);

      // See if we need a map from local sideset position to global sideset position
      // Only needed if there are sideeset variables (or dist factors which are currently ignored).
      // Assume all files have same number of variables...
      bool need_sideset_map = false;
      {
        int              num_vars;
        Excn::ExodusFile id(0);
        ex_get_variable_param(id, EX_SIDE_SET, &num_vars);
        need_sideset_map = num_vars > 0;
      }

      // Now get the sideset elements/sides. Ignoring DF for now...
      for (size_t ss = 0; ss < set_ids.size(); ss++) {
        // This is maximum possible size; will probably be reduced by duplicate elimination...
        using ElemSideMap = std::vector<std::pair<INT, INT>>;
        ElemSideMap elem_side(glob_ssets[ss].sideCount);
        int         ss_id  = glob_ssets[ss].id;
        size_t      offset = 0;

        for (size_t p = 0; p < part_count; p++) {
          for (size_t lss = 0; lss < sets[p].size(); lss++) {
            if (sets[p][lss].position_ == ss) {
              Excn::ExodusFile id(p);
              sets[p][lss].elems.resize(sets[p][lss].sideCount);
              sets[p][lss].sides.resize(sets[p][lss].sideCount);
              ex_get_set(id, EX_SIDE_SET, ss_id, Data(sets[p][lss].elems),
                         Data(sets[p][lss].sides));

              // Add these to the elem_side vector...
              for (size_t i = 0; i < sets[p][lss].sideCount; i++) {
                size_t global_elem =
                    local_mesh[p].localElementToGlobal[sets[p][lss].elems[i] - 1] + 1;
                elem_side[offset + i] =
                    std::make_pair((INT)global_elem, (INT)sets[p][lss].sides[i]);
              }
              offset += sets[p][lss].sideCount;
              break;
            }
          }
        }

        uniquify(elem_side);

        // Set the output sideset definition...
        glob_ssets[ss].sideCount = elem_side.size();
        glob_ssets[ss].dfCount   = 0;
        glob_ssets[ss].elems.resize(elem_side.size());
        glob_ssets[ss].sides.resize(elem_side.size());

        // Populate the global sideset elements and sides...
        for (size_t i = 0; i < elem_side.size(); i++) {
          std::tie(glob_ssets[ss].elems[i], glob_ssets[ss].sides[i]) = elem_side[i];
        }

        if (need_sideset_map) {
          // For this sideset in each part, figure out the mapping for
          // its (elem, side, variable) position into the corresponding
          // global sideset position...

          // Try the lower_bound searching of elem_side for now.  If
          // inefficient, fix later...
          for (size_t p = 0; p < part_count; p++) {
            for (size_t lss = 0; lss < sets[p].size(); lss++) {
              if (sets[p][lss].position_ == ss) {
                sets[p][lss].elemOrderMap.resize(sets[p][lss].sideCount);
                for (size_t i = 0; i < sets[p][lss].sideCount; i++) {
                  size_t global_elem =
                      local_mesh[p].localElementToGlobal[sets[p][lss].elems[i] - 1] + 1;
                  std::pair<INT, INT> es =
                      std::make_pair((INT)global_elem, (INT)sets[p][lss].sides[i]);

                  auto   iter = std::lower_bound(elem_side.begin(), elem_side.end(), es);
                  size_t pos  = iter - elem_side.begin();
                  sets[p][lss].elemOrderMap[i] = pos;
                }
                break;
              }
            }
          }
        }
      }

      // Calculate sideset offset
      for (size_t ss = 0; ss < glob_ssets.size(); ss++) {
        size_t sum = 0;
        for (size_t p = 0; p < part_count; p++) {
          for (size_t lss = 0; lss < sets[p].size(); lss++) {
            if (sets[p][lss].position_ == ss) {
              sets[p][lss].offset_ = sum;
              sum += sets[p][lss].entity_count();

              if (debug_level & 16) {
                fmt::print("Part {} ", p + 1);
                sets[p][lss].dump();
              }
              break;
            }
          }
        }
      }

      // Free some memory which is no longer needed...
      // (Could move up into sideset loop above)
      for (size_t p = 0; p < part_count; p++) {
        for (auto &elem : sets[p]) {
          clear(elem.elems);
          clear(elem.sides);
          clear(elem.distFactors);
        }
      }
    }
  }

  template <typename INT> void get_put_sidesets(std::vector<Excn::SideSet<INT>> &glob_ssets)
  {
    int exoid = Excn::ExodusFile::output(); // output file identifier
    for (auto &glob_sset : glob_ssets) {
      ex_put_set(exoid, EX_SIDE_SET, glob_sset.id, const_cast<INT *>(Data(glob_sset.elems)),
                 const_cast<INT *>(&glob_sset.sides[0]));
      if (glob_sset.dfCount > 0) {
        ex_put_set_dist_fact(exoid, EX_SIDE_SET, glob_sset.id,
                             reinterpret_cast<void *>(Data(glob_sset.distFactors)));
      }
    }

    for (auto &glob_sset : glob_ssets) {
      clear(glob_sset.elems);
      clear(glob_sset.sides);
      clear(glob_sset.distFactors);
    }
  }

  template <typename T, typename INT>
  void add_status_variable(int id_out, const Excn::Mesh<INT> &global,
                           const std::vector<Excn::Block> &blocks,
                           const std::vector<Excn::Block> &glob_blocks,
                           const std::vector<INT> &local_element_to_global, int step, int variable,
                           T alive, int combined_variable_index)
  {
    SMART_ASSERT(sizeof(T) == Excn::ExodusFile::io_word_size());
    std::vector<T> status;
    SMART_ASSERT(alive == 0.0 || alive == 1.0)(alive);

    for (size_t b = 0; b < global.count(Excn::ObjectType::EBLK); b++) {
      status.resize(glob_blocks[b].entity_count());
      std::fill(status.begin(), status.end(), (1.0 - alive));
      size_t boffset       = blocks[b].offset_;
      size_t goffset       = glob_blocks[b].offset_;
      size_t element_count = blocks[b].entity_count();
      for (size_t e = 0; e < element_count; e++) {
        size_t global_block_pos  = local_element_to_global[(e + boffset)] - goffset;
        status[global_block_pos] = alive;
      }

      // If combining this status variable with a mesh status
      // variable, do that now...
      if (combined_variable_index > 0) {
        std::vector<T> mesh_status(glob_blocks[b].entity_count());
        ex_get_var(id_out, step, EX_ELEM_BLOCK, combined_variable_index, glob_blocks[b].id,
                   glob_blocks[b].entity_count(), Data(mesh_status));

        if (alive == 1.0) {
          for (size_t i = 0; i < glob_blocks[b].entity_count(); i++) {
            status[i] = status[i] * mesh_status[i];
          }
        }
        else {
          SMART_ASSERT(alive == 0.0);
          for (size_t i = 0; i < glob_blocks[b].entity_count(); i++) {
            status[i] = 1.0 - ((1.0 - status[i]) * (1.0 - mesh_status[i]));
          }
        }
      }
      ex_put_var(id_out, step, EX_ELEM_BLOCK, variable, glob_blocks[b].id,
                 glob_blocks[b].entity_count(), Data(status));
    }
  }

  StringVector get_exodus_variable_names(int id, ex_entity_type elType, size_t var_count)
  {
    // Allocate space for variable names...
    char **name_list = get_name_array(var_count, Excn::ExodusFile::max_name_length());
    int    error     = ex_get_variable_names(id, elType, var_count, name_list);
    if (error != EX_NOERR) {
      fmt::print(stderr, "ERROR: Cannot get variable names\n");
      exit(EXIT_FAILURE);
    }

    StringVector names(var_count);
    for (size_t j = 0; j < var_count; j++) {
      compress_white_space(name_list[j]);
      names[j] = std::string(name_list[j]);
    }

    free_name_array(name_list, var_count);
    return names;
  }

  template <typename T, typename INT>
  void filter_truth_table(int id, Excn::Mesh<INT> &global, std::vector<T> &glob_blocks,
                          Excn::Variables &vars, const StringIdVector &variable_names)
  {
    // This routine checks the 'variable_names' list to see if the user
    // has restricted the output of certain variables to certain element
    // blocks. If so, then the truth table is modified to match the
    // users request.
    if (variable_names.empty()) {
      return;
    }

    // Check for a non-zero id entry in the variable_names list which
    // will signify a user-specified element block.
    bool found_it = false;
    for (size_t i = 0; i < variable_names.size() && !found_it; i++) {
      if (variable_names[i].second > 0) {
        found_it = true;
      }
    }
    if (!found_it) {
      return;
    }

    // At this point, know that there is at least one block-restricted
    // variable output specification. For each one, modify the global
    // truth table to match the specification.
    StringVector exo_names = get_exodus_variable_names(id, vars.type(), vars.count());

    std::string var_name;
    int         out_position = -1;
    for (auto [v_name, v_blkid] : variable_names) {
      if (v_blkid > 0 && var_name != v_name) {
        var_name = v_name;
        // Find which exodus variable matches this name
        out_position = -1;
        for (size_t j = 0; j < exo_names.size(); j++) {
          if (case_compare(exo_names[j], var_name)) {
            out_position = vars.index_[j] - 1;
            break;
          }
        }
        if (out_position < 0) {
          fmt::print(stderr, "ERROR: Variable '{}' does not exist on any block in this database.\n",
                     v_name);
          exit(EXIT_FAILURE);
        }

        // Set all truth table entries for this variable to negative
        // of current value and then iterate over specified blocks and
        // set those positive.  This way can make sure that the
        // variable truly exists for the block that the user specified.
        found_it = false;
        for (size_t b = 0; b < global.count(vars.objectType); b++) {
          if (glob_blocks[b].id == v_blkid) {
            if (glob_blocks[b].truthTable[out_position] == 0) {
              fmt::print(stderr, "ERROR: Variable '{}' does not exist on block {}.\n", v_name,
                         v_blkid);
              exit(EXIT_FAILURE);
            }
            else {
              found_it = true;
              break;
            }
          }
          else {
            // User does not want this variable output on these blocks...
            glob_blocks[b].truthTable[out_position] = 0;
          }
        }

        if (!found_it) {
          fmt::print(stderr,
                     "ERROR: User-specified block id of {} for variable '{}' does not exist.\n",
                     v_blkid, v_name);
          exit(EXIT_FAILURE);
        }
      }
    }
  }

  template <typename T, typename INT>
  void get_truth_table(Excn::Mesh<INT> &global, std::vector<std::vector<T>> &blocks,
                       std::vector<T> &glob_blocks, Excn::Variables &vars, int debug)
  {
    // read truth table - sum across all parts since many will
    // have values set to zero for zero length blocks the element
    // variable truth table organized as a 2D array:
    // [global.count(Excn::EBLK)][num_elem_vars]

    Excn::ObjectType object_type = vars.objectType;

    if (vars.count(Excn::InOut::OUT_) > 0) {
      // For each input exodus file, get it's truth table and fill
      // in the location in the output truth table...

      size_t part_count = blocks.size();
      for (size_t p = 0; p < part_count; p++) {
        Excn::ExodusFile id(p);

        for (size_t b = 0; b < global.count(object_type); b++) {
          size_t gb = 0;
          for (; gb < global.count(object_type); gb++) {
            if (glob_blocks[gb].id == blocks[p][b].id) {
              break;
            }
          }
          SMART_ASSERT(glob_blocks[gb].id == blocks[p][b].id);

          if (p == 0) {
            glob_blocks[gb].truthTable.resize(vars.count(Excn::InOut::OUT_));
          }

          if (vars.count(Excn::InOut::IN_) > 0) {
            blocks[p][b].truthTable.resize(vars.count(Excn::InOut::IN_));
            ex_get_object_truth_vector(id, vars.type(), blocks[p][b].id,
                                       vars.count(Excn::InOut::IN_), Data(blocks[p][b].truthTable));

            // Find global block corresponding to this block. (Ids match)
            for (int j = 0; j < vars.count(Excn::InOut::IN_); j++) {
              if (vars.index_[j] > 0) {
                glob_blocks[gb].truthTable[vars.index_[j] - 1] += blocks[p][b].truthTable[j];
              }
            }
          }
          if (vars.addStatus) {
            glob_blocks[gb].truthTable[vars.count(Excn::InOut::OUT_) - 1] = 1;
          }
        }
      }

      // reset truth table values that may be greater than 1
      for (size_t b = 0; b < global.count(object_type); b++) {
        for (int j = 0; j < vars.count(Excn::InOut::OUT_); j++) {
          if (glob_blocks[b].truthTable[j] > 0) {
            glob_blocks[b].truthTable[j] = 1;
          }
        }
      }

      if (debug_level & debug) {
        fmt::print("Truth table for {}\t{} variables\t{} sets\n", vars.label(),
                   vars.count(Excn::InOut::OUT_), global.count(object_type));
        for (size_t b = 0; b < global.count(object_type); b++) {
          for (int j = 0; j < vars.count(Excn::InOut::OUT_); j++) {
            fmt::print("{}", glob_blocks[b].truthTable[j]);
          }
          fmt::print("\n");
        }
      }
    }
  }

  bool case_compare(const std::string &s1, const std::string &s2)
  {
    return (s1.size() == s2.size()) &&
           std::equal(s1.begin(), s1.end(), s2.begin(),
                      [](char a, char b) { return std::tolower(a) == std::tolower(b); });
  }

  void add_info_record(char *info_record, int size)
  {
    // Add 'uname' output to the passed in character string.
    // Maximum size of string is 'size' (not including terminating nullptr)
    // This is used as information data in the concatenated results file
    // to help in tracking when/where/... the file was created
    auto info = sys_info("CONJOIN");
    copy_string(info_record, info, size + 1);
  }

  inline bool is_whitespace(char c)
  {
    static char white_space[] = {' ', '\t', '\n', '\r', ',', '\0'};
    return std::strchr(white_space, c) != nullptr;
  }

  void compress_white_space(char *str)
  {
    char *ibuf = str;
    char *obuf = str;

    int i   = 0;
    int cnt = 0;

    // Don't process an empty string.
    if (str == nullptr) {
      return;
    }

    // Skip leading...
    while (*ibuf != 0 && is_whitespace(*ibuf)) {
      ++ibuf;
    }

    while (*ibuf != 0) {
      if (is_whitespace(*ibuf) && cnt > 0) {
        ibuf++;
      }
      else {
        if (!is_whitespace(*ibuf)) {
          cnt = 0;
        }
        else {
          *ibuf = ' ';
          cnt   = 1;
        }
        obuf[i++] = *ibuf++;
      }
    }
    obuf[i--] = '\0';

    // Skip trailing whitespace
    while (i > 0 && is_whitespace(obuf[i])) {
      obuf[i--] = '\0';
    }
  }

  int get_width(size_t max_value)
  {
    // Returns the field width which will accommodate the
    // largest value.
    int width = 0;
    if (max_value >= 10) {
      width = int(log10(static_cast<double>(max_value)));
    }
    return width + 1;
  }

  template <typename T, typename INT>
  void map_element_vars(size_t loffset, size_t goffset, size_t entity_count, std::vector<T> &values,
                        std::vector<T> &global_values, INT *part_loc_elem_to_global)
  {
    // copy values to master element value information
    for (size_t j = 0; j < entity_count; j++) {
      size_t global_block_pos         = part_loc_elem_to_global[(j + loffset)] - goffset;
      global_values[global_block_pos] = values[j];
    }
  }

  template <typename T, typename U>
  void map_sideset_vars(U & /*unused*/, size_t /*unused*/, std::vector<T> & /*unused*/,
                        std::vector<T> & /*unused*/)
  {
    throw std::runtime_error("Internal Error.");
  }

  template <typename INT>
  void map_sideset_vars(Excn::SideSet<INT> &local_set, size_t entity_count,
                        std::vector<double> &values, std::vector<double> &global_values)
  {
    // copy values to master nodeset value information
    for (size_t j = 0; j < entity_count; j++) {
      size_t global_loc = local_set.elemOrderMap[j];
      SMART_ASSERT(global_loc < global_values.size());
      global_values[global_loc] = values[j];
    }
  }

  template <typename INT>
  void map_sideset_vars(Excn::SideSet<INT> &local_set, size_t entity_count,
                        std::vector<float> &values, std::vector<float> &global_values)
  {
    // copy values to master nodeset value information
    for (size_t j = 0; j < entity_count; j++) {
      size_t global_loc = local_set.elemOrderMap[j];
      SMART_ASSERT(global_loc < global_values.size());
      global_values[global_loc] = values[j];
    }
  }

  template <typename T, typename U>
  void map_nodeset_vars(U & /*unused*/, size_t /*unused*/, std::vector<T> & /*unused*/,
                        std::vector<T> & /*unused*/)
  {
    throw std::runtime_error("Internal Error.");
  }

  template <typename INT>
  void map_nodeset_vars(Excn::NodeSet<INT> &local_set, size_t entity_count,
                        std::vector<double> &values, std::vector<double> &global_values)
  {
    // copy values to master nodeset value information
    for (size_t j = 0; j < entity_count; j++) {
      size_t global_loc = local_set.nodeOrderMap[j];
      SMART_ASSERT(global_loc < global_values.size());
      global_values[global_loc] = values[j];
    }
  }

  template <typename INT>
  void map_nodeset_vars(Excn::NodeSet<INT> &local_set, size_t entity_count,
                        std::vector<float> &values, std::vector<float> &global_values)
  {
    // copy values to master nodeset value information
    for (size_t j = 0; j < entity_count; j++) {
      size_t global_loc = local_set.nodeOrderMap[j];
      SMART_ASSERT(global_loc < global_values.size());
      global_values[global_loc] = values[j];
    }
  }

  template <typename T, typename U, typename INT>
  int read_write_master_values(Excn::Variables &vars, const Excn::Mesh<INT> &global,
                               std::vector<U>               &global_sets,
                               std::vector<Excn::Mesh<INT>> &local_mesh,
                               std::vector<std::vector<U>> &local_sets, std::vector<T> &values,
                               size_t p, Excn::ExodusFile &id, int time_step, int time_step_out)
  {
    int error  = 0;
    int id_out = Excn::ExodusFile::output(); // output file identifier

    size_t max_size = 0;
    for (size_t b = 0; b < global.count(vars.objectType); b++) {
      if (max_size < global_sets[b].entity_count()) {
        max_size = global_sets[b].entity_count();
      }
    }
    std::vector<T> master_values(max_size);

    // Only needed for element, but haven't cleaned this up yet...
    INT *part_loc_elem_to_global = Data(local_mesh[p].localElementToGlobal);

    for (int i = 0; i < vars.count(Excn::InOut::IN_); i++) {
      if (vars.index_[i] > 0) {
        int ivar = vars.index_[i] - 1;

        for (size_t b = 0; b < global.count(vars.objectType); b++) {
          size_t lb = 0;
          for (; lb < global.count(vars.objectType); lb++) {
            if (global_sets[b].id == local_sets[p][lb].id) {
              break;
            }
          }
          SMART_ASSERT(global_sets[b].id == local_sets[p][lb].id);

          if (global_sets[b].truthTable[ivar] && local_sets[p][lb].entity_count() > 0) {

            int entity_count = local_sets[p][lb].entity_count();

            if (local_sets[p][lb].truthTable[i] > 0) {
              error += ex_get_var(id, time_step + 1, exodus_object_type(vars.objectType), i + 1,
                                  local_sets[p][lb].id, entity_count, Data(values));

              switch (vars.objectType) {
              case Excn::ObjectType::EBLK:
                map_element_vars(local_sets[p][lb].offset_, global_sets[b].offset_, entity_count,
                                 values, master_values, part_loc_elem_to_global);
                break;

              case Excn::ObjectType::SSET:
                map_sideset_vars(local_sets[p][lb], entity_count, values, master_values);
                break;

              case Excn::ObjectType::NSET:
                map_nodeset_vars(local_sets[p][lb], entity_count, values, master_values);
                break;
              default: break;
              }
            }
            ex_put_var(id_out, time_step_out, exodus_object_type(vars.objectType), ivar + 1,
                       global_sets[b].id, global_sets[b].entity_count(), Data(master_values));
          }
        }
      }
    }
    return error;
  }

  template <typename U>
  void create_output_truth_table(std::vector<U> &global_sets, Excn::Variables &vars,
                                 std::vector<int> &truth_table)
  {
    truth_table.resize(global_sets.size() * vars.count(Excn::InOut::OUT_));
    for (auto &global_set : global_sets) {
      int bout = global_set.position_;
      SMART_ASSERT(bout >= 0);
      for (int j = 0; j < vars.count(Excn::InOut::OUT_); j++) {
        int out_ttable_loc          = (bout * vars.count(Excn::InOut::OUT_)) + j;
        truth_table[out_ttable_loc] = global_set.truthTable[j];
      }
    }
  }

  template <typename INT>
  size_t find_max_entity_count(size_t part_count, std::vector<Excn::Mesh<INT>> &local_mesh,
                               const Excn::Mesh<INT>                        &global,
                               std::vector<std::vector<Excn::Block>>        &blocks,
                               std::vector<std::vector<Excn::NodeSet<INT>>> &nodesets,
                               std::vector<std::vector<Excn::SideSet<INT>>> &sidesets)
  {
    size_t max_ent = local_mesh[0].count(Excn::ObjectType::NODE);
    for (size_t p = 1; p < part_count; p++) {
      if ((size_t)local_mesh[p].count(Excn::ObjectType::NODE) > max_ent) {
        max_ent = local_mesh[p].count(Excn::ObjectType::NODE);
      }
    }

    for (size_t p = 0; p < part_count; p++) {
      for (size_t b = 0; b < global.count(Excn::ObjectType::EBLK); b++) {
        if (blocks[p][b].entity_count() > max_ent) {
          max_ent = blocks[p][b].entity_count();
        }
      }
    }

    // Nodesets...
    for (size_t p = 0; p < part_count; p++) {
      for (size_t b = 0; b < global.count(Excn::ObjectType::NSET); b++) {
        if (nodesets[p][b].entity_count() > max_ent) {
          max_ent = nodesets[p][b].entity_count();
        }
      }
    }

    // Sidesets...
    for (size_t p = 0; p < part_count; p++) {
      for (size_t b = 0; b < global.count(Excn::ObjectType::SSET); b++) {
        if (sidesets[p][b].entity_count() > max_ent) {
          max_ent = sidesets[p][b].entity_count();
        }
      }
    }
    return max_ent;
  }

  void sort_file_times(StringVector &input_files)
  {
    // Sort files based on minimum timestep time
    std::vector<std::pair<double, std::string>> file_time_name;
    file_time_name.reserve(input_files.size());
    for (auto &filename : input_files) {
      float version       = 0.0;
      int   cpu_word_size = sizeof(float);
      int   io_wrd_size   = 0;
      int   exoid = ex_open(filename.c_str(), EX_READ, &cpu_word_size, &io_wrd_size, &version);
      if (exoid < 0) {
        fmt::print(stderr, "ERROR: Cannot open file '{}'\n", filename);
        exit(EXIT_FAILURE);
      }

      int    nts  = ex_inquire_int(exoid, EX_INQ_TIME);
      double time = 0.0;
      if (nts > 0) {
        ex_get_time(exoid, 1, &time);
      }
      file_time_name.emplace_back(time, filename);
      ex_close(exoid);
    }

    std::sort(file_time_name.begin(), file_time_name.end());
    input_files.clear();
    input_files.reserve(file_time_name.size());

    for (const auto &entry : file_time_name) {
      input_files.push_back(entry.second);
    }
  }

} // namespace
