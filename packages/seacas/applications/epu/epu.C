/*
 * Copyright(C) 1999-2024 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
// concatenates EXODUS/GENESIS output from parallel processors to a single file

#include <algorithm>
#include <array>
#include <cfloat>
#include <climits>
#include <cmath>
#include <copy_string_cpp.h>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <fmt/chrono.h>
#include <fmt/ostream.h>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "format_time.h"
#include "open_file_limit.h"
#include "sys_info.h"
#include "time_stamp.h"
#include "vector_data.h"

#define USE_STD_SORT 1
#if !USE_STD_SORT
#include "pdqsort.h"
#endif
#include "hwm.h"

// Enable SMART_ASSERT even in Release mode...
#define SMART_ASSERT_DEBUG_MODE 1
#include "smart_assert.h"

#include <exodusII.h>

using StringVector = std::vector<std::string>;

#include "EP_ExodusEntity.h"
#include "EP_ExodusFile.h"
#include "EP_Internals.h"
#include "EP_ObjectType.h"
#include "EP_SystemInterface.h"
#include "EP_Variables.h"
#include "EP_Version.h"

#if EX_API_VERS_NODOT <= 467
#error "Requires exodusII version 4.68 or later"
#endif

#include "add_to_log.h"

// The main program templated to permit float/double transfer.
template <typename T, typename INT>
int epu(Excn::SystemInterface &interFace, int start_part, int part_count, int cycle, T /* dummy */,
        INT int_size_dummy);

class mpi
{
public:
  mpi(int argc, char *argv[])
  {
#if ENABLE_PARALLEL_EPU
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);           // CHECK: ALLOW MPI_COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &epu_proc_count); // CHECK: ALLOW MPI_COMM_WORLD
#else
    (void)(argc);
    (void)(argv);
#endif
  }

  ~mpi()
  {
#if ENABLE_PARALLEL_EPU
    MPI_Finalize();
#endif
  }

  int rank{0};
  int epu_proc_count{1};
};

using ExodusIdVector                          = std::vector<ex_entity_id>;
template <typename T> using MasterValueVector = std::vector<T>;

extern double seacas_timer();
namespace {
  unsigned int debug_level = 0;
  const double FILL_VALUE  = FLT_MAX;
  int          rank        = 0;
  std::string  tsFormat    = "[{:%H:%M:%S}] ";

  int get_width(int max_value);

  void LOG(const char *message)
  {
    if ((debug_level & 1) != 0u) {
      fmt::print("{}", time_stamp(tsFormat));
    }
    if (rank == 0) {
      fmt::print("{}", message);
    }
  }

  [[noreturn]] void exodus_error(int lineno)
  {
    std::ostringstream errmsg;
    fmt::print(errmsg,
               "Exodus error ({}) {} at line {} in file epu.C. Please report to gdsjaar@sandia.gov "
               "if you need help.",
               exerrval, ex_strerror(exerrval), lineno);

    ex_err(nullptr, nullptr, EX_PRTLASTMSG);
    throw std::runtime_error(errmsg.str());
  }

  template <typename T> void clear(std::vector<T> &vec)
  {
    vec.clear();
    vec.shrink_to_fit();
    SMART_ASSERT(vec.capacity() == 0);
  }

  ex_entity_type exodus_object_type(const Excn::ObjectType &epu_type)
  {
    switch (epu_type) {
    case Excn::ObjectType::EBLK: return EX_ELEM_BLOCK;
    case Excn::ObjectType::SSET: return EX_SIDE_SET;
    case Excn::ObjectType::NSET: return EX_NODE_SET;
    case Excn::ObjectType::EDBLK: return EX_EDGE_BLOCK;
    case Excn::ObjectType::FABLK: return EX_FACE_BLOCK;
    default:
      throw std::runtime_error("Invalid Object Type in exodus_object_type: " +
                               std::to_string(static_cast<int>(epu_type)));
    }
  }

  char **get_name_array(int size, int length)
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

  bool case_compare(const std::string &s1, const std::string &s2);

  template <typename INT>
  void get_id_map(int exoid, ex_entity_type type, ex_inquiry inq_type, std::vector<INT> &ids)
  {
    // Check whether there is a "original_global_id_map" map on
    // the database. If so, use it instead of the "elem_num_map".
    bool map_read  = false;
    int  map_count = ex_inquire_int(exoid, inq_type);
    if (map_count > 0) {
      char **names = get_name_array(map_count, Excn::ExodusFile::max_name_length());
      int    error = ex_get_names(exoid, type, names);
      if (error < 0) {
        exodus_error(__LINE__);
      }

      for (int i = 0; i < map_count; i++) {
        if (case_compare(names[i], "original_global_id_map")) {
          error = ex_get_num_map(exoid, type, i + 1, Data(ids));
          if (error < 0) {
            exodus_error(__LINE__);
          }
          map_read = true;
          break;
        }
      }
      free_name_array(names, map_count);
    }

    if (!map_read) {
      int error = ex_get_id_map(exoid, type, Data(ids));
      if (error < 0) {
        exodus_error(__LINE__);
      }
    }
  }

  template <typename INT> bool is_sequential(std::vector<INT> &map)
  {
    for (size_t i = 0; i < map.size(); i++) {
      if (map[i] != (INT)i + 1) {
        return false;
      }
    }
    return true;
  }

  // SEE: http://lemire.me/blog/2017/04/10/removing-duplicates-from-lists-quickly
  template <typename T> size_t unique(std::vector<T> &out)
  {
    if (out.empty())
      return 0;
    size_t i    = 1;
    size_t pos  = 1;
    T      oldv = out[0];
    for (; i < out.size(); ++i) {
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

  void compress_white_space(char *str);
  void add_info_record(char *info_record, int size);
  void put_global_info(const Excn::Mesh &global);
  void get_put_qa(int id, int id_out);
  void get_put_coordinate_names(int in, int out, int dimensionality);
  void get_put_assemblies(int in, int out, Excn::Mesh &global);

  template <typename T> void get_put_coordinate_frames(int id, int id_out, T float_or_double);

  template <typename T, typename INT>
  void get_put_coordinates(Excn::Mesh &global, int part_count, std::vector<Excn::Mesh> &local_mesh,
                           const std::vector<std::vector<INT>> &local_node_to_global,
                           T                                    float_or_double);

  template <typename T, typename INT>
  void get_coordinates(int id, int dimensionality, size_t num_nodes,
                       const std::vector<std::vector<INT>> &local_node_to_global, int proc,
                       std::vector<T> &x, std::vector<T> &y, std::vector<T> &z);

  template <typename INT>
  void get_put_nodal_communication_map(int                                  part_count,
                                       const std::vector<std::vector<INT>> &local_node_to_global,
                                       const std::vector<int> &processor_map, int output_processor);

  StringVector get_exodus_variable_names(int id, ex_entity_type elType, int var_count);

  template <typename T>
  void filter_truth_table(int id, Excn::Mesh &global, std::vector<T> &glob_blocks,
                          Excn::Variables &vars, const Excn::StringIdVector &variable_names);

  template <typename T>
  void get_truth_table(Excn::Mesh &global, std::vector<T> &glob_blocks,
                       std::vector<Excn::Mesh> &local, Excn::Variables &vars, int debug);

  template <typename U>
  void create_output_truth_table(const Excn::Mesh &global, std::vector<U> &global_sets,
                                 Excn::Variables &vars, std::vector<int> &truth_table);

  template <typename T, typename U, typename INT>
  void read_write_master_values(Excn::Variables &vars, const Excn::Mesh &global,
                                std::vector<U> &global_sets, std::vector<Excn::Mesh> &local_mesh,
                                std::vector<std::vector<U>> &local_sets,
                                MasterValueVector<T> &master_values, std::vector<T> &values,
                                int part_count, int time_step, int time_step_out,
                                const std::vector<std::vector<INT>> &local_entity_to_global);

  void get_variable_params(int id, Excn::Variables &vars,
                           const Excn::StringIdVector &variable_list);

  void get_put_variable_names(int id, int out, Excn::Variables &vars,
                              Excn::SystemInterface &interFace);

  template <typename INT>
  void build_reverse_element_map(std::vector<std::vector<INT>>         &local_element_to_global,
                                 const std::vector<Excn::Mesh>         &local_mesh,
                                 std::vector<std::vector<Excn::Block>> &blocks,
                                 std::vector<Excn::Block> &glob_blocks, Excn::Mesh *global,
                                 int part_count, std::vector<INT> &global_element_map,
                                 bool map_ids);

  template <typename T, typename INT>
  void get_nodesets(int part_count, size_t total_node_count,
                    const std::vector<std::vector<INT>>          &local_node_to_global,
                    std::vector<std::vector<Excn::NodeSet<INT>>> &nodesets,
                    std::vector<Excn::NodeSet<INT>> &glob_sets, T float_or_double);

  template <typename INT>
  void build_reverse_node_map(std::vector<std::vector<INT>> &local_node_to_global,
                              const std::vector<Excn::Mesh> &local_mesh, Excn::Mesh *global,
                              int part_count, std::vector<INT> &global_node_map);

  void get_element_blocks(int part_count, const std::vector<Excn::Mesh> &local_mesh,
                          const Excn::Mesh &global, std::vector<std::vector<Excn::Block>> &blocks,
                          std::vector<Excn::Block> &glob_blocks);

  template <typename T, typename INT>
  void put_element_blocks(int part_count, int start_part,
                          std::vector<std::vector<Excn::Block>> &blocks,
                          std::vector<Excn::Block>              &glob_blocks,
                          const std::vector<std::vector<INT>>   &local_node_to_global,
                          const std::vector<std::vector<INT>>   &local_element_to_global,
                          T                                      float_or_double);

  template <typename T, typename INT>
  void put_element_blocks(int part_count, int start_part,
                          std::vector<std::vector<Excn::Block>> &blocks,
                          std::vector<Excn::Block>              &glob_blocks,
                          const std::vector<std::vector<INT>>   &local_node_to_global,
                          T /* float_or_double */);

  template <typename INT>
  void get_edgeblocks(int part_count, const std::vector<Excn::Mesh> &local_mesh,
                      const Excn::Mesh                               &global,
                      std::vector<std::vector<Excn::EdgeBlock<INT>>> &edgeblocks,
                      std::vector<Excn::EdgeBlock<INT>>              &glob_edgeblocks);
  template <typename T, typename INT>
  void put_edgeblocks(int part_count, int start_part,
                      std::vector<std::vector<Excn::EdgeBlock<INT>>> &edgeblocks,
                      std::vector<Excn::EdgeBlock<INT>>              &glob_edgeblocks,
                      const std::vector<std::vector<INT>>            &local_node_to_global,
                      const std::vector<std::vector<INT>>            &local_element_to_global,
                      T /* float_or_double */);
  template <typename T, typename INT>
  void put_edgeblocks(int part_count, int start_part,
                      std::vector<std::vector<Excn::EdgeBlock<INT>>> &edgeblocks,
                      std::vector<Excn::EdgeBlock<INT>>              &glob_edgeblocks,
                      const std::vector<std::vector<INT>>            &local_node_to_global,
                      T /* float_or_double */);

  template <typename INT>
  void build_reverse_edge_map(std::vector<std::vector<INT>>                  &local_edge_to_global,
                              const std::vector<Excn::Mesh>                  &local_mesh,
                              std::vector<std::vector<Excn::EdgeBlock<INT>>> &edgeblocks,
                              std::vector<Excn::EdgeBlock<INT>>              &glob_edgeblocks,
                              Excn::Mesh *global, int part_count, std::vector<INT> &global_edge_map,
                              bool map_ids);

  template <typename INT>
  void get_faceblocks(int part_count, const std::vector<Excn::Mesh> &local_mesh,
                      const Excn::Mesh                               &global,
                      std::vector<std::vector<Excn::FaceBlock<INT>>> &faceblocks,
                      std::vector<Excn::FaceBlock<INT>>              &glob_faceblocks);
  template <typename T, typename INT>
  void put_faceblocks(int part_count, int start_part,
                      std::vector<std::vector<Excn::FaceBlock<INT>>> &faceblocks,
                      std::vector<Excn::FaceBlock<INT>>              &glob_faceblocks,
                      const std::vector<std::vector<INT>>            &local_node_to_global,
                      const std::vector<std::vector<INT>>            &local_element_to_global,
                      T /* float_or_double */);
  template <typename T, typename INT>
  void put_faceblocks(int part_count, int start_part,
                      std::vector<std::vector<Excn::FaceBlock<INT>>> &faceblocks,
                      std::vector<Excn::FaceBlock<INT>>              &glob_faceblocks,
                      const std::vector<std::vector<INT>>            &local_node_to_global,
                      T /* float_or_double */);

  template <typename INT>
  void build_reverse_face_map(std::vector<std::vector<INT>>                  &local_face_to_global,
                              const std::vector<Excn::Mesh>                  &local_mesh,
                              std::vector<std::vector<Excn::FaceBlock<INT>>> &faceblocks,
                              std::vector<Excn::FaceBlock<INT>>              &glob_faceblocks,
                              Excn::Mesh *global, int part_count, std::vector<INT> &global_face_map,
                              bool map_ids);

  template <typename INT> void put_nodesets(std::vector<Excn::NodeSet<INT>> &glob_sets);

  template <typename INT>
  void get_sideset_metadata(int part_count, std::vector<std::vector<Excn::SideSet<INT>>> &sets,
                            std::vector<Excn::SideSet<INT>> &glob_ssets);

  template <typename INT>
  void
  get_put_sidesets(int part_count, const std::vector<std::vector<INT>> &local_element_to_global,
                   std::vector<std::vector<Excn::SideSet<INT>>> &sets,
                   std::vector<Excn::SideSet<INT>> &glob_ssets, Excn::SystemInterface &interFace);

  template <typename INT>
  void add_processor_map(int id_out, int part_count, int start_part, const Excn::Mesh &global,
                         std::vector<std::vector<Excn::Block>> &blocks,
                         const std::vector<Excn::Block>        &glob_blocks,
                         const std::vector<std::vector<INT>>   &local_element_to_global);

  template <typename T, typename INT>
  void add_processor_variable(int id_out, int part_count, int start_part, const Excn::Mesh &global,
                              std::vector<std::vector<Excn::Block>> &blocks,
                              const std::vector<Excn::Block>        &glob_blocks,
                              const std::vector<std::vector<INT>>   &local_element_to_global,
                              int step, int variable, std::vector<T> &proc);

  template <typename INT>
  size_t find_max_entity_count(int part_count, std::vector<Excn::Mesh> &local_mesh,
                               const Excn::Mesh                               &global,
                               std::vector<std::vector<Excn::Block>>          &blocks,
                               std::vector<std::vector<Excn::NodeSet<INT>>>   &nodesets,
                               std::vector<std::vector<Excn::SideSet<INT>>>   &sidesets,
                               std::vector<std::vector<Excn::EdgeBlock<INT>>> &edgeblocks,
                               std::vector<std::vector<Excn::FaceBlock<INT>>> &faceblocks);

  template <typename INT>
  size_t find_max_global_entity_count(const Excn::Mesh &global, std::vector<Excn::Block> &blocks,
                                      std::vector<Excn::NodeSet<INT>>   &nodesets,
                                      std::vector<Excn::SideSet<INT>>   &sidesets,
                                      std::vector<Excn::EdgeBlock<INT>> &edgeblocks,
                                      std::vector<Excn::FaceBlock<INT>> &faceblocks);

} // namespace

using namespace Excn;

int main(int argc, char *argv[])
{
  mpi my_mpi(argc, argv);
  rank = my_mpi.rank;

  try {
    time_t begin_time = std::time(nullptr);
    SystemInterface::show_version(rank);
    if (rank == 0) {
#if ENABLE_PARALLEL_EPU
      fmt::print("\tParallel Capability Enabled.\n");
#else
      fmt::print("\tParallel Capability Not Enabled.\n");
#endif
    }

    SystemInterface interFace(rank);
    bool            execute = interFace.parse_options(argc, argv);

    if (!execute) {
      return EXIT_SUCCESS;
    }

    // Debug Options: (can be or'd together)
    //   1 -- time stamp
    //   2 -- check nodal variable consistency
    //   4 -- Element Blocks
    //   8 -- Nodes
    //  16 -- Sidesets
    //  32 -- Nodesets
    //  64 -- Edge Blocks
    // 128 -- Face Blocks
    // 256 -- exodus verbose.
    // 512 -- Check consistent global field values between processors
    debug_level = interFace.debug();

    if ((debug_level & 256) != 0U) {
      ex_opts(EX_VERBOSE | EX_DEBUG);
    }
    else {
      ex_opts(0);
    }

    int start_part      = interFace.start_part();
    int processor_count = interFace.processor_count();

    int part_count = interFace.part_count();
    if (part_count <= 1) {
      fmt::print("INFO: Only one processor or part, no concatenation needed.\n");
      return EXIT_SUCCESS;
    }

    int error = 0;
    if (my_mpi.epu_proc_count > 1) {
      interFace.subcycle(my_mpi.epu_proc_count);

      int per_proc = processor_count / my_mpi.epu_proc_count;
      int extra    = processor_count % my_mpi.epu_proc_count;

      part_count = per_proc + (rank < extra ? 1 : 0);

      if (rank < extra) {
        start_part = (per_proc + 1) * rank;
      }
      else {
        start_part = (per_proc + 1) * extra + per_proc * (rank - extra);
      }

      SMART_ASSERT(start_part + part_count <= processor_count);

      if (!ExodusFile::initialize(interFace, start_part, part_count, rank, false)) {
        throw std::runtime_error("ERROR: (EPU) Problem initializing input and/or output files.\n");
      }

      if (ExodusFile::io_word_size() == 4) { // Reals are floats
        if (interFace.int64()) {
          error = epu(interFace, start_part, part_count, rank, static_cast<float>(0.0),
                      static_cast<int64_t>(0));
        }
        else {
          error = epu(interFace, start_part, part_count, rank, static_cast<float>(0.0), 0);
        }
      }
      else { // Reals are doubles
        if (interFace.int64()) {
          error = epu(interFace, start_part, part_count, rank, 0.0, static_cast<int64_t>(0));
        }
        else {
          error = epu(interFace, start_part, part_count, rank, 0.0, 0);
        }
      }

      ExodusFile::close_all();
#if ENABLE_PARALLEL_EPU
      MPI_Barrier(MPI_COMM_WORLD); // CHECK: ALLOW MPI_COMM_WORLD
#endif
    }
    else {
      int max_open_file = open_file_limit() - 1; // -1 for output exodus file.

      // Only used to test the auto subcycle without requiring thousands of files...
      if (interFace.max_open_files() > 0) {
        max_open_file = interFace.max_open_files();
      }

      if (interFace.is_auto() && interFace.subcycle() < 0 && processor_count > max_open_file &&
          part_count == processor_count && interFace.cycle() == -1) {
        // Rule of thumb -- number of subcycles = cube_root(processor_count);
        // if that value > max_open_file, then use square root.
        // if that is still too large, just do no subcycles... and implement
        // a recursive subcycling capability at some point...
        int sub_cycle_count = (int)(std::pow(processor_count, 1.0 / 3) + 0.9);
        if (((processor_count + sub_cycle_count - 1) / sub_cycle_count) > max_open_file) {
          sub_cycle_count = (int)std::sqrt(processor_count);
        }

        if (((processor_count + sub_cycle_count - 1) / sub_cycle_count) < max_open_file) {
          interFace.subcycle(sub_cycle_count);
          if (rank == 0) {
            fmt::print("\tAutomatically activating subcyle mode\n\tNumber of processors ({}) "
                       "exceeds open file limit ({}).\n"
                       "\tUsing --subcycle={}\n\n",
                       processor_count, max_open_file, sub_cycle_count);
          }
          interFace.subcycle_join(true);
        }
      }

      int cycle = interFace.cycle();
      if (interFace.subcycle() >= 0) {
        start_part = 0;
        int cycles = interFace.subcycle();
        if (cycles > 0) {
          // use the specified number of cycles...
          part_count = (processor_count + cycles - 1) / cycles;
          if (cycle >= 0) {
            start_part = cycle * part_count;
          }
        }

        // Sanity check...
        if (part_count < 1) {
          throw std::runtime_error(
              "ERROR: (EPU) The subcycle specification results in less than 1 part per "
              "cycle which is not allowed.\n");
        }
        interFace.subcycle((processor_count + part_count - 1) / part_count);

        if (start_part + part_count > processor_count) {
          part_count = processor_count - start_part;
        }
      }

      if (cycle < 0) {
        cycle = 0;
      }
      while (start_part < processor_count) {

        if (start_part + part_count > processor_count) {
          part_count = processor_count - start_part;
        }

        SMART_ASSERT(part_count > 0);
        SMART_ASSERT(start_part + part_count <= processor_count);

        if (!ExodusFile::initialize(interFace, start_part, part_count, cycle, false)) {
          throw std::runtime_error(
              "ERROR: (EPU) Problem initializing input and/or output files.\n");
        }

        if (ExodusFile::io_word_size() == 4) { // Reals are floats
          if (interFace.int64()) {
            error = epu(interFace, start_part, part_count, cycle++, static_cast<float>(0.0),
                        static_cast<int64_t>(0));
          }
          else {
            error = epu(interFace, start_part, part_count, cycle++, static_cast<float>(0.0), 0);
          }
        }
        else { // Reals are doubles
          if (interFace.int64()) {
            error = epu(interFace, start_part, part_count, cycle++, 0.0, static_cast<int64_t>(0));
          }
          else {
            error = epu(interFace, start_part, part_count, cycle++, 0.0, 0);
          }
        }

        start_part += part_count;
        ExodusFile::close_all();
        if (interFace.subcycle() < 0 || (interFace.subcycle() > 0 && interFace.cycle() >= 0)) {
          break;
        }
      }
    }

    if (interFace.subcycle() > 0 && interFace.cycle() < 0 && interFace.subcycle_join() &&
        rank == 0) {
      // Now, join the subcycled parts into a single file...
      start_part = 0;
      part_count = interFace.subcycle();
      interFace.subcycle(0);
      interFace.processor_count(part_count);
      interFace.step_min(1);
      interFace.step_max(INT_MAX);
      interFace.step_interval(1);

      if (!ExodusFile::initialize(interFace, start_part, part_count, 0, true)) {
        throw std::runtime_error("ERROR: (EPU) Problem initializing input and/or output files.\n");
      }

      if (ExodusFile::io_word_size() == 4) { // Reals are floats
        if (interFace.int64()) {
          error = epu(interFace, start_part, part_count, 0, static_cast<float>(0.0),
                      static_cast<int64_t>(0));
        }
        else {
          error = epu(interFace, start_part, part_count, 0, static_cast<float>(0.0), 0);
        }
      }
      else { // Reals are doubles
        if (interFace.int64()) {
          error = epu(interFace, start_part, part_count, 0, 0.0, static_cast<int64_t>(0));
        }
        else {
          error = epu(interFace, start_part, part_count, 0, 0.0, 0);
        }
      }

      bool delete_temp = error == 0 && !interFace.keep_temporary();
      ExodusFile::handle_temporary_files(delete_temp);
    }

    if (error == 0 && interFace.remove_file_per_rank_files()) {
      ExodusFile::unlink_input_files();
    }

    time_t end_time = std::time(nullptr);
    if (rank == 0) {
      add_to_log(argv[0], static_cast<int>(end_time - begin_time));
    }
    return error;
  }
  catch (std::exception &e) {
    fmt::print(stderr, "{}\n", e.what());
    return EXIT_FAILURE;
  }
}

template <typename T, typename INT>
int epu(SystemInterface &interFace, int start_part, int part_count, int cycle, T float_or_double,
        INT /*unused*/)
{
  double execution_time = seacas_timer();
  SMART_ASSERT(sizeof(T) == ExodusFile::io_word_size());

  if (rank == 0) {
    fmt::print("\nIO Word sizes: {} bytes floating point and {} bytes integer.\n", sizeof(T),
               sizeof(INT));
  }
  int p; // file counter p=0..part_count-1

  auto *mytitle = new char[MAX_LINE_LENGTH + 1];
  memset(mytitle, '\0', MAX_LINE_LENGTH + 1);

  Mesh global;

  // contains the global node information from each file/processor
  std::vector<std::vector<INT>> local_node_to_global(part_count);

  // contains the global element information from each file/processor
  std::vector<std::vector<INT>> local_element_to_global(part_count);

  // contains the global edge information from each file/processor
  std::vector<std::vector<INT>> local_edge_to_global(part_count);

  // contains the global face information from each file/processor
  std::vector<std::vector<INT>> local_face_to_global(part_count);

  std::vector<Mesh> local_mesh(part_count);

  // ******************************************************************
  // 1. Read global info

  int error = 0;

  LOG("\n**** READ LOCAL (GLOBAL) INFO ****\n");
  std::string title0;

  // EPU assumes IDS are always passed through the API as 64-bit ints.
  SMART_ASSERT(ex_int64_status(ExodusFile(0)) & EX_IDS_INT64_API);

  if (sizeof(INT) == 8) {
    SMART_ASSERT((ex_int64_status(ExodusFile(0)) & EX_BULK_INT64_API) != 0);
    SMART_ASSERT((ex_int64_status(ExodusFile(0)) & EX_MAPS_INT64_API) != 0);
  }
  else {
    SMART_ASSERT(sizeof(INT) == 4);
    SMART_ASSERT((ex_int64_status(ExodusFile(0)) & EX_BULK_INT64_API) == 0);
    SMART_ASSERT((ex_int64_status(ExodusFile(0)) & EX_MAPS_INT64_API) == 0);
  }

  // If there are any processors with zero nodes, then
  // that node won't have any nodal variables defined.  We need to
  // find the first processor which has a non-zero node count to use
  // when we look for the nodal variable count.
  int64_t non_zero_node_count = -1;
  for (p = 0; p < part_count; p++) {
    ex_init_params exodus{};
    error = ex_get_init_ext(ExodusFile(p), &exodus);
    if (error < 0) {
      exodus_error(__LINE__);
    }

    local_mesh[p].dimensionality = exodus.num_dim;
    local_mesh[p].nodeCount      = exodus.num_nodes;
    local_mesh[p].elementCount   = exodus.num_elem;
    local_mesh[p].edgeCount      = exodus.num_edge;
    local_mesh[p].faceCount      = exodus.num_face;
    local_mesh[p].blockCount     = exodus.num_elem_blk;
    local_mesh[p].nodesetCount   = exodus.num_node_sets;
    local_mesh[p].sidesetCount   = exodus.num_side_sets;
    local_mesh[p].assemblyCount  = exodus.num_assembly;
    local_mesh[p].edgeBlockCount = exodus.num_edge_blk;
    local_mesh[p].faceBlockCount = exodus.num_face_blk;
    local_mesh[p].title          = exodus.title;

    if (local_mesh[p].nodeCount > 0 && non_zero_node_count == -1) {
      non_zero_node_count = p;
    }

    if (p == 0) {
      global.title          = mytitle;
      global.dimensionality = local_mesh[p].dimensionality;
      global.blockCount     = local_mesh[p].count(Excn::ObjectType::EBLK);
      global.nodesetCount   = local_mesh[p].count(Excn::ObjectType::NSET);
      global.sidesetCount   = local_mesh[p].count(Excn::ObjectType::SSET);
      global.assemblyCount  = local_mesh[p].count(Excn::ObjectType::ASSM);
      global.edgeBlockCount = local_mesh[p].count(Excn::ObjectType::EDBLK);
      global.faceBlockCount = local_mesh[p].count(Excn::ObjectType::FABLK);
    }
    else {
      SMART_ASSERT(global.dimensionality == local_mesh[p].dimensionality);
      SMART_ASSERT(global.count(Excn::ObjectType::EBLK) ==
                   local_mesh[p].count(Excn::ObjectType::EBLK));
      SMART_ASSERT(global.count(Excn::ObjectType::ASSM) ==
                   local_mesh[p].count(Excn::ObjectType::ASSM));
      if (!interFace.omit_nodesets()) {
        SMART_ASSERT(global.count(Excn::ObjectType::NSET) ==
                     local_mesh[p].count(Excn::ObjectType::NSET));
      }
      if (!interFace.omit_sidesets()) {
        SMART_ASSERT(global.count(Excn::ObjectType::SSET) ==
                     local_mesh[p].count(Excn::ObjectType::SSET));
      }
      if (!interFace.omit_edgeblocks()) {
        SMART_ASSERT(global.count(Excn::ObjectType::EDBLK) ==
                     local_mesh[p].count(Excn::ObjectType::EDBLK));
      }
      if (!interFace.omit_faceblocks()) {
        SMART_ASSERT(global.count(Excn::ObjectType::FABLK) ==
                     local_mesh[p].count(Excn::ObjectType::FABLK));
      }
    }

    local_node_to_global[p].resize(local_mesh[p].nodeCount);
    local_element_to_global[p].resize(local_mesh[p].elementCount);
    local_edge_to_global[p].resize(local_mesh[p].edgeCount);
    local_face_to_global[p].resize(local_mesh[p].faceCount);

    // sum required data
    // note that num_blocks is the same for every processor
    global.elementCount += local_mesh[p].elementCount;
    global.edgeCount += local_mesh[p].edgeCount;
    global.faceCount += local_mesh[p].faceCount;

  } // end for (p=0..part_count)

  if (non_zero_node_count == -1) {
    non_zero_node_count = 0; // No nodes on entire model...
  }

  delete[] mytitle;

  if (interFace.omit_edgeblocks()) {
    global.edgeBlockCount = 0;
  }

  if (interFace.omit_faceblocks()) {
    global.faceBlockCount = 0;
  }

  if (interFace.omit_nodesets()) {
    global.nodesetCount = 0;
  }

  if (interFace.omit_sidesets()) {
    global.sidesetCount = 0;
  }

  // Need these throughout run, so declare outside of this block...
  std::vector<Block>              glob_blocks(global.count(Excn::ObjectType::EBLK));
  std::vector<std::vector<Block>> blocks(part_count);

  std::vector<EdgeBlock<INT>>              glob_edgeblocks(global.count(Excn::ObjectType::EDBLK));
  std::vector<std::vector<EdgeBlock<INT>>> edgeblocks(part_count);

  std::vector<FaceBlock<INT>>              glob_faceblocks(global.count(Excn::ObjectType::FABLK));
  std::vector<std::vector<FaceBlock<INT>>> faceblocks(part_count);

  std::vector<SideSet<INT>>              glob_ssets;
  std::vector<std::vector<SideSet<INT>>> sidesets(part_count);

  std::vector<NodeSet<INT>>              glob_nsets;
  std::vector<std::vector<NodeSet<INT>>> nodesets(part_count);
  {
    // Now, build the reverse global node map which permits access of the
    // local id given the global id.
    std::vector<INT> global_node_map;
    build_reverse_node_map(local_node_to_global, local_mesh, &global, part_count, global_node_map);

    LOG("Finished reading/writing Global Info\n");
    if (interFace.output_shared_nodes()) {
      // Get list of all shared nodes...
      std::vector<std::vector<INT>> shared(part_count);
      std::vector<int>              num_shared(global_node_map.size());

      for (auto &part : shared) {
        part.resize(global_node_map.size());
      }

      for (int pc = 0; pc < part_count; pc++) {
        size_t node_count = local_mesh[pc].nodeCount;
        for (size_t i = 0; i < node_count; i++) {
          INT gloc         = local_node_to_global[pc][i];
          shared[pc][gloc] = i + 1;
          num_shared[gloc]++;
        }
      }

      if (rank == 0) {
        fmt::print("Node Sharing information: (Part:Local Node Id)\n");
        for (size_t i = 0; i < global_node_map.size(); i++) {
          if (num_shared[i] > 1) {
            fmt::print("Global Node {}:", fmt::group_digits(i + 1));
            for (int pc = 0; pc < part_count; pc++) {
              if (shared[pc][i] >= 1) {
                fmt::print("\t{}:{}", pc, fmt::group_digits(shared[pc][i]));
              }
            }
            fmt::print("\n");
          }
        }
      }
    }

    // ****************************************************************************
    // Get Block information including element attributes
    // must check for zero length blocks
    get_element_blocks(part_count, local_mesh, global, blocks, glob_blocks);

    bool map_element_ids = interFace.map_element_ids();
    if (interFace.subcycle() >= 0) {
      map_element_ids = false;
    }
    std::vector<INT> global_element_map(global.elementCount);
    build_reverse_element_map(local_element_to_global, local_mesh, blocks, glob_blocks, &global,
                              part_count, global_element_map, map_element_ids);

    get_edgeblocks(part_count, local_mesh, global, edgeblocks, glob_edgeblocks);
    if (glob_edgeblocks.size() > 0) {
      bool map_edge_ids = interFace.map_edge_ids();
      if (interFace.subcycle() >= 0) {
        map_edge_ids = false;
      }
      std::vector<INT> global_edge_map(global.edgeCount);
      build_reverse_edge_map(local_edge_to_global, local_mesh, edgeblocks, glob_edgeblocks, &global,
                             part_count, global_edge_map, map_edge_ids);
    }

    get_faceblocks(part_count, local_mesh, global, faceblocks, glob_faceblocks);
    if (glob_faceblocks.size() > 0) {
      bool map_face_ids = interFace.map_face_ids();
      if (interFace.subcycle() >= 0) {
        map_face_ids = false;
      }
      std::vector<INT> global_face_map(global.faceCount);
      build_reverse_face_map(local_face_to_global, local_mesh, faceblocks, glob_faceblocks, &global,
                             part_count, global_face_map, map_face_ids);
    }

    //
    //    NOTE:  Node set/side set information can be different for each processor
    /************************************************************************/
    // Get Side sets
    if (!interFace.omit_sidesets()) {
      LOG("\n**** GET SIDE SETS *****\n");
      get_sideset_metadata(part_count, sidesets, glob_ssets);
      if (global.count(Excn::ObjectType::SSET) != glob_ssets.size()) {
        fmt::print("\nWARNING: Invalid sidesets will not be written to output database.\n");
        global.sidesetCount = glob_ssets.size();
      }
    }

    /************************************************************************/
    // Get Node sets
    if (!interFace.omit_nodesets()) {
      LOG("\n**** GET NODE SETS *****\n");
      get_nodesets(part_count, global.nodeCount, local_node_to_global, nodesets, glob_nsets,
                   float_or_double);
      if (global.count(Excn::ObjectType::NSET) != glob_nsets.size()) {
        fmt::print("\nWARNING: Invalid nodesets will not be written to output database.\n");
        global.nodesetCount = glob_nsets.size();
      }
    }

    /************************************************************************/
    // Start writing the output file...

    LOG("\n**** BEGIN WRITING OUTPUT FILE *****\n");
    CommunicationMetaData comm_data;

    if (!interFace.int64()) {
      int64_t twoBill = 1;
      twoBill <<= 31;
      int64_t fourBill = 1;
      fourBill <<= 32;
      // Check whether output mesh requires 64-bit integers...
      if (global.nodeCount >= twoBill || global.elementCount >= twoBill) {
        throw std::runtime_error("\n\nERROR: (EPU) Output file requires 64-bit integers. You must "
                                 "rerun epu with the -64 option.\n\n");
      }

      if (!interFace.use_netcdf4() && !interFace.use_netcdf5()) {
        // Check size required to store coordinates and connectivity
        if (global.nodeCount * 8 >= fourBill) {
          fmt::print(stderr, "\nINFO: Output file requires NetCDF-4 format or NetCDF-5. Setting "
                             "NetCDF-4 automatically.\n\n");
          interFace.set_use_netcdf4();
        }

        for (const auto &block : glob_blocks) {
          int64_t element_count = block.entity_count();
          int64_t nnpe          = block.nodesPerElement;
          if (element_count * nnpe * 4 >= fourBill) {
            fmt::print(stderr, "\nINFO: Output file requires NetCDF-4 format or NetCDF-5. Setting "
                               "NetCDF-4 automatically.\n\n");
            interFace.set_use_netcdf4();
            break;
          }
        }
      }
    }
    // Create the output file...
    if (!ExodusFile::create_output(interFace, cycle)) {
      throw std::runtime_error("ERROR: (EPU) Problem creating output file.\n");
    }

    // EPU assumes IDS are always passed through the API as 64-bit ints.
    SMART_ASSERT(ex_int64_status(ExodusFile::output()) & EX_IDS_INT64_API);

    if (sizeof(INT) == 8) {
      SMART_ASSERT((ex_int64_status(ExodusFile::output()) & EX_BULK_INT64_API) != 0);
      SMART_ASSERT((ex_int64_status(ExodusFile::output()) & EX_MAPS_INT64_API) != 0);
    }
    else {
      SMART_ASSERT(sizeof(INT) == 4);
      SMART_ASSERT((ex_int64_status(ExodusFile::output()) & EX_BULK_INT64_API) == 0);
      SMART_ASSERT((ex_int64_status(ExodusFile::output()) & EX_MAPS_INT64_API) == 0);
    }

    // Define metadata for model....
    put_global_info(global);
    get_put_coordinate_frames(ExodusFile(0), ExodusFile::output(), float_or_double);

    Internals<INT> exodus(ExodusFile::output(), ExodusFile::max_name_length());

    if (interFace.append()) {
      bool matches = exodus.check_meta_data(global, glob_blocks, glob_nsets, glob_ssets,
                                            glob_edgeblocks, glob_faceblocks, comm_data);
      if (!matches) {
        throw std::runtime_error("\n\nERROR: (EPU) Current mesh dimensions do not match "
                                 "the mesh dimensions in the file being appended to.\n\n");
      }
    }
    else {
      global.needNodeMap    = !is_sequential(global_node_map);
      global.needElementMap = !is_sequential(global_element_map);

      exodus.write_meta_data(global, glob_blocks, glob_nsets, glob_ssets, glob_edgeblocks,
                             glob_faceblocks, comm_data);

      get_put_assemblies(ExodusFile(0), ExodusFile::output(), global);

      if (interFace.add_processor_id_map()) {
        add_processor_map(ExodusFile::output(), part_count, start_part, global, blocks, glob_blocks,
                          local_element_to_global);
      }

      // Output bulk mesh data....
      put_nodesets(glob_nsets);

      // c.2.  Write Global Node Number Map
      if (global.needNodeMap) {
        LOG("Writing global node number map...\n");
        error = ex_put_id_map(ExodusFile::output(), EX_NODE_MAP, Data(global_node_map));
        if (error < 0) {
          exodus_error(__LINE__);
        }
      }

      if (global.needElementMap) {
        LOG("Writing out master global elements information...\n");
        if (!global_element_map.empty()) {
          error = ex_put_id_map(ExodusFile::output(), EX_ELEM_MAP, Data(global_element_map));
          if (error < 0) {
            exodus_error(__LINE__);
          }
        }
      }

      // Needed on glory writing to Lustre or we end up with empty maps...
      ex_update(ExodusFile::output());

      if (interFace.map_element_ids()) {
        put_element_blocks(part_count, start_part, blocks, glob_blocks, local_node_to_global,
                           local_element_to_global, float_or_double);
        put_edgeblocks(part_count, start_part, edgeblocks, glob_edgeblocks, local_node_to_global,
                       local_edge_to_global, float_or_double);
        put_faceblocks(part_count, start_part, faceblocks, glob_faceblocks, local_node_to_global,
                       local_face_to_global, float_or_double);
      }
      else {
        put_element_blocks(part_count, start_part, blocks, glob_blocks, local_node_to_global,
                           float_or_double);
        put_edgeblocks(part_count, start_part, edgeblocks, glob_edgeblocks, local_node_to_global,
                       float_or_double);
        put_faceblocks(part_count, start_part, faceblocks, glob_faceblocks, local_node_to_global,
                       float_or_double);
      }
    }

    get_put_sidesets(part_count, local_element_to_global, sidesets, glob_ssets, interFace);
  }
  // ************************************************************************
  // 2. Get Coordinate Info.
  if (!interFace.append()) {
    LOG("\n\n**** GET COORDINATE INFO ****\n");
    get_put_coordinates(global, part_count, local_mesh, local_node_to_global, (T)0.0);
    LOG("Wrote coordinate information...\n");
  }

  if (interFace.add_nodal_communication_map() && interFace.subcycle() > 0) {
    LOG("\n\n**** GET NODAL COMMUNICATION MAP INFO ****\n");
    // Need a mapping from processors in the original mesh parts to
    // processors in the subcycle output...
    auto proc_count = interFace.processor_count();

    std::vector<int> processor_map(proc_count);
    auto             orig_part_count = interFace.part_count();
    for (int i = 0; i < proc_count; i++) {
      processor_map[i] = i / orig_part_count;
    }
    get_put_nodal_communication_map(part_count, local_node_to_global, processor_map, cycle);
    LOG("Wrote nodal communication map information...\n");
  }

  // ####################TRANSIENT DATA SECTION###########################
  // ***********************************************************************
  // 9. Get Variable Information and names
  LOG("\n**** GET VARIABLE INFORMATION AND NAMES ****\n");

  //  I. read number of variables for each type.
  //  NOTE: it is assumed that every processor has the same global, nodal,
  //        and element lists

  Variables global_vars(Excn::ObjectType::GLOBAL);
  Variables nodal_vars(Excn::ObjectType::NODE);
  Variables element_vars(Excn::ObjectType::EBLK);
  Variables nodeset_vars(Excn::ObjectType::NSET);
  Variables sideset_vars(Excn::ObjectType::SSET);
  Variables edgeblock_vars(Excn::ObjectType::EDBLK);
  Variables faceblock_vars(Excn::ObjectType::FABLK);

  element_vars.addProcessorId = interFace.add_processor_id_field();

  {
    ExodusFile id(non_zero_node_count);

    get_variable_params(id, global_vars, interFace.global_var_names());
    get_variable_params(id, nodal_vars, interFace.node_var_names());
    get_variable_params(id, element_vars, interFace.elem_var_names());
    if (!interFace.omit_nodesets()) {
      get_variable_params(id, nodeset_vars, interFace.nset_var_names());
    }
    if (!interFace.omit_sidesets()) {
      get_variable_params(id, sideset_vars, interFace.sset_var_names());
    }

    if (!interFace.omit_edgeblocks()) {
      get_variable_params(id, edgeblock_vars, interFace.edblk_var_names());
    }

    if (!interFace.omit_faceblocks()) {
      get_variable_params(id, faceblock_vars, interFace.fablk_var_names());
    }

    get_truth_table(global, glob_blocks, local_mesh, element_vars, 4);
    filter_truth_table(id, global, glob_blocks, element_vars, interFace.elem_var_names());

    if (!interFace.omit_nodesets()) {
      get_truth_table(global, glob_nsets, local_mesh, nodeset_vars, 32);
      filter_truth_table(id, global, glob_nsets, nodeset_vars, interFace.nset_var_names());
    }

    if (!interFace.omit_sidesets()) {
      get_truth_table(global, glob_ssets, local_mesh, sideset_vars, 16);
      filter_truth_table(id, global, glob_ssets, sideset_vars, interFace.sset_var_names());
    }

    if (!interFace.omit_edgeblocks()) {
      get_truth_table(global, glob_edgeblocks, local_mesh, edgeblock_vars, 64);
      filter_truth_table(id, global, glob_edgeblocks, edgeblock_vars, interFace.edblk_var_names());
    }

    if (!interFace.omit_faceblocks()) {
      get_truth_table(global, glob_faceblocks, local_mesh, faceblock_vars, 128);
      filter_truth_table(id, global, glob_faceblocks, faceblock_vars, interFace.fablk_var_names());
    }
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
  if (global_vars.count(InOut::OUT) + nodal_vars.count(InOut::OUT) +
          element_vars.count(InOut::OUT) + nodeset_vars.count(InOut::OUT) +
          sideset_vars.count(InOut::OUT) + edgeblock_vars.count(InOut::OUT) +
          faceblock_vars.count(InOut::OUT) >
      0) {

    std::vector<int> elem_truth_table(
        global.truthTable[static_cast<int>(Excn::ObjectType::EBLK)].size());
    create_output_truth_table(global, glob_blocks, element_vars, elem_truth_table);

    if (!interFace.append()) {
      error = ex_put_all_var_param(
          ExodusFile::output(), global_vars.count(InOut::OUT), nodal_vars.count(InOut::OUT),
          element_vars.count(InOut::OUT), Data(elem_truth_table), nodeset_vars.count(InOut::OUT),
          Data(global.truthTable[static_cast<int>(Excn::ObjectType::NSET)]),
          sideset_vars.count(InOut::OUT),
          Data(global.truthTable[static_cast<int>(Excn::ObjectType::SSET)]));
      if (error < 0) {
        exodus_error(__LINE__);
      }

      if (edgeblock_vars.count(InOut::OUT) > 0) {
        error = ex_put_variable_param(ExodusFile::output(), EX_EDGE_BLOCK,
                                      edgeblock_vars.count(InOut::OUT));
        if (error < 0) {
          exodus_error(__LINE__);
        }
        error =
            ex_put_truth_table(ExodusFile::output(), EX_EDGE_BLOCK, glob_edgeblocks.size(),
                               edgeblock_vars.count(InOut::OUT),
                               Data(global.truthTable[static_cast<int>(Excn::ObjectType::EDBLK)]));
        if (error < 0) {
          exodus_error(__LINE__);
        }
      }

      if (faceblock_vars.count(InOut::OUT) > 0) {
        error = ex_put_variable_param(ExodusFile::output(), EX_FACE_BLOCK,
                                      faceblock_vars.count(InOut::OUT));
        if (error < 0) {
          exodus_error(__LINE__);
        }
        error =
            ex_put_truth_table(ExodusFile::output(), EX_FACE_BLOCK, glob_faceblocks.size(),
                               faceblock_vars.count(InOut::OUT),
                               Data(global.truthTable[static_cast<int>(Excn::ObjectType::FABLK)]));
        if (error < 0) {
          exodus_error(__LINE__);
        }
      }
    }
  }

  // II. read/write the variable names
  {
    ExodusFile id(non_zero_node_count);
    get_put_variable_names(id, ExodusFile::output(), global_vars, interFace);
    get_put_variable_names(id, ExodusFile::output(), nodal_vars, interFace);
    get_put_variable_names(id, ExodusFile::output(), element_vars, interFace);
    if (!interFace.omit_nodesets()) {
      get_put_variable_names(id, ExodusFile::output(), nodeset_vars, interFace);
    }
    if (!interFace.omit_sidesets()) {
      get_put_variable_names(id, ExodusFile::output(), sideset_vars, interFace);
    }
    if (!interFace.omit_edgeblocks()) {
      get_put_variable_names(id, ExodusFile::output(), edgeblock_vars, interFace);
    }
    if (!interFace.omit_faceblocks()) {
      get_put_variable_names(id, ExodusFile::output(), faceblock_vars, interFace);
    }
  }
  if (!interFace.append()) {
    ex_update(ExodusFile::output());
  }

  /**********************************************************************/
  // 10. Get Transient Data
  //     This routine reads in a time dump from an EXODUSII file

  int time_step;
  int num_time_steps = 0;

  LOG("\n**** GET TRANSIENT NODAL, GLOBAL, AND ELEMENT DATA VALUES ****\n");
  // Stage I: Get the number_of_time_steps information

  bool differ = false;
  for (p = 0; p < part_count; p++) {
    ExodusFile id(p);

    int nts = ex_inquire_int(id, EX_INQ_TIME);
    if (p == 0) {
      num_time_steps = nts;
    }
    else {
      if (nts != num_time_steps) {
        differ = true;
      }
      num_time_steps = num_time_steps < nts ? num_time_steps : nts;
    }
  }
  if (differ) {
    fmt::print(stderr,
               "\nWARNING: The number of time steps is not the same on all input databases.\n"
               "         Using minimum count of {}\n\n",
               num_time_steps);
  }
  else {
    if (rank == 0) {
      fmt::print("\nNumber of time steps on input databases = {}\n\n", num_time_steps);
    }
  }

  std::vector<T> global_values(global_vars.count(InOut::IN));
  std::vector<T> output_global_values(global_vars.count(InOut::OUT));

  // TODO(gdsjaar): Handle variables via a class instead of 3-D array.
  // Determine maximum number of entities on any processor...
  size_t max_ent = find_max_entity_count(part_count, local_mesh, global, blocks, nodesets, sidesets,
                                         edgeblocks, faceblocks);
  std::vector<T> values(max_ent);

  size_t max_global_ent = find_max_global_entity_count(global, glob_blocks, glob_nsets, glob_ssets,
                                                       glob_edgeblocks, glob_faceblocks);
  std::vector<T> master_values(max_global_ent);

  // Stage II.  Extracting transient variable data.
  //            loop over time steps

  if (num_time_steps == 0 && element_vars.add_processor_id()) {
    // Add a fake timestep with just the processor id information.
    // If adding the processor_id field, do it here...
    T time_val = 0.0;
    error      = ex_put_time(ExodusFile::output(), 1, &time_val);
    if (error < 0) {
      exodus_error(__LINE__);
    }

    std::vector<T> proc;
    add_processor_variable(ExodusFile::output(), part_count, start_part, global, blocks,
                           glob_blocks, local_element_to_global, 1,
                           element_vars.index_[element_vars.count()], proc);
  }

  // Determine if user wants a subset of timesteps transferred to the output file.
  int ts_min  = interFace.step_min();
  int ts_max  = interFace.step_max();
  int ts_step = interFace.step_interval();

  if (ts_min < 0) {
    ts_min = num_time_steps + 1 + ts_min;
  }
  if (ts_max < 0) {
    ts_max = num_time_steps + 1 + ts_max;
  }

  // Time steps for output file
  int time_step_out     = 0;
  T   sentinel          = static_cast<T>(-FLT_MAX);
  T   min_time_to_write = sentinel;

  if (interFace.append()) {
    // See how many steps already exist on the output database
    // and the corresponding time.
    int nstep = ex_inquire_int(ExodusFile::output(), EX_INQ_TIME);

    // Get the time corresponding to this step...
    error = ex_get_time(ExodusFile::output(), nstep, &min_time_to_write);
    if (error < 0) {
      exodus_error(__LINE__);
    }

    time_step_out = nstep;
  }

  ts_max = ts_max < num_time_steps ? ts_max : num_time_steps;
  if (ts_min <= ts_max) {
    if (debug_level & 1) {
      fmt::print("{}", time_stamp(tsFormat));
    }
    if (rank == 0) {
      fmt::print("\tTransferring step {} to step {} by {}\n", ts_min, ts_max, ts_step);
    }
  }

  // Determine how many steps will be written...
  int output_steps = (ts_max - ts_min) / ts_step + 1;
  int subcycles    = interFace.subcycle();

  double start_time = seacas_timer();

  for (time_step = ts_min - 1; time_step < ts_max; time_step += ts_step) {
    time_step_out++;

    T time_val = -std::numeric_limits<T>::max();
    {
      // read in and write out the time step information
      ExodusFile id(0);

      error = ex_get_time(id, time_step + 1, &time_val);
      if (error < 0) {
        exodus_error(__LINE__);
      }
      if (time_val <= min_time_to_write) {
        continue;
      }

      if (min_time_to_write != sentinel) {
        if (rank == 0) {
          fmt::print("\tAppend Mode: Skipping {} input steps to align times with already written "
                     "steps on output file.\n\n",
                     time_step - (ts_min - 1));
        }
        min_time_to_write = sentinel;
      }

      error = ex_put_time(ExodusFile::output(), time_step_out, &time_val);
      if (error < 0) {
        exodus_error(__LINE__);
      }

      for (p = 1; p < part_count; p++) {
        ExodusFile idp(p);
        T          proc_time_val = 0.0;
        error                    = ex_get_time(idp, time_step + 1, &proc_time_val);
        if (error < 0) {
          exodus_error(__LINE__);
        }
        if (proc_time_val != time_val) {
          fmt::print(stderr,
                     "WARNING: (EPU) At step {}, the times on processors {} and {} do not match:\n"
                     "         {:.8} vs {:.8} (absolute diff: {:.8})\n"
                     "         This may indicate a corrupt database.\n",
                     time_step + 1, start_part, p + start_part, time_val, proc_time_val,
                     std::abs(time_val - proc_time_val));
        }
      }

      // NOTE: Assuming that each processor has the exact same global
      // information
      if (global_vars.count(InOut::OUT) > 0) {
        if (debug_level & 1) {
          if (rank == 0) {
            fmt::print("{}Global Variables...\n", time_stamp(tsFormat));
          }
        }
        error = ex_get_var(id, time_step + 1, EX_GLOBAL, 0, 0, global_vars.count(),
                           Data(global_values));
        if (error < 0) {
          exodus_error(__LINE__);
        }
        // Map ...
        for (int ig = 0; ig < global_vars.count(InOut::IN); ig++) {
          if (global_vars.index_[ig] > 0) {
            SMART_ASSERT(ig < (int)global_values.size());
            output_global_values[global_vars.index_[ig] - 1] = global_values[ig];
          }
        }
        error = ex_put_var(ExodusFile::output(), time_step_out, EX_GLOBAL, 1, 0,
                           global_vars.count(InOut::OUT), Data(output_global_values));
        if (error < 0) {
          exodus_error(__LINE__);
        }

        // Check global variable consistency...
        if (debug_level & 512) {
          std::vector<T> proc_global_values(global_vars.count(InOut::IN));
          for (p = 1; p < part_count; p++) {
            ExodusFile idp(p);
            error = ex_get_var(idp, time_step + 1, EX_GLOBAL, 0, 0, global_vars.count(InOut::IN),
                               Data(proc_global_values));
            if (error < 0) {
              exodus_error(__LINE__);
            }
            for (int ig = 0; ig < global_vars.count(InOut::IN); ig++) {
              if (proc_global_values[ig] != global_values[ig]) {
                fmt::print(stderr,
                           fmt::runtime("At step {:{}}, Global Variable {:{}}, P{:0{}} = {:15.8g}, P{:0{}} = "
					"{:15.8g}\n"),
                           time_step + 1, ts_max + 1, ig + 1,
                           get_width(global_vars.count(InOut::IN)), start_part,
                           get_width(interFace.processor_count()), start_part + p,
                           get_width(interFace.processor_count()), proc_global_values[ig]);
              }
            }
          }
        }
      }
    }

    // Original: Total Execution Time = 20.55 seconds, Maximum memory = 3,826 MiBytes.
    // Refactor: Total Execution Time = 18.72 seconds, Maximum memory = 2,252 MiBytes. (nodal
    // variables only)
    //           Total Execution Time = 18.25 seconds, Maximum memory =   869 MiBytes. (all
    //           variables) Total Execution Time = 18.07 seconds, Maximum memory =   863 MiBytes.
    //           (all variables, share master_value array)
    // ========================================================================
    // Nodal Values...
    if (debug_level & 1) {
      if (rank == 0) {
        fmt::print("{}Nodal Variables...\n", time_stamp(tsFormat));
      }
    }

    if (nodal_vars.count(InOut::OUT) > 0) {
      T fill_val = -std::numeric_limits<T>::max();
      for (int i = 0; i < nodal_vars.count(InOut::IN); i++) {
        if (debug_level & 2) {
          std::fill(master_values.begin(), master_values.end(), fill_val);
        }

        for (p = 0; p < part_count; p++) {
          ExodusFile id(p);

          size_t node_count = local_mesh[p].nodeCount;
          if (nodal_vars.index_[i] > 0) {
            error = ex_get_var(id, time_step + 1, EX_NODAL, i + 1, 0, node_count, Data(values));
            if (error < 0) {
              exodus_error(__LINE__);
            }

            if (debug_level & 2) {
              for (size_t j = 0; j < node_count; j++) {
                size_t nodal_value = local_node_to_global[p][j];
                if (master_values[nodal_value] != fill_val &&
                    master_values[nodal_value] != values[j]) {
                  fmt::print(stderr, "Variable {}, Node {}, old = {}, new = {}\n", i + 1,
                             fmt::group_digits(nodal_value), master_values[nodal_value], values[j]);
                }
              }
            }

            if (interFace.sum_shared_nodes()) {
              // sum values into master nodal value information. Note
              // that for non-shared nodes, this will be the same as a
              // copy; for shared nodes, it will be a true sum.
              for (size_t j = 0; j < node_count; j++) {
                // Map local nodal value to global location...
                size_t nodal_value = local_node_to_global[p][j];
                master_values[nodal_value] += values[j];
              }
            }
            else {
              // copy values to master nodal value information
              for (size_t j = 0; j < node_count; j++) {
                // Map local nodal value to global location...
                size_t nodal_value         = local_node_to_global[p][j];
                master_values[nodal_value] = values[j];
              }
            }
          }
        }
        // output nodal variable info. for specified time step
        int i_out = nodal_vars.index_[i];
        SMART_ASSERT(i_out <= nodal_vars.count(InOut::OUT));
        error = ex_put_var(ExodusFile::output(), time_step_out, EX_NODAL, i_out, 0,
                           global.nodeCount, Data(master_values));
        if (error < 0) {
          exodus_error(__LINE__);
        }
      }
    }

    // ========================================================================
    // Extracting element transient variable data
    if (debug_level & 1) {
      if (rank == 0) {
        fmt::print("{}Element Variables...\n", time_stamp(tsFormat));
      }
    }

    if (element_vars.count(InOut::IN) > 0) {
      read_write_master_values(element_vars, global, glob_blocks, local_mesh, blocks, master_values,
                               values, part_count, time_step, time_step_out,
                               local_element_to_global);
    }

    // If adding the processor_id field, do it here...
    // Use the output time step for writing data
    if (element_vars.add_processor_id()) {
      std::vector<T> proc;
      add_processor_variable(ExodusFile::output(), part_count, start_part, global, blocks,
                             glob_blocks, local_element_to_global, time_step_out,
                             element_vars.index_[element_vars.count(InOut::IN)], proc);
    }

    // ========================================================================
    // Extracting sideset transient variable data
    if (!interFace.omit_sidesets()) {
      if (debug_level & 1) {
        if (rank == 0) {
          fmt::print("{}Sideset Variables...\n", time_stamp(tsFormat));
        }
      }

      if (sideset_vars.count(InOut::IN) > 0) {
        read_write_master_values(sideset_vars, global, glob_ssets, local_mesh, sidesets,
                                 master_values, values, part_count, time_step, time_step_out,
                                 local_element_to_global);
      }
    }

    if (!interFace.omit_nodesets()) {
      // ========================================================================
      // Extracting nodeset transient variable data
      if (debug_level & 1) {
        if (rank == 0) {
          fmt::print("{}Nodeset Variables...\n", time_stamp(tsFormat));
        }
      }

      if (nodeset_vars.count(InOut::IN) > 0) {
        read_write_master_values(nodeset_vars, global, glob_nsets, local_mesh, nodesets,
                                 master_values, values, part_count, time_step, time_step_out,
                                 local_element_to_global);
      }
    }

    if (!interFace.omit_edgeblocks()) {
      // ========================================================================
      // Extracting edgeblock transient variable data
      if (debug_level & 1) {
        if (rank == 0) {
          fmt::print("{}Edgeblock Variables...\n", time_stamp(tsFormat));
        }
      }

      if (edgeblock_vars.count(InOut::IN) > 0) {
        read_write_master_values(edgeblock_vars, global, glob_edgeblocks, local_mesh, edgeblocks,
                                 master_values, values, part_count, time_step, time_step_out,
                                 local_edge_to_global);
      }
    }

    if (!interFace.omit_faceblocks()) {
      // ========================================================================
      // Extracting faceblock transient variable data
      if (debug_level & 1) {
        if (rank == 0) {
          fmt::print("{}Faceblock Variables...\n", time_stamp(tsFormat));
        }
      }

      if (faceblock_vars.count(InOut::IN) > 0) {
        read_write_master_values(faceblock_vars, global, glob_faceblocks, local_mesh, faceblocks,
                                 master_values, values, part_count, time_step, time_step_out,
                                 local_face_to_global);
      }
    }

    // ========================================================================
    if (debug_level & 1) {
      fmt::print("{}", time_stamp(tsFormat));
    }

    if (subcycles > 2) {
      if (rank == 0) {
        fmt::print("{}/{} ", cycle + 1, subcycles);
      }
    }

    double cur_time            = seacas_timer();
    double elapsed             = cur_time - start_time;
    double time_per_step       = elapsed / time_step_out;
    double percentage_done     = (time_step_out * 100.0) / output_steps;
    double estimated_remaining = time_per_step * (output_steps - time_step_out);
    fmt::print("Wrote step {:6}, time {:8.4e}\t\t[{:5.1f}%, Elapsed={}, ETA={}]    \r",
               fmt::group_digits(time_step + 1), time_val, percentage_done, format_time(elapsed),
               format_time(estimated_remaining));
    if (debug_level & 1) {
      fmt::print("\n");
    }
  }
  /*************************************************************************/
  // FINALIZE program
  if (debug_level & 1) {
    fmt::print("{}", time_stamp(tsFormat));
  }
  if (subcycles > 2) {
    fmt::print("{}/{} ", cycle + 1, subcycles);
  }
  fmt::print("\n\nTotal Execution Time = {:.2f} seconds, Maximum memory = {} MiBytes.\n******* "
             "END *******\n",
             seacas_timer() - execution_time,
             fmt::group_digits((get_hwm_memory_info() + 1024 * 1024 - 1) / (1024 * 1024)));
  return 0;
}

namespace {
  void get_put_assemblies(int in, int out, Excn::Mesh &global)
  {
    if (global.assemblyCount > 0) {
      std::vector<ex_assembly> assemblies(global.assemblyCount);
      for (int i = 0; i < global.assemblyCount; i++) {
        assemblies[i].name        = nullptr;
        assemblies[i].entity_list = nullptr;
      }
      ex_get_assemblies(in, Data(assemblies));

      for (int i = 0; i < global.assemblyCount; i++) {
        assemblies[i].entity_list = new int64_t[assemblies[i].entity_count];
      }

      // Now get the assembly entity lists...
      ex_get_assemblies(in, Data(assemblies));
      ex_put_assemblies(out, assemblies.size(), Data(assemblies));

      for (int i = 0; i < global.assemblyCount; i++) {
        delete[] assemblies[i].entity_list;
        // This is `calloc`d in `ex_get_assemblies` call.
        free(assemblies[i].name);
      }
    }
  }

  template <typename T> void get_put_coordinate_frames(int id, int id_out, T /* float_or_double */)
  {
    int num_frames = ex_inquire_int(id, EX_INQ_COORD_FRAMES);
    if (num_frames <= 0) {
      return;
    }

    ExodusIdVector    ids(num_frames);
    std::vector<T>    coordinates(9 * num_frames);
    std::vector<char> tags(num_frames);

    int error = ex_get_coordinate_frames(id, &num_frames, Data(ids), Data(coordinates), Data(tags));
    if (error < 0) {
      exodus_error(__LINE__);
    }

    // Now output to the combined file...
    error = ex_put_coordinate_frames(id_out, num_frames, Data(ids), Data(coordinates), Data(tags));
    if (error < 0) {
      exodus_error(__LINE__);
    }
  }

  void get_put_qa(int id, int id_out)
  {
    // NOTE: Assuming info and QA records for all processors
    int error = 0;

    // I. Get and store info strings, if they exist
    int   num_info_records = ex_inquire_int(id, EX_INQ_INFO);
    auto *info_records     = new char *[num_info_records + 1];
    int   info_string_len  = MAX_LINE_LENGTH;

    {
      for (int i = 0; i < num_info_records + 1; i++) {
        info_records[i] = new char[info_string_len + 1];
        memset(info_records[i], '\0', info_string_len + 1);
      }
    }

    if (num_info_records > 0) {
      error = ex_get_info(id, info_records);
      if (error < 0) {
        exodus_error(__LINE__);
      }
    }

    // Add an info record for EPU
    add_info_record(info_records[num_info_records], MAX_LINE_LENGTH);

    error = ex_put_info(id_out, num_info_records + 1, info_records);
    if (error < 0) {
      exodus_error(__LINE__);
    }

    {
      for (int i = 0; i < num_info_records + 1; i++) {
        delete[] info_records[i];
      }
      delete[] info_records;
    }

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
    if (num_qa_records != 0) {
      error = ex_get_qa(id, qaRecord[0].qa_record);
      if (error < 0) {
        exodus_error(__LINE__);
      }
    }

    std::string buffer;

    copy_string(qaRecord[num_qa_records].qa_record[0][0], qainfo[0], MAX_STR_LENGTH + 1); // Code
    copy_string(qaRecord[num_qa_records].qa_record[0][1], qainfo[2], MAX_STR_LENGTH + 1); // Version

    time_t date_time = std::time(nullptr);
#if defined __NVCC__
    auto *lt = std::localtime(&date_time);
    buffer   = fmt::format("{:%Y/%m/%d}", *lt);
#else
    auto const lt = fmt::localtime(date_time);
    buffer        = fmt::format("{:%Y/%m/%d}", lt);
#endif
    copy_string(qaRecord[num_qa_records].qa_record[0][2], buffer, MAX_STR_LENGTH + 1);

#if defined __NVCC__
    buffer = fmt::format("{:%H:%M:%S}", *lt);
#else
    buffer = fmt::format("{:%H:%M:%S}", lt);
#endif
    copy_string(qaRecord[num_qa_records].qa_record[0][3], buffer, MAX_STR_LENGTH + 1);

    error = ex_put_qa(id_out, num_qa_records + 1, qaRecord[0].qa_record);
    if (error < 0) {
      exodus_error(__LINE__);
    }

    for (int i = 0; i < num_qa_records + 1; i++) {
      for (int j = 0; j < 4; j++) {
        delete[] qaRecord[i].qa_record[0][j];
      }
    }
  }

  template <typename T, typename INT>
  void get_put_coordinates(Mesh &global, int part_count, std::vector<Mesh> &local_mesh,
                           const std::vector<std::vector<INT>> &local_node_to_global,
                           T /* float_or_double */)
  {
    T FillValue = static_cast<T>(FILL_VALUE);
    SMART_ASSERT(sizeof(T) == ExodusFile::io_word_size());
    std::vector<T> x(global.nodeCount);
    std::vector<T> y(global.nodeCount);
    std::vector<T> z(global.nodeCount);

    if (debug_level & 8) {
      std::fill(x.begin(), x.end(), FillValue);
      std::fill(y.begin(), y.end(), FillValue);
      std::fill(z.begin(), z.end(), FillValue);
    }

    int error = 0;
    for (int p = 0; p < part_count; p++) {
      get_coordinates(ExodusFile(p), global.dimensionality, local_mesh[p].nodeCount,
                      local_node_to_global, p, x, y, z);
    } // end for p=0..part_count

    // Get Coordinate Names
    // NOTE: Assuming coordinate names should be
    // the same for all files/processors therefore, only one
    // file/processor needs to be loaded
    get_put_coordinate_names(ExodusFile(0), ExodusFile::output(), global.dimensionality);

    // Write out coordinate information
    error = ex_put_coord(ExodusFile::output(), Data(x), Data(y), Data(z));
    if (error < 0) {
      exodus_error(__LINE__);
    }
  }

  void get_put_coordinate_names(int in, int out, int dimensionality)
  {
    int    error            = 0;
    char **coordinate_names = get_name_array(dimensionality, ExodusFile::max_name_length());

    error = ex_get_coord_names(in, coordinate_names);
    if (error < 0) {
      exodus_error(__LINE__);
    }
    error = ex_put_coord_names(out, coordinate_names);
    if (error < 0) {
      exodus_error(__LINE__);
    }
    LOG("Wrote coordinate names...\n");

    free_name_array(coordinate_names, dimensionality);
  }

  template <typename INT>
  void get_put_nodal_communication_map(int                                  part_count,
                                       const std::vector<std::vector<INT>> &local_node_to_global,
                                       const std::vector<int> &processor_map, int output_processor)
  {
    int                              error = 0;
    std::vector<std::pair<INT, int>> node_cmap;
    for (int p = 0; p < part_count; p++) {
      INT num_int_nodes, num_bor_nodes, num_ext_nodes, num_int_elems, num_bor_elems;
      INT num_node_cmaps, num_elem_cmaps;
      ex_get_loadbal_param(ExodusFile(p), &num_int_nodes, &num_bor_nodes, &num_ext_nodes,
                           &num_int_elems, &num_bor_elems, &num_node_cmaps, &num_elem_cmaps, p);

      std::vector<INT> nodeCmapIds(num_node_cmaps);
      std::vector<INT> nodeCmapNodeCnts(num_node_cmaps);
      ex_get_cmap_params(ExodusFile(p), Data(nodeCmapIds), Data(nodeCmapNodeCnts), nullptr, nullptr,
                         p);

      int64_t my_node_count =
          std::accumulate(nodeCmapNodeCnts.begin(), nodeCmapNodeCnts.end(), int64_t(0));
      std::vector<INT> nodes(my_node_count);
      std::vector<INT> procs(my_node_count);

      int64_t cm_offset = 0;
      for (INT i = 0; i < num_node_cmaps; i++) {
        ex_get_node_cmap(ExodusFile(p), nodeCmapIds[i], &nodes[cm_offset], &procs[cm_offset], p);
        cm_offset += nodeCmapNodeCnts[i];
      }
      for (size_t i = 0; i < nodes.size(); i++) {
        auto output_proc = processor_map[procs[i]];
        if (output_proc != output_processor) {
          node_cmap.emplace_back(local_node_to_global[p][nodes[i] - 1] + 1,
                                 processor_map[procs[i]]);
        }
      }
    } // end for p=0..part_count
    uniquify(node_cmap);

    std::vector<INT> gnodes;
    std::vector<INT> gprocs;
    gnodes.reserve(node_cmap.size());
    gprocs.reserve(node_cmap.size());

    for (auto &[node, proc] : node_cmap) {
      gnodes.push_back(node);
      gprocs.push_back(proc);
    }

    // NOTE: This is inefficient in general since going in and out of define mode.
    //       Would be better to do this at file-creation time, but this encapsulates all
    //       code related to node communication map output to this routine...

    // Write out nodal communication map information
    error = ex_put_init_info(ExodusFile::output(), processor_map.back() + 1, 1, (char *)"P");
    if (error < 0) {
      exodus_error(__LINE__);
    }

    error = ex_put_loadbal_param(ExodusFile::output(), 0, 0, 0, 0, 0, 1, 0, output_processor);
    if (error < 0) {
      exodus_error(__LINE__);
    }

    std::array<INT, 1> ids{1};
    std::array<INT, 1> cnts{(INT)gnodes.size()};
    error = ex_put_cmap_params(ExodusFile::output(), Data(ids), Data(cnts), nullptr, nullptr,
                               output_processor);
    if (error < 0) {
      exodus_error(__LINE__);
    }

    error = ex_put_node_cmap(ExodusFile::output(), 1, Data(gnodes), Data(gprocs), output_processor);
    if (error < 0) {
      exodus_error(__LINE__);
    }
  }

  template <typename T, typename INT>
  void get_coordinates(int id, int dimensionality, size_t num_nodes,
                       const std::vector<std::vector<INT>> &local_node_to_global, int proc,
                       std::vector<T> &x, std::vector<T> &y, std::vector<T> &z)
  {
    SMART_ASSERT(sizeof(T) == ExodusFile::io_word_size());
    T              FillValue = static_cast<T>(FILL_VALUE);
    int            error     = 0;
    std::vector<T> local_x(num_nodes);
    std::vector<T> local_y(num_nodes);
    std::vector<T> local_z(num_nodes);

    error = ex_get_coord(id, Data(local_x), Data(local_y), Data(local_z));
    if (error < 0) {
      exodus_error(__LINE__);
    }

    // Check for 2D or 3D coordinates

    if (dimensionality == 3) {
      if (debug_level & 8) {
        for (size_t i = 0; i < num_nodes; i++) {
          size_t node = local_node_to_global[proc][i];
          if (x[node] != FillValue && y[node] != FillValue && z[node] != FillValue) {
            if (x[node] != local_x[i] || y[node] != local_y[i] || z[node] != local_z[i]) {
              fmt::print(stderr,
                         "\nWARNING: Node {} has different coordinates in at least two files.\n"
                         "         cur value = {:14.6e} {:14.6e} {:14.6e}\n"
                         "         new value = {:14.6e} {:14.6e} {:14.6e} from processor {}\n",
                         fmt::group_digits(node + 1), x[node], y[node], z[node], local_x[i],
                         local_y[i], local_z[i], proc);
            }
          }
        }
      }

      for (size_t i = 0; i < num_nodes; i++) {

        // The following WILL overwrite x[node],y[node],z[node] if node is the
        // same for different processors
        size_t node = local_node_to_global[proc][i];
        x[node]     = local_x[i];
        y[node]     = local_y[i];
        z[node]     = local_z[i];
      }
    }
    else if (dimensionality == 2) {
      if (debug_level & 8) {
        for (size_t i = 0; i < num_nodes; i++) {
          size_t node = local_node_to_global[proc][i];
          if (x[node] != FillValue && y[node] != FillValue) {
            if (x[node] != local_x[i] || y[node] != local_y[i]) {
              fmt::print(stderr,
                         "\nWARNING: Node {} has different coordinates in at least two files.\n"
                         "         cur value = {:14.6e} {:14.6e}\n"
                         "         new value = {:14.6e} {:14.6e} from processor {}\n",
                         fmt::group_digits(node + 1), x[node], y[node], local_x[i], local_y[i],
                         proc);
            }
          }
        }
      }
      for (size_t i = 0; i < num_nodes; i++) {

        // The following WILL overwrite x[node],y[node] if node is the same for
        // different processors
        size_t node = local_node_to_global[proc][i];

        x[node] = local_x[i];
        y[node] = local_y[i];
      }
    }
    else { // dimensionality == 1
      if (debug_level & 8) {
        for (size_t i = 0; i < num_nodes; i++) {
          size_t node = local_node_to_global[proc][i];
          if (x[node] != FillValue && y[node] != FillValue) {
            if (x[node] != local_x[i]) {
              fmt::print(stderr,
                         "\nWARNING: Node {} has different coordinates in at least two files.\n"
                         "         cur value = {:14.6e}\tnew value = {:14.6e} from processor {}\n",
                         fmt::group_digits(node + 1), x[node], local_x[i], proc);
            }
          }
        }
      }
      for (size_t i = 0; i < num_nodes; i++) {

        // The following WILL overwrite x[node] if node is the same for
        // different processors
        size_t node = local_node_to_global[proc][i];
        x[node]     = local_x[i];
      }
    }
  }

  void get_element_blocks(int part_count, const std::vector<Mesh> &local_mesh, const Mesh &global,
                          std::vector<std::vector<Block>> &blocks, std::vector<Block> &glob_blocks)
  {
    LOG("\n\n**** GET BLOCK INFORMATION (INCL. ELEMENT ATTRIBUTES) ****\n");

    for (int ip = 0; ip < part_count; ip++) {
      blocks[ip].resize(local_mesh[ip].count(Excn::ObjectType::EBLK));
    }

    if (rank == 0) {
      fmt::print("Global block count = {}\n", global.count(Excn::ObjectType::EBLK));
    }

    ExodusIdVector block_id(global.count(Excn::ObjectType::EBLK));

    int error = 0;
    for (int p = 0; p < part_count; p++) {
      ExodusFile id(p);

      error = ex_get_ids(id, EX_ELEM_BLOCK, Data(block_id));
      if (error < 0) {
        exodus_error(__LINE__);
      }

      // Check that the block id ordering is consistent among files...
      if (p > 0) {
        for (size_t b = 0; b < global.count(Excn::ObjectType::EBLK); b++) {
          if (blocks[0][b].id != block_id[b]) {
            std::ostringstream errmsg;
            fmt::print(errmsg,
                       "ERROR: (EPU) The internal element block id ordering for part {}\n"
                       "       is not consistent with the ordering for part 0.\n",
                       p);
            throw std::runtime_error(errmsg.str());
          }
        }
      }

      if ((debug_level & 4) != 0U) {
        fmt::print("\nGetting element block info for processor {}...\n", p);
      }
      else {
        if (p == 0) {
          LOG("\nGetting element block info.\n");
        }
      }

      for (size_t b = 0; b < global.count(Excn::ObjectType::EBLK); b++) {
        if ((debug_level & 4) != 0U) {
          fmt::print("Block {}, Id = {}", b, block_id[b]);
        }

        ex_block temp_block{};
        temp_block.id   = block_id[b];
        temp_block.type = EX_ELEM_BLOCK;
        error           = ex_get_block_param(id, &temp_block);
        if (error < 0) {
          exodus_error(__LINE__);
        }

        std::vector<char> name(Excn::ExodusFile::max_name_length() + 1);
        error = ex_get_name(id, EX_ELEM_BLOCK, block_id[b], Data(name));
        if (error < 0) {
          exodus_error(__LINE__);
        }

        blocks[p][b].id = block_id[b];
        if (name[0] != '\0') {
          blocks[p][b].name_ = Data(name);
        }
        if (p == 0) {
          glob_blocks[b].id = block_id[b];
          if (name[0] != '\0') {
            glob_blocks[b].name_ = Data(name);
          }
        }

        if (temp_block.num_entry != 0) {
          blocks[p][b].elementCount    = temp_block.num_entry;
          blocks[p][b].nodesPerElement = temp_block.num_nodes_per_entry;
          blocks[p][b].attributeCount  = temp_block.num_attribute;
          blocks[p][b].offset_         = temp_block.num_entry;
          blocks[p][b].position_       = b;
          copy_string(blocks[p][b].elType, temp_block.topology);

          glob_blocks[b].elementCount += temp_block.num_entry;
          glob_blocks[b].nodesPerElement = temp_block.num_nodes_per_entry;
          glob_blocks[b].attributeCount  = temp_block.num_attribute;
          glob_blocks[b].position_       = (int)b;
          copy_string(glob_blocks[b].elType, temp_block.topology);
        }

        if (temp_block.num_attribute > 0 && glob_blocks[b].attributeNames.empty()) {
          // Get attribute names.  Assume the same on all processors
          // on which the block exists.
          char **names = get_name_array(temp_block.num_attribute, ExodusFile::max_name_length());

          error = ex_get_attr_names(id, EX_ELEM_BLOCK, block_id[b], names);
          if (error < 0) {
            exodus_error(__LINE__);
          }
          for (int i = 0; i < temp_block.num_attribute; i++) {
            glob_blocks[b].attributeNames.emplace_back(names[i]);
          }
          free_name_array(names, temp_block.num_attribute);
        }
        if ((debug_level & 4) != 0U) {
          fmt::print(", Name = '{}', Elements = {:12}, Nodes/element = {}, Attributes = {}\n",
                     blocks[p][b].name_, fmt::group_digits(blocks[p][b].entity_count()),
                     blocks[p][b].nodesPerElement, blocks[p][b].attributeCount);
        }
      }
    } // end for p=0..part_count

    // Convert block_offset from elements/block/processor to true offset
    for (int p = 0; p < part_count; p++) {
      size_t sum = 0;
      for (size_t b = 0; b < global.count(Excn::ObjectType::EBLK); b++) {
        size_t save          = blocks[p][b].offset_;
        blocks[p][b].offset_ = sum;
        sum += save;
      }
    }
  }

  template <typename T, typename INT>
  void put_element_blocks(int part_count, int start_part, std::vector<std::vector<Block>> &blocks,
                          std::vector<Block>                  &glob_blocks,
                          const std::vector<std::vector<INT>> &local_node_to_global,
                          const std::vector<std::vector<INT>> &local_element_to_global,
                          T /* float_or_double */)
  {
    SMART_ASSERT(sizeof(T) == ExodusFile::io_word_size());
    int global_num_blocks = glob_blocks.size();

    LOG("\nReading and Writing element connectivity & attributes\n");

    std::vector<INT> linkage;
    std::vector<T>   attributes;
    std::vector<INT> local_linkage;
    std::vector<T>   local_attr;

    for (int b = 0; b < global_num_blocks; b++) {

      if (debug_level & 4) {
        fmt::print("\nOutput element block info for...\n"
                   "Block {}, Id = {}, Name = '{}', Elements = {:12}, Nodes/element = {}, "
                   "Attributes = {}\n"
                   "B{}:\t",
                   b, glob_blocks[b].id, glob_blocks[b].name_,
                   fmt::group_digits(glob_blocks[b].entity_count()), glob_blocks[b].nodesPerElement,
                   glob_blocks[b].attributeCount, b);
      }

      size_t max_nodes = glob_blocks[b].entity_count() * glob_blocks[b].nodesPerElement;
      linkage.resize(max_nodes);
      attributes.resize(static_cast<size_t>(glob_blocks[b].attributeCount) *
                        glob_blocks[b].entity_count());

      int error = 0;
      for (int p = 0; p < part_count; p++) {
        ExodusFile id(p);

        size_t global_pos;
        size_t global_block_pos;

        if (blocks[p][b].entity_count() > 0) { // non-zero length block

          if (debug_level & 4) {
            fmt::print("#");
          }
          size_t maximum_nodes = blocks[p][b].entity_count();
          maximum_nodes *= blocks[p][b].nodesPerElement;
          local_linkage.resize(maximum_nodes);

          ex_entity_id bid = blocks[p][b].id;
          error = ex_get_conn(id, EX_ELEM_BLOCK, bid, Data(local_linkage), nullptr, nullptr);
          if (error < 0) {
            fmt::print(
                stderr,
                "ERROR: (EPU) Cannot get element block connectivity for block {} on part {}.\n",
                bid, p + start_part);
            exodus_error(__LINE__);
          }
          size_t                  pos                     = 0;
          size_t                  goffset                 = glob_blocks[b].offset_;
          size_t                  element_count           = blocks[p][b].entity_count();
          size_t                  boffset                 = blocks[p][b].offset_;
          size_t                  npe                     = blocks[p][b].nodesPerElement;
          const std::vector<INT> &proc_loc_elem_to_global = local_element_to_global[p];
          const std::vector<INT> &proc_loc_node_to_global = local_node_to_global[p];

          for (size_t e = 0; e < element_count; e++) {
            global_block_pos = proc_loc_elem_to_global[(e + boffset)] - goffset;
            global_pos       = global_block_pos * npe;

            for (size_t n = 0; n < npe; n++) {
              size_t node           = proc_loc_node_to_global[local_linkage[pos++] - 1];
              linkage[global_pos++] = node + 1;
            }
          }

          // Get attributes list,  if it exists
          if (blocks[p][b].attributeCount > 0) {
            size_t max_attr = blocks[p][b].entity_count() * blocks[p][b].attributeCount;
            local_attr.resize(max_attr);

            error = ex_get_attr(id, EX_ELEM_BLOCK, blocks[p][b].id, Data(local_attr));
            if (error < 0) {
              exodus_error(__LINE__);
            }

            pos = 0;

            size_t att_count = blocks[p][b].attributeCount;
            for (size_t e = 0; e < element_count; e++) {
              // global_pos is global position within this element block...
              global_block_pos = local_element_to_global[p][(e + boffset)] - goffset;
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

      // Write out block info
      int id_out = ExodusFile::output(); // output file identifier

      if (!linkage.empty()) {
        error =
            ex_put_conn(id_out, EX_ELEM_BLOCK, glob_blocks[b].id, Data(linkage), nullptr, nullptr);
        if (error < 0) {
          exodus_error(__LINE__);
        }
      }

      // Write out attributes list if it exists
      if (glob_blocks[b].attributeCount > 0) {
        error = ex_put_attr(id_out, EX_ELEM_BLOCK, glob_blocks[b].id, Data(attributes));
        if (error < 0) {
          exodus_error(__LINE__);
        }
      } // end for b=0..global_num_blocks-1
      if (debug_level & 4) {
        fmt::print("\n");
      }
    }
    fmt::print("\n");
  }

  template <typename T, typename INT>
  void put_element_blocks(int part_count, int start_part, std::vector<std::vector<Block>> &blocks,
                          std::vector<Block>                  &glob_blocks,
                          const std::vector<std::vector<INT>> &local_node_to_global,
                          T /* float_or_double */)
  {
    // This variant of `put_element_blocks` is used in the case of
    // `nomap` and uses much less memory (but may be slower).  It
    // relies on the fact that with `nomap`, a the elements for a
    // block in a parts are all contiguous in the output file, so we
    // can read them, map the nodes, and then output them using a
    // partial write function and they don't need to be mapped into a
    // global array (which also does not need to be allocated).  Phase
    // 1 of a "low-memory" option for epu.

    SMART_ASSERT(sizeof(T) == ExodusFile::io_word_size());
    int global_num_blocks = glob_blocks.size();

    LOG("\nReading and Writing element connectivity & attributes (Low Memory Method)\n");

    for (int b = 0; b < global_num_blocks; b++) {

      if (debug_level & 4) {
        fmt::print("\nOutput element block info for...\n"
                   "Block {}, Id = {}, Name = '{}', Elements = {:12}, Nodes/element = {}, "
                   "Attributes = {}\n",
                   b, glob_blocks[b].id, glob_blocks[b].name_,
                   fmt::group_digits(glob_blocks[b].entity_count()), glob_blocks[b].nodesPerElement,
                   glob_blocks[b].attributeCount);
      }

      int id_out = ExodusFile::output(); // output file identifier

      size_t part_block_offset = 1;
      for (int p = 0; p < part_count; p++) {
        ExodusFile id(p);

        if (blocks[p][b].entity_count() > 0) { // non-zero length block
          size_t           node_count = blocks[p][b].entity_count() * blocks[p][b].nodesPerElement;
          std::vector<INT> local_linkage(node_count);

          ex_entity_id bid = blocks[p][b].id;
          int error = ex_get_conn(id, EX_ELEM_BLOCK, bid, Data(local_linkage), nullptr, nullptr);
          if (error < 0) {
            fmt::print(
                stderr,
                "ERROR: (EPU) Cannot get element block connectivity for block {} on part {}.\n",
                bid, p + start_part);
            exodus_error(__LINE__);
          }

          size_t                  element_count           = blocks[p][b].entity_count();
          size_t                  npe                     = blocks[p][b].nodesPerElement;
          const std::vector<INT> &proc_loc_node_to_global = local_node_to_global[p];

          size_t pos = 0;
          for (size_t e = 0; e < element_count; e++) {
            for (size_t n = 0; n < npe; n++) {
              size_t node          = proc_loc_node_to_global[local_linkage[pos] - 1];
              local_linkage[pos++] = node + 1;
            }
          }

          if (debug_level & 4) {
            fmt::print(stderr, "part, block, offset, count = {} {} {} {}\n", p, bid,
                       part_block_offset, fmt::group_digits(element_count));
          }
          error = ex_put_partial_conn(id_out, EX_ELEM_BLOCK, bid, part_block_offset, element_count,
                                      Data(local_linkage), nullptr, nullptr);
          if (error < 0) {
            fmt::print(stderr,
                       "ERROR: (EPU) Cannot output element block connectivity for block {} on part "
                       "{} (offset = {}, count = {}).\n",
                       bid, p + start_part, part_block_offset, fmt::group_digits(element_count));
            exodus_error(__LINE__);
          }

          // Get attributes list,  if it exists
          if (blocks[p][b].attributeCount > 0) {
            size_t         max_attr = blocks[p][b].entity_count() * blocks[p][b].attributeCount;
            std::vector<T> local_attr(max_attr);

            error = ex_get_attr(id, EX_ELEM_BLOCK, blocks[p][b].id, Data(local_attr));
            if (error < 0) {
              exodus_error(__LINE__);
            }

            error = ex_put_partial_attr(id_out, EX_ELEM_BLOCK, bid, part_block_offset,
                                        element_count, Data(local_attr));
            if (error < 0) {
              exodus_error(__LINE__);
            }
          }
          part_block_offset += element_count;

        } // end if blocks[p][b].entity_count() (non-zero length block)
      } // end for p=0..part_count-1
    }
  }

  template <typename INT>
  void build_reverse_element_map(std::vector<std::vector<INT>>   &local_element_to_global,
                                 const std::vector<Mesh>         &local_mesh,
                                 std::vector<std::vector<Block>> &blocks,
                                 std::vector<Block> &glob_blocks, Mesh *global, int part_count,
                                 std::vector<INT> &global_element_map, bool map_ids)
  {
    // Global element map and count.
    std::vector<std::vector<INT>> global_element_numbers(part_count);

    size_t tot_size = 0;
    for (int p = 0; p < part_count; p++) {
      tot_size += local_mesh[p].elementCount;
      global_element_numbers[p].resize(local_mesh[p].elementCount);
    }
    global_element_map.resize(tot_size);

    {
      size_t offset = 0;
      for (int p = 0; p < part_count; p++) {
        ExodusFile id(p);
        get_id_map(id, EX_ELEM_MAP, EX_INQ_ELEM_MAP, global_element_numbers[p]);
        std::copy(global_element_numbers[p].begin(), global_element_numbers[p].end(),
                  &global_element_map[offset]);
        offset += local_mesh[p].elementCount;
      }
    }

    // Now, sort the global_element_map array.
#if USE_STD_SORT
    std::sort(global_element_map.begin(), global_element_map.end());
#else
    pdqsort(global_element_map.begin(), global_element_map.end());
#endif

    global->elementCount = global_element_map.size();

    // See if any duplicates...
    for (int64_t i = 1; i < global->elementCount; i++) {
      if (global_element_map[i - 1] == global_element_map[i]) {
        // Duplicates in the element id list...
        // This is not yet handled.  Notify the user and continue for now...
        fmt::print(stderr, "\n!!!! POSSIBLE ERROR: (EPU) There were at least 2"
                           " elements with duplicated ids detected.\n"
                           "!!!!\tThis may cause problems in the output file.\n\n");
        break;
      }
    }

    // See whether the element numbers are contiguous.  If so, we can map
    // the elements back to their original location. Since the elements are
    // sorted and there are no duplicates, we just need to see if the id
    // at global_element_map.size() == global_element_map.size();
    bool is_contiguous = global_element_map.empty() ||
                         ((size_t)global_element_map.back() == global_element_map.size());
    if (rank == 0) {
      fmt::print("Element id map {} contiguous.\n", (is_contiguous ? "is" : "is not"));
    }

  // Create the map that maps from a local processor element to the
  // global map. This combines the mapping local processor element to
  // 'global id' and then 'global id' to global position. The
  // mapping is now a direct lookup instead of a lookup followed by
  // a reverse map.
  //
  // If the map is contiguous, then the global_id to global_position map is 1->1
  REMAP:
    if (is_contiguous && map_ids) {
      auto   cur_pos = global_element_map.begin();
      size_t element_value;
      for (int p = 0; p < part_count; p++) {
        size_t element_count = local_mesh[p].elementCount;
        for (size_t i = 0; i < element_count; i++) {
          INT global_element = global_element_numbers[p][i];

          if (cur_pos == global_element_map.end() || *cur_pos != global_element) {
            auto iter = std::lower_bound(global_element_map.begin(), global_element_map.end(),
                                         global_element);
            SMART_ASSERT(iter != global_element_map.end());
            cur_pos = iter;
          }
          element_value                 = cur_pos - global_element_map.begin();
          local_element_to_global[p][i] = element_value;
          ++cur_pos;
        }
      }
    }
    else {
      IntVector proc_off(part_count);
      std::fill(proc_off.begin(), proc_off.end(), 0);

      size_t gpos = 0;
      for (size_t b = 0; b < glob_blocks.size(); b++) {
        for (int p = 0; p < part_count; p++) {
          size_t poff          = proc_off[p];
          size_t element_count = blocks[p][b].entity_count();
          for (size_t e = 0; e < element_count; e++) {
            local_element_to_global[p][e + poff] = gpos++;
          }
          proc_off[p] += element_count;
        }
      }

      std::vector<INT> element_map;
      for (int p = 0; p < part_count; p++) {
        size_t element_count = local_mesh[p].elementCount;
        element_map.resize(element_count);
        get_id_map(ExodusFile(p), EX_ELEM_MAP, EX_INQ_ELEM_MAP, element_map);

        for (size_t e = 0; e < element_count; e++) {
          gpos                     = local_element_to_global[p][e];
          global_element_map[gpos] = element_map[e];
        }
      }
    }

    // Now, need to set up the element block offsets.  Within an element
    // block, the ids are contiguous, so to map a global element to its
    // position within the element block, use the formula:
    //
    //    pos_in_block = global_pos - block_offset
    //
    // Since there is the possibility that the element blocks in this
    // file are in a different order than the element blocks in original
    // file, we need to find out the minimum element id in each element
    // block and that is the offset (0-based ids).  Will also get the
    // maximum id to use as a sanity check.
    // Note that ids are contiguous from 0..element_count-1

    for (size_t b = 0; b < glob_blocks.size(); b++) {
      size_t min_id = global_element_map.size();
      size_t max_id = 0;

      for (int p = 0; p < part_count; p++) {
        size_t offset        = blocks[p][b].offset_;
        size_t element_count = blocks[p][b].entity_count();
        for (size_t e = 0; e < element_count; e++) {
          size_t id = local_element_to_global[p][e + offset];
          min_id    = (id < min_id) ? id : min_id;
          max_id    = (id > max_id) ? id : max_id;
        }
      }

      if (glob_blocks[b].entity_count() == 0) {
        min_id = 0;
        max_id = 0;
      }
      else {
        if (max_id - min_id + 1 != glob_blocks[b].entity_count()) {
          if (map_ids) {
            map_ids = false;
            fmt::print(stderr,
                       "WARNING: The element ids are globally contiguous,\n"
                       "\tbut they are not consistent for element block {}.\n"
                       "\tRetrying with element id mapping turned off.\n",
                       glob_blocks[b].id);
            goto REMAP;
          }
          else {
            std::ostringstream errmsg;
            fmt::print(errmsg,
                       "ERROR: (EPU) The element ids for element block {} are not consistent.\n"
                       "Block {}, Id = {} min/max id = {}/{} size = {}.\n",
                       glob_blocks[b].id, b, glob_blocks[b].id, min_id + 1, max_id + 1,
                       fmt::group_digits(glob_blocks[b].entity_count()));
            throw std::runtime_error(errmsg.str());
          }
        }
      }
      glob_blocks[b].offset_ = min_id;
      if (debug_level & 4) {
        fmt::print("Block {}, Id = {} min/max id = {}/{} offset = {}\n", b, glob_blocks[b].id,
                   min_id + 1, max_id + 1, glob_blocks[b].offset_);
      }
    }
  }

  template <typename INT>
  void build_reverse_edge_map(std::vector<std::vector<INT>>            &local_edge_to_global,
                              const std::vector<Mesh>                  &local_mesh,
                              std::vector<std::vector<EdgeBlock<INT>>> &edgeblocks,
                              std::vector<EdgeBlock<INT>> &glob_edgeblocks, Mesh *global,
                              int part_count, std::vector<INT> &global_edge_map, bool map_ids)
  {
    // Global edge map and count
    std::vector<std::vector<INT>> global_edge_numbers(part_count);

    size_t tot_size = 0;
    for (int p = 0; p < part_count; p++) {
      tot_size += local_mesh[p].edgeCount;
      global_edge_numbers[p].resize(local_mesh[p].edgeCount);
    }
    global_edge_map.resize(tot_size);

    {
      size_t offset = 0;
      for (int p = 0; p < part_count; p++) {
        ExodusFile id(p);
        get_id_map(id, EX_EDGE_MAP, EX_INQ_EDGE_MAP, global_edge_numbers[p]);
        std::copy(global_edge_numbers[p].begin(), global_edge_numbers[p].end(),
                  &global_edge_map[offset]);
        offset += local_mesh[p].edgeCount;
      }
    }

    // Now, sort the global_edge_map array.
#if USE_STD_SORT
    std::sort(global_edge_map.begin(), global_edge_map.end());
#else
    pdqsort(global_edge_map.begin(), global_edge_map.end());
#endif

    global->edgeCount = global_edge_map.size();

    // See if any duplicates...
    for (int64_t i = 1; i < global->edgeCount; i++) {
      if (global_edge_map[i - 1] == global_edge_map[i]) {
        // Duplicates in the edge id list...
        // This is not yet handled.  Notify the user and continue for now...
        fmt::print(stderr, "\n!!!! POSSIBLE ERROR: (EPU) There were at least 2"
                           " edges with duplicated ids detected.\n"
                           "!!!!\tThis may cause problems in the output file.\n\n");
        break;
      }
    }

    // See whether the edge numbers are contiguous.  If so, we can map
    // the edges back to their original location. Since the edges are
    // sorted and there are no duplicates, we just need to see if the id
    // at global_edge_map.size() == global_edge_map.size();
    bool is_contiguous =
        global_edge_map.empty() || ((size_t)global_edge_map.back() == global_edge_map.size());
    if (rank == 0) {
      fmt::print("Edge id map {} contiguous.\n", (is_contiguous ? "is" : "is not"));
    }

  // Create the map that maps from a local processor edge to the
  // global map. This combines the mapping local processor edge to
  // 'global id' and then 'global id' to global position. The
  // mapping is now a direct lookup instead of a lookup followed by
  // a reverse map.
  //
  // If the map is contiguous, then the global_id to global_position map is 1->1
  REMAP:
    if (is_contiguous && map_ids) {
      auto   cur_pos = global_edge_map.begin();
      size_t edge_value;
      for (int p = 0; p < part_count; p++) {
        size_t edge_count = local_mesh[p].edgeCount;
        for (size_t i = 0; i < edge_count; i++) {
          INT global_edge = global_edge_numbers[p][i];

          if (cur_pos == global_edge_map.end() || *cur_pos != global_edge) {
            auto iter =
                std::lower_bound(global_edge_map.begin(), global_edge_map.end(), global_edge);
            SMART_ASSERT(iter != global_edge_map.end());
            cur_pos = iter;
          }
          edge_value                 = cur_pos - global_edge_map.begin();
          local_edge_to_global[p][i] = edge_value;
          ++cur_pos;
        }
      }
    }
    else {
      IntVector proc_off(part_count);
      std::fill(proc_off.begin(), proc_off.end(), 0);

      size_t gpos = 0;
      for (size_t b = 0; b < glob_edgeblocks.size(); b++) {
        for (int p = 0; p < part_count; p++) {
          size_t poff       = proc_off[p];
          size_t edge_count = edgeblocks[p][b].entity_count();
          for (size_t e = 0; e < edge_count; e++) {
            local_edge_to_global[p][e + poff] = gpos++;
          }
          proc_off[p] += edge_count;
        }
      }

      std::vector<INT> edge_map;
      for (int p = 0; p < part_count; p++) {
        size_t edge_count = local_mesh[p].edgeCount;
        edge_map.resize(edge_count);
        get_id_map(ExodusFile(p), EX_EDGE_MAP, EX_INQ_EDGE_MAP, edge_map);

        for (size_t e = 0; e < edge_count; e++) {
          gpos                  = local_edge_to_global[p][e];
          global_edge_map[gpos] = edge_map[e];
        }
      }
    }

    // Now, need to set up the edge block offsets.  Within an edge
    // block, the ids are contiguous, so to map a global edge to its
    // position within the edge block, use the formula:
    //
    //    pos_in_block = global_pos - block_offset
    //
    // Since there is the possibility that the edge blocks in this
    // file are in a different order than the edge blocks in original
    // file, we need to find out the minimum edge id in each edge
    // block and that is the offset (0-based ids).  Will also get the
    // maximum id to use as a sanity check.
    // Note that ids are contiguous from 0..edge_count-1

    for (size_t b = 0; b < glob_edgeblocks.size(); b++) {
      size_t min_id = global_edge_map.size();
      size_t max_id = 0;

      for (int p = 0; p < part_count; p++) {
        size_t offset     = edgeblocks[p][b].offset_;
        size_t edge_count = edgeblocks[p][b].entity_count();
        for (size_t e = 0; e < edge_count; e++) {
          size_t id = local_edge_to_global[p][e + offset];
          min_id    = (id < min_id) ? id : min_id;
          max_id    = (id > max_id) ? id : max_id;
        }
      }

      if (glob_edgeblocks[b].entity_count() == 0) {
        min_id = 0;
        max_id = 0;
      }
      else {
        if (max_id - min_id + 1 != glob_edgeblocks[b].entity_count()) {
          if (map_ids) {
            map_ids = false;
            fmt::print(stderr,
                       "WARNING: The edge ids are globally contiguous,\n"
                       "\tbut they are not consistent for edge block {}.\n"
                       "\tRetrying with edge id mapping turned off.\n",
                       glob_edgeblocks[b].id);
            goto REMAP;
          }
          else {
            std::ostringstream errmsg;
            fmt::print(errmsg,
                       "ERROR: (EPU) The edge ids for edge block {} are not consistent.\n"
                       "Block {}, Id = {} min/max id = {}/{} size = {}.\n",
                       glob_edgeblocks[b].id, b, glob_edgeblocks[b].id, min_id + 1, max_id + 1,
                       fmt::group_digits(glob_edgeblocks[b].entity_count()));
            throw std::runtime_error(errmsg.str());
          }
        }
      }
      glob_edgeblocks[b].offset_ = min_id;
      if (debug_level & 4) {
        fmt::print("Block {}, Id = {} min/max id = {}/{} offset = {}\n", b, glob_edgeblocks[b].id,
                   min_id + 1, max_id + 1, glob_edgeblocks[b].offset_);
      }
    }
  }

  template <typename INT>
  void build_reverse_face_map(std::vector<std::vector<INT>>            &local_face_to_global,
                              const std::vector<Mesh>                  &local_mesh,
                              std::vector<std::vector<FaceBlock<INT>>> &faceblocks,
                              std::vector<FaceBlock<INT>> &glob_faceblocks, Mesh *global,
                              int part_count, std::vector<INT> &global_face_map, bool map_ids)
  {
    // Global element map and count.
    std::vector<std::vector<INT>> global_face_numbers(part_count);

    size_t tot_size = 0;
    for (int p = 0; p < part_count; p++) {
      tot_size += local_mesh[p].faceCount;
      global_face_numbers[p].resize(local_mesh[p].faceCount);
    }
    global_face_map.resize(tot_size);

    {
      size_t offset = 0;
      for (int p = 0; p < part_count; p++) {
        ExodusFile id(p);
        get_id_map(id, EX_FACE_MAP, EX_INQ_FACE_MAP, global_face_numbers[p]);
        std::copy(global_face_numbers[p].begin(), global_face_numbers[p].end(),
                  &global_face_map[offset]);
        offset += local_mesh[p].faceCount;
      }
    }

    // Now, sort the global_face_map array.
#if USE_STD_SORT
    std::sort(global_face_map.begin(), global_face_map.end());
#else
    pdqsort(global_face_map.begin(), global_face_map.end());
#endif

    global->faceCount = global_face_map.size();

    // See if any duplicates...
    for (int64_t i = 1; i < global->faceCount; i++) {
      if (global_face_map[i - 1] == global_face_map[i]) {
        // Duplicates in the face id list...
        // This is not yet handled.  Notify the user and continue for now...
        fmt::print(stderr, "\n!!!! POSSIBLE ERROR: (EPU) There were at least 2"
                           " faces with duplicated ids detected.\n"
                           "!!!!\tThis may cause problems in the output file.\n\n");
        break;
      }
    }

    // See whether the face numbers are contiguous.  If so, we can map
    // the faces back to their original location. Since the faces are
    // sorted and there are no duplicates, we just need to see if the id
    // at global_face_map.size() == global_face_map.size();
    bool is_contiguous =
        global_face_map.empty() || ((size_t)global_face_map.back() == global_face_map.size());
    if (rank == 0) {
      fmt::print("Face id map {} contiguous.\n", (is_contiguous ? "is" : "is not"));
    }

  // Create the map that maps from a local processor face to the
  // global map. This combines the mapping local processor face to
  // 'global id' and then 'global id' to global position. The
  // mapping is now a direct lookup instead of a lookup followed by
  // a reverse map.
  //
  // If the map is contiguous, then the global_id to global_position map is 1->1
  REMAP:
    if (is_contiguous && map_ids) {
      auto   cur_pos = global_face_map.begin();
      size_t face_value;
      for (int p = 0; p < part_count; p++) {
        size_t face_count = local_mesh[p].faceCount;
        for (size_t i = 0; i < face_count; i++) {
          INT global_face = global_face_numbers[p][i];

          if (cur_pos == global_face_map.end() || *cur_pos != global_face) {
            auto iter =
                std::lower_bound(global_face_map.begin(), global_face_map.end(), global_face);
            SMART_ASSERT(iter != global_face_map.end());
            cur_pos = iter;
          }
          face_value                 = cur_pos - global_face_map.begin();
          local_face_to_global[p][i] = face_value;
          ++cur_pos;
        }
      }
    }
    else {
      IntVector proc_off(part_count);
      std::fill(proc_off.begin(), proc_off.end(), 0);

      size_t gpos = 0;
      for (size_t b = 0; b < glob_faceblocks.size(); b++) {
        for (int p = 0; p < part_count; p++) {
          size_t poff       = proc_off[p];
          size_t face_count = faceblocks[p][b].entity_count();
          for (size_t e = 0; e < face_count; e++) {
            local_face_to_global[p][e + poff] = gpos++;
          }
          proc_off[p] += face_count;
        }
      }

      std::vector<INT> face_map;
      for (int p = 0; p < part_count; p++) {
        size_t face_count = local_mesh[p].faceCount;
        face_map.resize(face_count);
        get_id_map(ExodusFile(p), EX_ELEM_MAP, EX_INQ_ELEM_MAP, face_map);

        for (size_t e = 0; e < face_count; e++) {
          gpos                  = local_face_to_global[p][e];
          global_face_map[gpos] = face_map[e];
        }
      }
    }

    // Now, need to set up the face block offsets.  Within an face
    // block, the ids are contiguous, so to map a global face to its
    // position within the face block, use the formula:
    //
    //    pos_in_block = global_pos - block_offset
    //
    // Since there is the possibility that the face blocks in this
    // file are in a different order than the face blocks in original
    // file, we need to find out the minimum face id in each face
    // block and that is the offset (0-based ids).  Will also get the
    // maximum id to use as a sanity check.
    // Note that ids are contiguous from 0..face_count-1

    for (size_t b = 0; b < glob_faceblocks.size(); b++) {
      size_t min_id = global_face_map.size();
      size_t max_id = 0;

      for (int p = 0; p < part_count; p++) {
        size_t offset     = faceblocks[p][b].offset_;
        size_t face_count = faceblocks[p][b].entity_count();
        for (size_t e = 0; e < face_count; e++) {
          size_t id = local_face_to_global[p][e + offset];
          min_id    = (id < min_id) ? id : min_id;
          max_id    = (id > max_id) ? id : max_id;
        }
      }

      if (glob_faceblocks[b].entity_count() == 0) {
        min_id = 0;
        max_id = 0;
      }
      else {
        if (max_id - min_id + 1 != glob_faceblocks[b].entity_count()) {
          if (map_ids) {
            map_ids = false;
            fmt::print(stderr,
                       "WARNING: The face ids are globally contiguous,\n"
                       "\tbut they are not consistent for face block {}.\n"
                       "\tRetrying with face id mapping turned off.\n",
                       glob_faceblocks[b].id);
            goto REMAP;
          }
          else {
            std::ostringstream errmsg;
            fmt::print(errmsg,
                       "ERROR: (EPU) The face ids for face block {} are not consistent.\n"
                       "Block {}, Id = {} min/max id = {}/{} size = {}.\n",
                       glob_faceblocks[b].id, b, glob_faceblocks[b].id, min_id + 1, max_id + 1,
                       fmt::group_digits(glob_faceblocks[b].entity_count()));
            throw std::runtime_error(errmsg.str());
          }
        }
      }
      glob_faceblocks[b].offset_ = min_id;
      if (debug_level & 4) {
        fmt::print("Block {}, Id = {} min/max id = {}/{} offset = {}\n", b, glob_faceblocks[b].id,
                   min_id + 1, max_id + 1, glob_faceblocks[b].offset_);
      }
    }
  }

  template <typename INT>
  void build_reverse_node_map(std::vector<std::vector<INT>> &local_node_to_global,
                              const std::vector<Mesh> &local_mesh, Mesh *global, int part_count,
                              std::vector<INT> &global_node_map)
  {
    // Instead of using <set> and <map>, consider using a sorted vector...
    // Append all local node maps to the global node map.
    // Sort the global node map
    // Remove duplicates.
    // Position within map is now the map...
    // When building the local-proc node to global id, use binary_search...

    // Global node map and count.
    std::vector<std::vector<INT>> global_node_numbers(part_count);

    size_t tot_size = 0;
    for (int p = 0; p < part_count; p++) {
      tot_size += local_mesh[p].nodeCount;
      global_node_numbers[p].resize(local_mesh[p].nodeCount);
    }
    global_node_map.resize(tot_size);

    size_t offset = 0;
    for (int p = 0; p < part_count; p++) {
      ExodusFile id(p);
      get_id_map(id, EX_NODE_MAP, EX_INQ_NODE_MAP, global_node_numbers[p]);
      std::copy(global_node_numbers[p].begin(), global_node_numbers[p].end(),
                &global_node_map[offset]);
      offset += local_mesh[p].nodeCount;
    }

    // Now, sort the global_node_map array and remove duplicates...
#if USE_STD_SORT
    std::sort(global_node_map.begin(), global_node_map.end());
#else
    pdqsort(global_node_map.begin(), global_node_map.end());
#endif
    global_node_map.resize(unique(global_node_map));
    global_node_map.shrink_to_fit();

    size_t total_num_nodes = global_node_map.size();
    global->nodeCount      = total_num_nodes;

    // See whether the node numbers are contiguous.  If so, we can map
    // the nodes back to their original location. Since the nodes are
    // sorted and there are no duplicates, we just need to see if the id
    // at global_node_map.size() == global_node_map.size();
    bool is_contiguous =
        global_node_map.empty() || ((size_t)global_node_map.back() == global_node_map.size());
    if (rank == 0) {
      fmt::print("Node map {} contiguous.\n", (is_contiguous ? "is" : "is not"));
    }

    // Create the map the maps from a local processor node to the
    // global map. This combines the mapping local processor node to
    // 'global id' and then 'global id' to global position. The
    // mapping is now a direct lookup instead of a lookup followed by
    // a reverse map.
    auto cur_pos = global_node_map.begin();
    INT  nodal_value;
    for (int p = 0; p < part_count; p++) {
      size_t node_count = local_mesh[p].nodeCount;
      for (size_t i = 0; i < node_count; i++) {
        INT global_node = global_node_numbers[p][i];

        if (cur_pos == global_node_map.end() || *cur_pos != global_node) {
          auto iter = std::lower_bound(global_node_map.begin(), global_node_map.end(), global_node);
          SMART_ASSERT(iter != global_node_map.end());
          cur_pos = iter;
        }
        nodal_value                = cur_pos - global_node_map.begin();
        local_node_to_global[p][i] = nodal_value;
        ++cur_pos;
      }
    }
  }

  void get_put_variable_names(int id, int out, Variables &vars, Excn::SystemInterface &interFace)
  {
    if (vars.count(InOut::OUT) > 0) {

      char **output_name_list =
          get_name_array(vars.count(InOut::OUT), ExodusFile::max_name_length());

      int extra          = vars.add_processor_id() ? 1 : 0;
      int num_input_vars = vars.index_.size();

      char **input_name_list = get_name_array(num_input_vars, ExodusFile::max_name_length());
      int error = ex_get_variable_names(id, vars.type(), num_input_vars - extra, input_name_list);
      if (error < 0) {
        exodus_error(__LINE__);
      }

      if (vars.add_processor_id()) {
        // Check that input list of names does not already contain 'processor_id'...
        // If found, create an 'epu'-specific name; don't redo the check; assume ok...
        bool found = false;
        for (int i = 0; i < num_input_vars && !found; i++) {
          if (case_compare(input_name_list[i], "processor_id")) {
            found = true;
          }
        }
        if (found) {
          fmt::print(stderr, "\nWARNING: Variable 'processor_id' already exists on database.\n"
                             "         Adding 'processor_id_epu' instead.\n\n");
          copy_string(input_name_list[num_input_vars - 1], "processor_id_epu",
                      ExodusFile::max_name_length() + 1);
        }
        else {
          copy_string(input_name_list[num_input_vars - 1], "processor_id",
                      ExodusFile::max_name_length() + 1);
        }
      }

      // Iterate through the 'var_index' and transfer
      // Assume that the number of pointers is limited to
      // the number of results variables, plus one
      // extra pointer for optional add_processor_id
      size_t maxlen = 0;
      for (int i = 0; i < num_input_vars; i++) {
        if (vars.index_[i] > 0) {
          copy_string(output_name_list[vars.index_[i] - 1], input_name_list[i],
                      ExodusFile::max_name_length() + 1);
          if (strlen(input_name_list[i]) > maxlen) {
            maxlen = strlen(input_name_list[i]);
          }
        }
      }
      maxlen += 2;
      int width = interFace.screen_width();
      // Assume 8 characters for initial tab...
      int nfield = (width - 8) / maxlen;
      if (nfield < 1) {
        nfield = 1;
      }

      if (rank == 0) {
        fmt::print("Found {} {} variables.\n\t", vars.count(InOut::OUT), vars.label());
        int i    = 0;
        int ifld = 1;
        while (i < vars.count(InOut::OUT)) {
          fmt::print("{:<{}}", output_name_list[i++], maxlen);
          if (++ifld > nfield && i < vars.count(InOut::OUT)) {
            fmt::print("\n\t");
            ifld = 1;
          }
        }
        fmt::print("\n\n");
      }

      if (!interFace.append()) {
        error = ex_put_variable_names(out, vars.type(), vars.count(InOut::OUT), output_name_list);
        if (error < 0) {
          exodus_error(__LINE__);
        }
      }

      free_name_array(output_name_list, vars.count(InOut::OUT));
      free_name_array(input_name_list, num_input_vars);
    }
  }

  void get_variable_params(int id, Variables &vars, const StringIdVector &variable_list)
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

    // If 'type' is ELEMENT and addProcessorId is true, then reserve
    // space for an additional variable.
    int extra = 0;
    if (vars.type() == EX_ELEM_BLOCK && vars.addProcessorId) {
      extra = 1;
    }

    int num_vars;
    int error = ex_get_variable_param(id, vars.type(), &num_vars);
    if (error < 0) {
      exodus_error(__LINE__);
    }

    vars.index_.resize(num_vars + extra);

    // Create initial index which defaults to no output...
    std::fill(vars.index_.begin(), vars.index_.end(), 0);
    // ...Except for the processor_id (if any)
    if (vars.addProcessorId) {
      vars.index_[num_vars] = 1;
    }

    // If 'variable_list' is empty or specified 'ALL', then all
    // variables are to be output
    if (variable_list.empty() ||
        (variable_list.size() == 1 && case_compare(variable_list[0].first, "all"))) {
      for (size_t i = 0; i < vars.index_.size(); i++) {
        vars.index_[i] = i + 1;
      }
      vars.outputCount = num_vars + extra;
      return;
    }

    // Another possibility is user specifies "NONE" for the variable
    // list so no variables will be written.  Just return 0 or 1 based
    // on the 'add_processor_id' setting. Default var_index is ok.
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
          std::ostringstream errmsg;
          fmt::print(errmsg, "ERROR: (EPU) Variable '{}' is not valid.\n", elem.first);
          throw std::runtime_error(errmsg.str());
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

      if (vars.addProcessorId) {
        vars.index_[num_vars] = nz_count; // Already counted above...
      }
      vars.outputCount = nz_count;
      return;
    }
  }

  void put_global_info(const Mesh &global)
  {
    // Write out Global info
    if (rank == 0) {
      fmt::print(" Title: {}\n\n"
                 " Number of coordinates per node  = {:15}\n"
                 " Number of nodes                 = {:15}\n"
                 " Number of elements              = {:15}\n"
                 " Number of element blocks        = {:15}\n"
                 " Number of assemblies            = {:15}\n\n"
                 " Number of nodal point sets      = {:15}\n"
                 " Number of element side sets     = {:15}\n\n"
                 " Number of edge blocks           = {:15}\n"
                 " Number of face blocks           = {:15}\n\n",
                 global.title, fmt::group_digits(global.dimensionality),
                 fmt::group_digits(global.nodeCount), fmt::group_digits(global.elementCount),
                 fmt::group_digits(global.count(Excn::ObjectType::EBLK)),
                 fmt::group_digits(global.count(Excn::ObjectType::ASSM)),
                 fmt::group_digits(global.count(Excn::ObjectType::NSET)),
                 fmt::group_digits(global.count(Excn::ObjectType::SSET)),
                 fmt::group_digits(global.count(Excn::ObjectType::EDBLK)),
                 fmt::group_digits(global.count(Excn::ObjectType::FABLK)));
    }
    int id_out = ExodusFile::output();
    get_put_qa(ExodusFile(0), id_out);
  }

  template <typename T, typename INT>
  void get_nodesets(int part_count, size_t total_node_count,
                    const std::vector<std::vector<INT>>    &local_node_to_global,
                    std::vector<std::vector<NodeSet<INT>>> &nodesets,
                    std::vector<NodeSet<INT>>              &glob_sets, T /* float_or_double */)
  {
    // Find number of nodesets in the global model...
    std::set<ex_entity_id> set_ids;
    ExodusIdVector         ids;

    int bad_ns = 0;
    for (int p = 0; p < part_count; p++) {
      ExodusFile id(p);
      int        ns_count = ex_inquire_int(id, EX_INQ_NODE_SETS);

      // Get the ids for these
      ids.resize(ns_count);
      int error = ex_get_ids(id, EX_NODE_SET, Data(ids));
      if (error < 0) {
        exodus_error(__LINE__);
      }

      for (int iset = 0; iset < ns_count; iset++) {
        if (ids[iset] != 0) {
          set_ids.insert(ids[iset]);
        }
        else {
          bad_ns++;
        }
      }
    }

    if (bad_ns != 0) {
      fmt::print(
          stderr,
          "ERROR: (EPU) There were {} nodesets (counting all files) which had an id equal to "
          "0 which is not allowed.\n",
          bad_ns);
    }

    if (set_ids.empty()) {
      return;
    }

    // set_ids now contains the union of all nodeset ids...
    for (int p = 0; p < part_count; p++) {
      ExodusFile id(p);

      nodesets[p].resize(set_ids.size());
      auto I  = set_ids.begin();
      auto IE = set_ids.end();

      // Get the ids again so we can map current order back to file order...
      int error = ex_get_ids(id, EX_NODE_SET, Data(ids));
      if (error < 0) {
        exodus_error(__LINE__);
      }

      std::vector<ex_set> sets(set_ids.size());
      int                 i = 0;
      while (I != IE) {
        sets[i].id                       = *I++;
        sets[i].type                     = EX_NODE_SET;
        sets[i].entry_list               = nullptr;
        sets[i].extra_list               = nullptr;
        sets[i].distribution_factor_list = nullptr;
        i++;
      }

      error = ex_get_sets(id, set_ids.size(), Data(sets));
      if (error < 0) {
        exodus_error(__LINE__);
      }

      for (size_t iset = 0; iset < set_ids.size(); iset++) {
        nodesets[p][iset].id        = sets[iset].id;
        nodesets[p][iset].nodeCount = sets[iset].num_entry;
        nodesets[p][iset].dfCount   = sets[iset].num_distribution_factor;

        std::vector<char> name(Excn::ExodusFile::max_name_length() + 1);
        error = ex_get_name(id, EX_NODE_SET, nodesets[p][iset].id, Data(name));
        if (error < 0) {
          exodus_error(__LINE__);
        }
        if (name[0] != '\0') {
          nodesets[p][iset].name_ = Data(name);
        }

        if (nodesets[p][iset].dfCount != 0 &&
            nodesets[p][iset].dfCount != nodesets[p][iset].nodeCount) {
          fmt::print(stderr,
                     "WARNING: Nodeset {} with id {} on processor {} has an invalid distribution "
                     "factor count ({})."
                     " The distribution factors will be set to 1.0 on this nodeset.\n",
                     nodesets[p][iset].name_, nodesets[p][iset].id, p,
                     fmt::group_digits(nodesets[p][iset].dfCount));
          nodesets[p][iset].dfCount = 0;
        }

        nodesets[p][iset].position_ = -1;
        for (size_t j = 0; j < ids.size(); j++) {
          if (nodesets[p][iset].id == ids[j]) {
            nodesets[p][iset].position_ = j;
            break;
          }
        }

        if (debug_level & 32) {
          fmt::print("Processor {} ", p);
          nodesets[p][iset].dump();
        }
      }
    }

    glob_sets.resize(set_ids.size());
    {
      // Now get the nodeset nodes and df.

      // This is inefficient since the processor loop is on
      // the inside...  The other ordering would use more memory...

      std::vector<INT> ns_nodes;
      std::vector<T>   ns_df;

      for (size_t ns = 0; ns < set_ids.size(); ns++) {

        std::vector<INT> glob_ns_nodes(total_node_count + 1);
        std::vector<T>   glob_ns_df(total_node_count + 1);

        for (int p = 0; p < part_count; p++) {
          ExodusFile id(p);
          if (p == 0) {
            glob_sets[ns].name_ = nodesets[p][ns].name_;
            glob_sets[ns].id    = nodesets[p][ns].id;
          }

          if (nodesets[p][ns].position_ != -1) {
            if (glob_sets[ns].position_ != -1) {
              SMART_ASSERT(glob_sets[ns].position_ == nodesets[p][ns].position_);
            }
            else {
              glob_sets[ns].position_ = nodesets[p][ns].position_;
            }
          }

          size_t size = nodesets[p][ns].entity_count();
          if (size > 0) {
            ns_nodes.resize(size);
            ns_df.resize(size);

            int error = ex_get_set(id, EX_NODE_SET, nodesets[p][ns].id, Data(ns_nodes), nullptr);
            if (error < 0) {
              exodus_error(__LINE__);
            }
            if (nodesets[p][ns].dfCount > 0) {
              ex_get_set_dist_fact(id, EX_NODE_SET, nodesets[p][ns].id, Data(ns_df));
            }
            else {
              std::fill(ns_df.begin(), ns_df.end(), T(1.0));
            }

            // The node ids are in local space -- map to global; bring df along (if any).
            for (size_t iset = 0; iset < size; iset++) {
              SMART_ASSERT(ns_nodes[iset] > 0)(p)(ns)(iset)(ns_nodes[iset]);
              size_t global_node         = local_node_to_global[p][ns_nodes[iset] - 1] + 1;
              glob_ns_nodes[global_node] = 1;
              glob_ns_df[global_node]    = ns_df[iset];
            }
          }
        }
        // Count number of nonzero entries and transfer to the
        // output nodeset
        // NOTE: global_node above is 1-based.
        glob_sets[ns].nodeCount =
            std::accumulate(glob_ns_nodes.begin(), glob_ns_nodes.end(), (INT)0);
        glob_sets[ns].nodeSetNodes.resize(glob_sets[ns].entity_count());
        glob_sets[ns].dfCount = glob_sets[ns].nodeCount;

        // distFactors is a vector of 'char' to allow storage of either float or double.
        glob_sets[ns].distFactors.resize(glob_sets[ns].dfCount * ExodusFile::io_word_size());

        T     *glob_df = (T *)(Data(glob_sets[ns].distFactors));
        size_t j       = 0;
        for (size_t i = 1; i <= total_node_count; i++) {
          if (glob_ns_nodes[i] == 1) {
            glob_df[j]                      = glob_ns_df[i];
            glob_sets[ns].nodeSetNodes[j++] = i;
            glob_ns_nodes[i]                = j;
          }
        }

        SMART_ASSERT(j == glob_sets[ns].entity_count())(j)(glob_sets[ns].entity_count())(ns);

        // See if we need a map from local nodeset position to global nodeset position
        // Only needed if there are nodeset variables
        // Assume all files have same number of variables...
        int num_vars;
        {
          ExodusFile id(0);
          int        error = ex_get_variable_param(id, EX_NODE_SET, &num_vars);
          if (error < 0) {
            exodus_error(__LINE__);
          }
        }
        if (num_vars > 0) {
          for (int p = 0; p < part_count; p++) {
            ExodusFile id(p);
            // Get the nodelist, but store it in nodeOrderMap.
            // global_pos = nodeOrderMap[i]
            NodeSet<INT> &nset   = nodesets[p][ns];
            size_t        nnodes = nset.entity_count();
            nset.nodeOrderMap.resize(nnodes);
            int error = ex_get_set(id, EX_NODE_SET, nset.id, Data(nset.nodeOrderMap), nullptr);
            if (error < 0) {
              exodus_error(__LINE__);
            }

            for (size_t i = 0; i < nnodes; i++) {
              size_t local_node    = nset.nodeOrderMap[i];                        // 1-based
              size_t global_node   = local_node_to_global[p][local_node - 1] + 1; // 1-based
              size_t global_pos    = glob_ns_nodes[global_node];                  // 1-based
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

  template <typename INT> void put_nodesets(std::vector<NodeSet<INT>> &glob_sets)
  {
    int exoid = ExodusFile::output();

    if (debug_level & 32) {
      fmt::print("\nOutput NodeSets:\n");
    }
    for (auto &glob_set : glob_sets) {
      int error = ex_put_set(exoid, EX_NODE_SET, glob_set.id, Data(glob_set.nodeSetNodes), nullptr);
      if (error < 0) {
        exodus_error(__LINE__);
      }
      if (glob_set.dfCount > 0) {
        error = ex_put_set_dist_fact(exoid, EX_NODE_SET, glob_set.id, Data(glob_set.distFactors));
        if (error < 0) {
          exodus_error(__LINE__);
        }
      }

      // Done with the memory; clear out the vector containing the bulk data nodes and distFactors.
      clear(glob_set.nodeSetNodes);
      clear(glob_set.distFactors);

      if (debug_level & 32) {
        glob_set.dump();
      }
    }
  }

  template <typename INT>
  void get_sideset_metadata(int part_count, std::vector<std::vector<SideSet<INT>>> &sets,
                            std::vector<SideSet<INT>> &glob_ssets)
  {
    // Find number of sidesets in the global model...
    std::set<ex_entity_id> set_ids;
    ExodusIdVector         ids;

    int bad_ss = 0;
    {
      for (int p = 0; p < part_count; p++) {
        ExodusFile id(p);
        int        ss_count = ex_inquire_int(id, EX_INQ_SIDE_SETS);

        // Get the ids for these
        ids.resize(ss_count);
        int error = ex_get_ids(id, EX_SIDE_SET, Data(ids));
        if (error < 0) {
          exodus_error(__LINE__);
        }
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
                 "ERROR: (EPU) There were {} sidesets (counting all files) which had an id equal "
                 "to 0 which is not allowed.\n",
                 bad_ss);
    }

    if (set_ids.empty()) {
      return;
    }

    // set_ids now contains the union of all sideset ids...
    glob_ssets.resize(set_ids.size());

    {
      for (int p = 0; p < part_count; p++) {
        ExodusFile id(p);

        sets[p].resize(set_ids.size());
        auto I  = set_ids.begin();
        auto IE = set_ids.end();

        // Get the ids again so we can map current order back to file order...
        int error = ex_get_ids(id, EX_SIDE_SET, Data(ids));
        if (error < 0) {
          exodus_error(__LINE__);
        }

        std::vector<ex_set> exosets(set_ids.size());
        int                 j = 0;
        while (I != IE) {
          exosets[j].id                       = *I++;
          exosets[j].type                     = EX_SIDE_SET;
          exosets[j].entry_list               = nullptr;
          exosets[j].extra_list               = nullptr;
          exosets[j].distribution_factor_list = nullptr;
          j++;
        }

        error = ex_get_sets(id, set_ids.size(), Data(exosets));
        if (error < 0) {
          exodus_error(__LINE__);
        }

        for (size_t i = 0; i < set_ids.size(); i++) {
          sets[p][i].id        = exosets[i].id;
          sets[p][i].sideCount = exosets[i].num_entry;
          sets[p][i].dfCount   = exosets[i].num_distribution_factor;

          glob_ssets[i].id = sets[p][i].id;
          glob_ssets[i].sideCount += sets[p][i].entity_count();
          glob_ssets[i].dfCount += sets[p][i].dfCount;

          std::vector<char> name(Excn::ExodusFile::max_name_length() + 1);
          error = ex_get_name(id, EX_SIDE_SET, sets[p][i].id, Data(name));
          if (error < 0) {
            exodus_error(__LINE__);
          }
          if (name[0] != '\0') {
            sets[p][i].name_ = Data(name);
            if (p == 0) {
              glob_ssets[i].name_ = Data(name);
            }
          }

          sets[p][i].position_ = -1;
          for (size_t jj = 0; jj < ids.size(); jj++) {
            if (sets[p][i].id == ids[jj]) {
              sets[p][i].position_ = jj;
              break;
            }
          }

          if (sets[p][i].position_ != -1) {
            if (glob_ssets[i].position_ != -1) {
              SMART_ASSERT(glob_ssets[i].position_ == sets[p][i].position_);
            }
            else {
              glob_ssets[i].position_ = sets[p][i].position_;
            }
          }
        }
      }

      // Calculate sideset offset
      for (size_t b = 0; b < glob_ssets.size(); b++) {
        size_t sum = 0;
        for (int p = 0; p < part_count; p++) {
          sets[p][b].offset_ = sum;
          sum += sets[p][b].entity_count();

          if (debug_level & 16) {
            fmt::print("Processor {} ", p);
            sets[p][b].dump();
          }
        }
      }
    }
  }

  template <typename INT>
  void get_put_sidesets(int                                     part_count,
                        const std::vector<std::vector<INT>>    &local_element_to_global,
                        std::vector<std::vector<SideSet<INT>>> &sets,
                        std::vector<SideSet<INT>> &glob_ssets, Excn::SystemInterface &interFace)
  {
    // TODO(gdsjaar): See what work is really needed if in append mode...

    // Get a temporary vector to maintain the current
    // offset into the glob_ssets for storing sides
    Int64Vector offset(glob_ssets.size());
    std::fill(offset.begin(), offset.end(), int64_t(0));

    Int64Vector df_offset(glob_ssets.size());
    std::fill(df_offset.begin(), df_offset.end(), int64_t(0));

    {
      for (auto &glob_sset : glob_ssets) {
        glob_sset.elems.resize(glob_sset.entity_count());
        glob_sset.sides.resize(glob_sset.entity_count());
        glob_sset.distFactors.resize(glob_sset.dfCount * ExodusFile::io_word_size());
      }
    }

    // Now get the sideset elements, sides and df.
    for (int p = 0; p < part_count; p++) {
      ExodusFile id(p);
      for (size_t ss = 0; ss < glob_ssets.size(); ss++) {

        size_t size = sets[p][ss].entity_count();
        if (size > 0) {
          size_t off   = offset[ss];
          int    error = ex_get_set(id, EX_SIDE_SET, sets[p][ss].id, &glob_ssets[ss].elems[off],
                                    &glob_ssets[ss].sides[off]);
          if (error < 0) {
            exodus_error(__LINE__);
          }

          // The element ids are in local space -- map to global
          for (size_t i = 0; i < size; i++) {
            size_t local_elem = glob_ssets[ss].elems[off + i];
            SMART_ASSERT(local_elem > 0)(p)(ss)(i)(local_elem);
            SMART_ASSERT(glob_ssets[ss].sides[off + i] > 0 && glob_ssets[ss].sides[off + i] <= 6);
            size_t global_elem            = local_element_to_global[p][local_elem - 1];
            glob_ssets[ss].elems[off + i] = global_elem + 1;
          }
          offset[ss] += size;
        }

        // Distribution factors...
        size_t df_size = sets[p][ss].dfCount;
        if (df_size > 0) {
          size_t df_off = df_offset[ss] * ExodusFile::io_word_size();
          int    error  = ex_get_set_dist_fact(id, EX_SIDE_SET, sets[p][ss].id,
                                               &glob_ssets[ss].distFactors[df_off]);
          if (error < 0) {
            exodus_error(__LINE__);
          }
          df_offset[ss] += df_size;
        }
      }
    }

    if (debug_level & 16) {
      fmt::print("\nOutput SideSets:\n");
      for (auto &glob_sset : glob_ssets) {
        glob_sset.dump();
      }
    }

    if (!interFace.append()) {
      // Now write the actual sideset data...
      int exoid = ExodusFile::output(); // output file identifier
      for (auto &glob_sset : glob_ssets) {
        int error =
            ex_put_set(exoid, EX_SIDE_SET, glob_sset.id, const_cast<INT *>(Data(glob_sset.elems)),
                       const_cast<INT *>(Data(glob_sset.sides)));
        if (error < 0) {
          exodus_error(__LINE__);
        }
        if (glob_sset.dfCount > 0) {
          error = ex_put_set_dist_fact(exoid, EX_SIDE_SET, glob_sset.id,
                                       reinterpret_cast<void *>(Data(glob_sset.distFactors)));
          if (error < 0) {
            exodus_error(__LINE__);
          }
        }
      }
    }

    for (auto &glob_sset : glob_ssets) {
      clear(glob_sset.elems);
      clear(glob_sset.sides);
      clear(glob_sset.distFactors);
    }
  }

  template <typename INT>
  void get_edgeblocks(int part_count, const std::vector<Mesh> &local_mesh, const Mesh &global,
                      std::vector<std::vector<EdgeBlock<INT>>> &edgeblocks,
                      std::vector<EdgeBlock<INT>>              &glob_edgeblocks)
  {
    LOG("\n\n**** GET EDGE BLOCK INFORMATION (INCL. ELEMENT ATTRIBUTES) ****\n");

    for (int ip = 0; ip < part_count; ip++) {
      edgeblocks[ip].resize(local_mesh[ip].count(Excn::ObjectType::EDBLK));
    }

    if (rank == 0) {
      fmt::print("Global edge block count = {}\n",
                 fmt::group_digits(global.count(Excn::ObjectType::EDBLK)));
    }

    if (global.count(Excn::ObjectType::EDBLK) == 0) {
      return;
    }

    ExodusIdVector edgeblock_id(global.count(Excn::ObjectType::EDBLK));

    int error = 0;
    for (int p = 0; p < part_count; p++) {
      ExodusFile id(p);

      error = ex_get_ids(id, EX_EDGE_BLOCK, Data(edgeblock_id));
      if (error < 0) {
        exodus_error(__LINE__);
      }

      // Check that the block id ordering is consistent among files...
      if (p > 0) {
        for (size_t b = 0; b < global.count(Excn::ObjectType::EDBLK); b++) {
          if (edgeblocks[0][b].id != edgeblock_id[b]) {
            std::ostringstream errmsg;
            fmt::print(errmsg,
                       "ERROR: (EPU) The internal edge block id ordering for part {}\n"
                       "       is not consistent with the ordering for part 0.\n",
                       p);
            throw std::runtime_error(errmsg.str());
          }
        }
      }

      if ((debug_level & 64) != 0U) {
        fmt::print("\nGetting edge block info for processor {}...\n", p);
      }
      else {
        if (p == 0) {
          LOG("\nGetting edge block info.\n");
        }
      }

      for (size_t b = 0; b < global.count(Excn::ObjectType::EDBLK); b++) {
        if ((debug_level & 64) != 0U) {
          fmt::print("Edge Block {}, Id = {}", b, edgeblock_id[b]);
        }

        ex_block temp_block{};
        temp_block.id   = edgeblock_id[b];
        temp_block.type = EX_EDGE_BLOCK;
        error           = ex_get_block_param(id, &temp_block);
        if (error < 0) {
          exodus_error(__LINE__);
        }

        std::vector<char> name(Excn::ExodusFile::max_name_length() + 1);
        error = ex_get_name(id, EX_EDGE_BLOCK, edgeblock_id[b], Data(name));
        if (error < 0) {
          exodus_error(__LINE__);
        }

        edgeblocks[p][b].id = edgeblock_id[b];
        if (name[0] != '\0') {
          edgeblocks[p][b].name_ = Data(name);
        }
        if (p == 0) {
          glob_edgeblocks[b].id = edgeblock_id[b];
          if (name[0] != '\0') {
            glob_edgeblocks[b].name_ = Data(name);
          }
        }

        if (temp_block.num_entry != 0) {
          edgeblocks[p][b].edgeCount      = temp_block.num_entry;
          edgeblocks[p][b].nodesPerEdge   = temp_block.num_nodes_per_entry;
          edgeblocks[p][b].attributeCount = temp_block.num_attribute;
          edgeblocks[p][b].offset_        = temp_block.num_entry;
          edgeblocks[p][b].position_      = b;
          copy_string(edgeblocks[p][b].elType, temp_block.topology);

          glob_edgeblocks[b].edgeCount += temp_block.num_entry;
          glob_edgeblocks[b].nodesPerEdge   = temp_block.num_nodes_per_entry;
          glob_edgeblocks[b].attributeCount = temp_block.num_attribute;
          glob_edgeblocks[b].position_      = (int)b;
          copy_string(glob_edgeblocks[b].elType, temp_block.topology);
        }

        if (temp_block.num_attribute > 0 && glob_edgeblocks[b].attributeNames.empty()) {
          // Get attribute names.  Assume the same on all processors
          // on which the block exists.
          char **names = get_name_array(temp_block.num_attribute, ExodusFile::max_name_length());

          error = ex_get_attr_names(id, EX_EDGE_BLOCK, edgeblock_id[b], names);
          if (error < 0) {
            exodus_error(__LINE__);
          }
          for (int i = 0; i < temp_block.num_attribute; i++) {
            glob_edgeblocks[b].attributeNames.emplace_back(names[i]);
          }
          free_name_array(names, temp_block.num_attribute);
        }
        if ((debug_level & 4) != 0U) {
          fmt::print(", Name = '{}', Edges = {:12}, Nodes/Edge= {}, Attributes = {}\n",
                     edgeblocks[p][b].name_, fmt::group_digits(edgeblocks[p][b].entity_count()),
                     edgeblocks[p][b].nodesPerEdge, edgeblocks[p][b].attributeCount);
        }
      }
    } // end for p=0..part_count

    // Convert block_offset from edges/block/processor to true offset
    for (int p = 0; p < part_count; p++) {
      size_t sum = 0;
      for (size_t b = 0; b < global.count(Excn::ObjectType::EDBLK); b++) {
        size_t save              = edgeblocks[p][b].offset_;
        edgeblocks[p][b].offset_ = sum;
        sum += save;
      }
    }
  }

  template <typename T, typename INT>
  void put_edgeblocks(int part_count, int start_part,
                      std::vector<std::vector<EdgeBlock<INT>>> &edgeblocks,
                      std::vector<EdgeBlock<INT>>              &glob_edgeblocks,
                      const std::vector<std::vector<INT>>      &local_node_to_global,
                      const std::vector<std::vector<INT>>      &local_edge_to_global,
                      T /* float_or_double */)
  {
    SMART_ASSERT(sizeof(T) == ExodusFile::io_word_size());
    int global_num_edgeblocks = glob_edgeblocks.size();

    auto linkage    = new INT *[global_num_edgeblocks];
    auto attributes = new T *[global_num_edgeblocks];

    LOG("\nReading and Writing edge connectivity & attributes\n");

    for (int b = 0; b < global_num_edgeblocks; b++) {

      if (debug_level & 64) {
        fmt::print("\nOutput edge block info for...\n"
                   "Edge Block {}, Id = {}, Name = '{}', Edges = {:12}, Nodes/Edge = {}, "
                   "Attributes = {}\n"
                   "B{}:\t",
                   b, glob_edgeblocks[b].id, glob_edgeblocks[b].name_,
                   fmt::group_digits(glob_edgeblocks[b].entity_count()),
                   glob_edgeblocks[b].nodesPerEdge, glob_edgeblocks[b].attributeCount, b);
      }

      size_t max_nodes = glob_edgeblocks[b].entity_count();
      max_nodes *= glob_edgeblocks[b].nodesPerEdge;

      if (max_nodes > 0) {
        linkage[b] = new INT[max_nodes];
      }
      else {
        linkage[b] = nullptr;
      }
      INT *block_linkage = linkage[b];

      // Initialize attributes list, if it exists
      if (glob_edgeblocks[b].attributeCount > 0) {
        attributes[b] = new T[static_cast<size_t>(glob_edgeblocks[b].attributeCount) *
                              glob_edgeblocks[b].entity_count()];
      }

      int error = 0;
      for (int p = 0; p < part_count; p++) {
        ExodusFile id(p);

        size_t global_pos;
        size_t global_block_pos;

        if (edgeblocks[p][b].entity_count() > 0) { // non-zero length block

          if (debug_level & 64) {
            fmt::print("#");
          }
          size_t maximum_nodes = edgeblocks[p][b].entity_count();
          maximum_nodes *= edgeblocks[p][b].nodesPerEdge;
          std::vector<INT> local_linkage(maximum_nodes);
          SMART_ASSERT(block_linkage != nullptr);

          ex_entity_id bid = edgeblocks[p][b].id;
          error = ex_get_conn(id, EX_EDGE_BLOCK, bid, Data(local_linkage), nullptr, nullptr);
          if (error < 0) {
            fmt::print(stderr,
                       "ERROR: (EPU) Cannot get edge block connectivity for block {} on part {}.\n",
                       bid, p + start_part);
            exodus_error(__LINE__);
          }
          size_t                  pos                     = 0;
          size_t                  goffset                 = glob_edgeblocks[b].offset_;
          size_t                  edge_count              = edgeblocks[p][b].entity_count();
          size_t                  boffset                 = edgeblocks[p][b].offset_;
          size_t                  npe                     = edgeblocks[p][b].nodesPerEdge;
          const std::vector<INT> &proc_loc_edge_to_global = local_edge_to_global[p];
          const std::vector<INT> &proc_loc_node_to_global = local_node_to_global[p];

          for (size_t e = 0; e < edge_count; e++) {
            global_block_pos = proc_loc_edge_to_global[(e + boffset)] - goffset;
            global_pos       = global_block_pos * npe;

            for (size_t n = 0; n < npe; n++) {
              size_t node                 = proc_loc_node_to_global[local_linkage[pos++] - 1];
              block_linkage[global_pos++] = node + 1;
            }
          }

          // Get attributes list,  if it exists
          if (edgeblocks[p][b].attributeCount > 0) {
            size_t max_attr = edgeblocks[p][b].entity_count() * edgeblocks[p][b].attributeCount;
            std::vector<T> local_attr(max_attr);

            error = ex_get_attr(id, EX_EDGE_BLOCK, edgeblocks[p][b].id, Data(local_attr));
            if (error < 0) {
              exodus_error(__LINE__);
            }

            pos = 0;

            size_t att_count = edgeblocks[p][b].attributeCount;
            for (size_t e = 0; e < edge_count; e++) {
              // global_pos is global position within this element block...
              global_block_pos = local_edge_to_global[p][(e + boffset)] - goffset;
              global_pos       = global_block_pos * att_count;
              for (size_t n = 0; n < att_count; n++) {
                attributes[b][global_pos++] = local_attr[pos++];
              }
            }
          }

        } // end if edgeblocks[p][b].entity_count() (non-zero length block)
        else if (debug_level & 64) {
          fmt::print(".");
        }
      } // end for p=0..part_count-1

      // Write out block info
      int id_out = ExodusFile::output(); // output file identifier

      if (linkage[b] != nullptr) {
        error =
            ex_put_conn(id_out, EX_EDGE_BLOCK, glob_edgeblocks[b].id, linkage[b], nullptr, nullptr);
        if (error < 0) {
          exodus_error(__LINE__);
        }
        delete[] linkage[b];
      }

      // Write out attributes list if it exists
      if (glob_edgeblocks[b].attributeCount > 0) {
        error = ex_put_attr(id_out, EX_EDGE_BLOCK, glob_edgeblocks[b].id, attributes[b]);
        if (error < 0) {
          exodus_error(__LINE__);
        }
        delete[] attributes[b];
      } // end for b=0..global_num_edgeblocks-1
      if (debug_level & 64) {
        fmt::print("\n");
      }
    }
    fmt::print("\n");
    delete[] linkage;
    delete[] attributes;
  }

  template <typename T, typename INT>
  void put_edgeblocks(int part_count, int start_part,
                      std::vector<std::vector<EdgeBlock<INT>>> &edgeblocks,
                      std::vector<EdgeBlock<INT>>              &glob_edgeblocks,
                      const std::vector<std::vector<INT>>      &local_node_to_global,
                      T /* float_or_double */)
  {
    // This variant of `put_edgeblocks` is used in the case of
    // `nomap` and uses much less memory (but may be slower).  It
    // relies on the fact that with `nomap`, a the elements for a
    // block in a parts are all contiguous in the output file, so we
    // can read them, map the nodes, and then output them using a
    // partial write function and they don't need to be mapped into a
    // global array (which also does not need to be allocated).  Phase
    // 1 of a "low-memory" option for epu.

    SMART_ASSERT(sizeof(T) == ExodusFile::io_word_size());
    int global_num_edgeblocks = glob_edgeblocks.size();

    LOG("\nReading and Writing edge connectivity & attributes (Low Memory Method)\n");

    for (int b = 0; b < global_num_edgeblocks; b++) {

      if (debug_level & 64) {
        fmt::print("\nOutput edge block info for...\n"
                   "Edge Block {}, Id = {}, Name = '{}', Edges = {:12}, Nodes/Edge = {}, "
                   "Attributes = {}\n",
                   b, glob_edgeblocks[b].id, glob_edgeblocks[b].name_,
                   fmt::group_digits(glob_edgeblocks[b].entity_count()),
                   glob_edgeblocks[b].nodesPerEdge, glob_edgeblocks[b].attributeCount);
      }

      int id_out = ExodusFile::output(); // output file identifier

      size_t part_block_offset = 1;
      for (int p = 0; p < part_count; p++) {
        ExodusFile id(p);

        if (edgeblocks[p][b].entity_count() > 0) { // non-zero length block
          size_t node_count = edgeblocks[p][b].entity_count() * edgeblocks[p][b].nodesPerEdge;
          std::vector<INT> local_linkage(node_count);

          ex_entity_id bid = edgeblocks[p][b].id;
          int error = ex_get_conn(id, EX_EDGE_BLOCK, bid, Data(local_linkage), nullptr, nullptr);
          if (error < 0) {
            fmt::print(stderr,
                       "ERROR: (EPU) Cannot get edge block connectivity for block {} on part {}.\n",
                       bid, p + start_part);
            exodus_error(__LINE__);
          }

          size_t                  edge_count              = edgeblocks[p][b].entity_count();
          size_t                  npe                     = edgeblocks[p][b].nodesPerEdge;
          const std::vector<INT> &proc_loc_node_to_global = local_node_to_global[p];

          size_t pos = 0;
          for (size_t e = 0; e < edge_count; e++) {
            for (size_t n = 0; n < npe; n++) {
              size_t node          = proc_loc_node_to_global[local_linkage[pos] - 1];
              local_linkage[pos++] = node + 1;
            }
          }

          if (debug_level & 64) {
            fmt::print(stderr, "part, block, offset, count = {} {} {} {}\n", p, bid,
                       part_block_offset, fmt::group_digits(edge_count));
          }
          error = ex_put_partial_conn(id_out, EX_EDGE_BLOCK, bid, part_block_offset, edge_count,
                                      Data(local_linkage), nullptr, nullptr);
          if (error < 0) {
            fmt::print(stderr,
                       "ERROR: (EPU) Cannot output edge block connectivity for block {} on part "
                       "{} (offset = {}, count = {}).\n",
                       bid, p + start_part, part_block_offset, fmt::group_digits(edge_count));
            exodus_error(__LINE__);
          }

          // Get attributes list,  if it exists
          if (edgeblocks[p][b].attributeCount > 0) {
            size_t max_attr = edgeblocks[p][b].entity_count() * edgeblocks[p][b].attributeCount;
            std::vector<T> local_attr(max_attr);

            error = ex_get_attr(id, EX_EDGE_BLOCK, edgeblocks[p][b].id, Data(local_attr));
            if (error < 0) {
              exodus_error(__LINE__);
            }

            error = ex_put_partial_attr(id_out, EX_EDGE_BLOCK, bid, part_block_offset, edge_count,
                                        Data(local_attr));
            if (error < 0) {
              exodus_error(__LINE__);
            }
          }
          part_block_offset += edge_count;

        } // end if edgeblocks[p][b].entity_count() (non-zero length block)
      } // end for p=0..part_count-1
    }
  }

  template <typename INT>
  void get_faceblocks(int part_count, const std::vector<Mesh> &local_mesh, const Mesh &global,
                      std::vector<std::vector<FaceBlock<INT>>> &faceblocks,
                      std::vector<FaceBlock<INT>>              &glob_faceblocks)
  {
    LOG("\n\n**** GET FACE BLOCK INFORMATION (INCL. ELEMENT ATTRIBUTES) ****\n");

    for (int ip = 0; ip < part_count; ip++) {
      faceblocks[ip].resize(local_mesh[ip].count(Excn::ObjectType::FABLK));
    }

    if (rank == 0) {
      fmt::print("Global face block count = {}\n", global.count(Excn::ObjectType::FABLK));
    }

    if (global.count(Excn::ObjectType::FABLK) == 0) {
      return;
    }

    ExodusIdVector faceblock_id(global.count(Excn::ObjectType::FABLK));

    int error = 0;
    for (int p = 0; p < part_count; p++) {
      ExodusFile id(p);

      error = ex_get_ids(id, EX_FACE_BLOCK, Data(faceblock_id));
      if (error < 0) {
        exodus_error(__LINE__);
      }

      // Check that the block id ordering is consistent among files...
      if (p > 0) {
        for (size_t b = 0; b < global.count(Excn::ObjectType::FABLK); b++) {
          if (faceblocks[0][b].id != faceblock_id[b]) {
            std::ostringstream errmsg;
            fmt::print(errmsg,
                       "ERROR: (EPU) The internal face block id ordering for part {}\n"
                       "       is not consistent with the ordering for part 0.\n",
                       p);
            throw std::runtime_error(errmsg.str());
          }
        }
      }

      if ((debug_level & 64) != 0U) {
        fmt::print("\nGetting face block info for processor {}...\n", p);
      }
      else {
        if (p == 0) {
          LOG("\nGetting face block info.\n");
        }
      }

      for (size_t b = 0; b < global.count(Excn::ObjectType::FABLK); b++) {
        if ((debug_level & 64) != 0U) {
          fmt::print("Face Block {}, Id = {}", b, faceblock_id[b]);
        }

        ex_block temp_block{};
        temp_block.id   = faceblock_id[b];
        temp_block.type = EX_FACE_BLOCK;
        error           = ex_get_block_param(id, &temp_block);
        if (error < 0) {
          exodus_error(__LINE__);
        }

        std::vector<char> name(Excn::ExodusFile::max_name_length() + 1);
        error = ex_get_name(id, EX_FACE_BLOCK, faceblock_id[b], Data(name));
        if (error < 0) {
          exodus_error(__LINE__);
        }

        faceblocks[p][b].id = faceblock_id[b];
        if (name[0] != '\0') {
          faceblocks[p][b].name_ = Data(name);
        }
        if (p == 0) {
          glob_faceblocks[b].id = faceblock_id[b];
          if (name[0] != '\0') {
            glob_faceblocks[b].name_ = Data(name);
          }
        }

        if (temp_block.num_entry != 0) {
          faceblocks[p][b].faceCount      = temp_block.num_entry;
          faceblocks[p][b].nodesPerFace   = temp_block.num_nodes_per_entry;
          faceblocks[p][b].attributeCount = temp_block.num_attribute;
          faceblocks[p][b].offset_        = temp_block.num_entry;
          faceblocks[p][b].position_      = b;
          copy_string(faceblocks[p][b].elType, temp_block.topology);

          glob_faceblocks[b].faceCount += temp_block.num_entry;
          glob_faceblocks[b].nodesPerFace   = temp_block.num_nodes_per_entry;
          glob_faceblocks[b].attributeCount = temp_block.num_attribute;
          glob_faceblocks[b].position_      = (int)b;
          copy_string(glob_faceblocks[b].elType, temp_block.topology);
        }

        if (temp_block.num_attribute > 0 && glob_faceblocks[b].attributeNames.empty()) {
          // Get attribute names.  Assume the same on all processors
          // on which the block exists.
          char **names = get_name_array(temp_block.num_attribute, ExodusFile::max_name_length());

          error = ex_get_attr_names(id, EX_FACE_BLOCK, faceblock_id[b], names);
          if (error < 0) {
            exodus_error(__LINE__);
          }
          for (int i = 0; i < temp_block.num_attribute; i++) {
            glob_faceblocks[b].attributeNames.emplace_back(names[i]);
          }
          free_name_array(names, temp_block.num_attribute);
        }
        if ((debug_level & 4) != 0U) {
          fmt::print(", Name = '{}', Faces = {:12}, Nodes/Face= {}, Attributes = {}\n",
                     faceblocks[p][b].name_, fmt::group_digits(faceblocks[p][b].entity_count()),
                     faceblocks[p][b].nodesPerFace, faceblocks[p][b].attributeCount);
        }
      }
    } // end for p=0..part_count

    // Convert block_offset from faces/block/processor to true offset
    for (int p = 0; p < part_count; p++) {
      size_t sum = 0;
      for (size_t b = 0; b < global.count(Excn::ObjectType::FABLK); b++) {
        size_t save              = faceblocks[p][b].offset_;
        faceblocks[p][b].offset_ = sum;
        sum += save;
      }
    }
  }

  template <typename T, typename INT>
  void put_faceblocks(int part_count, int start_part,
                      std::vector<std::vector<FaceBlock<INT>>> &faceblocks,
                      std::vector<FaceBlock<INT>>              &glob_faceblocks,
                      const std::vector<std::vector<INT>>      &local_node_to_global,
                      const std::vector<std::vector<INT>>      &local_face_to_global,
                      T /* float_or_double */)
  {
    SMART_ASSERT(sizeof(T) == ExodusFile::io_word_size());
    int global_num_faceblocks = glob_faceblocks.size();

    auto linkage    = new INT *[global_num_faceblocks];
    auto attributes = new T *[global_num_faceblocks];

    LOG("\nReading and Writing face connectivity & attributes\n");

    for (int b = 0; b < global_num_faceblocks; b++) {

      if (debug_level & 64) {
        fmt::print("\nOutput face block info for...\n"
                   "Face Block {}, Id = {}, Name = '{}', Faces = {:12}, Nodes/Face = {}, "
                   "Attributes = {}\n"
                   "B{}:\t",
                   b, glob_faceblocks[b].id, glob_faceblocks[b].name_,
                   fmt::group_digits(glob_faceblocks[b].entity_count()),
                   glob_faceblocks[b].nodesPerFace, glob_faceblocks[b].attributeCount, b);
      }

      size_t max_nodes = glob_faceblocks[b].entity_count();
      max_nodes *= glob_faceblocks[b].nodesPerFace;

      if (max_nodes > 0) {
        linkage[b] = new INT[max_nodes];
      }
      else {
        linkage[b] = nullptr;
      }
      INT *block_linkage = linkage[b];

      // Initialize attributes list, if it exists
      if (glob_faceblocks[b].attributeCount > 0) {
        attributes[b] = new T[static_cast<size_t>(glob_faceblocks[b].attributeCount) *
                              glob_faceblocks[b].entity_count()];
      }

      int error = 0;
      for (int p = 0; p < part_count; p++) {
        ExodusFile id(p);

        size_t global_pos;
        size_t global_block_pos;

        if (faceblocks[p][b].entity_count() > 0) { // non-zero length block

          if (debug_level & 64) {
            fmt::print("#");
          }
          size_t maximum_nodes = faceblocks[p][b].entity_count();
          maximum_nodes *= faceblocks[p][b].nodesPerFace;
          std::vector<INT> local_linkage(maximum_nodes);
          SMART_ASSERT(block_linkage != nullptr);

          ex_entity_id bid = faceblocks[p][b].id;
          error = ex_get_conn(id, EX_FACE_BLOCK, bid, Data(local_linkage), nullptr, nullptr);
          if (error < 0) {
            fmt::print(stderr,
                       "ERROR: (EPU) Cannot get face block connectivity for block {} on part {}.\n",
                       bid, p + start_part);
            exodus_error(__LINE__);
          }
          size_t      pos                     = 0;
          size_t      goffset                 = glob_faceblocks[b].offset_;
          size_t      face_count              = faceblocks[p][b].entity_count();
          size_t      boffset                 = faceblocks[p][b].offset_;
          size_t      npe                     = faceblocks[p][b].nodesPerFace;
          const auto &proc_loc_face_to_global = local_face_to_global[p];
          const auto &proc_loc_node_to_global = local_node_to_global[p];

          for (size_t e = 0; e < face_count; e++) {
            global_block_pos = proc_loc_face_to_global[(e + boffset)] - goffset;
            global_pos       = global_block_pos * npe;

            for (size_t n = 0; n < npe; n++) {
              size_t node                 = proc_loc_node_to_global[local_linkage[pos++] - 1];
              block_linkage[global_pos++] = node + 1;
            }
          }

          // Get attributes list,  if it exists
          if (faceblocks[p][b].attributeCount > 0) {
            size_t max_attr = faceblocks[p][b].entity_count() * faceblocks[p][b].attributeCount;
            std::vector<T> local_attr(max_attr);

            error = ex_get_attr(id, EX_FACE_BLOCK, faceblocks[p][b].id, Data(local_attr));
            if (error < 0) {
              exodus_error(__LINE__);
            }

            pos = 0;

            size_t att_count = faceblocks[p][b].attributeCount;
            for (size_t e = 0; e < face_count; e++) {
              // global_pos is global position within this element block...
              global_block_pos = local_face_to_global[p][(e + boffset)] - goffset;
              global_pos       = global_block_pos * att_count;
              for (size_t n = 0; n < att_count; n++) {
                attributes[b][global_pos++] = local_attr[pos++];
              }
            }
          }

        } // end if faceblocks[p][b].entity_count() (non-zero length block)
        else if (debug_level & 64) {
          fmt::print(".");
        }
      } // end for p=0..part_count-1

      // Write out block info
      int id_out = ExodusFile::output(); // output file identifier

      if (linkage[b] != nullptr) {
        fmt::print(stderr, "(EPU) id_out({}) glob_faceblocks[{}].id({})\n", id_out, b,
                   glob_faceblocks[b].id);
        error =
            ex_put_conn(id_out, EX_FACE_BLOCK, glob_faceblocks[b].id, linkage[b], nullptr, nullptr);
        if (error < 0) {
          exodus_error(__LINE__);
        }
        delete[] linkage[b];
      }

      // Write out attributes list if it exists
      if (glob_faceblocks[b].attributeCount > 0) {
        error = ex_put_attr(id_out, EX_FACE_BLOCK, glob_faceblocks[b].id, attributes[b]);
        if (error < 0) {
          exodus_error(__LINE__);
        }
        delete[] attributes[b];
      } // end for b=0..global_num_faceblocks-1
      if (debug_level & 64) {
        fmt::print("\n");
      }
    }
    fmt::print("\n");
    delete[] linkage;
    delete[] attributes;
  }

  template <typename T, typename INT>
  void put_faceblocks(int part_count, int start_part,
                      std::vector<std::vector<FaceBlock<INT>>> &faceblocks,
                      std::vector<FaceBlock<INT>>              &glob_faceblocks,
                      const std::vector<std::vector<INT>>      &local_node_to_global,
                      T /* float_or_double */)
  {
    // This variant of `put_faceblocks` is used in the case of
    // `nomap` and uses much less memory (but may be slower).  It
    // relies on the fact that with `nomap`, a the elements for a
    // block in a parts are all contiguous in the output file, so we
    // can read them, map the nodes, and then output them using a
    // partial write function and they don't need to be mapped into a
    // global array (which also does not need to be allocated).  Phase
    // 1 of a "low-memory" option for epu.

    SMART_ASSERT(sizeof(T) == ExodusFile::io_word_size());
    int global_num_faceblocks = glob_faceblocks.size();

    LOG("\nReading and Writing face connectivity & attributes (Low Memory Method)\n");

    for (int b = 0; b < global_num_faceblocks; b++) {

      if (debug_level & 64) {
        fmt::print("\nOutput face block info for...\n"
                   "Face Block {}, Id = {}, Name = '{}', Faces = {:12}, Nodes/Face = {}, "
                   "Attributes = {}\n",
                   b, glob_faceblocks[b].id, glob_faceblocks[b].name_,
                   fmt::group_digits(glob_faceblocks[b].entity_count()),
                   glob_faceblocks[b].nodesPerFace, glob_faceblocks[b].attributeCount);
      }

      int id_out = ExodusFile::output(); // output file identifier

      size_t part_block_offset = 1;
      for (int p = 0; p < part_count; p++) {
        ExodusFile id(p);

        if (faceblocks[p][b].entity_count() > 0) { // non-zero length block
          size_t node_count = faceblocks[p][b].entity_count() * faceblocks[p][b].nodesPerFace;
          std::vector<INT> local_linkage(node_count);

          ex_entity_id bid = faceblocks[p][b].id;
          int error = ex_get_conn(id, EX_FACE_BLOCK, bid, Data(local_linkage), nullptr, nullptr);
          if (error < 0) {
            fmt::print(stderr,
                       "ERROR: (EPU) Cannot get face block connectivity for block {} on part {}.\n",
                       bid, p + start_part);
            exodus_error(__LINE__);
          }

          size_t                  face_count              = faceblocks[p][b].entity_count();
          size_t                  npe                     = faceblocks[p][b].nodesPerFace;
          const std::vector<INT> &proc_loc_node_to_global = local_node_to_global[p];

          size_t pos = 0;
          for (size_t e = 0; e < face_count; e++) {
            for (size_t n = 0; n < npe; n++) {
              size_t node          = proc_loc_node_to_global[local_linkage[pos] - 1];
              local_linkage[pos++] = node + 1;
            }
          }

          if (debug_level & 64) {
            fmt::print(stderr, "part, block, offset, count = {} {} {} {}\n", p, bid,
                       part_block_offset, face_count);
          }
          error = ex_put_partial_conn(id_out, EX_FACE_BLOCK, bid, part_block_offset, face_count,
                                      Data(local_linkage), nullptr, nullptr);
          if (error < 0) {
            fmt::print(stderr,
                       "ERROR: (EPU) Cannot output face block connectivity for block {} on part "
                       "{} (offset = {}, count = {}).\n",
                       bid, p + start_part, part_block_offset, face_count);
            exodus_error(__LINE__);
          }

          // Get attributes list,  if it exists
          if (faceblocks[p][b].attributeCount > 0) {
            size_t max_attr = faceblocks[p][b].entity_count() * faceblocks[p][b].attributeCount;
            std::vector<T> local_attr(max_attr);

            error = ex_get_attr(id, EX_FACE_BLOCK, faceblocks[p][b].id, Data(local_attr));
            if (error < 0) {
              exodus_error(__LINE__);
            }

            error = ex_put_partial_attr(id_out, EX_FACE_BLOCK, bid, part_block_offset, face_count,
                                        Data(local_attr));
            if (error < 0) {
              exodus_error(__LINE__);
            }
          }
          part_block_offset += face_count;

        } // end if faceblocks[p][b].entity_count() (non-zero length block)
      } // end for p=0..part_count-1
    }
  }

  template <typename T, typename INT>
  void add_processor_variable(int id_out, int part_count, int start_part, const Mesh &global,
                              std::vector<std::vector<Block>>     &blocks,
                              const std::vector<Block>            &glob_blocks,
                              const std::vector<std::vector<INT>> &local_element_to_global,
                              int step, int variable, std::vector<T> &proc)
  {
    SMART_ASSERT(sizeof(T) == ExodusFile::io_word_size());
    for (size_t b = 0; b < global.count(Excn::ObjectType::EBLK); b++) {
      proc.resize(glob_blocks[b].entity_count());
      for (int p = 0; p < part_count; p++) {
        size_t boffset       = blocks[p][b].offset_;
        size_t goffset       = glob_blocks[b].offset_;
        size_t element_count = blocks[p][b].entity_count();
        for (size_t e = 0; e < element_count; e++) {
          size_t global_block_pos = local_element_to_global[p][(e + boffset)] - goffset;
          proc[global_block_pos]  = p + start_part;
        }
      }
      int error = ex_put_var(id_out, step, EX_ELEM_BLOCK, variable, glob_blocks[b].id,
                             glob_blocks[b].entity_count(), Data(proc));
      if (error < 0) {
        exodus_error(__LINE__);
      }
    }
  }

  template <typename INT>
  void add_processor_map(int id_out, int part_count, int start_part, const Mesh &global,
                         std::vector<std::vector<Block>>     &blocks,
                         const std::vector<Block>            &glob_blocks,
                         const std::vector<std::vector<INT>> &local_element_to_global)
  {
    std::vector<INT> proc(global.elementCount);

    for (size_t b = 0; b < global.count(Excn::ObjectType::EBLK); b++) {
      proc.resize(glob_blocks[b].entity_count());
      for (int p = 0; p < part_count; p++) {
        size_t boffset       = blocks[p][b].offset_;
        size_t element_count = blocks[p][b].entity_count();
        for (size_t e = 0; e < element_count; e++) {
          size_t global_elem_pos = local_element_to_global[p][(e + boffset)];
          proc[global_elem_pos]  = p + start_part;
        }
      }
    }

    if (ex_put_map_param(id_out, 0, 1) < 0) {
      exodus_error(__LINE__);
    }
    if (ex_put_num_map(id_out, EX_ELEM_MAP, 1, Data(proc)) < 0) {
      exodus_error(__LINE__);
    }
    if (ex_put_name(id_out, EX_ELEM_MAP, 1, "processor_id") < 0) {
      exodus_error(__LINE__);
    }
  }

  StringVector get_exodus_variable_names(int id, ex_entity_type elType, int var_count)
  {
    // Allocate space for variable names...
    char **name_list = get_name_array(var_count, ExodusFile::max_name_length());
    int    error     = ex_get_variable_names(id, elType, var_count, name_list);
    if (error < 0) {
      exodus_error(__LINE__);
    }

    StringVector names(var_count);
    for (int j = 0; j < var_count; j++) {
      compress_white_space(name_list[j]);
      names[j] = std::string(name_list[j]);
    }

    free_name_array(name_list, var_count);
    return names;
  }

  template <typename T>
  void filter_truth_table(int id, Mesh &global, std::vector<T> &glob_blocks, Variables &vars,
                          const StringIdVector &variable_names)
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
    for (const auto &[variable_name, variable_block] : variable_names) {
      if (variable_block > 0) {
        if (var_name != variable_name) {
          var_name = variable_name;
          // Find which exodus variable matches this name
          out_position = -1;
          for (size_t j = 0; j < exo_names.size(); j++) {
            if (case_compare(exo_names[j], var_name)) {
              out_position = vars.index_[j] - 1;
              break;
            }
          }
          SMART_ASSERT(out_position >= 0);

          // Set all truth table entries for this variable to negative
          // of current value and then iterate over specified blocks and
          // set those positive.  This way can make sure that the
          // variable truly exists for the block that the user specified.
          for (size_t b = 0; b < global.count(vars.objectType); b++) {
            int truth_table_loc = (b * vars.count(InOut::OUT)) + out_position;
            global.truthTable[static_cast<int>(vars.objectType)][truth_table_loc] *= -1;
          }
        }
        // Find out which block corresponds to the specified id.
        int block = -1;
        for (size_t b = 0; b < global.count(vars.objectType); b++) {
          if (glob_blocks[b].id == variable_block) {
            block = b;
            break;
          }
        }

        if (block == -1) {
          std::ostringstream errmsg;
          fmt::print(
              errmsg,
              "ERROR: (EPU) User-specified block id of {} for variable '{}' does not exist.\n",
              variable_block, variable_name);
          throw std::runtime_error(errmsg.str());
        }

        int truth_table_loc = block * vars.count(InOut::OUT) + out_position;
        if (global.truthTable[static_cast<int>(vars.objectType)][truth_table_loc] == 0) {
          std::ostringstream errmsg;
          fmt::print(errmsg, "ERROR: (EPU) Variable '{}' does not exist on block {}.\n",
                     variable_name, variable_block);
          throw std::runtime_error(errmsg.str());
        }
        global.truthTable[static_cast<int>(vars.objectType)][truth_table_loc] = 1;
      }
    }

    // reset truth table values that may be negative
    int output_truth_table_length = vars.count(InOut::OUT) * global.count(vars.objectType);
    for (int j = 0; j < output_truth_table_length; j++) {
      if (global.truthTable[static_cast<int>(vars.objectType)][j] < 0) {
        global.truthTable[static_cast<int>(vars.objectType)][j] = 0;
      }
    }
  }

  template <typename T>
  void get_truth_table(Mesh &global, std::vector<T> &glob_blocks, std::vector<Mesh> &local,
                       Variables &vars, int debug)
  {
    // read truth table - sum across all processors since many will
    // have values set to zero for zero length blocks the element
    // variable truth table organized as a 2D array:
    // [global.count(EBLK)][num_elem_vars]

    ObjectType object_type               = vars.objectType;
    int        input_truth_table_length  = vars.count(InOut::IN) * global.count(object_type);
    int        output_truth_table_length = vars.count(InOut::OUT) * global.count(object_type);

    if (output_truth_table_length) {

      global.truthTable[static_cast<int>(object_type)].resize(output_truth_table_length);
      std::fill(global.truthTable[static_cast<int>(object_type)].begin(),
                global.truthTable[static_cast<int>(object_type)].end(), 0);

      // For each input exodus file, get it's truth table and fill
      // in the location in the output truth table...

      bool is_sidenodeset =
          vars.objectType == Excn::ObjectType::NSET || vars.objectType == Excn::ObjectType::SSET;
      int part_count = local.size();
      for (int p = 0; p < part_count; p++) {
        ExodusFile id(p);

        if (vars.count(InOut::IN) > 0) { // Could be zero if add_processor_id
          // is the only variable...
          local[p].truthTable[static_cast<int>(object_type)].resize(input_truth_table_length);
          int error =
              ex_get_truth_table(id, vars.type(), global.count(object_type), vars.count(InOut::IN),
                                 Data(local[p].truthTable[static_cast<int>(object_type)]));
          if (error < 0) {
            exodus_error(__LINE__);
          }
        }
        for (size_t b = 0; b < global.count(object_type); b++) {
          size_t bin = b;
          if (is_sidenodeset) {
            bin = glob_blocks[b].position_;
          }

          for (int j = 0; j < vars.count(InOut::IN); j++) {
            if (vars.index_[j] > 0) {
              int ki = (bin * vars.count(InOut::IN)) + j;
              int ko = (b * vars.count(InOut::OUT)) + vars.index_[j] - 1;
              SMART_ASSERT(ko < output_truth_table_length);
              SMART_ASSERT(ki < input_truth_table_length);
              global.truthTable[static_cast<int>(object_type)][ko] +=
                  local[p].truthTable[static_cast<int>(object_type)][ki];
            }
          }
          if (vars.addProcessorId) {
            int ko = (b * vars.count(InOut::OUT)) + vars.count(InOut::OUT) - 1;
            SMART_ASSERT(ko < output_truth_table_length);
            global.truthTable[static_cast<int>(object_type)][ko] = 1;
          }
        }
      }

      // reset truth table values that may be greater than 1
      for (int j = 0; j < output_truth_table_length; j++) {
        SMART_ASSERT(global.truthTable[static_cast<int>(object_type)][j] >= 0)
        (global.truthTable[static_cast<int>(object_type)][j]);
        if (global.truthTable[static_cast<int>(object_type)][j] > 0) {
          global.truthTable[static_cast<int>(object_type)][j] = 1;
        }
      }

      if (debug_level & debug) {
        fmt::print("Truth table for {}\n", vars.label());
        int k = 0;
        for (size_t b = 0; b < global.count(object_type); b++) {
          for (int j = 0; j < vars.count(InOut::OUT); j++) {
            fmt::print("{}", global.truthTable[static_cast<int>(object_type)][k++]);
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
    std::string info = sys_info("EPU");
    copy_string(info_record, info, size + 1);
  }

  inline bool is_whitespace(char c)
  {
    static char white_space[] = {' ', '\t', '\n', '\r', ',', '\0'};
    return strchr(white_space, c) != nullptr;
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

  int get_width(int max_value)
  {
    // Returns the field width which will accommodate the
    // largest value.
    int width = 1;
    if (max_value >= 10) {
      width = int(std::log10(static_cast<double>(max_value)));
    }
    return width + 1;
  }

  template <typename T, typename INT>
  void map_element_vars(size_t loffset, size_t goffset, size_t entity_count, std::vector<T> &values,
                        T *global_values, const std::vector<INT> &proc_loc_elem_to_global)
  {
    // copy values to master element value information
    T *local_values = Data(values);
    for (size_t j = 0; j < entity_count; j++) {
      size_t global_block_pos         = proc_loc_elem_to_global[(j + loffset)] - goffset;
      global_values[global_block_pos] = local_values[j];
    }
  }

  template <typename T>
  void map_sideset_vars(size_t loffset, size_t entity_count, std::vector<T> &values,
                        T *global_values)
  {
    // copy values to master sideset value information
    T *local_values = Data(values);
    for (size_t j = 0; j < entity_count; j++) {
      global_values[j + loffset] = local_values[j];
    }
  }

  template <typename T, typename U>
  void map_nodeset_vars(U & /*unused*/, size_t /*unused*/, size_t /*unused*/,
                        std::vector<T> & /*unused*/, T * /*unused*/)
  {
    throw std::runtime_error("Internal Error!");
  }

  template <typename INT>
  void map_nodeset_vars(Excn::NodeSet<INT> &local_set, size_t entity_count,
                        size_t glob_entity_count, std::vector<double> &values,
                        double *global_values)
  {
    // copy values to master nodeset value information
    double *local_values = Data(values);
    for (size_t j = 0; j < entity_count; j++) {
      size_t global_loc = local_set.nodeOrderMap[j];
      SMART_ASSERT(global_loc < glob_entity_count);
      global_values[global_loc] = local_values[j];
    }
  }

  template <typename INT>
  void map_nodeset_vars(Excn::NodeSet<INT> &local_set, size_t entity_count,
                        size_t glob_entity_count, std::vector<float> &values, float *global_values)
  {
    // copy values to master nodeset value information
    float *local_values = Data(values);
    for (size_t j = 0; j < entity_count; j++) {
      size_t global_loc = local_set.nodeOrderMap[j];
      SMART_ASSERT(global_loc < glob_entity_count);
      global_values[global_loc] = local_values[j];
    }
  }

  template <typename T, typename INT>
  void map_edgeblock_vars(size_t loffset, size_t goffset, size_t entity_count,
                          std::vector<T> &values, T *global_values,
                          const std::vector<INT> &proc_loc_elem_to_global)
  {
    // copy values to master edgeblock value information
    T *local_values = Data(values);
    for (size_t j = 0; j < entity_count; j++) {
      size_t global_block_pos         = proc_loc_elem_to_global[(j + loffset)] - goffset;
      global_values[global_block_pos] = local_values[j];
    }
  }

  template <typename T, typename U, typename INT>
  void read_write_master_values(Excn::Variables &vars, const Excn::Mesh &global,
                                std::vector<U> &global_sets, std::vector<Excn::Mesh> &local_mesh,
                                std::vector<std::vector<U>> &local_sets,
                                std::vector<T> &master_values, std::vector<T> &values,
                                int part_count, int time_step, int time_step_out,
                                const std::vector<std::vector<INT>> &local_entity_to_global)
  {
    bool is_sidenodeset =
        vars.objectType == Excn::ObjectType::NSET || vars.objectType == Excn::ObjectType::SSET;

    int id_out = ExodusFile::output(); // output file identifier

    for (int i = 0; i < vars.count(InOut::IN); i++) {
      if (vars.index_[i] > 0) {
        int ivar = vars.index_[i] - 1;

        for (size_t b = 0; b < global.count(vars.objectType); b++) {
          size_t bin = b;
          if (is_sidenodeset) {
            bin = global_sets[b].position_;
          }
          int output_truth_table_loc = (b * vars.count(InOut::OUT)) + ivar;
          int input_truth_table_loc  = (bin * vars.count(InOut::IN)) + i;
          if (debug_level & 4) {
            std::fill(master_values.begin(), master_values.end(), T(0.0));
          }

          for (int p = 0; p < part_count; p++) {
            ExodusFile id(p);

            const std::vector<INT> &proc_loc_elem_to_global = local_entity_to_global[p];
            const std::vector<INT> &proc_loc_edge_to_global = local_entity_to_global[p];
            const std::vector<INT> &proc_loc_face_to_global = local_entity_to_global[p];

            if (global.truthTable[static_cast<int>(vars.objectType)][output_truth_table_loc] &&
                local_sets[p][b].entity_count() > 0) {

              T     *iv_block_mev = Data(master_values);
              size_t entity_count = local_sets[p][b].entity_count();

              if (local_mesh[p]
                      .truthTable[static_cast<int>(vars.objectType)][input_truth_table_loc] > 0) {
                int error = ex_get_var(id, time_step + 1, exodus_object_type(vars.objectType),
                                       i + 1, local_sets[p][b].id, entity_count, Data(values));
                if (error < 0) {
                  exodus_error(__LINE__);
                }

                switch (vars.objectType) {
                case Excn::ObjectType::EBLK:
                  map_element_vars(local_sets[p][b].offset_, global_sets[b].offset_, entity_count,
                                   values, iv_block_mev, proc_loc_elem_to_global);
                  break;

                case Excn::ObjectType::SSET:
                  map_sideset_vars(local_sets[p][b].offset_, entity_count, values, iv_block_mev);
                  break;

                case Excn::ObjectType::NSET:
                  map_nodeset_vars(local_sets[p][b], entity_count, global_sets[b].entity_count(),
                                   values, iv_block_mev);
                  break;

                case Excn::ObjectType::EDBLK:
                  map_edgeblock_vars(local_sets[p][b].offset_, global_sets[b].offset_, entity_count,
                                     values, iv_block_mev, proc_loc_edge_to_global);
                  break;

                case Excn::ObjectType::FABLK:
                  map_edgeblock_vars(local_sets[p][b].offset_, global_sets[b].offset_, entity_count,
                                     values, iv_block_mev, proc_loc_face_to_global);
                  break;
                default: break;
                }
              }
            }
          }
          int truth_table_loc = (b * vars.count(InOut::OUT)) + ivar;

          if (global.truthTable[static_cast<int>(vars.objectType)][truth_table_loc]) {
            int error =
                ex_put_var(id_out, time_step_out, exodus_object_type(vars.objectType), ivar + 1,
                           global_sets[b].id, global_sets[b].entity_count(), Data(master_values));
            if (error < 0) {
              exodus_error(__LINE__);
            }
          }
        }
      }
    }
  }

  template <typename U>
  void create_output_truth_table(const Excn::Mesh &global, std::vector<U> &global_sets,
                                 Excn::Variables &vars, std::vector<int> &truth_table)
  {
    for (size_t b = 0; b < global.count(vars.objectType); b++) {
      int bout = global_sets[b].position_;
      SMART_ASSERT(bout >= 0);
      for (int j = 0; j < vars.count(InOut::OUT); j++) {
        int inp_ttable_loc = (b * vars.count(InOut::OUT)) + j;
        int out_ttable_loc = (bout * vars.count(InOut::OUT)) + j;
        truth_table[out_ttable_loc] =
            global.truthTable[static_cast<int>(vars.objectType)][inp_ttable_loc];
      }
    }
  }

  template <typename INT>
  size_t find_max_entity_count(int part_count, std::vector<Excn::Mesh> &local_mesh,
                               const Excn::Mesh &global, std::vector<std::vector<Block>> &blocks,
                               std::vector<std::vector<NodeSet<INT>>>   &nodesets,
                               std::vector<std::vector<SideSet<INT>>>   &sidesets,
                               std::vector<std::vector<EdgeBlock<INT>>> &edgeblocks,
                               std::vector<std::vector<FaceBlock<INT>>> &faceblocks)
  {
    size_t max_ent = local_mesh[0].nodeCount;
    for (int p = 1; p < part_count; p++) {
      if (static_cast<size_t>(local_mesh[p].nodeCount) > max_ent) {
        max_ent = local_mesh[p].nodeCount;
      }
    }

    for (int p = 0; p < part_count; p++) {
      for (size_t b = 0; b < global.count(Excn::ObjectType::EBLK); b++) {
        if (blocks[p][b].entity_count() > max_ent) {
          max_ent = blocks[p][b].entity_count();
        }
      }
    }

    // Nodesets...
    for (int p = 0; p < part_count; p++) {
      for (size_t b = 0; b < global.count(Excn::ObjectType::NSET); b++) {
        if (nodesets[p][b].entity_count() > max_ent) {
          max_ent = nodesets[p][b].entity_count();
        }
      }
    }

    // Sidesets...
    for (int p = 0; p < part_count; p++) {
      for (size_t b = 0; b < global.count(Excn::ObjectType::SSET); b++) {
        if (sidesets[p][b].entity_count() > max_ent) {
          max_ent = sidesets[p][b].entity_count();
        }
      }
    }

    // Edgeblocks...
    for (int p = 0; p < part_count; p++) {
      for (size_t b = 0; b < global.count(Excn::ObjectType::EDBLK); b++) {
        if (edgeblocks[p][b].entity_count() > max_ent) {
          max_ent = edgeblocks[p][b].entity_count();
        }
      }
    }

    // Faceblocks...
    for (int p = 0; p < part_count; p++) {
      for (size_t b = 0; b < global.count(Excn::ObjectType::FABLK); b++) {
        if (faceblocks[p][b].entity_count() > max_ent) {
          max_ent = faceblocks[p][b].entity_count();
        }
      }
    }
    return max_ent;
  }

  template <typename INT>
  size_t find_max_global_entity_count(const Excn::Mesh &global, std::vector<Block> &blocks,
                                      std::vector<NodeSet<INT>>   &nodesets,
                                      std::vector<SideSet<INT>>   &sidesets,
                                      std::vector<EdgeBlock<INT>> &edgeblocks,
                                      std::vector<FaceBlock<INT>> &faceblocks)
  {
    size_t max_ent = global.nodeCount;

    for (size_t b = 0; b < global.count(Excn::ObjectType::EBLK); b++) {
      if (blocks[b].entity_count() > max_ent) {
        max_ent = blocks[b].entity_count();
      }
    }

    // Nodesets...
    for (size_t b = 0; b < global.count(Excn::ObjectType::NSET); b++) {
      if (nodesets[b].entity_count() > max_ent) {
        max_ent = nodesets[b].entity_count();
      }
    }

    // Sidesets...
    for (size_t b = 0; b < global.count(Excn::ObjectType::SSET); b++) {
      if (sidesets[b].entity_count() > max_ent) {
        max_ent = sidesets[b].entity_count();
      }
    }

    // Edge Blocks...
    for (size_t b = 0; b < global.count(Excn::ObjectType::EDBLK); b++) {
      if (edgeblocks[b].entity_count() > max_ent) {
        max_ent = edgeblocks[b].entity_count();
      }
    }

    // Face Blocks...
    for (size_t b = 0; b < global.count(Excn::ObjectType::FABLK); b++) {
      if (faceblocks[b].entity_count() > max_ent) {
        max_ent = faceblocks[b].entity_count();
      }
    }
    return max_ent;
  }
} // namespace
