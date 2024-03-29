// Copyright(C) 1999-2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fmt/chrono.h>
#include <fmt/ostream.h>
#include <fstream>
#include <iostream>
#include <numeric>

#include "ED_SystemInterface.h"
#include "ED_Version.h"
#include "FileInfo.h"
#include "MinMaxData.h"
#include "Norm.h"
#include "Tolerance.h"
#include "edge_block.h"
#include "exoII_read.h"
#include "exo_block.h"
#include "exodiff.h"
#include "exodusII.h"
#include "face_block.h"
#include "map.h"
#include "node_set.h"
#include "side_set.h"
#include "smart_assert.h"
#include "stringx.h"
#include "util.h"

#include "add_to_log.h"

SystemInterface interFace;

struct TimeInterp
{
  TimeInterp() = default;

  int step1{-1}; // step at beginning of interval. -1 if time prior to time at step1
  int step2{-1}; // step at end of interval. -1 if time after time at step2

  double time{0.0}; // Time being interpolated to.

  // If t1 = time at step1 and t2 = time at step2,
  // then proportion = (time-t1)/(t2-t1)
  // Or, value at time = (1.0-proportion)*v1 + proportion*v2
  double proportion{0.0};
};

std::string Date()
{
  time_t calendar_time = time(nullptr);
#if defined __NVCC__
  char       tbuf[32];
  struct tm *local_time = localtime(&calendar_time);
  strftime(tbuf, 32, "%Y/%m/%d   %H:%M:%S %Z", local_time);
  std::string time_string(tbuf);
#else
  auto const local_time  = fmt::localtime(calendar_time);
  auto       time_string = fmt::format("{:%Y/%m/%d   %H:%M:%S %Z}", local_time);
#endif
  return time_string;
}

bool Invalid_Values(const double *values, size_t count);
bool Equal_Values(const double *values, size_t count, double *value);

void Print_Banner(const char *prefix)
{
  fmt::print("\n"
             "{0}  *****************************************************************\n"
             "{0}             ",
             prefix);
  SystemInterface::show_version();
  fmt::print("{0}             Authors:  Richard Drake, rrdrake@sandia.gov           \n"
             "{0}                       Greg Sjaardema, gdsjaar@sandia.gov          \n"
             "{0}             Run on    {1}\n"
             "{0}  *****************************************************************\n\n",
             prefix, Date());
}

// Issues: - When mapping element numbers, blocks are irrelevant.  Problem is
//           the variables that are determined to be stored in each file are
//           NOT independent of blocks .. in fact, that is how it determines
//           if the two files have the same element variable stored.  The
//           mapping will still run ok, just don't expect it to work if the
//           blocks don't line up and different variables are stored in
//           different blocks.

template <typename INT>
extern void Build_Variable_Names(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, bool *diff_found);

template <typename INT> extern bool Check_Global(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2);

template <typename INT>
extern void Check_Compatible_Meshes(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, bool check_only,
                                    const std::vector<INT> &node_map,
                                    const std::vector<INT> &elmt_map, const INT *node_id_map);

template <typename INT>
int Create_File(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, const std::string &diffile_name,
                bool *diff_found);

double To_Double(const std::string &str_val);

double FileDiff(double v1, double v2, ToleranceMode type);

void Die_TS(double ts);

template <typename INT> size_t global_elmt_num(ExoII_Read<INT> &file, size_t b_idx, size_t e_idx);

template <typename INT> double Find_Min_Coord_Sep(ExoII_Read<INT> &file);

int timeStepIsExcluded(int ts);

template <typename INT>
const double *get_nodal_values(ExoII_Read<INT> &filen, int time_step, size_t idx, size_t fno,
                               const std::string &name, bool *diff_flag);
template <typename INT>
const double *get_nodal_values(ExoII_Read<INT> &filen, const TimeInterp &t, size_t idx, size_t fno,
                               const std::string &name, bool *diff_flag);

template <typename INT>
void do_summaries(ExoII_Read<INT> &file, int time_step, std::vector<MinMaxData> &mm_glob,
                  std::vector<MinMaxData> &mm_node, std::vector<MinMaxData> &mm_elmt,
                  std::vector<MinMaxData> &mm_ns, std::vector<MinMaxData> &mm_ss,
                  std::vector<MinMaxData> &mm_eb, std::vector<MinMaxData> &mm_fb,
                  const std::vector<INT> &elmt_map, bool *diff_flag);

template <typename INT>
void do_diffs(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, int time_step1, const TimeInterp &t2,
              int out_file_id, int output_step, const std::vector<INT> &node_map,
              const INT *node_id_map, const std::vector<INT> &elmt_map, const INT *elem_id_map,
              Exo_Block<INT> **blocks2, std::vector<double> &var_vals, bool *diff_flag);

template <typename INT>
bool summarize_globals(ExoII_Read<INT> &file, int step, std::vector<MinMaxData> &mm_glob);
template <typename INT>
bool summarize_nodals(ExoII_Read<INT> &file, int step, std::vector<MinMaxData> &mm_node);
template <typename INT>
bool summarize_element(ExoII_Read<INT> &file, int step, const std::vector<INT> &elmt_map,
                       std::vector<MinMaxData> &mm_elmt);
template <typename INT>
bool summarize_nodeset(ExoII_Read<INT> &file, int step, std::vector<MinMaxData> &mm_ns);
template <typename INT>
bool summarize_sideset(ExoII_Read<INT> &file, int step, std::vector<MinMaxData> &mm_ss);
template <typename INT>
bool summarize_edgeblock(ExoII_Read<INT> &file, int step, std::vector<MinMaxData> &mm_eb);
template <typename INT>
bool summarize_faceblock(ExoII_Read<INT> &file, int step, std::vector<MinMaxData> &mm_fb);

template <typename INT>
bool diff_globals(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, int step1, const TimeInterp &t2,
                  int out_file_id, int output_step, std::vector<double> &gvals);
template <typename INT>
bool diff_nodals(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, int step1, const TimeInterp &t2,
                 int out_file_id, int output_step, const std::vector<INT> &node_map,
                 const INT *id_map, std::vector<double> &nvals);
template <typename INT>
bool diff_element(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, int step1, const TimeInterp &t2,
                  int out_file_id, int output_step, const std::vector<INT> &elmt_map,
                  const INT *id_map, Exo_Block<INT> **blocks2, std::vector<double> &evals);

template <typename INT>
bool diff_element_attributes(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2,
                             const std::vector<INT> &elmt_map, const INT *id_map,
                             Exo_Block<INT> **blocks2);

template <typename INT>
bool diff_nodeset(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, int step1, const TimeInterp &t2,
                  int out_file_id, int output_step, const INT *id_map, std::vector<double> &vals);

template <typename INT>
bool diff_sideset(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, int step1, const TimeInterp &t2,
                  int out_file_id, int output_step, const INT *id_map, std::vector<double> &vals);

template <typename INT>
bool diff_sideset_df(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, const INT *id_map);

template <typename INT>
void output_summary(ExoII_Read<INT> &file1, MinMaxData &mm_time, std::vector<MinMaxData> &mm_glob,
                    std::vector<MinMaxData> &mm_node, std::vector<MinMaxData> &mm_elmt,
                    std::vector<MinMaxData> &mm_ns, std::vector<MinMaxData> &mm_ss,
                    std::vector<MinMaxData> &mm_eb, std::vector<MinMaxData> &mm_fb,
                    const INT *node_id_map, const INT *elem_id_map);

#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(_MSC_VER) ||                \
    defined(__MINGW32__) || defined(_WIN64) || defined(__MINGW64__)
#define __ED_WINDOWS__ 1
#endif

#if !defined(__ED_WINDOWS__)
#include <csignal>
// bit of a hack to get GNU's functions to enable floating point error trapping
#ifdef LINUX
#ifdef __USE_GNU
#include <fenv.h>
#else
#define __USE_GNU
#include <fenv.h>
#undef __USE_GNU
#endif
#endif

struct sigaction sigact; // the signal handler & blocked signals

#endif

bool checking_invalid = false;
bool invalid_data     = false;
extern "C" {
void floating_point_exception_handler(int signo)
{
  if (!checking_invalid) {
    Error(fmt::format("caught floating point exception ({}) bad data?\n", signo));
  }
  else {
    invalid_data = true;
  }
}
}

namespace {
  int get_int_size(const std::string &file_name)
  {
    if (file_name.empty()) {
      return 0;
    }

    int   ws      = 0;
    int   comp_ws = 8;
    float dumb    = 0.0;
    int   exoid   = ex_open(file_name.c_str(), EX_READ, &comp_ws, &ws, &dumb);
    if (exoid < 0) {
      Error(fmt::format("Couldn't open file \"{}\".\n", file_name));
    }
    int size = (ex_int64_status(exoid) & EX_ALL_INT64_DB) != 0 ? 8 : 4;
    ex_close(exoid);
    return size;
  }

  template <typename INT> TimeInterp get_surrounding_times(double time, ExoII_Read<INT> &file)
  {
    TimeInterp tprop;
    tprop.time = time;

    int num_times = file.Num_Times();
    if (num_times == 0 || time < file.Time(1)) {
      tprop.step2 = 0;
      return tprop;
    }

    if (time > file.Time(num_times)) {
      tprop.step1 = 0;
      return tprop;
    }

    int tbef = 1;
    for (int i = 2; i <= num_times; i++) {
      if (file.Time(i) <= time) {
        tbef = i;
      }
      else if (interFace.time_tol.type != ToleranceMode::IGNORE_ &&
               !interFace.time_tol.Diff(time, file.Time(i))) {
        tbef = i;
      }
      else {
        break;
      }
    }

    if (!interFace.time_tol.Diff(time, file.Time(tbef))) {
      tprop.step1 = tprop.step2 = tbef;
      return tprop;
    }

    SMART_ASSERT(tbef + 1 <= num_times)(tbef + 1)(num_times);
    tprop.step1 = tbef;
    tprop.step2 = tbef + 1;

    // Calculate proprtion...
    double t1        = file.Time(tbef);
    double t2        = file.Time(tbef + 1);
    tprop.proportion = (time - t1) / (t2 - t1);
    return tprop;
  }

  template <typename INT> void output_init(ExoII_Read<INT> &file, int count, const char *prefix)
  {
    FileInfo fi(file.File_Name());
    fmt::print(
        "{0}  FILE {19}: {1}\n"
        "{0}   Title: {2}\n"
        "{0}          Dim = {3}, Nodes = {5}, Elements = {6}, Faces = {20}, Edges = {21}\n"
        "{0}          Element Blocks = {4}, Face Blocks = {10}, Edge Blocks = {9}, Nodesets = {7}, "
        "Sidesets = {8}, Assemblies = {22}\n"
        "{0}    Vars: Global = {11}, Nodal = {12}, Element = {13}, Face = {17}, Edge = {18}, "
        "Nodeset = {14}, Sideset = {15}, Times = {16}\n\n",
        prefix, fi.realpath(), file.Title(), file.Dimension(), file.Num_Element_Blocks(),
        file.Num_Nodes(), file.Num_Elements(), file.Num_Node_Sets(), file.Num_Side_Sets(),
        file.Num_Edge_Blocks(), file.Num_Face_Blocks(), file.Num_Global_Vars(),
        file.Num_Nodal_Vars(), file.Num_Element_Vars(), file.Num_NS_Vars(), file.Num_SS_Vars(),
        file.Num_Times(), file.Num_FB_Vars(), file.Num_EB_Vars(), count, file.Num_Faces(),
        file.Num_Edges(), file.Num_Assembly());
  }

  void initialize(std::vector<MinMaxData> &mm_entity, size_t size, const ToleranceType &ttype)
  {
    mm_entity.resize(size);
    for (auto &mm : mm_entity) {
      mm.type = ttype;
    }
  }

  template <typename INT> bool exodiff(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2);
} // namespace

int main(int argc, char *argv[])
{
  bool ok = interFace.parse_options(argc, argv);

  if (!ok) {
    exit(1);
  }

  checking_invalid = false;
  invalid_data     = false;

#if !defined(__ED_WINDOWS__)
  sigfillset(&(sigact.sa_mask));
  sigact.sa_handler = floating_point_exception_handler;
  if (sigaction(SIGFPE, &sigact, nullptr) == -1) {
    perror("sigaction failed");
  }
#endif

#if defined(LINUX) && defined(GNU) && !defined(__ED_WINDOWS__)
  // for GNU, this seems to be needed to turn on trapping
  feenableexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID);
#endif

  std::string file1_name = interFace.file1;
  std::string file2_name = interFace.file2;

  if (interFace.summary_flag && file1_name.empty()) {
    Error(fmt::format("Summary option specified but an exodus "
                      "file was not specified.\n"));
  }

  if (interFace.summary_flag) {
    file2_name                     = "";
    interFace.glob_var_do_all_flag = true;
    interFace.node_var_do_all_flag = true;
    interFace.elmt_var_do_all_flag = true;
    interFace.elmt_att_do_all_flag = true;
    interFace.ns_var_do_all_flag   = true;
    interFace.ss_var_do_all_flag   = true;
    interFace.eb_var_do_all_flag   = true;
    interFace.fb_var_do_all_flag   = true;
    interFace.map_flag             = MapType::FILE_ORDER;
    interFace.quiet_flag           = false;
    Print_Banner("#");
  }

  if (!interFace.quiet_flag && !interFace.summary_flag) {
    Print_Banner(" ");
  }

  // Check integer sizes in input file(s)...
  int int_size = 4;
  if (interFace.ints_64_bits) {
    int_size = 8;
  }
  else if (get_int_size(file1_name) == 8) {
    int_size = 8;
  }
  else if (!interFace.summary_flag && get_int_size(file2_name) == 8) {
    int_size = 8;
  }

  bool diff_flag = true;
  if (int_size == 4) {
    // Open input files.
    ExoII_Read<int> file1(file1_name);
    file1.modify_time_values(interFace.time_value_scale, interFace.time_value_offset);

    ExoII_Read<int> file2(file2_name);
    diff_flag = exodiff(file1, file2);
  }
  else {
    // Open input files.
    ExoII_Read<int64_t> file1(file1_name);
    ExoII_Read<int64_t> file2(file2_name);
    diff_flag = exodiff(file1, file2);
  }
#if 0
    add_to_log(argv[0], 0);
#else
  // Temporarily differentiate this version from previous version in logs.
  std::string code = "exodiff-" + version;
  add_to_log(code.c_str(), 0);
#endif

  if (interFace.exit_status_switch && diff_flag) {
    return 2;
  }

  return 0;
}

namespace {
  template <typename INT> bool exodiff(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2)
  {
    if (!interFace.quiet_flag && !interFace.summary_flag) {
      fmt::print("Reading first file ... \n");
    }
    std::string serr = file1.Open_File();
    if (!serr.empty()) {
      Error(fmt::format("{}\n", serr));
    }
    if (!interFace.summary_flag) {
      if (!interFace.quiet_flag) {
        fmt::print("Reading second file ... \n");
      }
      serr = file2.Open_File();
      if (!serr.empty()) {
        Error(fmt::format("{}\n", serr));
      }
    }

    if (interFace.summary_flag) {
      output_init(file1, 1, "#");
    }
    else {
      if (!interFace.quiet_flag) {
        output_init(file1, 1, "");
        output_init(file2, 2, "");
        if (interFace.pedantic) {
          fmt::print("  Pedantic Checking Enabled\n");
        }
        if (!interFace.command_file.empty()) {
          FileInfo fi(interFace.command_file);
          fmt::print("  COMMAND FILE: {}\n\n", fi.realpath());
        }
      }
    }

    if (!interFace.summary_flag) {
      bool is_same = Check_Global(file1, file2);
      if (!is_same) {
        file1.Close_File();
        file2.Close_File();
        DIFF_OUT("\nexodiff: Files are different\n");
        return interFace.exit_status_switch;
      }
    }

    // When mapping is on ("-m"), node_map maps indexes from file1 to indexes
    // into file2.  Similarly with elmt_map.
    std::vector<INT> node_map;
    std::vector<INT> elmt_map;
    if (interFace.map_flag == MapType::DISTANCE) {
      Compute_Maps(node_map, elmt_map, file1, file2);
    }
    else if (interFace.map_flag == MapType::PARTIAL) {
      // Same as distance, but ok if not all nodes/elements are matched
      Compute_Partial_Maps(node_map, elmt_map, file1, file2);
    }
    else if (interFace.map_flag == MapType::USE_FILE_IDS) {
      if (!interFace.ignore_maps) {
        // Node/element X in file 1 matches node/element X in file 2 no matter what order they are
        // in
        Compute_FileId_Maps(node_map, elmt_map, file1, file2);
      }
      else {
        size_t num_nodes = file1.Num_Nodes();
        node_map.resize(num_nodes);
        std::iota(node_map.begin(), node_map.end(), 0);

        size_t num_elem = file1.Num_Elements();
        elmt_map.resize(num_elem);
        std::iota(elmt_map.begin(), elmt_map.end(), 0);
      }
    }
    else if (interFace.map_flag == MapType::FILE_ORDER) {
      // Match by implicit ordering... IDs in that ordering must match (checked later)
      size_t num_nodes = file1.Num_Nodes();
      node_map.resize(num_nodes);
      std::iota(node_map.begin(), node_map.end(), 0);

      size_t num_elem = file1.Num_Elements();
      elmt_map.resize(num_elem);
      std::iota(elmt_map.begin(), elmt_map.end(), 0);
    }
    else {
      Error("Invalid map option.\n");
    }

    if (interFace.dump_mapping) {
      Dump_Maps(node_map, elmt_map, file1);
    }

    bool diff_flag = false; // Set to 'true' to indicate files contain diffs
    // Call this before checking for compatible meshes since it sets which variables
    // are going to be compared.  If no variables of a specific type, then not an error
    // if the meshes are different in that type.
    Build_Variable_Names(file1, file2, &diff_flag);

    // Get node and element number maps which map internal implicit ids into
    // global ids...
    const INT *node_id_map = nullptr;
    const INT *elem_id_map = nullptr;
    if (!interFace.ignore_maps) {
      file1.Load_Node_Map();
      file1.Load_Element_Map();
      node_id_map = file1.Get_Node_Map();
      elem_id_map = file1.Get_Element_Map();
      if (!interFace.summary_flag) {
        bool diff =
            Compare_Maps(file1, file2, node_map, elmt_map, interFace.map_flag == MapType::PARTIAL);
        if (diff && (interFace.map_flag == MapType::FILE_ORDER)) {
          fmt::print(stderr,
                     "exodiff: Exiting due to node/element mismatch with `-match_file_order` "
                     "option enabled.\n");
          if (interFace.exit_status_switch) {
            exit(2);
          }
          else {
            exit(1);
          }
        }
      }
    }
    else {
      // Ignoring the maps from the file, so create a dummy
      // map to make logic later in program consistent.
      node_id_map  = new INT[file1.Num_Nodes()];
      INT *tmp_map = const_cast<INT *>(node_id_map);
      std::iota(tmp_map, tmp_map + file1.Num_Nodes(), 1);

      elem_id_map = new INT[file1.Num_Elements()];
      tmp_map     = const_cast<INT *>(elem_id_map);
      std::iota(tmp_map, tmp_map + file1.Num_Elements(), 1);
    }

    int out_file_id = -1;
    if (!interFace.summary_flag) {
      std::string diffile_name = interFace.diff_file;
      Check_Compatible_Meshes(file1, file2, (diffile_name.empty()), node_map, elmt_map,
                              node_id_map);
      // Doesn't return if meshes are not compatible...

      out_file_id = Create_File(file1, file2, diffile_name, &diff_flag);
    }

    SMART_ASSERT(!(interFace.summary_flag && out_file_id >= 0));

    if (!interFace.quiet_flag || interFace.summary_flag) {
      fmt::print("\n{0} ==============================================================\n"
                 "{0}  NOTE: All node and element ids are reported as {1} ids.\n\n",
                 interFace.summary_flag ? "#" : " ", interFace.ignore_maps ? "local" : "global");
      if (interFace.interpolating) {
        fmt::print("{}  NOTE: Interpolation mode is enabled.\n\n",
                   interFace.summary_flag ? "#" : " ");
      }
    }

    std::vector<double> var_vals;
    if (out_file_id >= 0) {
      size_t max_ent = interFace.glob_var_names.size();
      if (file1.Num_Nodes() > max_ent) {
        max_ent = file1.Num_Nodes();
      }
      if (file1.Num_Elements() > max_ent) {
        max_ent = file1.Num_Elements();
      }
      if (file1.Num_Faces() > max_ent) {
        max_ent = file1.Num_Faces();
      }
      if (file1.Num_Edges() > max_ent) {
        max_ent = file1.Num_Edges();
      }

      var_vals.resize(max_ent);
    }

    // When mapping is in effect, it is efficient to grab pointers to all blocks.
    Exo_Block<INT> **blocks2 = nullptr;
    if (!elmt_map.empty()) {
      blocks2 = new Exo_Block<INT> *[file2.Num_Element_Blocks()];
      for (size_t b = 0; b < file2.Num_Element_Blocks(); ++b) {
        blocks2[b] = file2.Get_Element_Block_by_Index(b);
      }
    }

    // Diff attributes...
    if (!interFace.ignore_attributes && elmt_map.empty() && !interFace.summary_flag) {
      if (diff_element_attributes(file1, file2, elmt_map, elem_id_map, blocks2)) {
        diff_flag = true;
      }
    }

    // Diff sideset distribution factors...
    if (!interFace.ignore_sideset_df && !interFace.summary_flag) {
      if (diff_sideset_df(file1, file2, elem_id_map)) {
        diff_flag = true;
      }
    }

    int min_num_times = file1.Num_Times();
    int output_step   = 1;

    MinMaxData mm_time;
    mm_time.type = ToleranceType::mm_time;
    std::vector<MinMaxData> mm_glob;
    std::vector<MinMaxData> mm_node;
    std::vector<MinMaxData> mm_elmt;
    std::vector<MinMaxData> mm_ns;
    std::vector<MinMaxData> mm_ss;
    std::vector<MinMaxData> mm_eb;
    std::vector<MinMaxData> mm_fb;

    if (interFace.summary_flag) {
      initialize(mm_glob, interFace.glob_var_names.size(), ToleranceType::mm_global);
      initialize(mm_node, interFace.node_var_names.size(), ToleranceType::mm_nodal);
      initialize(mm_elmt, interFace.elmt_var_names.size(), ToleranceType::mm_element);
      initialize(mm_ns, interFace.ns_var_names.size(), ToleranceType::mm_nodeset);
      initialize(mm_ss, interFace.ss_var_names.size(), ToleranceType::mm_sideset);
      initialize(mm_eb, interFace.eb_var_names.size(), ToleranceType::mm_edgeblock);
      initialize(mm_fb, interFace.fb_var_names.size(), ToleranceType::mm_faceblock);
    }
    else {
      min_num_times =
          (file1.Num_Times() < file2.Num_Times() ? file1.Num_Times() : file2.Num_Times());

      if (interFace.interpolating) {
        min_num_times = file1.Num_Times();
      }

      if (interFace.time_step_stop > 0 && interFace.time_step_stop < min_num_times) {
        min_num_times = interFace.time_step_stop;
      }
    }

    // If explicit times are set, then only want to diff a single time at those
    // specified times....
    if (interFace.explicit_steps.first != 0 && interFace.explicit_steps.second != 0) {
      int ts1 = interFace.explicit_steps.first;
      if (ts1 == -1) {
        ts1 = file1.Num_Times();
      }
      int ts2 = interFace.explicit_steps.second;
      if (ts2 == -1) {
        ts2 = file2.Num_Times();
      }
      TimeInterp t2;
      t2.step1      = ts2;
      t2.step2      = ts2;
      t2.time       = file2.Time(ts2);
      t2.proportion = 0.0;

      if (!interFace.quiet_flag) {
        if (out_file_id >= 0) {
          fmt::print("Processing explicit time steps. File 1 step = {}  File 2 step = {}\n", ts1,
                     ts2);
        }
        else {
          std::string buf =
              fmt::format("  --------- Explicit Time step File 1: {}, {:13.7e} ~ File 2: {}, "
                          "{:13.7e} ---------",
                          ts1, file1.Time(ts1), ts2, t2.time);
          DIFF_OUT(buf, fmt::color::green);
        }
      }

      if (interFace.summary_flag) {
        do_summaries(file1, ts1, mm_glob, mm_node, mm_elmt, mm_ns, mm_ss, mm_eb, mm_fb, elmt_map,
                     &diff_flag);
      }
      else {
        do_diffs(file1, file2, ts1, t2, out_file_id, output_step, node_map, node_id_map, elmt_map,
                 elem_id_map, blocks2, var_vals, &diff_flag);
      }
    }
    else {

      // If time_step_offset == -1, then determine the offset automatically.
      // Assumes file1 has more steps than file2 and that the last step(s)
      // on file2 match the last step(s) on file1.
      if (interFace.time_step_offset == -1) {
        interFace.time_step_offset = file1.Num_Times() - file2.Num_Times();
        if (interFace.time_step_offset < 0) {
          Error("Second database must have less timesteps than "
                "first database.\n");
        }
      }

      // If time_step_offset == -2, then determine the offset automatically.
      // Find the closest time on file1 to the first time on file2.
      // Assumes file1 has more steps than file2.
      if (interFace.time_step_offset == -2) {
        if (file1.Num_Times() < file2.Num_Times()) {
          Error("Second database must have less timesteps than "
                "first database.\n");
        }

        double t2      = file2.Time(1);
        double mindiff = fabs(t2 - file1.Time(1));
        int    step    = 1;
        for (int i = 2; i < file1.Num_Times(); i++) {
          double t1   = file1.Time(i);
          double diff = fabs(t2 - t1);
          if (diff < mindiff) {
            step    = i;
            mindiff = diff;
          }
        }
        interFace.time_step_offset = step - 1;
      }

      if (interFace.time_step_offset > 0) {
        if (interFace.time_step_start > 0) {
          fmt::print(
              "The first {} timesteps in the first database will be skipped because of time step "
              "offset and time step start settings.\n\n",
              interFace.time_step_offset + interFace.time_step_start - 1);
        }
        else {
          fmt::print(
              "The first {} timesteps in the first database will be skipped because of time step "
              "offset setting.\n\n",
              interFace.time_step_offset);
        }
      }

      if (interFace.time_step_start == -1) {
        // Want to compare the last timestep on both databases...
        int time_step1             = file1.Num_Times();
        int time_step2             = file2.Num_Times();
        interFace.time_step_start  = time_step2;
        interFace.time_step_offset = time_step1 - time_step2;
        min_num_times              = interFace.time_step_start;
        fmt::print("Comparing only the final step (step {} on first, step {}"
                   " on second) on each database.\n\n",
                   time_step1, time_step2);
      }
      else if (interFace.time_step_start < 0) {
        interFace.time_step_start = min_num_times;
      }
      else if (interFace.time_step_start < 1) {
        interFace.time_step_start = 1;
      }

      if (interFace.time_step_start > min_num_times && min_num_times > 0) {
        Warning("Time step options resulted in no timesteps being compared.\n");
        diff_flag = true;
      }

      for (int time_step = interFace.time_step_start; time_step <= min_num_times;
           time_step += interFace.time_step_increment) {
        if (timeStepIsExcluded(time_step) || interFace.ignore_steps) {
          continue;
        }

        int time_step1 = time_step + interFace.time_step_offset;
        int time_step2 = time_step;
        SMART_ASSERT(time_step1 <= file1.Num_Times());

        TimeInterp t2;
        if (!interFace.summary_flag) {
          t2 = get_surrounding_times(file1.Time(time_step1), file2);
          if (!interFace.interpolating) {
            t2.step1      = time_step2;
            t2.step2      = time_step2;
            t2.time       = file2.Time(time_step2);
            t2.proportion = 0.0;
          }
          SMART_ASSERT(t2.step1 <= file2.Num_Times());
          SMART_ASSERT(t2.step2 <= file2.Num_Times());
        }

        if (interFace.summary_flag) {
          double t = file1.Time(time_step1);
          mm_time.spec_min_max(t, time_step1);
        }
        else if (out_file_id >= 0 && !interFace.quiet_flag) {
          fmt::print("Processing time step {}  (Difference in time values = {})\n", time_step1,
                     (file1.Time(time_step1) - file2.Time(time_step2)));
        }
        else if (out_file_id < 0) {
          if (!interFace.quiet_flag) {
            std::string buf;
            if (interFace.interpolating) {
              if (t2.step1 == -1) {
                buf = fmt::format(
                    "  --------- Time step {}, {:13.7e} ~ Skipping - Before all times on "
                    "file2 (INTERPOLATING)",
                    time_step1, file1.Time(time_step1));
              }
              else if (t2.step2 == -1) {
                buf = fmt::format(
                    "  --------- Time step {}, {:13.7e} ~ Skipping - After all times on "
                    "file2 (INTERPOLATING)",
                    time_step1, file1.Time(time_step1));
              }
              else if (t2.step1 == t2.step2) {
                buf = fmt::format(
                    "  --------- Time step {}, {:13.7e} ~ Matches step {}, {:13.7e} on file2 "
                    "{} diff: {:12.5e}",
                    time_step1, file1.Time(time_step1), t2.step1, file2.Time(t2.step1),
                    interFace.time_tol.abrstr(),
                    FileDiff(file1.Time(time_step1), file2.Time(t2.step1),
                             interFace.time_tol.type));
              }
              else {
                buf = fmt::format(
                    "  --------- Time step {}, {:13.7e} ~ Interpolating step {}, {:13.7e} and "
                    "step {}, {:13.7e}, proportion {:10.4e} on file2",
                    time_step1, file1.Time(time_step1), t2.step1, file2.Time(t2.step1), t2.step2,
                    file2.Time(t2.step2), t2.proportion);
              }
            }
            else {
              buf = fmt::format("  --------- Time step {}, {:13.7e} ~ {:13.7e}, {} diff: {:12.5e}",
                                time_step1, file1.Time(time_step1), file2.Time(time_step2),
                                interFace.time_tol.abrstr(),
                                FileDiff(file1.Time(time_step1), file2.Time(time_step2),
                                         interFace.time_tol.type));
            }
            fmt::print("{}", buf);
          }

          if (!interFace.interpolating &&
              interFace.time_tol.Diff(file1.Time(time_step1), file2.Time(time_step2))) {
            diff_flag = true;
            if (interFace.quiet_flag) {
              Die_TS(time_step1);
            }
            else {
              DIFF_OUT(" (FAILED) \n");
            }
          }
          else if (!interFace.quiet_flag) {
            fmt::print(" ---------\n");
          }
          if (interFace.interpolating && time_step == min_num_times) {
            // last time.  Check if final database times match within specified tolerance...
            int final2 = file2.Num_Times();
            if (interFace.final_time_tol.Diff(file1.Time(time_step1), file2.Time(final2))) {
              diff_flag = true;
              std::ostringstream diff;
              fmt::print(diff,
                         "\tFinal database times differ by {}  which is not within specified {}"
                         " tolerance of {} (FAILED)",
                         FileDiff(file1.Time(time_step1), file2.Time(final2),
                                  interFace.final_time_tol.type),
                         interFace.final_time_tol.typestr(), interFace.final_time_tol.value);
              DIFF_OUT(diff);
            }
          }
        }

        if (out_file_id >= 0) {
          double t = file1.Time(time_step1);
          ex_put_time(out_file_id, output_step, &t);
        }

        if (interFace.interpolating && (t2.step1 == -1 || t2.step2 == -1)) {
          continue;
        }

        if (interFace.summary_flag) {
          do_summaries(file1, time_step1, mm_glob, mm_node, mm_elmt, mm_ns, mm_ss, mm_eb, mm_fb,
                       elmt_map, &diff_flag);
        }
        else {
          do_diffs(file1, file2, time_step1, t2, out_file_id, output_step, node_map, node_id_map,
                   elmt_map, elem_id_map, blocks2, var_vals, &diff_flag);
        }

        output_step++;
      } // End of time step loop.

      // Make sure there is an operation to perform (compare times, variables, ...)
      if (!interFace.ignore_steps) {
        if ((min_num_times == 0 && interFace.coord_tol.type == ToleranceMode::IGNORE_) ||
            (min_num_times > 0 && interFace.time_tol.type == ToleranceMode::IGNORE_ &&
             interFace.glob_var_names.empty() && interFace.node_var_names.empty() &&
             interFace.elmt_var_names.empty() && interFace.elmt_att_names.empty() &&
             interFace.ns_var_names.empty() && interFace.ss_var_names.empty() &&
             interFace.eb_var_names.empty() && interFace.fb_var_names.empty())) {
          DIFF_OUT("\nWARNING: No comparisons were performed during this execution.");
          diff_flag = true;
        }
      }
    }

    if (interFace.summary_flag) {
      output_summary(file1, mm_time, mm_glob, mm_node, mm_elmt, mm_ns, mm_ss, mm_eb, mm_fb,
                     node_id_map, elem_id_map);
    }
    else if (out_file_id >= 0) {
      ex_close(out_file_id);
    }
    else if (diff_flag) {
      DIFF_OUT("\nexodiff: Files are different\n");
    }
    else if (interFace.ignore_steps && (file1.Num_Times() != 0 || file2.Num_Times() != 0)) {
      DIFF_OUT("\nexodiff: Files are the same, but all transient data was ignored due to "
               "-ignore_steps option",
               fmt::color::green);
    }
    else if (file1.Num_Times() != file2.Num_Times()) {
      if ((file1.Num_Times() - interFace.time_step_offset == file2.Num_Times()) ||
          (interFace.time_step_stop > 0) ||
          (interFace.explicit_steps.first != 0 && interFace.explicit_steps.second != 0) ||
          (interFace.interpolating)) {
        std::ostringstream diff;
        fmt::print(diff, "\nexodiff: Files are the same\n"
                         "         The number of timesteps are different but "
                         "the timesteps that were compared are the same.\n");
        DIFF_OUT(diff);
      }
      else {
        DIFF_OUT("\nexodiff: Files are different (# time steps differ)");
        diff_flag = true;
      }
    }
    else if (interFace.map_flag == MapType::PARTIAL) {
      DIFF_OUT("\nexodiff: Files are the same (partial match selected)\n", fmt::color::green);
    }
    else {
      DIFF_OUT("\nexodiff: Files are the same\n", fmt::color::green);
    }

    if (!interFace.ignore_maps) {
      file1.Free_Node_Map();
      file1.Free_Element_Map();
    }
    else {
      delete[] node_id_map;
      delete[] elem_id_map;
    }

    delete[] blocks2;

    file1.Close_File();
    if (!interFace.summary_flag) {
      file2.Close_File();
    }

    return diff_flag;
  }
} // namespace
double FileDiff(double v1, double v2, ToleranceMode type)
{
  if (type == ToleranceMode::IGNORE_) { // ignore
    return 0.0;
  }
  if (type == ToleranceMode::RELATIVE_) { // relative diff
    if (v1 == 0.0 && v2 == 0.0) {
      return 0.0;
    }
    double max = fabs(v1) < fabs(v2) ? fabs(v2) : fabs(v1);
    return (v1 - v2) / max;
  }
  if (type == ToleranceMode::COMBINED_) {
    // if (Abs(x - y) <= Max(absTol, relTol * Max(Abs(x), Abs(y))))
    // In the current implementation, absTol == relTol;
    // In summary, use abs tolerance if both values are less than 1.0;
    // else use relative tolerance.

    double max = fabs(v1) < fabs(v2) ? fabs(v2) : fabs(v1);
    double tol = 1.0 < max ? max : 1.0;
    return fabs(v1 - v2) / tol;
  }
  if (type == ToleranceMode::ABSOLUTE_) {
    return (v1 - v2);
  }
  if (type == ToleranceMode::EIGEN_REL_) { // relative diff
    if (v1 == 0.0 && v2 == 0.0) {
      return 0.0;
    }
    double max = fabs(v1) < fabs(v2) ? fabs(v2) : fabs(v1);
    return (fabs(v1) - fabs(v2)) / max;
  }
  if (type == ToleranceMode::EIGEN_COM_) {
    // if (Abs(x - y) <= Max(absTol, relTol * Max(Abs(x), Abs(y))))
    // In the current implementation, absTol == relTol;
    // In summary, use abs tolerance if both values are less than 1.0;
    // else use relative tolerance.

    double max = fabs(v1) < fabs(v2) ? fabs(v2) : fabs(v1);
    double tol = 1.0 < max ? max : 1.0;
    return fabs(fabs(v1) - fabs(v2)) / tol;
  }
  if (type == ToleranceMode::EIGEN_ABS_) {
    return (fabs(v1) - fabs(v2));
  }
  return 0.0;
}

void Die_TS(double ts)
{
  std::ostringstream diff;
  fmt::print(diff, "exodiff: Files are different (time step {})", ts);
  DIFF_OUT(diff);
  if (interFace.exit_status_switch) {
    exit(2);
  }
  else {
    exit(1);
  }
}

template <typename INT> size_t global_elmt_num(ExoII_Read<INT> &file, size_t b_idx, size_t e_idx)
{
  SMART_ASSERT(b_idx < file.Num_Element_Blocks());

  size_t g = 0;
  for (size_t b = 0; b < file.Num_Element_Blocks(); ++b) {
    if (b_idx == b) {
      return g + e_idx + 1;
    }

    SMART_ASSERT(file.Get_Element_Block_by_Index(b) != 0);
    g += file.Get_Element_Block_by_Index(b)->Size();
  }
  SMART_ASSERT(0);
  return 0;
}

bool Invalid_Values(const double *values, size_t count)
{
  bool valid = true;
  if (!interFace.ignore_nans) {
    checking_invalid = true;
    invalid_data     = false;

    SMART_ASSERT(values != nullptr);

    for (size_t i = 0; i < count; i++) {
#if defined(interix)
      if (values[i] != values[i])
#else
      if (std::isnan(values[i]))
#endif
      {
        valid = false;
        break;
      }
      if (invalid_data) { // may get set by SIGFPE handler
        valid = false;
        break;
      }
    }

    checking_invalid = false;
    invalid_data     = false;
  }
  return !valid;
}

bool Equal_Values(const double *values, size_t count, double *value)
{
  SMART_ASSERT(values != nullptr);
  *value = values[0];
  return (std::adjacent_find(values, values + count, std::not_equal_to<>()) == values + count);
}

template <typename INT>
const double *get_nodal_values(ExoII_Read<INT> &filen, int time_step, size_t idx, int fno,
                               const std::string &name, bool *diff_flag)
{
  const double *vals = nullptr;
  if (fno == 1 || !interFace.summary_flag) {
    filen.Load_Nodal_Results(time_step, idx);
    vals = filen.Get_Nodal_Results(idx);

    if (vals != nullptr) {
      if (Invalid_Values(vals, filen.Num_Nodes())) {
        Warning(fmt::format("NaN found for nodal variable '{}' in file {}\n", name, fno));
        *diff_flag = true;
      }
    }
  }
  return vals;
}

template <typename INT>
const double *get_nodal_values(ExoII_Read<INT> &filen, const TimeInterp &t, size_t idx, int fno,
                               const std::string &name, bool *diff_flag)
{
  const double *vals = nullptr;
  if (fno == 1 || !interFace.summary_flag) {
    vals = filen.Get_Nodal_Results(t.step1, t.step2, t.proportion, idx);

    if (vals != nullptr) {
      if (Invalid_Values(vals, filen.Num_Nodes())) {
        Warning(fmt::format("NaN found for nodal variable '{}' in file {}\n", name, fno));
        *diff_flag = true;
      }
    }
  }
  return vals;
}

template <typename INT>
bool summarize_globals(ExoII_Read<INT> &file, int step, std::vector<MinMaxData> &mm_glob)
{
  bool diff_flag = false;
  if (interFace.glob_var_names.empty()) {
    return diff_flag;
  }

  // Global variables.
  file.Load_Global_Results(step);
  const double *vals = file.Get_Global_Results();
  if (vals == nullptr) {
    Error("Could not find global variables on file 1.\n");
  }

  for (unsigned out_idx = 0; out_idx < interFace.glob_var_names.size(); ++out_idx) {
    const std::string &name = (interFace.glob_var_names)[out_idx];
    int                idx = find_string(file.Global_Var_Names(), name, interFace.nocase_var_names);
    if (idx < 0) {
      Error(fmt::format("Unable to find global variable named '{}' on database.\n", name));
    }
    mm_glob[out_idx].spec_min_max(vals[idx], step);
  }
  return diff_flag;
}

template <typename INT>
bool summarize_nodals(ExoII_Read<INT> &file, int step, std::vector<MinMaxData> &mm_node)
{
  bool diff_flag = false;
  for (unsigned n_idx = 0; n_idx < interFace.node_var_names.size(); ++n_idx) {
    const std::string &name = (interFace.node_var_names)[n_idx];
    int                idx  = find_string(file.Nodal_Var_Names(), name, interFace.nocase_var_names);
    if (idx < 0) {
      Error(fmt::format("Unable to find nodal variable named '{}' on database.\n", name));
    }
    const double *vals = get_nodal_values(file, step, idx, 1, name, &diff_flag);

    if (vals == nullptr) {
      Error("Could not find nodal variables on file 1\n");
    }

    size_t ncount = file.Num_Nodes();
    for (size_t n = 0; n < ncount; ++n) {
      mm_node[n_idx].spec_min_max(vals[n], step, n);
    }
    file.Free_Nodal_Results(idx);
  }
  file.Free_Nodal_Results();
  return diff_flag;
}

const double *get_validated_variable(Exo_Entity *entity, int step, int vidx,
                                     const std::string &name, bool *diff_flag)
{
  if (entity->Size() == 0) {
    return nullptr;
  }
  if (!entity->is_valid_var(vidx)) {
    return nullptr;
  }

  entity->Load_Results(step, vidx);
  const double *vals = entity->Get_Results(vidx);
  if (vals == nullptr) {
    Warning(fmt::format("Could not find variable '{}' in {} {}, file 1.\n", name,
                        entity->short_label(), entity->Id()));
    *diff_flag = true;
    return vals;
  }

  if (Invalid_Values(vals, entity->Size())) {
    Warning(fmt::format("NaN found for variable '{}' in {} {}, file 1\n", name,
                        entity->short_label(), entity->Id()));
    *diff_flag = true;
  }
  return vals;
}

const double *get_validated_variable(Exo_Entity *entity, const TimeInterp &t2, int vidx,
                                     const std::string &name, bool *diff_flag)
{
  if (entity == nullptr) {
    return nullptr;
  }
  if (entity->Size() == 0) {
    return nullptr;
  }
  if (!entity->is_valid_var(vidx)) {
    return nullptr;
  }

  entity->Load_Results(t2.step1, t2.step2, t2.proportion, vidx);
  const double *vals = entity->Get_Results(vidx);
  if (vals == nullptr) {
    Warning(fmt::format("Could not find variable '{}' in {} {}, file 2.\n", name,
                        entity->short_label(), entity->Id()));
    *diff_flag = true;
    return vals;
  }

  if (Invalid_Values(vals, entity->Size())) {
    Warning(fmt::format("NaN found for variable '{}' in {} {}, file 2.\n", name,
                        entity->short_label(), entity->Id()));
    *diff_flag = true;
  }
  return vals;
}

template <typename INT>
bool summarize_element(ExoII_Read<INT> &file, int step, const std::vector<INT> &elmt_map,
                       std::vector<MinMaxData> &mm_elmt)
{
  bool diff_flag = false;

  for (unsigned e_idx = 0; e_idx < interFace.elmt_var_names.size(); ++e_idx) {
    const std::string &name = (interFace.elmt_var_names)[e_idx];
    int vidx = find_string(file.Element_Var_Names(), name, interFace.nocase_var_names);
    if (vidx < 0) {
      Error(fmt::format("Unable to find element variable named '{}' on database.\n", name));
    }

    size_t global_elmt_index = 0;
    for (size_t b = 0; b < file.Num_Element_Blocks(); ++b) {
      Exo_Block<INT> *eblock = file.Get_Element_Block_by_Index(b);
      const double   *vals   = get_validated_variable(eblock, step, vidx, name, &diff_flag);
      if (vals == nullptr) {
        global_elmt_index += eblock->Size();
        continue;
      }

      size_t ecount = eblock->Size();
      for (size_t e = 0; e < ecount; ++e) {
        INT el_flag = 1;
        if (!elmt_map.empty()) {
          el_flag = elmt_map[global_elmt_index];
        }

        if (el_flag >= 0) {
          mm_elmt[e_idx].spec_min_max(vals[e], step, global_elmt_index, eblock->Id());
        }
        ++global_elmt_index;
      }

      eblock->Free_Results();
    }
  }
  return diff_flag;
}

template <typename INT>
bool summarize_nodeset(ExoII_Read<INT> &file, int step, std::vector<MinMaxData> &mm_ns)
{
  bool diff_flag = false;
  for (unsigned e_idx = 0; e_idx < interFace.ns_var_names.size(); ++e_idx) {
    const std::string &name = (interFace.ns_var_names)[e_idx];
    int                vidx = find_string(file.NS_Var_Names(), name, interFace.nocase_var_names);
    if (vidx < 0) {
      Error(fmt::format("Unable to find nodeset variable named '{}' on database.\n", name));
    }

    for (size_t b = 0; b < file.Num_Node_Sets(); ++b) {
      Node_Set<INT> *nset = file.Get_Node_Set_by_Index(b);

      const double *vals = get_validated_variable(nset, step, vidx, name, &diff_flag);
      if (vals == nullptr) {
        continue;
      }

      size_t ncount = nset->Size();
      for (size_t e = 0; e < ncount; ++e) {
        int idx = nset->Node_Index(e);
        mm_ns[e_idx].spec_min_max(vals[idx], step, e, nset->Id());
      }
      nset->Free_Results();
    }
  }
  return diff_flag;
}

template <typename INT>
bool summarize_sideset(ExoII_Read<INT> &file, int step, std::vector<MinMaxData> &mm_ss)
{
  bool diff_flag = false;
  for (unsigned e_idx = 0; e_idx < interFace.ss_var_names.size(); ++e_idx) {
    const std::string &name = (interFace.ss_var_names)[e_idx];
    int                vidx = find_string(file.SS_Var_Names(), name, interFace.nocase_var_names);
    if (vidx < 0) {
      Error(fmt::format("Unable to find sideset variable named '{}' on database.\n", name));
    }

    for (size_t b = 0; b < file.Num_Side_Sets(); ++b) {
      Side_Set<INT> *sset = file.Get_Side_Set_by_Index(b);

      const double *vals = get_validated_variable(sset, step, vidx, name, &diff_flag);
      if (vals == nullptr) {
        continue;
      }

      size_t ecount = sset->Size();
      for (size_t e = 0; e < ecount; ++e) {
        size_t ind = sset->Side_Index(e);
        mm_ss[e_idx].spec_min_max(vals[ind], step, e, sset->Id());
      }
      sset->Free_Results();
    }
  }
  return diff_flag;
}

template <typename INT>
bool summarize_edgeblock(ExoII_Read<INT> &file, int step, std::vector<MinMaxData> &mm_eb)
{
  bool diff_flag = false;
  for (unsigned e_idx = 0; e_idx < interFace.eb_var_names.size(); ++e_idx) {
    const std::string &name = (interFace.eb_var_names)[e_idx];
    int                vidx = find_string(file.EB_Var_Names(), name, interFace.nocase_var_names);
    if (vidx < 0) {
      Error(fmt::format("Unable to find edge block variable named '{}' on database.\n", name));
    }

    for (size_t b = 0; b < file.Num_Edge_Blocks(); ++b) {
      Edge_Block<INT> *eblock = file.Get_Edge_Block_by_Index(b);

      const double *vals = get_validated_variable(eblock, step, vidx, name, &diff_flag);
      if (vals == nullptr) {
        continue;
      }

      size_t ecount = eblock->Size();
      for (size_t e = 0; e < ecount; ++e) {
        size_t ind = eblock->Edge_Index(e);
        mm_eb[e_idx].spec_min_max(vals[ind], step, e, eblock->Id());
      }

      eblock->Free_Results();
    }
  }
  return diff_flag;
}

template <typename INT>
bool summarize_faceblock(ExoII_Read<INT> &file, int step, std::vector<MinMaxData> &mm_fb)
{
  bool diff_flag = false;
  for (unsigned f_idx = 0; f_idx < interFace.fb_var_names.size(); ++f_idx) {
    const std::string &name = (interFace.fb_var_names)[f_idx];
    int                vidx = find_string(file.FB_Var_Names(), name, interFace.nocase_var_names);
    if (vidx < 0) {
      Error(fmt::format("Unable to find face block variable named '{}' on database.\n", name));
    }

    for (size_t b = 0; b < file.Num_Face_Blocks(); ++b) {
      Face_Block<INT> *fblock = file.Get_Face_Block_by_Index(b);

      const double *vals = get_validated_variable(fblock, step, vidx, name, &diff_flag);
      if (vals == nullptr) {
        continue;
      }

      size_t fcount = fblock->Size();
      for (size_t f = 0; f < fcount; ++f) {
        size_t ind = fblock->Face_Index(f);
        mm_fb[f_idx].spec_min_max(vals[ind], step, f, fblock->Id());
      }

      fblock->Free_Results();
    }
  }
  return diff_flag;
}

template <typename INT>
void do_diffs(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, int time_step1, const TimeInterp &t2,
              int out_file_id, int output_step, const std::vector<INT> &node_map,
              const INT *node_id_map, const std::vector<INT> &elmt_map, const INT *elem_id_map,
              Exo_Block<INT> **blocks2, std::vector<double> &var_vals, bool *diff_flag)
{
  SMART_ASSERT(!interFace.summary_flag);
  if (diff_globals(file1, file2, time_step1, t2, out_file_id, output_step, var_vals)) {
    *diff_flag = true;
  }

  // Nodal variables.
  if (diff_nodals(file1, file2, time_step1, t2, out_file_id, output_step, node_map, node_id_map,
                  var_vals)) {
    *diff_flag = true;
  }

  // Element variables.
  if (diff_element(file1, file2, time_step1, t2, out_file_id, output_step, elmt_map, elem_id_map,
                   blocks2, var_vals)) {
    *diff_flag = true;
  }

  if (interFace.map_flag != MapType::PARTIAL) {
    // Nodeset variables.
    if (diff_nodeset(file1, file2, time_step1, t2, out_file_id, output_step, node_id_map,
                     var_vals)) {
      *diff_flag = true;
    }

    // Sideset variables.
    if (diff_sideset(file1, file2, time_step1, t2, out_file_id, output_step, elem_id_map,
                     var_vals)) {
      *diff_flag = true;
    }

    // Edge Block variables.
    if (diff_edgeblock(file1, file2, time_step1, t2, out_file_id, output_step, elem_id_map,
                       var_vals)) {
      *diff_flag = true;
    }

    // Face Block variables.
    if (diff_faceblock(file1, file2, time_step1, t2, out_file_id, output_step, elem_id_map,
                       var_vals)) {
      *diff_flag = true;
    }
  }
  else {
    if (!interFace.ns_var_names.empty() || !interFace.ss_var_names.empty() ||
        !interFace.eb_var_names.empty() || !interFace.fb_var_names.empty()) {
      fmt::print("WARNING: nodeset, sideset, edge block and face block variables not (yet) "
                 "compared for partial map\n");
    }
  }
}

template <typename INT>
void do_summaries(ExoII_Read<INT> &file, int time_step, std::vector<MinMaxData> &mm_glob,
                  std::vector<MinMaxData> &mm_node, std::vector<MinMaxData> &mm_elmt,
                  std::vector<MinMaxData> &mm_ns, std::vector<MinMaxData> &mm_ss,
                  std::vector<MinMaxData> &mm_eb, std::vector<MinMaxData> &mm_fb,
                  const std::vector<INT> &elmt_map, bool *diff_flag)
{
  SMART_ASSERT(interFace.summary_flag);
  if (summarize_globals(file, time_step, mm_glob)) {
    *diff_flag = true;
  }
  if (summarize_nodals(file, time_step, mm_node)) {
    *diff_flag = true;
  }
  if (summarize_element(file, time_step, elmt_map, mm_elmt)) {
    *diff_flag = true;
  }
  if (summarize_nodeset(file, time_step, mm_ns)) {
    *diff_flag = true;
  }
  if (summarize_sideset(file, time_step, mm_ss)) {
    *diff_flag = true;
  }
  if (summarize_edgeblock(file, time_step, mm_eb)) {
    *diff_flag = true;
  }
  if (summarize_faceblock(file, time_step, mm_fb)) {
    *diff_flag = true;
  }
}

void output_norms(Norm &norm, const std::string &name)
{
  if (interFace.doL1Norm && norm.diff(1) > 0.0) {
    std::string buf =
        fmt::format("   {:<{}} L1 norm of diff={:14.7e} ({:11.5e} ~ {:11.5e}) rel={:14.7e}", name,
                    name_length(), norm.diff(1), norm.left(1), norm.right(1), norm.relative(1));
    DIFF_OUT(buf, fmt::color::green);
  }
  if (interFace.doL2Norm && norm.diff(2) > 0.0) {
    std::string buf =
        fmt::format("   {:<{}} L2 norm of diff={:14.7e} ({:11.5e} ~ {:11.5e}) rel={:14.7e}", name,
                    name_length(), norm.diff(2), norm.left(2), norm.right(2), norm.relative(2));
    DIFF_OUT(buf, fmt::color::green);
  }
}

template <typename INT>
bool diff_globals(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, int step1, const TimeInterp &t2,
                  int out_file_id, int output_step, std::vector<double> &gvals)
{
  bool diff_flag = false;
  if (interFace.glob_var_names.empty()) {
    return diff_flag;
  }

  // Global variables.
  file1.Load_Global_Results(step1);
  const double *vals1 = file1.Get_Global_Results();
  if (vals1 == nullptr) {
    Error("Could not find global variables on file 1.\n");
  }

  file2.Load_Global_Results(t2.step1, t2.step2, t2.proportion);
  const double *vals2 = file2.Get_Global_Results();
  if (vals2 == nullptr) {
    Error("Could not find global variables on file 2.\n");
  }

  // ----------------------------------------------------------------------
  // Output file containing differences...
  if (out_file_id >= 0) {
    SMART_ASSERT(!gvals.empty());
    for (unsigned out_idx = 0; out_idx < interFace.glob_var_names.size(); ++out_idx) {
      const std::string &name = (interFace.glob_var_names)[out_idx];
      int idx1 = find_string(file1.Global_Var_Names(), name, interFace.nocase_var_names);
      int idx2 = find_string(file2.Global_Var_Names(), name, interFace.nocase_var_names);
      if (idx1 < 0 || idx2 < 0) {
        Error(fmt::format("Unable to find global variable named '{}' on database.\n", name));
      }
      gvals[out_idx] = FileDiff(vals1[idx1], vals2[idx2], interFace.output_type);
    }
    ex_put_var(out_file_id, output_step, EX_GLOBAL, 1, 0, interFace.glob_var_names.size(),
               Data(gvals));
    return diff_flag;
  }

  // -------------------------------------------------------------------
  // Determine if any diffs and output to terminal
  if (!interFace.quiet_flag && !interFace.glob_var_names.empty()) {
    fmt::print("Global variables:\n");
  }
  for (unsigned out_idx = 0; out_idx < interFace.glob_var_names.size(); ++out_idx) {
    const std::string &name = (interFace.glob_var_names)[out_idx];
    int idx1 = find_string(file1.Global_Var_Names(), name, interFace.nocase_var_names);
    int idx2 = find_string(file2.Global_Var_Names(), name, interFace.nocase_var_names);
    if (idx1 < 0 || idx2 < 0) {
      Error(fmt::format("Unable to find global variable named '{}' on database.\n", name));
    }

    if (Invalid_Values(&vals1[idx1], 1)) {
      Warning(fmt::format("NaN found for global variable '{}' in file 1\n", name));
      diff_flag = true;
    }

    if (Invalid_Values(&vals2[idx2], 1)) {
      Warning(fmt::format("NaN found for global variable '{}' in file 2\n", name));
      diff_flag = true;
    }

    if (interFace.glob_var[out_idx].Diff(vals1[idx1], vals2[idx2])) {
      diff_flag = true;

      if (!interFace.quiet_flag) {
        std::string buf =
            fmt::format("   {:<{}} {} diff: {:14.7e} ~ {:14.7e} ={:12.5e} (FAILED)", name,
                        name_length(), interFace.glob_var[out_idx].abrstr(), vals1[idx1],
                        vals2[idx2], interFace.glob_var[out_idx].Delta(vals1[idx1], vals2[idx2]));
        DIFF_OUT(buf);
      }
      else {
        Die_TS(step1);
      }
    }
  }
  return diff_flag;
}

template <typename INT>
bool diff_nodals(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, int step1, const TimeInterp &t2,
                 int out_file_id, int output_step, const std::vector<INT> &node_map,
                 const INT *id_map, std::vector<double> &nvals)
{
  bool diff_flag = false;

  // ---------------------------------------------------------------------
  // Output file containing differences...
  if (out_file_id >= 0) {
    SMART_ASSERT(!nvals.empty());
    int step2 = t2.step1;
    for (unsigned n_idx = 0; n_idx < interFace.node_var_names.size(); ++n_idx) {
      const std::string &name = (interFace.node_var_names)[n_idx];
      int idx1 = find_string(file1.Nodal_Var_Names(), name, interFace.nocase_var_names);
      int idx2 = find_string(file2.Nodal_Var_Names(), name, interFace.nocase_var_names);
      if (idx1 < 0 || idx2 < 0) {
        Error(fmt::format("Unable to find nodal variable named '{}' on database.\n", name));
      }

      const double *vals1 = get_nodal_values(file1, step1, idx1, 1, name, &diff_flag);
      const double *vals2 = get_nodal_values(file2, step2, idx2, 2, name, &diff_flag);

      if (vals1 == nullptr) {
        Error("Could not find nodal variables on file 1\n");
      }

      if (vals2 == nullptr) {
        Error("Could not find nodal variables on file 2\n");
      }

      size_t ncount = file1.Num_Nodes();
      for (size_t n = 0; n < ncount; ++n) {

        // Should this node be processed...
        if (node_map.empty() || node_map[n] >= 0) {
          INT n2   = node_map.empty() ? n : node_map[n];
          nvals[n] = FileDiff(vals1[n], vals2[n2], interFace.output_type);
        }
        else {
          nvals[n] = 0.;
        }
      } // End of node iteration...
      ex_put_var(out_file_id, output_step, EX_NODAL, n_idx + 1, 0, file1.Num_Nodes(), Data(nvals));
      file1.Free_Nodal_Results(idx1);
      file2.Free_Nodal_Results(idx2);
    }
    file1.Free_Nodal_Results();
    file2.Free_Nodal_Results();
    return diff_flag;
  }

  SMART_ASSERT(!interFace.summary_flag && out_file_id < 0);
  // ----------------------------------------------------------------------
  // Determine if any diffs and output to terminal
  if (!interFace.quiet_flag && !interFace.node_var_names.empty()) {
    fmt::print("Nodal variables:\n");
  }
  for (unsigned n_idx = 0; n_idx < interFace.node_var_names.size(); ++n_idx) {
    const std::string &name = (interFace.node_var_names)[n_idx];
    int idx1 = find_string(file1.Nodal_Var_Names(), name, interFace.nocase_var_names);
    int idx2 = find_string(file2.Nodal_Var_Names(), name, interFace.nocase_var_names);
    if (idx1 < 0 || idx2 < 0) {
      Error(fmt::format("Unable to find nodal variable named '{}' on database.\n", name));
    }

    const double *vals1 = get_nodal_values(file1, step1, idx1, 1, name, &diff_flag);
    const double *vals2 = get_nodal_values(file2, t2, idx2, 2, name, &diff_flag);

    if (vals1 == nullptr) {
      Warning(fmt::format("Could not find nodal variable '{}' on file 1.\n", name));
      diff_flag = true;
      continue;
    }

    if (vals2 == nullptr) {
      Warning(fmt::format("Could not find nodal variable '{}' on file 2.\n", name));
      diff_flag = true;
      continue;
    }

    DiffData max_diff;
    Norm     norm;

    size_t ncount = file1.Num_Nodes();
    for (size_t n = 0; n < ncount; ++n) {

      // Should this node be processed...
      if (node_map.empty() || node_map[n] >= 0) {
        INT    n2 = node_map.empty() ? n : node_map[n];
        double d  = interFace.node_var[n_idx].Delta(vals1[n], vals2[n2]);
        if (interFace.show_all_diffs) {
          if (d > interFace.node_var[n_idx].value) {
            diff_flag       = true;
            std::string buf = fmt::format(
                "   {:<{}} {} diff: {:14.7e} ~ {:14.7e} ={:12.5e} (node {})", name, name_length(),
                interFace.node_var[n_idx].abrstr(), vals1[n], vals2[n2], d, id_map[n]);
            DIFF_OUT(buf);
          }
        }
        else {
          max_diff.set_max(d, vals1[n], vals2[n2], n);
        }
        norm.add_value(vals1[n], vals2[n2]);
      }
    } // End of node iteration...

    output_norms(norm, name);

    if (max_diff.diff > interFace.node_var[n_idx].value) {
      diff_flag = true;
      if (!interFace.quiet_flag) {
        std::string buf =
            fmt::format("   {:<{}} {} diff: {:14.7e} ~ {:14.7e} ={:12.5e} (node {})", name,
                        name_length(), interFace.node_var[n_idx].abrstr(), max_diff.val1,
                        max_diff.val2, max_diff.diff, id_map[max_diff.id]);
        DIFF_OUT(buf);
      }
      else {
        Die_TS(step1);
      }
    }
    file1.Free_Nodal_Results(idx1);
    file2.Free_Nodal_Results(idx2);
  }
  file1.Free_Nodal_Results();
  file2.Free_Nodal_Results();
  return diff_flag;
}

template <typename INT>
bool diff_element(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, int step1, const TimeInterp &t2,
                  int out_file_id, int output_step, const std::vector<INT> &elmt_map,
                  const INT *id_map, Exo_Block<INT> **blocks2, std::vector<double> &evals)
{
  bool diff_flag = false;

  if (out_file_id >= 0) {
    SMART_ASSERT(!evals.empty());
  }

  if (out_file_id < 0 && !interFace.quiet_flag && !interFace.elmt_var_names.empty()) {
    fmt::print("Element variables:\n");
  }

  for (unsigned e_idx = 0; e_idx < interFace.elmt_var_names.size(); ++e_idx) {
    const std::string &name = (interFace.elmt_var_names)[e_idx];
    int vidx1 = find_string(file1.Element_Var_Names(), name, interFace.nocase_var_names);
    int vidx2 = find_string(file2.Element_Var_Names(), name, interFace.nocase_var_names);
    if (vidx1 < 0 || vidx2 < 0) {
      Error(fmt::format("Unable to find element variable named '{}' on database.\n", name));
    }

    Norm norm;

    if (!elmt_map.empty()) { // Load variable for all blocks in file 2.
      for (size_t b = 0; b < file2.Num_Element_Blocks(); ++b) {
        Exo_Block<INT> *block2 = file2.Get_Element_Block_by_Index(b);
        block2->Load_Results(t2.step1, t2.step2, t2.proportion, vidx2);
      }
    }

    size_t   global_elmt_index = 0;
    DiffData max_diff;
    for (size_t b = 0; b < file1.Num_Element_Blocks(); ++b) {
      Exo_Block<INT> *eblock1 = file1.Get_Element_Block_by_Index(b);
      if (!eblock1->is_valid_var(vidx1)) {
        global_elmt_index += eblock1->Size();
        continue;
      }
      if (eblock1->Size() == 0) {
        continue;
      }

      Exo_Block<INT> *eblock2 = nullptr;
      if (elmt_map.empty()) {
        if (interFace.by_name) {
          eblock2 = file2.Get_Element_Block_by_Name(eblock1->Name());
        }
        else {
          eblock2 = file2.Get_Element_Block_by_Id(eblock1->Id());
        }

        SMART_ASSERT(eblock2 != nullptr);
        if (!eblock2->is_valid_var(vidx2)) {
          continue;
        }
      }

      eblock1->Load_Results(step1, vidx1);
      const double *vals1 = eblock1->Get_Results(vidx1);
      if (vals1 == nullptr) {
        Warning(fmt::format("Could not find element variable '{}' in block {}, file 1.\n", name,
                            eblock1->Id()));
        diff_flag = true;
        continue;
      }

      if (Invalid_Values(vals1, eblock1->Size())) {
        Warning(fmt::format("NaN found for element variable '{}' in block {}, file 1\n", name,
                            eblock1->Id()));
        diff_flag = true;
      }

      double        v2    = 0;
      const double *vals2 = nullptr;

      if (elmt_map.empty()) {
        // Without mapping, get result for this block.
        size_t id = eblock1->Id();
        if (interFace.by_name) {
          eblock2 = file2.Get_Element_Block_by_Name(eblock1->Name());
        }
        else {
          eblock2 = file2.Get_Element_Block_by_Id(id);
        }
        vals2 = get_validated_variable(eblock2, t2, vidx2, name, &diff_flag);
        if (vals2 == nullptr) {
          continue;
        }
      }

      size_t ecount   = eblock1->Size();
      size_t block_id = eblock1->Id();
      for (size_t e = 0; e < ecount; ++e) {
        if (out_file_id >= 0) {
          evals[e] = 0.;
        }
        INT el_flag = 1;
        if (!elmt_map.empty()) {
          el_flag = elmt_map[global_elmt_index];
        }

        if (el_flag >= 0) {
          if (elmt_map.empty()) {
            if (vals2 != nullptr) {
              v2 = vals2[e];
            }
          }
          else {
            // With mapping, map global index from file 1 to global index
            // for file 2.  Then convert to block index and elmt index.
            auto bl_idx = file2.Global_to_Block_Local(elmt_map[global_elmt_index] + 1);
            SMART_ASSERT(blocks2[bl_idx.first] != nullptr);
            if (blocks2[bl_idx.first]->is_valid_var(vidx2)) {
              auto *tmp = blocks2[bl_idx.first]->Get_Results(vidx2);
              if (tmp != nullptr) {
                v2 = tmp[bl_idx.second]; // Get value from file 2.
              }
              else {
                v2 = vals1[e]; // Should never happen...
              }
            }
            else {
              // Easiest from logic standpoint to just set v2 equal to v1 at
              // this point and continue through rest of loop.
              v2 = vals1[e];
            }
          }

          if (out_file_id >= 0) {
            evals[e] = FileDiff(vals1[e], v2, interFace.output_type);
          }
          else if (interFace.show_all_diffs) {
            double d = interFace.elmt_var[e_idx].Delta(vals1[e], v2);
            if (d > interFace.elmt_var[e_idx].value) {
              diff_flag       = true;
              std::string buf = fmt::format(
                  "   {:<{}} {} diff: {:14.7e} ~ {:14.7e} ={:12.5e} (block {}, elmt {})", name,
                  name_length(), interFace.elmt_var[e_idx].abrstr(), vals1[e], v2, d, block_id,
                  id_map[global_elmt_index]);
              DIFF_OUT(buf);
            }
          }
          else {
            double d = interFace.elmt_var[e_idx].Delta(vals1[e], v2);
            max_diff.set_max(d, vals1[e], v2, global_elmt_index, block_id);
          }
          norm.add_value(vals1[e], v2);
        }
        ++global_elmt_index;
      }

      if (out_file_id >= 0) {
        ex_put_var(out_file_id, output_step, EX_ELEM_BLOCK, e_idx + 1, eblock1->Id(),
                   eblock1->Size(), Data(evals));
      }

      eblock1->Free_Results();
      if (elmt_map.empty() && eblock2 != nullptr) {
        eblock2->Free_Results();
      }

    } // End of element block loop.

    output_norms(norm, name);

    if (max_diff.diff > interFace.elmt_var[e_idx].value) {
      diff_flag = true;

      if (!interFace.quiet_flag) {
        std::string buf =
            fmt::format("   {:<{}} {} diff: {:14.7e} ~ {:14.7e} ={:12.5e} (block {}, elmt {})",
                        name, name_length(), interFace.elmt_var[e_idx].abrstr(), max_diff.val1,
                        max_diff.val2, max_diff.diff, max_diff.blk, id_map[max_diff.id]);
        DIFF_OUT(buf);
      }
      else {
        Die_TS(step1);
      }
    }

  } // End of element variable loop.
  return diff_flag;
}

template <typename INT>
bool diff_nodeset(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, int step1, const TimeInterp &t2,
                  int out_file_id, int output_step, const INT *id_map, std::vector<double> &vals)
{
  bool diff_flag = false;

  if (out_file_id >= 0) {
    SMART_ASSERT(!vals.empty());
  }

  if (out_file_id < 0 && !interFace.quiet_flag && !interFace.ns_var_names.empty()) {
    fmt::print("Nodeset variables:\n");
  }
  for (unsigned e_idx = 0; e_idx < interFace.ns_var_names.size(); ++e_idx) {
    const std::string &name  = (interFace.ns_var_names)[e_idx];
    int                vidx1 = find_string(file1.NS_Var_Names(), name, interFace.nocase_var_names);
    int                vidx2 = find_string(file2.NS_Var_Names(), name, interFace.nocase_var_names);
    if (vidx1 < 0 || vidx2 < 0) {
      Error(fmt::format("Unable to find nodeset variable named '{}' on database.\n", name));
    }

    DiffData max_diff;
    Norm     norm;

    for (size_t b = 0; b < file1.Num_Node_Sets(); ++b) {
      Node_Set<INT> *nset1 = file1.Get_Node_Set_by_Index(b);
      const double  *vals1 = get_validated_variable(nset1, step1, vidx1, name, &diff_flag);
      if (vals1 == nullptr) {
        continue;
      }

      Node_Set<INT> *nset2 = nullptr;
      size_t         id    = nset1->Id();
      if (interFace.by_name) {
        nset2 = file2.Get_Node_Set_by_Name(nset1->Name());
      }
      else {
        nset2 = file2.Get_Node_Set_by_Id(id);
      }
      const double *vals2 = get_validated_variable(nset2, t2, vidx2, name, &diff_flag);
      if (vals2 == nullptr) {
        continue;
      }

      size_t ncount = nset1->Size();
      if (nset2->Size() == ncount) {
        for (size_t e = 0; e < ncount; ++e) {
          int idx1 = nset1->Node_Index(e);
          int idx2 = nset2->Node_Index(e);

          if (out_file_id >= 0) {
            vals[idx1] = FileDiff(vals1[idx1], vals2[idx2], interFace.output_type);
          }
          else if (interFace.show_all_diffs) {
            double d = interFace.ns_var[e_idx].Delta(vals1[idx1], vals2[idx2]);
            if (d > interFace.ns_var[e_idx].value) {
              diff_flag = true;
              std::string buf =
                  fmt::format("   {:<{}} {} diff: {:14.7e} ~ {:14.7e} ={:12.5e} (set {}, node {})",
                              name, name_length(), interFace.ns_var[e_idx].abrstr(), vals1[idx1],
                              vals2[idx2], d, nset1->Id(), e);
              DIFF_OUT(buf);
            }
          }
          else {
            double d = interFace.ns_var[e_idx].Delta(vals1[idx1], vals2[idx2]);
            max_diff.set_max(d, vals1[idx1], vals2[idx2], e, nset1->Id());
          }
          norm.add_value(vals1[idx1], vals2[idx2]);
        }

        if (out_file_id >= 0) {
          ex_put_var(out_file_id, output_step, EX_NODE_SET, e_idx + 1, nset1->Id(), nset1->Size(),
                     Data(vals));
        }
      }
      else {
        std::string buf =
            fmt::format("   {:<{}}     diff: nodeset node counts differ for nodeset {}", name,
                        name_length(), nset1->Id());
        DIFF_OUT(buf);
        diff_flag = true;
      }

      nset1->Free_Results();
      nset2->Free_Results();
    } // End of nodeset loop.

    output_norms(norm, name);

    if (max_diff.diff > interFace.ns_var[e_idx].value) {
      diff_flag = true;

      if (!interFace.quiet_flag) {
        Node_Set<INT> *nset = file1.Get_Node_Set_by_Id(max_diff.blk);
        std::string    buf  = fmt::format(
            "   {:<{}} {} diff: {:14.7e} ~ {:14.7e} ={:12.5e} (set {}, node {})", name,
            name_length(), interFace.ns_var[e_idx].abrstr(), max_diff.val1, max_diff.val2,
            max_diff.diff, max_diff.blk, id_map[nset->Node_Id(max_diff.id) - 1]);
        DIFF_OUT(buf);
      }
      else {
        Die_TS(step1);
      }
    }
  } // End of nodeset variable loop.
  return diff_flag;
}

template <typename INT>
bool diff_sideset(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, int step1, const TimeInterp &t2,
                  int out_file_id, int output_step, const INT *id_map, std::vector<double> &vals)
{
  bool diff_flag = false;

  if (out_file_id >= 0) {
    SMART_ASSERT(!vals.empty());
  }

  if (out_file_id < 0 && !interFace.quiet_flag && !interFace.ss_var_names.empty()) {
    fmt::print("Sideset variables:\n");
  }

  for (unsigned e_idx = 0; e_idx < interFace.ss_var_names.size(); ++e_idx) {
    const std::string &name  = (interFace.ss_var_names)[e_idx];
    int                vidx1 = find_string(file1.SS_Var_Names(), name, interFace.nocase_var_names);
    int                vidx2 = find_string(file2.SS_Var_Names(), name, interFace.nocase_var_names);

    if (vidx1 < 0 || vidx2 < 0) {
      Error(fmt::format("Unable to find sideset variable named '{}' on database.\n", name));
    }

    DiffData max_diff;
    Norm     norm;

    for (size_t b = 0; b < file1.Num_Side_Sets(); ++b) {
      Side_Set<INT> *sset1 = file1.Get_Side_Set_by_Index(b);
      SMART_ASSERT(sset1 != nullptr);
      const double *vals1 = get_validated_variable(sset1, step1, vidx1, name, &diff_flag);
      if (vals1 == nullptr) {
        continue;
      }

      Side_Set<INT> *sset2 = nullptr;
      if (interFace.by_name) {
        sset2 = file2.Get_Side_Set_by_Name(sset1->Name());
      }
      else {
        sset2 = file2.Get_Side_Set_by_Id(sset1->Id());
      }
      const double *vals2 = get_validated_variable(sset2, t2, vidx2, name, &diff_flag);
      if (vals2 == nullptr) {
        continue;
      }

      size_t ecount = sset1->Size();
      if (sset2->Size() == ecount) {
        for (size_t e = 0; e < ecount; ++e) {
          size_t ind1 = sset1->Side_Index(e);
          size_t ind2 = sset2->Side_Index(e);

          if (out_file_id >= 0) {
            vals[ind1] = FileDiff(vals1[ind1], vals2[ind2], interFace.output_type);
          }
          else if (interFace.show_all_diffs) {
            double d = interFace.ss_var[e_idx].Delta(vals1[ind1], vals2[ind2]);
            if (d > interFace.ss_var[e_idx].value) {
              diff_flag       = true;
              std::string buf = fmt::format(
                  "   {:<{}} {} diff: {:14.7e} ~ {:14.7e} ={:12.5e} (set {}, side {}.{})", name,
                  name_length(), interFace.ss_var[e_idx].abrstr(), vals1[ind1], vals2[ind2], d,
                  sset1->Id(), id_map[sset1->Side_Id(e).first - 1], (int)sset1->Side_Id(e).second);
              DIFF_OUT(buf);
            }
          }
          else {
            double d = interFace.ss_var[e_idx].Delta(vals1[ind1], vals2[ind2]);
            max_diff.set_max(d, vals1[ind1], vals2[ind2], e, sset1->Id());
          }
          norm.add_value(vals1[ind1], vals2[ind2]);
        }
        if (out_file_id >= 0) {
          ex_put_var(out_file_id, output_step, EX_SIDE_SET, e_idx + 1, sset1->Id(), sset1->Size(),
                     Data(vals));
        }
      }
      else {
        std::string buf =
            fmt::format("   {:<{}}     diff: sideset side counts differ for sideset {}", name,
                        name_length(), sset1->Id());
        DIFF_OUT(buf);
        diff_flag = true;
      }

      sset1->Free_Results();
      sset2->Free_Results();
    } // End of sideset loop.

    if (max_diff.diff > interFace.ss_var[e_idx].value) {
      diff_flag = true;

      output_norms(norm, name);

      if (!interFace.quiet_flag) {
        Side_Set<INT> *sset = file1.Get_Side_Set_by_Id(max_diff.blk);
        std::string    buf  = fmt::format(
            "   {:<{}} {} diff: {:14.7e} ~ {:14.7e} ={:12.5e} (set {}, side {}.{})", name,
            name_length(), interFace.ss_var[e_idx].abrstr(), max_diff.val1, max_diff.val2,
            max_diff.diff, max_diff.blk, id_map[sset->Side_Id(max_diff.id).first - 1],
            (int)sset->Side_Id(max_diff.id).second);
        DIFF_OUT(buf);
      }
      else {
        Die_TS(step1);
      }
    }

  } // End of sideset variable loop.
  return diff_flag;
}

template <typename INT>
bool diff_sideset_df(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, const INT *id_map)
{
  bool diff_flag = false;

  std::string name        = "Distribution Factors";
  int         length_name = name.length();

  if (!interFace.quiet_flag && file1.Num_Side_Sets() > 0) {
    fmt::print("Sideset Distribution Factors:\n");
  }
  DiffData max_diff;
  for (size_t b = 0; b < file1.Num_Side_Sets(); ++b) {
    Side_Set<INT> *sset1 = file1.Get_Side_Set_by_Index(b);
    SMART_ASSERT(sset1 != nullptr);

    Side_Set<INT> *sset2 = nullptr;
    if (interFace.by_name) {
      sset2 = file2.Get_Side_Set_by_Name(sset1->Name());
    }
    else {
      sset2 = file2.Get_Side_Set_by_Id(sset1->Id());
    }
    if (sset2 == nullptr) {
      continue;
    }

    if (sset1->Distribution_Factor_Count() == 0 || sset2->Distribution_Factor_Count() == 0) {
      continue;
    }

    const double *vals1 = sset1->Distribution_Factors();

    if (vals1 == nullptr) {
      Warning(
          fmt::format("Could not read distribution factors in sideset {}, file 1.\n", sset1->Id()));
      diff_flag = true;
      continue;
    }

    double value1 = 0.0;
    double value2 = 0.0;
    bool   same1  = false;
    bool   same2  = false;

    size_t ecount = sset1->Size();

    {
      std::pair<INT, INT> range1 = sset1->Distribution_Factor_Range(ecount - 1);
      if (Invalid_Values(vals1, range1.second)) {
        Warning(fmt::format("NaN found for distribution factors in sideset {}, file 1.\n",
                            sset1->Id()));
        diff_flag = true;
      }

      // See if all df are the same value:
      same1 = Equal_Values(vals1, range1.second, &value1);
    }

    auto *vals2 = sset2->Distribution_Factors();

    if (vals2 == nullptr) {
      Warning(
          fmt::format("Could not read distribution factors in sideset {}, file 2.\n", sset2->Id()));
      diff_flag = true;
      continue;
    }

    {
      std::pair<INT, INT> range2 = sset2->Distribution_Factor_Range(sset2->Size() - 1);
      if (Invalid_Values(vals2, range2.second)) {
        Warning(fmt::format("NaN found for distribution factors in sideset {}, file 2.\n",
                            sset2->Id()));
        diff_flag = true;
      }

      // See if all df are the same value:
      same2 = Equal_Values(vals2, range2.second, &value2);
    }

    if (same1 && same2 && (value1 == value2)) {
      continue;
    }

    if (sset2->Size() == ecount) {
      for (size_t e = 0; e < ecount; ++e) {
        std::pair<INT, INT> range1 = sset1->Distribution_Factor_Range(e);
        std::pair<INT, INT> range2 = sset2->Distribution_Factor_Range(e);
        SMART_ASSERT(range1.second - range1.first == range2.second - range2.first);

        for (INT i = 0; i < range1.second - range1.first; i++) {
          double v1 = vals1[range1.first + i];
          double v2 = vals2[range2.first + i];

          if (interFace.show_all_diffs) {
            double d = interFace.ss_df_tol.Delta(v1, v2);
            if (d > interFace.ss_df_tol.value) {
              diff_flag       = true;
              std::string buf = fmt::format(
                  "   {:<{}} {} diff: {:14.7e} ~ {:14.7e} ={:12.5e} (set {}, side {}"
                  ".{}-{})",
                  name, length_name, interFace.ss_df_tol.abrstr(), v1, v2, d, sset1->Id(),
                  id_map[sset1->Side_Id(e).first - 1], (int)sset1->Side_Id(e).second, (int)i + 1);
              DIFF_OUT(buf);
            }
          }
          else {
            double d = interFace.ss_df_tol.Delta(v1, v2);
            max_diff.set_max(d, v1, v2, e, sset1->Id());
          }
        }
      }
    }
    else {
      std::string buf = fmt::format("   {:<{}}     diff: sideset side counts differ for sideset {}",
                                    name, length_name, sset1->Id());
      DIFF_OUT(buf);
      diff_flag = true;
    }

    sset1->Free_Distribution_Factors();
    sset2->Free_Distribution_Factors();
  } // End of sideset loop.

  if (max_diff.diff > interFace.ss_df_tol.value) {
    diff_flag = true;

    if (!interFace.quiet_flag) {
      Side_Set<INT> *sset = file1.Get_Side_Set_by_Id(max_diff.blk);
      std::string    buf =
          fmt::format("   {:<{}} {} diff: {:14.7e} ~ {:14.7e} ={:12.5e} (set {}, side {}.{})", name,
                      length_name, interFace.ss_df_tol.abrstr(), max_diff.val1, max_diff.val2,
                      max_diff.diff, max_diff.blk, id_map[sset->Side_Id(max_diff.id).first - 1],
                      (int)sset->Side_Id(max_diff.id).second);
      DIFF_OUT(buf);
    }
    else {
      Die_TS(-1);
    }
  }

  return diff_flag;
}

template <typename INT>
bool diff_edgeblock(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, int step1, const TimeInterp &t2,
                    int out_file_id, int output_step, const INT * /* id_map */,
                    std::vector<double> &vals)
{
  bool diff_flag = false;

  if (out_file_id >= 0) {
    SMART_ASSERT(!vals.empty());
  }

  if (out_file_id < 0 && !interFace.quiet_flag && !interFace.eb_var_names.empty()) {
    fmt::print("Edge Block variables:\n");
  }

  for (unsigned e_idx = 0; e_idx < interFace.eb_var_names.size(); ++e_idx) {
    const std::string &name  = (interFace.eb_var_names)[e_idx];
    int                vidx1 = find_string(file1.EB_Var_Names(), name, interFace.nocase_var_names);
    int                vidx2 = find_string(file2.EB_Var_Names(), name, interFace.nocase_var_names);
    if (vidx1 < 0 || vidx2 < 0) {
      Error(fmt::format("Unable to find edge block variable named '{}' on database.\n", name));
    }

    DiffData max_diff;
    Norm     norm;

    for (size_t b = 0; b < file1.Num_Edge_Blocks(); ++b) {
      Edge_Block<INT> *eblock1 = file1.Get_Edge_Block_by_Index(b);
      const double    *vals1   = get_validated_variable(eblock1, step1, vidx1, name, &diff_flag);
      if (vals1 == nullptr) {
        continue;
      }

      Edge_Block<INT> *eblock2 = nullptr;
      if (interFace.by_name) {
        eblock2 = file2.Get_Edge_Block_by_Name(eblock1->Name());
      }
      else {
        eblock2 = file2.Get_Edge_Block_by_Id(eblock1->Id());
      }
      const double *vals2 = get_validated_variable(eblock2, t2, vidx2, name, &diff_flag);
      if (vals2 == nullptr) {
        continue;
      }

      size_t ecount = eblock1->Size();
      if (eblock2->Size() == ecount) {
        for (size_t e = 0; e < ecount; ++e) {
          size_t ind1 = eblock1->Edge_Index(e);
          size_t ind2 = eblock2->Edge_Index(e);

          if (out_file_id >= 0) {
            vals[ind1] = FileDiff(vals1[ind1], vals2[ind2], interFace.output_type);
          }
          else if (interFace.show_all_diffs) {
            double d = interFace.eb_var[e_idx].Delta(vals1[ind1], vals2[ind2]);
            if (d > interFace.eb_var[e_idx].value) {
              diff_flag       = true;
              std::string buf = fmt::format(
                  "   {:<{}} {} diff: {:14.7e} ~ {:14.7e} ={:12.5e} (edge block {}, edge {})", name,
                  name_length(), interFace.eb_var[e_idx].abrstr(), vals1[ind1], vals2[ind2], d,
                  eblock1->Id(), eblock1->Edge_Index(e) + 1);
              DIFF_OUT(buf);
            }
          }
          else {
            double d = interFace.eb_var[e_idx].Delta(vals1[ind1], vals2[ind2]);
            max_diff.set_max(d, vals1[ind1], vals2[ind2], e, eblock1->Id());
          }
          norm.add_value(vals1[ind1], vals2[ind2]);
        }
        if (out_file_id >= 0) {
          ex_put_var(out_file_id, output_step, EX_EDGE_BLOCK, e_idx + 1, eblock1->Id(),
                     eblock1->Size(), Data(vals));
        }
      }
      else {
        std::string buf =
            fmt::format("   {:<{}}     diff: edge block edge counts differ for edge block {}", name,
                        name_length(), eblock1->Id());
        DIFF_OUT(buf);
        diff_flag = true;
      }

      eblock1->Free_Results();
      eblock2->Free_Results();
    } // End of edgeblock loop.

    if (max_diff.diff > interFace.eb_var[e_idx].value) {
      diff_flag = true;

      output_norms(norm, name);

      if (!interFace.quiet_flag) {
        Edge_Block<INT> *eblock = file1.Get_Edge_Block_by_Id(max_diff.blk);
        std::string      buf    = fmt::format(
            "   {:<{}} {} diff: {:14.7e} ~ {:14.7e} ={:12.5e} (edge block {}, edge {})", name,
            name_length(), interFace.eb_var[e_idx].abrstr(), max_diff.val1, max_diff.val2,
            max_diff.diff, max_diff.blk, eblock->Edge_Index(max_diff.id) + 1);
        DIFF_OUT(buf);
      }
      else {
        Die_TS(step1);
      }
    }

  } // End of edgeblock variable loop.
  return diff_flag;
}

template <typename INT>
bool diff_faceblock(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, int step1, const TimeInterp &t2,
                    int out_file_id, int output_step, const INT * /* id_map */,
                    std::vector<double> &vals)
{
  bool diff_flag = false;

  if (out_file_id >= 0) {
    SMART_ASSERT(!vals.empty());
  }

  if (out_file_id < 0 && !interFace.quiet_flag && !interFace.fb_var_names.empty()) {
    fmt::print("Face Block variables:\n");
  }

  for (unsigned f_idx = 0; f_idx < interFace.fb_var_names.size(); ++f_idx) {
    const std::string &name  = (interFace.fb_var_names)[f_idx];
    int                vidx1 = find_string(file1.FB_Var_Names(), name, interFace.nocase_var_names);
    int                vidx2 = find_string(file2.FB_Var_Names(), name, interFace.nocase_var_names);
    if (vidx1 < 0 || vidx2 < 0) {
      Error(fmt::format("Unable to find face block variable named '{}' on database.\n", name));
    }

    DiffData max_diff;
    Norm     norm;

    for (size_t b = 0; b < file1.Num_Face_Blocks(); ++b) {
      Face_Block<INT> *fblock1 = file1.Get_Face_Block_by_Index(b);
      const double    *vals1   = get_validated_variable(fblock1, step1, vidx1, name, &diff_flag);
      if (vals1 == nullptr) {
        continue;
      }

      Face_Block<INT> *fblock2 = nullptr;
      if (interFace.by_name) {
        fblock2 = file2.Get_Face_Block_by_Name(fblock1->Name());
      }
      else {
        fblock2 = file2.Get_Face_Block_by_Id(fblock1->Id());
      }
      const double *vals2 = get_validated_variable(fblock2, t2, vidx2, name, &diff_flag);
      if (vals2 == nullptr) {
        continue;
      }

      size_t fcount = fblock1->Size();
      if (fblock2->Size() == fcount) {
        for (size_t f = 0; f < fcount; ++f) {
          size_t ind1 = fblock1->Face_Index(f);
          size_t ind2 = fblock2->Face_Index(f);

          if (out_file_id >= 0) {
            vals[ind1] = FileDiff(vals1[ind1], vals2[ind2], interFace.output_type);
          }
          else if (interFace.show_all_diffs) {
            double d = interFace.fb_var[f_idx].Delta(vals1[ind1], vals2[ind2]);
            if (d > interFace.fb_var[f_idx].value) {
              diff_flag       = true;
              std::string buf = fmt::format(
                  "   {:<{}} {} diff: {:14.7e} ~ {:14.7e} ={:12.5e} (face block {}, face {})", name,
                  name_length(), interFace.fb_var[f_idx].abrstr(), vals1[ind1], vals2[ind2], d,
                  fblock1->Id(), fblock1->Face_Index(f) + 1);
              DIFF_OUT(buf);
            }
          }
          else {
            double d = interFace.fb_var[f_idx].Delta(vals1[ind1], vals2[ind2]);
            max_diff.set_max(d, vals1[ind1], vals2[ind2], f, fblock1->Id());
          }
          norm.add_value(vals1[ind1], vals2[ind2]);
        }
        if (out_file_id >= 0) {
          ex_put_var(out_file_id, output_step, EX_FACE_BLOCK, f_idx + 1, fblock1->Id(),
                     fblock1->Size(), Data(vals));
        }
      }
      else {
        std::string buf =
            fmt::format("   {:<{}}     diff: face block face counts differ for face block {}", name,
                        name_length(), fblock1->Id());
        DIFF_OUT(buf);
        diff_flag = true;
      }

      fblock1->Free_Results();
      fblock2->Free_Results();
    } // End of faceblock loop.

    if (max_diff.diff > interFace.fb_var[f_idx].value) {
      diff_flag = true;

      output_norms(norm, name);

      if (!interFace.quiet_flag) {
        std::string buf =
            fmt::format("   {:<{}} {} diff: {:14.7e} ~ {:14.7e} ={:12.5e} (face block {}, face {})",
                        name, name_length(), interFace.fb_var[f_idx].abrstr(), max_diff.val1,
                        max_diff.val2, max_diff.diff, max_diff.blk, max_diff.id + 1);
        DIFF_OUT(buf);
      }
      else {
        Die_TS(step1);
      }
    }

  } // End of faceblock variable loop.
  return diff_flag;
}

template <typename INT>
bool diff_element_attributes(ExoII_Read<INT> &file1, ExoII_Read<INT>          &file2,
                             const std::vector<INT> & /*elmt_map*/, const INT *id_map,
                             Exo_Block<INT> ** /*blocks2*/)
{
  if (interFace.summary_flag) {
    return false;
  }

  if (file1.Num_Elements() == 0 || file2.Num_Elements() == 0) {
    return false;
  }

  bool diff_was_output = false;
  bool diff_flag       = false;

  size_t global_elmt_offset = 0;
  for (size_t b = 0; b < file1.Num_Element_Blocks(); ++b) {
    Exo_Block<INT> *eblock1 = file1.Get_Element_Block_by_Index(b);
    SMART_ASSERT(eblock1 != nullptr);

    size_t block_id = eblock1->Id();

    Exo_Block<INT> *eblock2 = nullptr;
    if (interFace.by_name) {
      eblock2 = file2.Get_Element_Block_by_Name(eblock1->Name());
    }
    else {
      eblock2 = file2.Get_Element_Block_by_Id(block_id);
    }

    SMART_ASSERT(eblock2 != nullptr);

    if (!diff_was_output && (eblock1->attr_count() > 0 || eblock2->attr_count() > 0)) {
      diff_was_output = true;
      fmt::print("Element attributes:\n");
    }

    for (int idx1 = 0; idx1 < eblock1->attr_count(); idx1++) {
      size_t global_elmt_index = global_elmt_offset;

      DiffData           max_diff;
      const std::string &name = eblock1->Get_Attribute_Name(idx1);

      // Find same attribute in eblock2...
      int idx2 = eblock2->Find_Attribute_Index(name);
      if (idx2 < 0) {
        continue;
      }

      // Find name in interFace.elmt_att_names
      int tol_idx = -1;
      for (unsigned e_idx = 0; e_idx < interFace.elmt_att_names.size(); ++e_idx) {
        if (name == (interFace.elmt_att_names)[e_idx]) {
          tol_idx = e_idx;
          break;
        }
      }

      if (tol_idx == -1) {
        continue;
      }

      Norm norm;

      eblock1->Load_Attributes(idx1);
      const double *vals1 = eblock1->Get_Attributes(idx1);

      if (vals1 == nullptr) {
        Warning(fmt::format("Could not find element attribute '{}' in block {}, file 1.\n", name,
                            eblock1->Id()));
        diff_flag = true;
        continue;
      }

      if (Invalid_Values(vals1, eblock1->Size())) {
        Warning(fmt::format("NaN found for element attribute '{}' in block {}, file 1.\n", name,
                            eblock1->Id()));
        diff_flag = true;
      }

      // Without mapping, get result for this block.
      eblock2->Load_Attributes(idx2);
      const double *vals2 = eblock2->Get_Attributes(idx2);

      if (vals2 == nullptr) {
        Warning(fmt::format("Could not find element attribute '{}' in block {}, file 2.\n", name,
                            eblock2->Id()));
        diff_flag = true;
        continue;
      }

      if (Invalid_Values(vals2, eblock2->Size())) {
        Warning(fmt::format("NaN found for element attribute '{}' in block {}, file 2.\n", name,
                            eblock2->Id()));
        diff_flag = true;
      }

      size_t ecount = eblock1->Size();
      for (size_t e = 0; e < ecount; ++e) {

        if (interFace.show_all_diffs) {
          double d = interFace.elmt_att[tol_idx].Delta(vals1[e], vals2[e]);
          if (d > interFace.elmt_att[tol_idx].value) {
            diff_flag = true;
            std::string buf =
                fmt::format("   {:<{}} {} diff: {:14.7e} ~ {:14.7e} ={:12.5e} (block {}, elmt {})",
                            name, name_length(), interFace.elmt_att[tol_idx].abrstr(), vals1[e],
                            vals2[e], d, block_id, id_map[global_elmt_index]);
            DIFF_OUT(buf);
          }
        }
        else {
          double d = interFace.elmt_att[tol_idx].Delta(vals1[e], vals2[e]);
          max_diff.set_max(d, vals1[e], vals2[e], global_elmt_index, block_id);
        }
        norm.add_value(vals1[e], vals2[e]);
        ++global_elmt_index;
      }

      output_norms(norm, name);

      if (max_diff.diff > interFace.elmt_att[tol_idx].value) {
        diff_flag = true;

        if (!interFace.quiet_flag) {
          std::string buf =
              fmt::format("   {:<{}} {} diff: {:14.7e} ~ {:14.7e} ={:12.5e} (block {}, elmt {})",
                          name, name_length(), interFace.elmt_att[tol_idx].abrstr(), max_diff.val1,
                          max_diff.val2, max_diff.diff, max_diff.blk, id_map[max_diff.id]);
          DIFF_OUT(buf);
        }
        else {
          Die_TS(-1);
        }
      }
    } // End of attribute loop.
    eblock1->Free_Attributes();
    eblock2->Free_Attributes();

    global_elmt_offset += eblock1->Size();
  } // End of element block loop.
  return diff_flag;
}

template <typename INT>
void output_summary(ExoII_Read<INT> &file1, MinMaxData &mm_time, std::vector<MinMaxData> &mm_glob,
                    std::vector<MinMaxData> &mm_node, std::vector<MinMaxData> &mm_elmt,
                    std::vector<MinMaxData> &mm_ns, std::vector<MinMaxData> &mm_ss,
                    std::vector<MinMaxData> &mm_eb, std::vector<MinMaxData> &mm_fb,
                    const INT *node_id_map, const INT *elem_id_map)
{
  int i;
  int n;

  fmt::print("# NOTES:  - The min/max values are reporting the min/max in absolute value.\n"
             "#         - Time values (t) are 1-offset time step numbers.\n"
             "#         - Element block numbers are the block ids.\n"
             "#         - Node(n) and element(e) numbers are 1-offset.\n");

  if (interFace.coord_sep) {
    double min_separation = Find_Min_Coord_Sep(file1);
    fmt::print("\nCOORDINATES absolute 1.e-6    # min separation = {}\n", min_separation);
  }
  else {
    fmt::print("\nCOORDINATES absolute 1.e-6    # min separation not calculated\n");
  }

  if (file1.Num_Times() > 0) {
    fmt::print("\nTIME STEPS relative 1.e-6 floor 0.0     # min: ");
    fmt::print("{:15.8g} @ t{} max: {:15.8g} @ t{}\n", mm_time.min_val, mm_time.min_step,
               mm_time.max_val, mm_time.max_step);
  }
  else {
    fmt::print("\n# No TIME STEPS\n");
  }

  n = interFace.glob_var_names.size();
  if (n > 0) {
    fmt::print("GLOBAL VARIABLES relative 1.e-6 floor 0.0\n");
    for (i = 0; i < n; ++i) {
      fmt::print("\t{:<{}}  # min: {:15.8g} @ t{}\tmax: {:15.8g} @ t{}\n",
                 ((interFace.glob_var_names)[i]), name_length(), mm_glob[i].min_val,
                 mm_glob[i].min_step, mm_glob[i].max_val, mm_glob[i].max_step);
    }
  }
  else {
    fmt::print("\n# No GLOBAL VARIABLES\n");
  }

  n = interFace.node_var_names.size();
  if (n > 0 && file1.Num_Nodes() > 0) {
    fmt::print("\nNODAL VARIABLES relative 1.e-6 floor 0.0\n");
    for (i = 0; i < n; ++i) {
      fmt::print("\t{:<{}}  # min: {:15.8g} @ t{},n{}\tmax: {:15.8g} @ t{},n{}\n",
                 ((interFace.node_var_names)[i]), name_length(), mm_node[i].min_val,
                 mm_node[i].min_step, node_id_map[mm_node[i].min_id], mm_node[i].max_val,
                 mm_node[i].max_step, node_id_map[mm_node[i].max_id]);
    }
  }
  else {
    fmt::print("\n# No NODAL VARIABLES and/or NODES\n");
  }

  n = interFace.elmt_var_names.size();
  if (n > 0 && file1.Num_Elements() > 0) {
    fmt::print("\nELEMENT VARIABLES relative 1.e-6 floor 0.0\n");
    for (i = 0; i < n; ++i) {
      fmt::print("\t{:<{}}  # min: {:15.8g} @ t{},b{},e{}\tmax: {:15.8g} @ t{},b{}"
                 ",e{}\n",
                 ((interFace.elmt_var_names)[i]), name_length(), mm_elmt[i].min_val,
                 mm_elmt[i].min_step, mm_elmt[i].min_blk, elem_id_map[mm_elmt[i].min_id],
                 mm_elmt[i].max_val, mm_elmt[i].max_step, mm_elmt[i].max_blk,
                 elem_id_map[mm_elmt[i].max_id]);
    }
  }
  else {
    fmt::print("\n# No ELEMENT VARIABLES and/or ELEMENTS\n");
  }

  n = interFace.ns_var_names.size();
  if (n > 0) {
    fmt::print("\nNODESET VARIABLES relative 1.e-6 floor 0.0\n");
    for (i = 0; i < n; ++i) {
      Node_Set<INT> *nsmin = file1.Get_Node_Set_by_Id(mm_ns[i].min_blk);
      Node_Set<INT> *nsmax = file1.Get_Node_Set_by_Id(mm_ns[i].max_blk);
      fmt::print("\t{:<{}}  # min: {:15.8g} @ t{},s{},n{}\tmax: {:15.8g} @ t{},s{}"
                 ",n{}\n",
                 ((interFace.ns_var_names)[i]), name_length(), mm_ns[i].min_val, mm_ns[i].min_step,
                 mm_ns[i].min_blk, node_id_map[nsmin->Node_Id(mm_ns[i].min_id) - 1],
                 mm_ns[i].max_val, mm_ns[i].max_step, mm_ns[i].max_blk,
                 node_id_map[nsmax->Node_Id(mm_ns[i].max_id) - 1]);
    }
  }
  else {
    fmt::print("\n# No NODESET VARIABLES\n");
  }

  n = interFace.ss_var_names.size();
  if (n > 0) {
    fmt::print("\nSIDESET VARIABLES relative 1.e-6 floor 0.0\n");
    for (i = 0; i < n; ++i) {
      Side_Set<INT>      *ssmin    = file1.Get_Side_Set_by_Id(mm_ss[i].min_blk);
      Side_Set<INT>      *ssmax    = file1.Get_Side_Set_by_Id(mm_ss[i].max_blk);
      std::pair<INT, INT> min_side = ssmin->Side_Id(mm_ss[i].min_id);
      std::pair<INT, INT> max_side = ssmax->Side_Id(mm_ss[i].max_id);
      fmt::print("\t{:<{}}  # min: {:15.8g} @ t{},s{},f{}.{}\tmax: {:15.8g} @ t{},s{}"
                 ",f{}.{}\n",
                 ((interFace.ss_var_names)[i]), name_length(), mm_ss[i].min_val, mm_ss[i].min_step,
                 mm_ss[i].min_blk, elem_id_map[min_side.first - 1], min_side.second,
                 mm_ss[i].max_val, mm_ss[i].max_step, mm_ss[i].max_blk,
                 elem_id_map[max_side.first - 1], max_side.second);
    }
  }
  else {
    fmt::print("\n# No SIDESET VARIABLES\n");
  }

  n = interFace.eb_var_names.size();
  if (n > 0) {
    fmt::print("\nEDGE BLOCK VARIABLES relative 1.e-6 floor 0.0\n");
    for (i = 0; i < n; ++i) {
      fmt::print("\t{:<{}}  # min: {:15.8g} @ t{},b{},e{}\tmax: {:15.8g} @ t{},b{}"
                 ",e{}\n",
                 ((interFace.eb_var_names)[i]), name_length(), mm_eb[i].min_val, mm_eb[i].min_step,
                 mm_eb[i].min_blk, mm_eb[i].min_id + 1, mm_eb[i].max_val, mm_eb[i].max_step,
                 mm_eb[i].max_blk, mm_eb[i].max_id + 1);
    }
  }
  else {
    fmt::print("\n# No EDGE BLOCK VARIABLES\n");
  }

  n = interFace.fb_var_names.size();
  if (n > 0) {
    fmt::print("\nFACE BLOCK VARIABLES relative 1.e-6 floor 0.0\n");
    for (i = 0; i < n; ++i) {
      fmt::print("\t{:<{}}  # min: {:15.8g} @ t{},b{},f{}\tmax: {:15.8g} @ t{},b{}"
                 ",f{}\n",
                 ((interFace.fb_var_names)[i]), name_length(), mm_fb[i].min_val, mm_fb[i].min_step,
                 mm_fb[i].min_blk, mm_fb[i].min_id + 1, mm_fb[i].max_val, mm_fb[i].max_step,
                 mm_fb[i].max_blk, mm_fb[i].max_id + 1);
    }
  }
  else {
    fmt::print("\n# No FACE BLOCK VARIABLES\n");
  }
  fmt::print("\n");
}

int timeStepIsExcluded(int ts)
{
  for (const auto &elem : interFace.exclude_steps) {
    if (ts == elem) {
      return 1;
    }
  }
  return 0;
}
