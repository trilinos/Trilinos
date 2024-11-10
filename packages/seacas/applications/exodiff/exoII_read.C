// Copyright(C) 1999-2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "ED_SystemInterface.h" // for SystemInterface, etc
#include "assembly.h"
#include "edge_block.h" // for Edge_Block
#include "exoII_read.h"
#include "exo_block.h"  // for Exo_Block
#include "exodusII.h"   // for ex_init_params, ex_opts, etc
#include "face_block.h" // for Face_Block
#include "fmt/ostream.h"
#include "node_set.h"     // for Node_Set
#include "side_set.h"     // for Side_Set
#include "smart_assert.h" // for SMART_ASSERT, Assert, etc
#include "stringx.h"      // for chop_whitespace
#include "util.h"         // for free_name_array, etc
#include <algorithm>      // for copy
#include <cstdint>        // for int64_t
#include <cstdio>         // for fclose, FILE, fopen
#include <cstdlib>        // for exit
#include <cstring>        // for strlen
#include <iostream>
#include <set>    // for set
#include <string> // for string, char_traits, etc
#include <vector> // for vector

namespace {
  void read_vars(int file_id, EXOTYPE flag, const char *type, int num_vars,
                 std::vector<std::string> &varlist);
} // namespace

template <typename INT> ExoII_Read<INT>::ExoII_Read() = default;

template <typename INT> ExoII_Read<INT>::ExoII_Read(std::string fname) : file_name(std::move(fname))
{
}

template <typename INT> ExoII_Read<INT>::~ExoII_Read()
{
  try {
    SMART_ASSERT(Check_State());

    if (file_id >= 0) {
      std::string err = Close_File();
      if (!err.empty()) {
        Error(fmt::format("ExoII_Read destructor(): closing file: \"{}\"\n", err));
      }
    }

    delete[] eblocks;
    delete[] nsets;
    delete[] ssets;
    delete[] nodes;
    delete[] times;
    delete[] edge_blocks;
    delete[] face_blocks;
    delete[] assemblies;

    if (results) {
      for (unsigned i = 0; i < nodal_vars.size(); ++i) {
        delete[] results[i];
      }
      delete[] results;
    }
    delete[] global_vals;
    delete[] global_vals2;
    delete[] node_map;
    delete[] elmt_map;
    delete[] elmt_order;
  }
  catch (...) {
  }
}

template <typename INT> std::string ExoII_Read<INT>::Close_File()
{
  SMART_ASSERT(Check_State());

  if (file_id < 0) {
    return "exodiff: ERROR: File is not open!";
  }
  int err = ex_close(file_id);

  if (err < 0) {
    Error(fmt::format("ExoII_Read::Close_File(): {}: Unable to close file!  Aborting...\n", err));
  }
  if (err > 0) {
    return fmt::format("WARNING: {} issued upon close", err);
  }

  file_id = -1;

  return "";
}

template <typename INT> double ExoII_Read<INT>::Time(int time_num) const
{
  SMART_ASSERT(Check_State());
  SMART_ASSERT(time_num > 0 && time_num <= num_times)(time_num)(num_times);
  return times[time_num - 1];
}

template <typename INT> const std::string &ExoII_Read<INT>::Global_Var_Name(int index) const
{
  SMART_ASSERT(Check_State());
  SMART_ASSERT(index >= 0 && (unsigned)index < global_vars.size());
  return global_vars[index];
}

template <typename INT> const std::string &ExoII_Read<INT>::Nodal_Var_Name(int index) const
{
  SMART_ASSERT(Check_State());
  SMART_ASSERT(index >= 0 && (unsigned)index < nodal_vars.size());
  return nodal_vars[index];
}

template <typename INT> const std::string &ExoII_Read<INT>::Element_Var_Name(int index) const
{
  SMART_ASSERT(Check_State());
  SMART_ASSERT(index >= 0 && (unsigned)index < elmt_vars.size());
  return elmt_vars[index];
}

template <typename INT> const std::string &ExoII_Read<INT>::Element_Att_Name(int index) const
{
  SMART_ASSERT(Check_State());
  SMART_ASSERT(index >= 0 && (unsigned)index < elmt_atts.size());
  return elmt_atts[index];
}

template <typename INT> const std::string &ExoII_Read<INT>::NS_Var_Name(int index) const
{
  SMART_ASSERT(Check_State());
  SMART_ASSERT(index >= 0 && (unsigned)index < ns_vars.size());
  return ns_vars[index];
}

template <typename INT> const std::string &ExoII_Read<INT>::SS_Var_Name(int index) const
{
  SMART_ASSERT(Check_State());
  SMART_ASSERT(index >= 0 && (unsigned)index < ss_vars.size());
  return ss_vars[index];
}

template <typename INT> const std::string &ExoII_Read<INT>::EB_Var_Name(int index) const
{
  SMART_ASSERT(Check_State());
  SMART_ASSERT(index >= 0 && (unsigned)index < eb_vars.size());
  return eb_vars[index];
}

template <typename INT> const std::string &ExoII_Read<INT>::FB_Var_Name(int index) const
{
  SMART_ASSERT(Check_State());
  SMART_ASSERT(index >= 0 && (unsigned)index < fb_vars.size());
  return fb_vars[index];
}

template <typename INT>
Exo_Block<INT> *ExoII_Read<INT>::Get_Element_Block_by_Index(size_t block_index) const
{
  SMART_ASSERT(Check_State());
  SMART_ASSERT(block_index < num_elmt_blocks);
  return &eblocks[block_index];
}

template <typename INT> Exo_Block<INT> *ExoII_Read<INT>::Get_Element_Block_by_Id(size_t id) const
{
  SMART_ASSERT(Check_State());
  for (size_t i = 0; i < num_elmt_blocks; i++) {
    if (eblocks[i].Id() == id) {
      return &eblocks[i];
    }
  }
  return nullptr;
}

template <typename INT>
Exo_Block<INT> *ExoII_Read<INT>::Get_Element_Block_by_Name(const std::string &name) const
{
  SMART_ASSERT(Check_State());
  for (size_t i = 0; i < num_elmt_blocks; i++) {
    if (eblocks[i].Name() == name) {
      return &eblocks[i];
    }
  }
  return nullptr;
}

template <typename INT>
Assembly<INT> *ExoII_Read<INT>::Get_Assembly_by_Index(size_t block_index) const
{
  SMART_ASSERT(Check_State());
  SMART_ASSERT(block_index < num_assemblies);
  return &assemblies[block_index];
}

template <typename INT>
Assembly<INT> *ExoII_Read<INT>::Get_Assembly_by_Name(const std::string &name) const
{
  SMART_ASSERT(Check_State());
  for (size_t i = 0; i < num_assemblies; i++) {
    if (assemblies[i].Name() == name) {
      return &assemblies[i];
    }
  }
  return nullptr;
}

template <typename INT> Assembly<INT> *ExoII_Read<INT>::Get_Assembly_by_Id(size_t set_id) const
{
  SMART_ASSERT(Check_State());
  for (size_t i = 0; i < num_assemblies; i++) {
    if (assemblies[i].Id() == set_id) {
      return &assemblies[i];
    }
  }
  return nullptr;
}

template <typename INT>
Exo_Entity *ExoII_Read<INT>::Get_Entity_by_Index(EXOTYPE type, size_t block_index) const
{
  SMART_ASSERT(Check_State());

  switch (type) {
  case EX_ELEM_BLOCK: SMART_ASSERT(block_index < num_elmt_blocks); return &eblocks[block_index];
  case EX_NODE_SET: SMART_ASSERT(block_index < num_node_sets); return &nsets[block_index];
  case EX_SIDE_SET: SMART_ASSERT(block_index < num_side_sets); return &ssets[block_index];
  case EX_EDGE_BLOCK: SMART_ASSERT(block_index < num_edge_blocks); return &edge_blocks[block_index];
  case EX_FACE_BLOCK: SMART_ASSERT(block_index < num_face_blocks); return &face_blocks[block_index];
  case EX_ASSEMBLY: SMART_ASSERT(block_index < num_assemblies); return &assemblies[block_index];
  default: return nullptr;
  }
}

template <typename INT> Exo_Entity *ExoII_Read<INT>::Get_Entity_by_Id(EXOTYPE type, size_t id) const
{
  SMART_ASSERT(Check_State());
  switch (type) {
  case EX_ELEM_BLOCK:
    for (size_t i = 0; i < num_elmt_blocks; i++) {
      if (eblocks[i].Id() == id) {
        return &eblocks[i];
      }
    }
    break;
  case EX_NODE_SET:
    for (size_t i = 0; i < num_node_sets; i++) {
      if (nsets[i].Id() == id) {
        return &nsets[i];
      }
    }
    break;
  case EX_SIDE_SET:
    for (size_t i = 0; i < num_side_sets; i++) {
      if (ssets[i].Id() == id) {
        return &ssets[i];
      }
    }
    break;
  case EX_EDGE_BLOCK:
    for (size_t i = 0; i < num_edge_blocks; i++) {
      if (edge_blocks[i].Id() == id) {
        return &edge_blocks[i];
      }
    }
    break;
  case EX_FACE_BLOCK:
    for (size_t i = 0; i < num_face_blocks; i++) {
      if (face_blocks[i].Id() == id) {
        return &face_blocks[i];
      }
    }
    break;
  case EX_ASSEMBLY:
    for (size_t i = 0; i < num_assemblies; i++) {
      if (assemblies[i].Id() == id) {
        return &assemblies[i];
      }
    }
    break;
  default: return nullptr;
  }
  return nullptr;
}

template <typename INT>
Exo_Entity *ExoII_Read<INT>::Get_Entity_by_Name(EXOTYPE type, const std::string &name) const
{
  SMART_ASSERT(Check_State());
  switch (type) {
  case EX_ELEM_BLOCK:
    for (size_t i = 0; i < num_elmt_blocks; i++) {
      if (eblocks[i].Name() == name) {
        return &eblocks[i];
      }
    }
    break;
  case EX_NODE_SET:
    for (size_t i = 0; i < num_node_sets; i++) {
      if (nsets[i].Name() == name) {
        return &nsets[i];
      }
    }
    break;
  case EX_SIDE_SET:
    for (size_t i = 0; i < num_side_sets; i++) {
      if (ssets[i].Name() == name) {
        return &ssets[i];
      }
    }
    break;
  case EX_EDGE_BLOCK:
    for (size_t i = 0; i < num_edge_blocks; i++) {
      if (edge_blocks[i].Name() == name) {
        return &edge_blocks[i];
      }
    }
    break;
  case EX_FACE_BLOCK:
    for (size_t i = 0; i < num_face_blocks; i++) {
      if (face_blocks[i].Name() == name) {
        return &face_blocks[i];
      }
    }
    break;
  case EX_ASSEMBLY:
    for (size_t i = 0; i < num_assemblies; i++) {
      if (assemblies[i].Name() == name) {
        return &assemblies[i];
      }
    }
    break;
  default: return nullptr;
  }
  return nullptr;
}

template <typename INT> Node_Set<INT> *ExoII_Read<INT>::Get_Node_Set_by_Id(size_t set_id) const
{
  SMART_ASSERT(Check_State());
  for (size_t i = 0; i < num_node_sets; i++) {
    if (nsets[i].Id() == set_id) {
      return &nsets[i];
    }
  }
  return nullptr;
}

template <typename INT>
Node_Set<INT> *ExoII_Read<INT>::Get_Node_Set_by_Name(const std::string &name) const
{
  SMART_ASSERT(Check_State());
  for (size_t i = 0; i < num_node_sets; i++) {
    if (nsets[i].Name() == name) {
      return &nsets[i];
    }
  }
  return nullptr;
}

template <typename INT> Side_Set<INT> *ExoII_Read<INT>::Get_Side_Set_by_Id(size_t set_id) const
{
  SMART_ASSERT(Check_State());
  for (size_t i = 0; i < num_side_sets; i++) {
    if (ssets[i].Id() == set_id) {
      return &ssets[i];
    }
  }
  return nullptr;
}

template <typename INT>
Side_Set<INT> *ExoII_Read<INT>::Get_Side_Set_by_Name(const std::string &name) const
{
  SMART_ASSERT(Check_State());
  for (size_t i = 0; i < num_side_sets; i++) {
    if (ssets[i].Name() == name) {
      return &ssets[i];
    }
  }
  return nullptr;
}

template <typename INT>
Edge_Block<INT> *ExoII_Read<INT>::Get_Edge_Block_by_Id(size_t block_id) const
{
  SMART_ASSERT(Check_State());
  for (size_t i = 0; i < num_edge_blocks; i++) {
    if (edge_blocks[i].Id() == block_id) {
      return &edge_blocks[i];
    }
  }
  return nullptr;
}

template <typename INT>
Edge_Block<INT> *ExoII_Read<INT>::Get_Edge_Block_by_Name(const std::string &name) const
{
  SMART_ASSERT(Check_State());
  for (size_t i = 0; i < num_edge_blocks; i++) {
    if (edge_blocks[i].Name() == name) {
      return &edge_blocks[i];
    }
  }
  return nullptr;
}

template <typename INT>
Face_Block<INT> *ExoII_Read<INT>::Get_Face_Block_by_Id(size_t block_id) const
{
  SMART_ASSERT(Check_State());
  for (size_t i = 0; i < num_face_blocks; i++) {
    if (face_blocks[i].Id() == block_id) {
      return &face_blocks[i];
    }
  }
  return nullptr;
}

template <typename INT>
Face_Block<INT> *ExoII_Read<INT>::Get_Face_Block_by_Name(const std::string &name) const
{
  SMART_ASSERT(Check_State());
  for (size_t i = 0; i < num_face_blocks; i++) {
    if (face_blocks[i].Name() == name) {
      return &face_blocks[i];
    }
  }
  return nullptr;
}

template <typename INT>
std::string ExoII_Read<INT>::Load_Element_Block_Description(size_t block_index) const
{
  SMART_ASSERT(Check_State());
  if (!Open()) {
    return "exodiff: ERROR:  Must open file before loading blocks!";
  }
  SMART_ASSERT(block_index < num_elmt_blocks);

  eblocks[block_index].Load_Connectivity();
  //  eblocks[idx].Load_Connectivity();
  //  eblocks[idx].Load_Attributes();

  return "";
}

template <typename INT> std::string ExoII_Read<INT>::Load_Element_Block_Descriptions() const
{
  SMART_ASSERT(Check_State());
  if (!Open()) {
    return "exodiff: ERROR:  Must open file before loading blocks!";
  }
  for (size_t b = 0; b < num_elmt_blocks; ++b) {
    eblocks[b].Load_Connectivity();
  }

  return "";
}

template <typename INT> std::string ExoII_Read<INT>::Free_Element_Block(size_t block_index) const
{
  SMART_ASSERT(Check_State());

  SMART_ASSERT(block_index < num_elmt_blocks);

  eblocks[block_index].Free_Connectivity();
  eblocks[block_index].Free_Attributes();

  return "";
}

template <typename INT> std::string ExoII_Read<INT>::Free_Element_Blocks() const
{
  SMART_ASSERT(Check_State());

  for (size_t b = 0; b < num_elmt_blocks; ++b) {
    eblocks[b].Free_Connectivity();
    eblocks[b].Free_Attributes();
  }

  return "";
}

template <typename INT> size_t ExoII_Read<INT>::Block_Id(size_t block_index) const
{
  SMART_ASSERT(Check_State());
  SMART_ASSERT(block_index < num_elmt_blocks);
  return eblocks[block_index].Id();
}

template <typename INT> std::string ExoII_Read<INT>::Load_Node_Map()
{
  SMART_ASSERT(Check_State());

  if (!Open()) {
    return "WARNING:  File not open!";
  }
  delete[] node_map;
  node_map = nullptr;

  if (num_nodes == 0) {
    return "WARNING:  There are no nodes!";
  }
  node_map = new INT[num_nodes];
  SMART_ASSERT(node_map != nullptr);

  ex_opts(0); // Temporarily turn off error reporting in case map isn't stored.
  int err = ex_get_id_map(file_id, EX_NODE_MAP, node_map);
  ex_opts(EX_VERBOSE);

  if (err < 0) {
    Error(fmt::format("Unable to load node map; Exodus error = {}.  Aborting...\n", err));
  }
  else if (err > 0) {
    return "WARNING: Default node map being used.";
  }
  return "";
}

template <typename INT> std::string ExoII_Read<INT>::Free_Node_Map()
{
  SMART_ASSERT(Check_State());

  delete[] node_map;
  node_map = nullptr;

  return "";
}

template <typename INT> std::string ExoII_Read<INT>::Load_Element_Map()
{
  SMART_ASSERT(Check_State());

  if (!Open()) {
    return "WARNING:  File not open!";
  }
  delete[] elmt_map;
  elmt_map = nullptr;

  if (num_elmts == 0) {
    return "WARNING:  There are no elements!";
  }
  elmt_map = new INT[num_elmts];
  SMART_ASSERT(elmt_map != nullptr);

  ex_opts(0); // Temporarily turn off error reporting in case map isn't stored.
  int err = ex_get_id_map(file_id, EX_ELEM_MAP, elmt_map);
  ex_opts(EX_VERBOSE);

  if (err < 0) {
    Error(fmt::format("Unable to load element map; Exodus error = {}.  Aborting...\n", err));
  }
  else if (err > 0) {
    return "WARNING: Default element map being used.";
  }
  return "";
}

template <typename INT> std::string ExoII_Read<INT>::Free_Element_Map()
{
  SMART_ASSERT(Check_State());

  delete[] elmt_map;
  elmt_map = nullptr;

  return "";
}

template <typename INT> std::string ExoII_Read<INT>::Load_Nodal_Coordinates()
{
  SMART_ASSERT(Check_State());

  if (!Open()) {
    return "WARNING:  File not open!";
  }
  if (num_nodes) {
    size_t count = num_nodes * dimension;
    nodes        = new double[count];
    SMART_ASSERT(nodes != nullptr);
    double *x = nodes, *y = nodes, *z = nodes;
    if (dimension > 1) {
      y = nodes + num_nodes;
    }
    if (dimension > 2) {
      z = nodes + (2 * num_nodes);
    }

    int err = ex_get_coord(file_id, x, y, z);
    if (err < 0) {
      Error("Failed to get nodal coordinates!  Aborting...\n");
    }
    else if (err > 0) {
      delete[] nodes;
      nodes = nullptr;
      return fmt::format("exodiff: WARNING:  "
                         "Exodus issued warning \"{}\" on call to ex_get_coord()!"
                         "  I'm not going to keep what it gave me for coordinates.",
                         err);
    }
  }
  else {
    return "WARNING:  There are no nodes!";
  }
  return "";
}

template <typename INT> void ExoII_Read<INT>::Free_Nodal_Coordinates()
{
  SMART_ASSERT(Check_State());
  delete[] nodes;
  nodes = nullptr;
}

template <typename INT>
std::string ExoII_Read<INT>::Load_Nodal_Results(int time_step_num, int var_index)
{
  SMART_ASSERT(Check_State());
  SMART_ASSERT(time_step_num > 0 && time_step_num <= num_times);
  SMART_ASSERT(var_index >= 0 && (unsigned)var_index < nodal_vars.size());

  if (!Open()) {
    return "WARNING:  File not open!";
  }
  if (cur_time != time_step_num) {
    for (unsigned i = 0; i < nodal_vars.size(); ++i) {
      delete[] results[i];
      results[i] = nullptr;
    }
    cur_time = time_step_num;
  }

  if (num_nodes) {
    results[var_index] = new double[num_nodes];

    int err =
        ex_get_var(file_id, cur_time, EX_NODAL, var_index + 1, 0, num_nodes, results[var_index]);
    if (err < 0) {
      Error("ExoII_Read::Load_Nodal_Results(): Failed to get "
            "nodal variable values!  Aborting...\n");
    }
    else if (err > 0) {
      delete[] results[var_index];
      results[var_index] = nullptr;
      return fmt::format("ExoII_Read::Load_Nodal_Results(): WARNING:  "
                         "Exodus issued warning \"{}\" on call to ex_get_var()!"
                         "  I'm not going to keep what it gave me for values.",
                         err);
    }
  }
  else {
    return "WARNING:  There are no nodes!";
  }
  return "";
}

template <typename INT>
const double *ExoII_Read<INT>::Get_Nodal_Results(int t1, int t2, double proportion,
                                                 int var_index) const // Interpolated results.
{
  static std::vector<double> st_results;
  static std::vector<double> st_results2;

  SMART_ASSERT(Check_State());
  SMART_ASSERT(t1 > 0 && t1 <= num_times);
  SMART_ASSERT(t2 > 0 && t2 <= num_times);
  SMART_ASSERT(var_index >= 0 && (unsigned)var_index < nodal_vars.size());

  if (!Open()) {
    return nullptr;
  }

  if (st_results.empty()) {
    st_results.resize(num_nodes);
  }

  int err = ex_get_var(file_id, t1, EX_NODAL, var_index + 1, 0, num_nodes, st_results.data());
  if (err < 0) {
    Error("ExoII_Read::Get_Nodal_Results(): Failed to get "
          "nodal variable values!  Aborting...\n");
  }

  if (t1 != t2) {
    if (st_results2.empty()) {
      st_results2.resize(num_nodes);
    }

    err = ex_get_var(file_id, t2, EX_NODAL, var_index + 1, 0, num_nodes, st_results2.data());
    if (err < 0) {
      Error("ExoII_Read::Load_Nodal_Results(): Failed to get "
            "nodal variable values!  Aborting...\n");
    }

    // Interpolate the values...
    for (size_t i = 0; i < num_nodes; i++) {
      st_results[i] = (1.0 - proportion) * st_results[i] + proportion * st_results2[i];
    }
  }
  return st_results.data();
}

template <typename INT> void ExoII_Read<INT>::Free_Nodal_Results()
{
  SMART_ASSERT(Check_State());
  if (results) {
    for (unsigned i = 0; i < nodal_vars.size(); ++i) {
      delete[] results[i];
      results[i] = nullptr;
    }
  }
}

template <typename INT> void ExoII_Read<INT>::Free_Nodal_Results(int var_index)
{
  SMART_ASSERT(Check_State());
  if (results) {
    if (results[var_index]) {
      delete[] results[var_index];
      results[var_index] = nullptr;
    }
  }
}

template <typename INT> const double *ExoII_Read<INT>::Get_Nodal_Results(int var_index) const
{
  SMART_ASSERT(Check_State());
  if (cur_time == 0) {
    return nullptr;
  }
  SMART_ASSERT(var_index >= 0 && (unsigned)var_index < nodal_vars.size());

  return results[var_index];
}

template <typename INT> std::string ExoII_Read<INT>::Load_Global_Results(int time_step_num)
{
  SMART_ASSERT(Check_State());
  SMART_ASSERT(time_step_num > 0 && time_step_num <= num_times);

  if (!Open()) {
    return "WARNING:  File not open!";
  }
  if (global_vars.empty()) {
    return "WARNING:  No global variables! (doing nothing)";
  }

  if (global_vals == nullptr) {
    global_vals = new double[global_vars.size()];
    SMART_ASSERT(global_vals != nullptr);
  }

  for (unsigned j = 0; j < global_vars.size(); ++j) {
    global_vals[j] = 0.0;
  }

  int err = ex_get_var(file_id, time_step_num, EX_GLOBAL, 1, 1, global_vars.size(), global_vals);

  if (err < 0) {
    Error("ExoII_Read::Load_Global_Results(): Failed to get "
          "global variable values!  Aborting...\n");
  }
  else if (err > 0) {
    return fmt::format("ExoII_Read::Load_Global_Results(): WARNING:  "
                       "Exodus issued warning \"{}\" on call to ex_get_glob_vars()!",
                       err);
  }
  return "";
}

template <typename INT>
std::string ExoII_Read<INT>::Load_Global_Results(int t1, int t2, double proportion)
{
  SMART_ASSERT(Check_State());
  SMART_ASSERT(t1 > 0 && t1 <= num_times);
  SMART_ASSERT(t2 > 0 && t2 <= num_times);

  if (!Open()) {
    return "WARNING:  File not open!";
  }
  if (global_vars.empty()) {
    return "WARNING:  No global variables! (doing nothing)";
  }
  if (global_vals == nullptr) {
    global_vals = new double[global_vars.size()];
    SMART_ASSERT(global_vals != nullptr);
  }

  if (t2 != t1 && (global_vals2 == nullptr)) {
    global_vals2 = new double[global_vars.size()];
    SMART_ASSERT(global_vals2 != nullptr);
  }

  for (unsigned j = 0; j < global_vars.size(); ++j) {
    global_vals[j] = 0.0;
  }

  int err = ex_get_var(file_id, t1, EX_GLOBAL, 1, 1, global_vars.size(), global_vals);

  if (err < 0) {
    Error("ExoII_Read::Load_Global_Results(): Failed to get "
          "global variable values!  Aborting...\n");
  }

  if (t2 != t1) {
    err = ex_get_var(file_id, t2, EX_GLOBAL, 1, 1, global_vars.size(), global_vals2);
    if (err < 0) {
      Error("ExoII_Read::Load_Global_Results(): Failed to get "
            "global variable values!  Aborting...\n");
    }

    // Do the interpolation...
    for (size_t j = 0; j < global_vars.size(); j++) {
      global_vals[j] = (1.0 - proportion) * global_vals[j] + proportion * global_vals2[j];
    }
  }
  return "";
}

template <typename INT>
Side_Set<INT> *ExoII_Read<INT>::Get_Side_Set_by_Index(size_t side_set_index) const
{
  SMART_ASSERT(Check_State());

  if (side_set_index >= num_side_sets) {
    return nullptr;
  }

  return &ssets[side_set_index];
}

template <typename INT>
Node_Set<INT> *ExoII_Read<INT>::Get_Node_Set_by_Index(size_t set_index) const
{
  SMART_ASSERT(Check_State());

  if (set_index >= num_node_sets) {
    return nullptr;
  }

  return &nsets[set_index];
}

template <typename INT>
Edge_Block<INT> *ExoII_Read<INT>::Get_Edge_Block_by_Index(size_t edge_block_index) const
{
  SMART_ASSERT(Check_State());

  if (edge_block_index >= num_edge_blocks) {
    return nullptr;
  }

  return &edge_blocks[edge_block_index];
}

template <typename INT>
Face_Block<INT> *ExoII_Read<INT>::Get_Face_Block_by_Index(size_t face_block_index) const
{
  SMART_ASSERT(Check_State());

  if (face_block_index >= num_face_blocks) {
    return nullptr;
  }

  return &face_blocks[face_block_index];
}

// **********************  Misc functions  *************************** //

// This function converts an Exodus global element number (1-offset) into
// its block index (0-offset) and block element index (0-offset).
template <typename INT>
std::pair<int, size_t> ExoII_Read<INT>::Global_to_Block_Local(size_t global_elmt_num) const
{
  SMART_ASSERT(Check_State());

  if (!Open()) {
    Error("exodiff: ERROR:  File not open!");
  }
  if (global_elmt_num < 1 || global_elmt_num > num_elmts) {
    Error(fmt::format("exodiff: ERROR:  global_elmt_num = {} is out of bounds [1, {}]!",
                      fmt::group_digits(global_elmt_num), fmt::group_digits(num_elmts)));
  }

  int block_index = 0;

  size_t total = 0;
  while (total + eblocks[block_index].Size() < global_elmt_num) {
    total += eblocks[block_index++].Size();
  }

  return std::make_pair(block_index, global_elmt_num - total - 1);
}

template <typename INT> int ExoII_Read<INT>::Check_State() const
{
  SMART_ASSERT(file_id >= -1);
  SMART_ASSERT(db_version >= 0.0);
  SMART_ASSERT(api_version >= 0.0);
  SMART_ASSERT(io_word_size == 0 || io_word_size == 4 || io_word_size == 8);

  SMART_ASSERT(!(file_id >= 0 && io_word_size == 0));
  SMART_ASSERT(!(file_id >= 0 && file_name.empty()));

  SMART_ASSERT(!(num_elmt_blocks > 0 && !eblocks));
  SMART_ASSERT(!(num_node_sets > 0 && !nsets));
  SMART_ASSERT(!(num_side_sets > 0 && !ssets));

  SMART_ASSERT(!(num_nodes == 0 && nodes));

  SMART_ASSERT(num_times >= 0);
  SMART_ASSERT(!(num_times > 0 && !times));

  SMART_ASSERT(cur_time >= 0 && cur_time <= num_times);
  SMART_ASSERT(!(!nodal_vars.empty() && !results));
  SMART_ASSERT(!(nodal_vars.empty() && results));

  return 1;
}

template <typename INT> std::string ExoII_Read<INT>::File_Name(const char *fname)
{
  SMART_ASSERT(Check_State());

  if (Open()) {
    return "exodiff: ERROR: File is already open!";
  }
  if ((fname == nullptr) || std::strlen(fname) == 0) {
    return "exodiff: ERROR: File name is empty!";
  }
  file_name = fname;

  return "";
}

template <typename INT> std::string ExoII_Read<INT>::Open_File(const char *fname)
{
  SMART_ASSERT(Check_State());

  if (Open()) {
    return "exodiff: ERROR: File already open!";
  }
  if ((fname != nullptr) && std::strlen(fname) > 0) {
    file_name = fname;
  }
  else if (file_name.empty()) {
    return "No file name to open!";
  }
  int   ws = 0, comp_ws = 8;
  float dumb = 0.0;
  int   mode = EX_READ;
  if (sizeof(INT) == 8) {
    mode |= EX_ALL_INT64_API;
  }
  auto old_opt = ex_opts(EX_VERBOSE);
  int  err     = ex_open(file_name.c_str(), mode, &comp_ws, &ws, &dumb);
  ex_opts(old_opt);
  if (err < 0) {
    std::ostringstream oss;
    fmt::print(oss, "Couldn't open file \"{}\".", file_name);

    // ExodusII library could not open file.  See if a file (exodusII
    // or not) exists with the specified name.
    FILE *fid = fopen(file_name.c_str(), "r");
    if (fid != nullptr) {
      fmt::print(oss, " File exists, but library could not open.");
      fclose(fid);
    }
    else {
      fmt::print(oss, " File does not exist.");
    }
    return oss.str();
  }

  file_id      = err;
  io_word_size = ws;

  Get_Init_Data();

  return "";
}

template <typename INT> void ExoII_Read<INT>::Get_Init_Data()
{
  SMART_ASSERT(Check_State());
  SMART_ASSERT(file_id >= 0);

  // Determine max size of entity and variable names on the database
  int length_name = ex_inquire_int(file_id, EX_INQ_DB_MAX_USED_NAME_LENGTH);
  ex_set_max_name_length(file_id, length_name);

  ex_init_params info{};
  info.title[0] = '\0';

  int err = ex_get_init_ext(file_id, &info);
  if (err < 0) {
    Error(fmt::format("Failed to get init data!"
                      " Error number = {}.  Aborting...\n",
                      err));
  }

  dimension       = info.num_dim;
  num_nodes       = info.num_nodes;
  num_elmts       = info.num_elem;
  num_faces       = info.num_face;
  num_edges       = info.num_edge;
  num_elmt_blocks = info.num_elem_blk;
  num_node_sets   = info.num_node_sets;
  num_side_sets   = info.num_side_sets;
  num_edge_blocks = info.num_edge_blk;
  num_face_blocks = info.num_face_blk;
  num_assemblies  = info.num_assembly;
  title           = info.title;

  if (err > 0 && !interFace.quiet_flag) {
    fmt::print(stderr, "exodiff: WARNING: was issued, number = {}\n", err);
  }
  if (dimension < 1 || dimension > 3) {
    Error(fmt::format("Init data appears corrupt:\n"
                      "         dimension = {}\n"
                      "         num_nodes = {}\n"
                      "         num_elmts = {}\n"
                      "         num_elmt_blocks = {}\n"
                      "         num_node_sets = {}\n"
                      "         num_side_sets = {}\n"
                      "         num_edge_blocks = {}\n"
                      "         num_face_blocks = {}\n"
                      " ... Aborting...\n",
                      dimension, fmt::group_digits(num_nodes), fmt::group_digits(num_elmts),
                      num_elmt_blocks, num_node_sets, num_edge_blocks, num_face_blocks,
                      num_side_sets));
  }

  int num_qa   = ex_inquire_int(file_id, EX_INQ_QA);
  int num_info = ex_inquire_int(file_id, EX_INQ_INFO);

  if (num_qa < 0 || num_info < 0) {
    Error(fmt::format("inquire data appears corrupt:\n"
                      "         num_qa   = {}\n"
                      "         num_info = {}\n"
                      " ... Aborting...\n",
                      num_qa, num_info));
  }

  //                   Coordinate Names...

  char **coords = get_name_array(3, length_name);
  err           = ex_get_coord_names(file_id, coords);
  if (err < 0) {
    Error("Failed to get coordinate names!  Aborting...\n");
  }

  coord_names.clear();
  for (int i = 0; i < dimension; ++i) {
    coord_names.emplace_back(coords[i]);
  }
  free_name_array(coords, 3);

  // Assembly Data...
  delete[] assemblies;
  assemblies = nullptr;
  if (num_assemblies > 0) {
    assemblies = new Assembly<INT>[num_assemblies];
    SMART_ASSERT(assemblies != nullptr);
    std::vector<INT> ids(num_assemblies);

    err = ex_get_ids(file_id, EX_ASSEMBLY, Data(ids));

    if (err < 0) {
      Error("Failed to get assembly ids!  Aborting...\n");
    }

    for (size_t b = 0; b < num_assemblies; ++b) {
      if (ids[b] <= EX_INVALID_ID) {
        fmt::print(stderr,
                   "EXODIFF  WARNING:  Assembly Id "
                   "for assembly index {} is {} which is negative. This was returned by call to "
                   "ex_get_ids().\n",
                   b, ids[b]);
      }

      assemblies[b].initialize(file_id, ids[b]);
    }
  }
  //                 Element Block Data...

  delete[] eblocks;
  eblocks = nullptr;
  if (num_elmt_blocks > 0) {
    eblocks = new Exo_Block<INT>[num_elmt_blocks];
    SMART_ASSERT(eblocks != nullptr);
    std::vector<INT> ids(num_elmt_blocks);

    err = ex_get_ids(file_id, EX_ELEM_BLOCK, Data(ids));

    if (err < 0) {
      Error("Failed to get element block ids!  Aborting...\n");
    }

    size_t e_count = 0;
    for (size_t b = 0; b < num_elmt_blocks; ++b) {
      if (ids[b] <= EX_INVALID_ID) {
        fmt::print(stderr,
                   "EXODIFF  WARNING:  Element block Id "
                   "for block index {} is {} which is negative. This was returned by call to "
                   "ex_get_elem_blk_ids().\n",
                   b, ids[b]);
      }

      eblocks[b].initialize(file_id, ids[b]);
      eblocks[b].offset(e_count);
      e_count += eblocks[b].Size();
    }

    if (e_count != num_elmts && !interFace.quiet_flag) {
      fmt::print(stderr,
                 "exodiff: WARNING: Total number of elements {}"
                 " does not equal the sum of the number of elements "
                 "in each block {}\n",
                 fmt::group_digits(num_elmts), fmt::group_digits(e_count));
    }

    // Gather the attribute names (even though not all attributes are on all blocks)
    std::set<std::string> names;
    for (size_t b = 0; b < num_elmt_blocks; ++b) {
      for (int a = 0; a < eblocks[b].attr_count(); a++) {
        names.insert(eblocks[b].Get_Attribute_Name(a));
      }
    }
    elmt_atts.resize(names.size());
    std::copy(names.begin(), names.end(), elmt_atts.begin());
  }

  //                     Node & Side sets...

  delete[] nsets;
  nsets = nullptr;
  if (num_node_sets > 0) {
    nsets = new Node_Set<INT>[num_node_sets];
    SMART_ASSERT(nsets != nullptr);
    std::vector<INT> ids(num_node_sets);

    err = ex_get_ids(file_id, EX_NODE_SET, Data(ids));

    if (err < 0) {
      Error("Failed to get nodeset ids!  Aborting...\n");
    }

    for (size_t nset = 0; nset < num_node_sets; ++nset) {
      if (ids[nset] <= EX_INVALID_ID) {
        fmt::print(stderr,
                   "EXODIFF  WARNING: Nodeset Id "
                   "for nodeset index {} is {}"
                   " which is negative.  This was returned by call to ex_get_ids().\n",
                   nset, ids[nset]);
      }

      nsets[nset].initialize(file_id, ids[nset]);
    }
  }

  delete[] ssets;
  ssets = nullptr;
  if (num_side_sets) {
    ssets = new Side_Set<INT>[num_side_sets];
    SMART_ASSERT(ssets != nullptr);
    std::vector<INT> ids(num_side_sets);

    err = ex_get_ids(file_id, EX_SIDE_SET, Data(ids));

    if (err < 0) {
      Error("Failed to get sideset ids!  Aborting...\n");
    }

    for (size_t sset = 0; sset < num_side_sets; ++sset) {
      if (ids[sset] <= EX_INVALID_ID) {
        fmt::print(stderr,
                   "EXODIFF  WARNING: Sideset Id for sideset index {} is {}"
                   " which is negative.  This was returned by call to ex_get_ids().\n",
                   sset, ids[sset]);
      }
      ssets[sset].initialize(file_id, ids[sset]);
    }
  }

  //                     Edge & Face blocks...

  delete[] edge_blocks;
  edge_blocks = nullptr;
  if (num_edge_blocks > 0) {
    edge_blocks = new Edge_Block<INT>[num_edge_blocks];
    SMART_ASSERT(edge_blocks != nullptr);
    std::vector<INT> ids(num_edge_blocks);

    err = ex_get_ids(file_id, EX_EDGE_BLOCK, Data(ids));

    if (err < 0) {
      Error("Failed to get edgeblock ids!  Aborting...\n");
    }

    for (size_t edge_block = 0; edge_block < num_edge_blocks; ++edge_block) {
      if (ids[edge_block] <= EX_INVALID_ID) {
        fmt::print(stderr,
                   "EXODIFF  WARNING: Edgeblock Id "
                   "for edgeblock index {} is {}"
                   " which is negative.  This was returned by call to ex_get_ids().\n",
                   edge_block, ids[edge_block]);
      }

      edge_blocks[edge_block].initialize(file_id, ids[edge_block]);
    }
  }

  delete[] face_blocks;
  face_blocks = nullptr;
  if (num_face_blocks > 0) {
    face_blocks = new Face_Block<INT>[num_face_blocks];
    SMART_ASSERT(face_blocks != nullptr);
    std::vector<INT> ids(num_face_blocks);

    err = ex_get_ids(file_id, EX_FACE_BLOCK, Data(ids));

    if (err < 0) {
      Error("Failed to get faceblock ids!  Aborting...\n");
    }

    for (size_t face_block = 0; face_block < num_face_blocks; ++face_block) {
      if (ids[face_block] <= EX_INVALID_ID) {
        fmt::print(stderr,
                   "EXODIFF  WARNING: Faceblock Id "
                   "for faceblock index {} is {}"
                   " which is negative.  This was returned by call to ex_get_ids().\n",
                   face_block, ids[face_block]);
      }

      face_blocks[face_block].initialize(file_id, ids[face_block]);
    }
  }

  //  **************  RESULTS info  ***************  //

  int num_global_vars, num_nodal_vars, num_elmt_vars, num_ns_vars, num_ss_vars, num_edge_vars,
      num_face_vars;

  err = ex_get_variable_param(file_id, EX_GLOBAL, &num_global_vars);
  if (err < 0) {
    Error("Failed to get number of global variables!  Aborting...\n");
  }

  err = ex_get_variable_param(file_id, EX_NODAL, &num_nodal_vars);
  if (err < 0) {
    Error("Failed to get number of nodal variables!  Aborting...\n");
  }

  err = ex_get_variable_param(file_id, EX_ELEM_BLOCK, &num_elmt_vars);
  if (err < 0) {
    Error("Failed to get number of element variables!  Aborting...\n");
  }

  err = ex_get_variable_param(file_id, EX_NODE_SET, &num_ns_vars);
  if (err < 0) {
    Error("Failed to get number of nodeset variables!  Aborting...\n");
  }

  err = ex_get_variable_param(file_id, EX_SIDE_SET, &num_ss_vars);
  if (err < 0) {
    Error("Failed to get number of sideset variables!  Aborting...\n");
  }

  err = ex_get_variable_param(file_id, EX_EDGE_BLOCK, &num_edge_vars);
  if (err < 0) {
    Error("Failed to get number of edgeblock variables!  Aborting...\n");
  }

  err = ex_get_variable_param(file_id, EX_FACE_BLOCK, &num_face_vars);
  if (err < 0) {
    Error("Failed to get number of faceblock variables!  Aborting...\n");
  }

  if (num_global_vars < 0 || num_nodal_vars < 0 || num_elmt_vars < 0 || num_ns_vars < 0 ||
      num_ss_vars < 0 || num_edge_vars < 0 || num_face_vars < 0) {
    Error(fmt::format("Data appears corrupt for"
                      " number of variables !\n"
                      "\tnum global vars  = {}\n"
                      "\tnum nodal vars   = {}\n"
                      "\tnum element vars = {}\n"
                      "\tnum nodeset vars = {}\n"
                      "\tnum sideset vars = {}\n"
                      "\tnum edgeblock vars = {}\n"
                      "\tnum faceblock vars = {}\n"
                      " ... Aborting...\n",
                      num_global_vars, num_nodal_vars, num_elmt_vars, num_ns_vars, num_ss_vars,
                      num_edge_vars, num_face_vars));
  }

  read_vars(file_id, EX_GLOBAL, "Global", num_global_vars, global_vars);
  read_vars(file_id, EX_NODAL, "Nodal", num_nodal_vars, nodal_vars);
  read_vars(file_id, EX_ELEM_BLOCK, "Element", num_elmt_vars, elmt_vars);
  read_vars(file_id, EX_NODE_SET, "Nodeset", num_ns_vars, ns_vars);
  read_vars(file_id, EX_SIDE_SET, "Sideset", num_ss_vars, ss_vars);
  read_vars(file_id, EX_EDGE_BLOCK, "Edgeblock", num_edge_vars, eb_vars);
  read_vars(file_id, EX_FACE_BLOCK, "Faceblock", num_face_vars, fb_vars);

  // Times:
  num_times = ex_inquire_int(file_id, EX_INQ_TIME);
  if (num_times < 0) {
    Error(fmt::format("Number of time steps came back negative ({})!  Aborting...\n", num_times));
  }

  if ((num_global_vars > 0 || num_nodal_vars > 0 || num_elmt_vars > 0 || num_ns_vars > 0 ||
       num_ss_vars > 0 || num_edge_vars > 0 || num_face_vars > 0) &&
      num_times == 0) {
    Warning("Consistency error -- The database contains transient variables, but no "
            "timesteps!\n");
  }

  if (num_times) {
    times = new double[num_times];
    SMART_ASSERT(times != nullptr);
    ex_get_all_times(file_id, times);
    if (time_scale != 1.0 || time_offset != 0.0) {
      for (int i = 0; i < num_times; i++) {
        times[i] = time_scale * times[i] + time_offset;
      }
    }
  }

  if (num_nodal_vars != 0) {
    if (num_times == 0) {
      Warning(fmt::format("Consistency error--The database contains {}"
                          " nodal variables, but there are no time steps defined.\n",
                          num_nodal_vars));
    }
    if (num_times) {
      results = new double *[num_nodal_vars];
      for (int i = 0; i < num_nodal_vars; ++i) {
        results[i] = nullptr;
      }
    }
  }

} // End of EXODIFF

namespace {
  void read_vars(int file_id, EXOTYPE flag, const char *type, int num_vars,
                 std::vector<std::string> &varlist)
  {
    if (num_vars != 0) {
      int    name_size = ex_inquire_int(file_id, EX_INQ_MAX_READ_NAME_LENGTH);
      char **varnames  = get_name_array(num_vars, name_size);
      int    err       = ex_get_variable_names(file_id, flag, num_vars, varnames);

      if (err < 0) {
        Error(fmt::format("Failed to get {} variable names!  Aborting...\n", type));
      }
      else if (err > 0 && !interFace.quiet_flag) {
        fmt::print(
            stderr,
            "exodiff: WARNING: Exodus issued warning \"{}\" on call to ex_get_var_names()!\n", err);
      }
      for (int vg = 0; vg < num_vars; ++vg) {
        SMART_ASSERT(varnames[vg] != nullptr);
        if (std::strlen(varnames[vg]) == 0 ||
            static_cast<int>(std::strlen(varnames[vg])) > name_size) {
          std::ostringstream out;
          fmt::print(out,
                     "exodiff: ERROR: {} variable names appear corrupt\n"
                     "                A length is 0 or greater than name_size({})\n"
                     "                Here are the names that I received from"
                     " a call to ex_get_var_names(...):\n",
                     type, name_size);
          for (int k = 1; k <= num_vars; ++k) {
            fmt::print(out, "\t\t{}) \"{}\"\n", k, varnames[k - 1]);
          }
          fmt::print(out, "                 Aborting...\n");
          Error(out);
        }

        std::string n(varnames[vg]);
        chop_whitespace(n);
        varlist.push_back(n);
      }
      free_name_array(varnames, num_vars);
    }
  }
} // namespace
template class ExoII_Read<int>;
template class ExoII_Read<int64_t>;
