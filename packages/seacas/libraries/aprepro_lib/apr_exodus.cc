// Copyright(C) 1999-2022, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#if defined(EXODUS_SUPPORT)
#include "aprepro.h"
#include "exodusII.h"
#include "exodusII_int.h"

#include "apr_symrec.h"
#include "apr_util.h"
#include "aprepro_parser.h"

#include "fmt/format.h"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>

namespace {
  std::string LowerCase(std::string name);

  bool matches_prefix(const char *pre, const char *str)
  {
    return strncmp(pre, str, strlen(pre)) == 0;
  }

  std::string entity_type_name(ex_entity_type ent_type)
  {
    switch (ent_type) {
    case EX_ASSEMBLY: return "assembly_";
    case EX_ELEM_BLOCK: return "block_";
    case EX_NODE_SET: return "nodeset_";
    case EX_SIDE_SET: return "sideset_";
    default: return "invalid_";
    }
  }

  void add_name(int exoid, ex_entity_type ent_type, int64_t id, std::vector<char> &name,
                std::string &names)
  {
    std::string str_name;
    ex_get_name(exoid, ent_type, id, name.data());
    if (name[0] == '\0') {
      str_name = entity_type_name(ent_type) + std::to_string(id);
    }
    else {
      str_name = name.data();
    }

    if (!names.empty()) {
      names += ",";
    }
    names += str_name;
  }

  void get_change_set_names(int exoid, SEAMS::Aprepro *aprepro)
  {
    int   idum;
    float rdum;

    // Get root of database...
    int rootid = exoid & EX_FILE_ID_MASK;

    std::string cs_names;

    int               group_name_length = ex_inquire_int(rootid, EX_INQ_GROUP_NAME_LEN);
    std::vector<char> group_name(group_name_length + 1, '\0');

    int num_children = ex_inquire_int(rootid, EX_INQ_NUM_CHILD_GROUPS);
    for (int i = 0; i < num_children; i++) {
      ex_inquire(rootid + 1 + i, EX_INQ_GROUP_NAME, &idum, &rdum, group_name.data());
      if (i > 0) {
        cs_names += ",";
      }
      cs_names += group_name.data();
    }
    aprepro->add_variable("ex_change_set_names", cs_names);
  }
} // namespace

namespace SEAMS {
  extern SEAMS::Aprepro *aprepro;

  int open_exodus_file(char *filename, int cs_idx, bool do_all_cs_defines)
  {
    int   cpu = sizeof(double);
    int   io  = 0;
    float version;

    int exo = ex_open(filename, EX_READ | EX_ALL_INT64_API, &cpu, &io, &version);
    if (exo < 0) {
      // If there is an include path specified, try opening file there
      std::string file_path(aprepro->ap_options.include_path);
      if (!file_path.empty()) {
        file_path += "/";
        file_path += filename;
        exo = ex_open(file_path.c_str(), EX_READ | EX_ALL_INT64_API, &cpu, &io, &version);
      }
      if (exo < 0) {
        yyerror(*aprepro, "Error opening exodusII file.");
      }
    }

    // See if the file contains change sets.  If it does, open the first one.
    int num_change_sets   = ex_inquire_int(exo, EX_INQ_NUM_CHILD_GROUPS);
    int active_change_set = 0;
    if (num_change_sets >= 1) {
      if (do_all_cs_defines) {
        if (cs_idx == 0) {
          cs_idx = 1;
          aprepro->warning(
              fmt::format("Input database contains {} change sets. Rreading from change set {}.",
                          num_change_sets, cs_idx));
        }
      }
      if (cs_idx <= num_change_sets) {
        active_change_set = cs_idx;
        exo += cs_idx;
        get_change_set_names(exo, aprepro);
      }
      else {
        yyerror(*aprepro, fmt::format("Specified change set index {} exceeds count {}", cs_idx,
                                      num_change_sets));
        return -1;
      }
    }
    else {
      if (cs_idx > 0) {
        aprepro->warning(
            fmt::format("Input database does not contain change sets, but a change set "
                        "index {} was specified.  Ignoring.",
                        cs_idx));
      }
    }

    aprepro->add_variable("ex_change_set_count", num_change_sets);
    if (do_all_cs_defines) {
      aprepro->add_variable("ex_active_change_set", active_change_set);
      aprepro->add_variable("ex_version", version);
    }
    return exo;
  }

  const char *do_exodus_query_change_sets(char *filename)
  {
    // This will open the file and determine whether there are any change sets on the file.
    // It will define the `ex_change_set_count`, and `ex_change_set_names` variables.

    // Open the specified exodusII file, read the info records
    // then parse them as input to aprepro.
    int exoid = open_exodus_file(filename, 0, false);
    if (exoid > 0) {
      ex_close(exoid);
    }
    return "";
  }

  const char *do_exodus_info(char *filename, char *prefix)
  {
    char *ret_string = nullptr;

    // Open the specified exodusII file, read the info records
    // then parse them as input to aprepro.
    int exoid = open_exodus_file(filename, 0, true);
    if (exoid <= 0) {
      return "";
    }

    int count = ex_inquire_int(exoid, EX_INQ_INFO);

    if (count > 0) {
      auto info = new char *[count];
      for (int i = 0; i < count; i++) {
        info[i] = new char[MAX_LINE_LENGTH + 1];
        memset(info[i], '\0', MAX_LINE_LENGTH + 1);
      }

      ex_get_info(exoid, info);

      std::string lines;
      size_t      prefix_len = strlen(prefix);
      for (int i = 0; i < count; i++) {
        if (matches_prefix(prefix, info[i])) {
          lines += std::string(info[i]).substr(prefix_len);
          lines += "\n";
        }
      }

      new_string(lines, &ret_string);

      for (int i = 0; i < count; i++) {
        delete[] info[i];
      }
      delete[] info;

      ex_close(exoid);
      return ret_string;
    }

    ex_close(exoid);
    return "";
  }

  const char *do_exodus_info_range(char *filename, char *beg, char *end)
  {
    char *ret_string = nullptr;

    // Open the specified exodusII file, read the info records
    // then parse them as input to aprepro.
    int exoid = open_exodus_file(filename, 0, true);
    if (exoid <= 0) {
      return "";
    }

    int count = ex_inquire_int(exoid, EX_INQ_INFO);

    if (count > 0) {
      auto info = new char *[count];
      for (int i = 0; i < count; i++) {
        info[i] = new char[MAX_LINE_LENGTH + 1];
        memset(info[i], '\0', MAX_LINE_LENGTH + 1);
      }

      ex_get_info(exoid, info);

      bool        in_range = false;
      std::string lines;
      for (int i = 0; i < count; i++) {
        if (in_range && strcmp(info[i], end) == 0) {
          in_range = false;
        }
        if (in_range) {
          lines += std::string(info[i]);
          lines += "\n";
        }
        if (!in_range && strcmp(info[i], beg) == 0) {
          in_range = true;
        }
      }

      new_string(lines, &ret_string);

      for (int i = 0; i < count; i++) {
        delete[] info[i];
      }
      delete[] info;

      ex_close(exoid);
      return ret_string;
    }

    ex_close(exoid);
    return "";
  }

  const char *do_exodus_meta_cd(char *filename, double cs_index)
  {
    // Open the specified exodusII file, read the metadata and set
    // variables for each item.
    // Examples include "node_count", "element_count", ...
    int exoid = open_exodus_file(filename, (int)cs_index, true);
    if (exoid <= 0) {
      return "";
    }

    // read database parameters
    ex_init_params info;
    ex_get_init_ext(exoid, &info);

    aprepro->add_variable("ex_title", info.title);
    aprepro->add_variable("ex_dimension", info.num_dim);
    aprepro->add_variable("ex_node_count", info.num_nodes);
    aprepro->add_variable("ex_element_count", info.num_elem);
    aprepro->add_variable("ex_block_count", info.num_elem_blk);
    aprepro->add_variable("ex_assembly_count", info.num_assembly);
    aprepro->add_variable("ex_nodeset_count", info.num_node_sets);
    aprepro->add_variable("ex_sideset_count", info.num_side_sets);

    { // Nemesis Information
      int  proc_count;
      int  proc_in_file;
      char file_type[MAX_STR_LENGTH + 1];

      ex_get_init_info(exoid, &proc_count, &proc_in_file, file_type);

      if (proc_count > 1) {
        int64_t global_nodes;
        int64_t global_elements;
        int64_t global_blocks;
        int64_t global_nsets;
        int64_t global_ssets;

        aprepro->add_variable("ex_processor_count", proc_count);

        ex_get_init_global(exoid, &global_nodes, &global_elements, &global_blocks, &global_nsets,
                           &global_ssets);

        aprepro->add_variable("ex_node_count_global", global_nodes);
        aprepro->add_variable("ex_element_count_global", global_elements);
      }
    }

    // Read The Element Blocks, Node Sets, and Side Sets and set variables for each of these.
    // The Scheme Is:
    // -- 'ex_block_ids' Is an array of the element block ids. (ex_block_count, 1)
    // -- 'ex_block_info' is an array of the element block info (id, num_elem, num_node_per_element,
    // num_attrib) for each block (ex_block_count,4)
    // -- 'ex_nodeset_ids'
    // -- 'ex_nodeset_info'
    // -- 'ex_sideset_ids'
    // -- 'ex_sideset_info'
    int max_name_length = ex_inquire_int(exoid, EX_INQ_DB_MAX_USED_NAME_LENGTH);
    ex_set_max_name_length(exoid, max_name_length);
    std::vector<char> name(max_name_length + 1);
    std::string       str_name;

    if (info.num_elem_blk > 0) {
      auto array_data       = aprepro->make_array(info.num_elem_blk, 1);
      auto array_block_info = aprepro->make_array(info.num_elem_blk, 4);

      std::vector<int64_t> ids(info.num_elem_blk);
      ex_get_ids(exoid, EX_ELEM_BLOCK, ids.data());

      char    type[MAX_STR_LENGTH + 1];
      int64_t nel;
      int64_t nnel;
      int64_t natr;

      std::string names;
      std::string topology;

      int64_t idx = 0;
      for (int64_t i = 0; i < info.num_elem_blk; i++) {
        ex_get_block(exoid, EX_ELEM_BLOCK, ids[i], type, &nel, &nnel, nullptr, nullptr, &natr);
        array_data->data[i]           = ids[i];
        array_block_info->data[idx++] = ids[i];
        array_block_info->data[idx++] = nel;
        array_block_info->data[idx++] = nnel;
        array_block_info->data[idx++] = natr;

        if (i > 0) {
          topology += ",";
        }
        topology += type;
        add_name(exoid, EX_ELEM_BLOCK, ids[i], name, names);
      }

      topology = LowerCase(topology);
      aprepro->add_variable("ex_block_topology", topology);
      aprepro->add_variable("ex_block_names", names);
      aprepro->add_variable("ex_block_ids", array_data);
      aprepro->add_variable("ex_block_info", array_block_info);
    }

    if (info.num_assembly > 0) {
      std::vector<int64_t> ids(info.num_assembly);
      ex_get_ids(exoid, EX_ASSEMBLY, ids.data());

      std::string names;
      std::string type;
      auto        array_data = aprepro->make_array(info.num_assembly, 1);
      auto        array_info = aprepro->make_array(info.num_assembly, 1);

      for (int64_t i = 0; i < info.num_assembly; i++) {
        ex_assembly assembly;
        assembly.id          = ids[i];
        assembly.name        = name.data();
        assembly.entity_list = nullptr;

        ex_get_assembly(exoid, &assembly);
        if (i > 0) {
          names += ",";
          type += ",";
        }
        array_data->data[i] = ids[i];
        array_info->data[i] = assembly.entity_count;
        names += assembly.name;
        type += ex_name_of_object(assembly.type);
      }
      aprepro->add_variable("ex_assembly_type", type);
      aprepro->add_variable("ex_assembly_names", names);
      aprepro->add_variable("ex_assembly_ids", array_data);
      aprepro->add_variable("ex_assembly_info", array_info);
    }

    // Nodesets...
    if (info.num_node_sets > 0) {
      auto array_data     = aprepro->make_array(info.num_node_sets, 1);
      auto array_set_info = aprepro->make_array(info.num_node_sets, 3);

      std::vector<int64_t> ids(info.num_node_sets);
      ex_get_ids(exoid, EX_NODE_SET, ids.data());

      std::string names;
      int64_t     idx = 0;
      for (int64_t i = 0; i < info.num_node_sets; i++) {
        int64_t num_entry;
        int64_t num_dist;
        ex_get_set_param(exoid, EX_NODE_SET, ids[i], &num_entry, &num_dist);
        array_data->data[i]         = ids[i];
        array_set_info->data[idx++] = ids[i];
        array_set_info->data[idx++] = num_entry;
        array_set_info->data[idx++] = num_dist;

        add_name(exoid, EX_NODE_SET, ids[i], name, names);
      }

      aprepro->add_variable("ex_nodeset_names", names);
      aprepro->add_variable("ex_nodeset_ids", array_data);
      aprepro->add_variable("ex_nodeset_info", array_set_info);
    }

    // Sidesets...
    if (info.num_side_sets > 0) {
      auto array_data     = aprepro->make_array(info.num_side_sets, 1);
      auto array_set_info = aprepro->make_array(info.num_side_sets, 3);

      std::vector<int64_t> ids(info.num_side_sets);
      ex_get_ids(exoid, EX_SIDE_SET, ids.data());

      std::string names;
      int64_t     idx = 0;
      for (int64_t i = 0; i < info.num_side_sets; i++) {
        int64_t num_entry;
        int64_t num_dist;
        ex_get_set_param(exoid, EX_SIDE_SET, ids[i], &num_entry, &num_dist);
        array_data->data[i]         = ids[i];
        array_set_info->data[idx++] = ids[i];
        array_set_info->data[idx++] = num_entry;
        array_set_info->data[idx++] = num_dist;

        add_name(exoid, EX_SIDE_SET, ids[i], name, names);
      }

      aprepro->add_variable("ex_sideset_names", names);
      aprepro->add_variable("ex_sideset_ids", array_data);
      aprepro->add_variable("ex_sideset_info", array_set_info);
    }

    // Get timestep count
    int64_t ts_count = ex_inquire_int(exoid, EX_INQ_TIME);
    aprepro->add_variable("ex_timestep_count", ts_count);

    if (ts_count > 0) {
      std::vector<double> timesteps(ts_count);
      ex_get_all_times(exoid, timesteps.data());

      auto ts_array_data = aprepro->make_array(ts_count, 1);
      for (int64_t i = 0; i < ts_count; i++) {
        ts_array_data->data[i] = timesteps[i];
      }
      aprepro->add_variable("ex_timestep_times", ts_array_data);
    }

    // See if any global variables on file...
    int num_global = 0;
    ex_get_variable_param(exoid, EX_GLOBAL, &num_global);
    if (num_global > 0) {
      std::string names;
      for (int i = 0; i < num_global; i++) {
        ex_get_variable_name(exoid, EX_GLOBAL, i + 1, name.data());
        if (i > 0) {
          names += ",";
        }
        names += name.data();
      }
      aprepro->add_variable("ex_global_var_names", names);

      auto                glo_array_data = aprepro->make_array(ts_count, num_global);
      std::vector<double> globals(num_global);
      int                 index = 0;
      for (int64_t i = 0; i < ts_count; i++) {
        ex_get_var(exoid, i + 1, EX_GLOBAL, 0, 0, num_global, globals.data());
        for (int j = 0; j < num_global; j++) {
          glo_array_data->data[index++] = globals[j];
        }
      }
      aprepro->add_variable("ex_global_var_value", glo_array_data);
    }

    ex_close(exoid);
    return "";
  }

  const char *do_exodus_meta(char *filename) { return do_exodus_meta_cd(filename, 0); }

} // namespace SEAMS

namespace {
  inline int  to_lower(int c) { return std::tolower(c); }
  std::string LowerCase(std::string name)
  {
    std::transform(name.begin(), name.end(), name.begin(), to_lower);
    return name;
  }
} // namespace
#endif
