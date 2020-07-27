// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "ED_SystemInterface.h" // for SystemInterface, interFace
#include "Tolerance.h"          // for Tolerance, etc
#include "exo_entity.h"         // for Exo_Entity, EXOTYPE
#include "exodusII.h"
#include "fmt/color.h"
#include "fmt/ostream.h"
#include "smart_assert.h" // for SMART_ASSERT
#include "stringx.h"      // for find_string, etc
#include "util.h"
#include <cstddef> // for size_t
#include <cstdio>  // for nullptr
#include <string>  // for string, char_traits, etc
#include <vector>  // for vector
template <typename INT> class ExoII_Read;

namespace {
  void build_variable_names(const char *type, std::vector<std::string> &names,
                            std::vector<Tolerance> &tols, const Tolerance &default_tol,
                            bool do_all_flag, const std::vector<std::string> &var_names1,
                            const std::vector<std::string> &var_names2, bool *diff_found);

  template <typename INT>
  void build_truth_table(EXOTYPE type, const char *label, std::vector<std::string> &names,
                         size_t num_entity, ExoII_Read<INT> &file1, ExoII_Read<INT> &file2,
                         const std::vector<std::string> &var_names1,
                         const std::vector<std::string> &var_names2, std::vector<int> &truth_tab,
                         bool quiet_flag, bool *diff_found);

  void output_exodus_names(int file_id, EXOTYPE type, const std::vector<std::string> &names);
  void output_diff_names(const char *type, const std::vector<std::string> &names);
  void output_compare_names(const char *type, const std::vector<std::string> &names,
                            const std::vector<Tolerance> &tol, int num_vars1, int num_vars2);
} // namespace

template <typename INT>
void Build_Variable_Names(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, bool *diff_found)
{
  // Build (and compare) global variable names.
  build_variable_names("global", interFace.glob_var_names, interFace.glob_var,
                       interFace.glob_var_default, interFace.glob_var_do_all_flag,
                       file1.Global_Var_Names(), file2.Global_Var_Names(), diff_found);

  // Build (and compare) nodal variable names.
  build_variable_names("nodal", interFace.node_var_names, interFace.node_var,
                       interFace.node_var_default, interFace.node_var_do_all_flag,
                       file1.Nodal_Var_Names(), file2.Nodal_Var_Names(), diff_found);

  // Build (and compare) element variable names.
  build_variable_names("element", interFace.elmt_var_names, interFace.elmt_var,
                       interFace.elmt_var_default, interFace.elmt_var_do_all_flag,
                       file1.Elmt_Var_Names(), file2.Elmt_Var_Names(), diff_found);

  // Build (and compare) element variable names.
  if (!interFace.ignore_attributes) {
    build_variable_names("element attribute", interFace.elmt_att_names, interFace.elmt_att,
                         interFace.elmt_att_default, interFace.elmt_att_do_all_flag,
                         file1.Elmt_Att_Names(), file2.Elmt_Att_Names(), diff_found);
  }

  // Build (and compare) nodeset variable names.
  build_variable_names("nodeset", interFace.ns_var_names, interFace.ns_var,
                       interFace.ns_var_default, interFace.ns_var_do_all_flag, file1.NS_Var_Names(),
                       file2.NS_Var_Names(), diff_found);

  // Build (and compare) sideset variable names.
  build_variable_names("sideset", interFace.ss_var_names, interFace.ss_var,
                       interFace.ss_var_default, interFace.ss_var_do_all_flag, file1.SS_Var_Names(),
                       file2.SS_Var_Names(), diff_found);
}

template <typename INT>
int Create_File(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, const std::string &diffile_name,
                bool *diff_found)
{
  // Multiple modes:
  // summary_flag == true   --> Single file, output summary and variable names, return
  // diffile_name == ""     --> Dual file, output summary, variable names, check compatibility,
  // diffile_name != ""     --> Three files (2 in, 1 out)
  //                            create output file which is diff of input.
  //                            output summary, variable names, check compatibility
  // quiet_flag == true     --> don't output summary information

  SMART_ASSERT(!interFace.summary_flag);
  //========================================================================
  // From here on down, have two input files and possibly 1 output file...
  // Create output file.

  int out_file_id = -1;
  if (!diffile_name.empty()) {

    // Take minimum word size for output file.
    int iows =
        file1.IO_Word_Size() < file2.IO_Word_Size() ? file1.IO_Word_Size() : file2.IO_Word_Size();
    int compws = sizeof(double);

    int mode = EX_CLOBBER;
    if (sizeof(INT) == 8) {
      mode |= EX_ALL_INT64_DB;
      mode |= EX_ALL_INT64_API;
    }
    out_file_id = ex_create(diffile_name.c_str(), mode, &compws, &iows);
    if (out_file_id < 0) {
      Error(fmt::format("Couldn't create output file \"{}\".\n", diffile_name));
      exit(1);
    }
    ex_copy(file1.File_ID(), out_file_id);
  }

  if (!interFace.quiet_flag) {
    if (out_file_id >= 0) { // The files are to be differenced .. just list names.
      if (interFace.coord_tol.type != IGNORE_) {
        fmt::print("Coordinates:  tol: {:8g} {}, floor: {:8g}\n", interFace.coord_tol.value,
                   interFace.coord_tol.typestr(), interFace.coord_tol.floor);
      }
      else {
        fmt::print("Locations of nodes will not be considered.\n");
      }
      if (interFace.time_tol.type != IGNORE_) {
        fmt::print("Time step values:  tol: {:8g} {}, floor: {:8g}\n", interFace.time_tol.value,
                   interFace.time_tol.typestr(), interFace.time_tol.floor);
      }
      else {
        fmt::print("Time step time values will not be differenced.\n");
      }
      output_diff_names("Global", interFace.glob_var_names);
      output_diff_names("Nodal", interFace.node_var_names);
      output_diff_names("Element", interFace.elmt_var_names);
      output_diff_names("Element Attribute", interFace.elmt_att_names);
      output_diff_names("Nodeset", interFace.ns_var_names);
      output_diff_names("Sideset", interFace.ss_var_names);
    }
    else { // The files are to be compared .. echo additional info.
      if (Tolerance::use_old_floor) {
        std::ostringstream info;
        fmt::print(info, "INFO: Using old definition of floor tolerance. |a-b|<floor.\n\n");
        DIFF_OUT(info, fmt::color::yellow);
      }
      if (interFace.coord_tol.type != IGNORE_) {
        fmt::print("\nNodal coordinates will be compared .. tol: {:8g} ({}), floor: {:8g}\n",
                   interFace.coord_tol.value, interFace.coord_tol.typestr(),
                   interFace.coord_tol.floor);
      }
      else {
        std::ostringstream info;
        fmt::print(info, "\nNodal coordinates will not be compared.\n");
        DIFF_OUT(info, fmt::color::yellow);
      }

      if (interFace.time_tol.type != IGNORE_) {
        fmt::print("Time step values will be compared  .. tol: {:8g} ({}), floor: {:8g}\n",
                   interFace.time_tol.value, interFace.time_tol.typestr(),
                   interFace.time_tol.floor);
      }
      else {
        std::ostringstream info;
        fmt::print(info, "Time step time values will not be compared.\n");
        DIFF_OUT(info, fmt::color::yellow);
      }

      output_compare_names("Global", interFace.glob_var_names, interFace.glob_var,
                           file1.Num_Global_Vars(), file2.Num_Global_Vars());

      output_compare_names("Nodal", interFace.node_var_names, interFace.node_var,
                           file1.Num_Nodal_Vars(), file2.Num_Nodal_Vars());

      output_compare_names("Element", interFace.elmt_var_names, interFace.elmt_var,
                           file1.Num_Elmt_Vars(), file2.Num_Elmt_Vars());

      output_compare_names("Element Attribute", interFace.elmt_att_names, interFace.elmt_att,
                           file1.Num_Elmt_Atts(), file2.Num_Elmt_Atts());

      output_compare_names("Nodeset", interFace.ns_var_names, interFace.ns_var, file1.Num_NS_Vars(),
                           file2.Num_NS_Vars());

      output_compare_names("Sideset", interFace.ss_var_names, interFace.ss_var, file1.Num_SS_Vars(),
                           file2.Num_SS_Vars());
      if (!interFace.ignore_sideset_df && interFace.ss_df_tol.type != IGNORE_ &&
          file1.Num_Side_Sets() > 0 && file2.Num_Side_Sets() > 0) {
        fmt::print(
            "Sideset Distribution Factors will be compared .. tol: {:8g} ({}), floor: {:8g}\n",
            interFace.ss_df_tol.value, interFace.ss_df_tol.typestr(), interFace.ss_df_tol.floor);
      }
      else {
        if (interFace.ignore_sideset_df || interFace.ss_df_tol.type == IGNORE_) {
          std::ostringstream info;
          fmt::print(info, "Sideset Distribution Factors will not be compared.\n");
          DIFF_OUT(info, fmt::color::yellow);
        }
        else {
          fmt::print("No Sideset Distribution Factors on either file.\n");
        }
      }
    }
  }

  std::vector<int> truth_tab;
  build_truth_table(EX_ELEM_BLOCK, "Element Block", interFace.elmt_var_names,
                    file1.Num_Elmt_Blocks(), file1, file2, file1.Elmt_Var_Names(),
                    file2.Elmt_Var_Names(), truth_tab, interFace.quiet_flag, diff_found);

  std::vector<int> ns_truth_tab;
  build_truth_table(EX_NODE_SET, "Nodeset", interFace.ns_var_names, file1.Num_Node_Sets(), file1,
                    file2, file1.NS_Var_Names(), file2.NS_Var_Names(), ns_truth_tab,
                    interFace.quiet_flag, diff_found);

  std::vector<int> ss_truth_tab;
  build_truth_table(EX_SIDE_SET, "Sideset", interFace.ss_var_names, file1.Num_Side_Sets(), file1,
                    file2, file1.SS_Var_Names(), file2.SS_Var_Names(), ss_truth_tab,
                    interFace.quiet_flag, diff_found);

  // Put out the concatenated variable parameters here and then
  // put out the names....
  if (out_file_id >= 0) {
    ex_put_all_var_param(out_file_id, interFace.glob_var_names.size(),
                         interFace.node_var_names.size(), interFace.elmt_var_names.size(),
                         truth_tab.data(), interFace.ns_var_names.size(), ns_truth_tab.data(),
                         interFace.ss_var_names.size(), ss_truth_tab.data());

    output_exodus_names(out_file_id, EX_GLOBAL, interFace.glob_var_names);
    output_exodus_names(out_file_id, EX_NODAL, interFace.node_var_names);
    output_exodus_names(out_file_id, EX_ELEM_BLOCK, interFace.elmt_var_names);
    output_exodus_names(out_file_id, EX_NODE_SET, interFace.ns_var_names);
    output_exodus_names(out_file_id, EX_SIDE_SET, interFace.ss_var_names);
  }
  return out_file_id;
}

namespace {
  void output_exodus_names(int file_id, EXOTYPE type, const std::vector<std::string> &names)
  {
    if (!names.empty()) {
      std::vector<char *> vars(names.size());
      for (unsigned i = 0; i < names.size(); ++i) {
        vars[i] = const_cast<char *>(names[i].c_str());
        SMART_ASSERT(vars[i] != nullptr);
      }
      ex_put_variable_names(file_id, type, names.size(), vars.data());
    }
  }

  void output_compare_names(const char *type, const std::vector<std::string> &names,
                            const std::vector<Tolerance> &tol, int num_vars1, int num_vars2)
  {
    if (!names.empty()) {
      fmt::print("{} variables to be compared:\n", type);
      for (unsigned v = 0; v < names.size(); ++v) {
        if (v == 0) {
          fmt::print("{:<32} tol: {:8g} ({}), floor: {:8g}\n", names[v], tol[v].value,
                     tol[v].typestr(), tol[v].floor);
        }
        else {
          fmt::print("{:<32}      {:8g} ({}),        {:8g}\n", names[v], tol[v].value,
                     tol[v].typestr(), tol[v].floor);
        }
      }
    }
    else if (num_vars1 == 0 && num_vars2 == 0) {
      fmt::print("No {} variables on either file.\n", type);
    }
    else {
      fmt::print("{} variables will not be compared.\n", type);
    }
  }

  void output_diff_names(const char *type, const std::vector<std::string> &names)
  {
    if (!names.empty()) {
      fmt::print("{} variables to be differenced:\n", type);
      for (auto &name : names) {
        fmt::print("\t{}\n", name);
      }
    }
    else {
      fmt::print("No {} variables will be differenced.\n", type);
    }
  }

  void build_variable_names(const char *type, std::vector<std::string> &names,
                            std::vector<Tolerance> &tols, const Tolerance &default_tol,
                            bool do_all_flag, const std::vector<std::string> &var_names1,
                            const std::vector<std::string> &var_names2, bool *diff_found)
  {
    std::vector<std::string> x_list; // exclusion list
    for (auto name : names) {
      chop_whitespace(name);
      SMART_ASSERT(!name.empty());
      if (name[0] == '!') {
        x_list.push_back(extract_token(name, "!")); // remove "!" & add
      }
    }

    if (do_all_flag) {
      int n;
      int name_length = var_names1.size();
      for (n = 0; n < name_length; ++n) {
        const std::string &name = var_names1[n];
        if (!interFace.summary_flag &&
            find_string(var_names2, name, interFace.nocase_var_names) < 0) {
          if (find_string(x_list, name, interFace.nocase_var_names) < 0) {
            if (interFace.allowNameMismatch) {
              x_list.push_back(name);
            }
            else {
              *diff_found = true;
              if (!interFace.quiet_flag) {
                std::ostringstream diff;
                fmt::print(diff,
                           "exodiff: DIFFERENCE .. The {} variable \"{}\" is in the first file but "
                           "not the second.\n",
                           type, name);
                DIFF_OUT(diff);
              }
              continue;
            }
          }
        }
        if (find_string(names, name, interFace.nocase_var_names) < 0 &&
            find_string(x_list, name, interFace.nocase_var_names) < 0) {
          int idx = names.size();
          names.push_back(name);
          tols[idx] = default_tol;
        }
      }

      if (!interFace.noSymmetricNameCheck) {
        name_length = var_names2.size();
        for (n = 0; n < name_length; ++n) {
          const std::string &name = var_names2[n];
          if (!interFace.summary_flag &&
              find_string(var_names1, name, interFace.nocase_var_names) < 0) {
            if (find_string(x_list, name, interFace.nocase_var_names) < 0) {
              *diff_found = true;
              if (!interFace.quiet_flag) {
                std::ostringstream diff;
                fmt::print(diff,
                           "exodiff: DIFFERENCE .. The {} variable \"{}\" is in the second file "
                           "but not the first.\n",
                           type, name);
                DIFF_OUT(diff);
              }
              continue;
            }
          }
          SMART_ASSERT(find_string(names, name, interFace.nocase_var_names) >= 0 ||
                       find_string(x_list, name, interFace.nocase_var_names) >= 0);
        }
      }
    }

    std::vector<std::string> tmp_list;
    for (unsigned n = 0; n < names.size(); ++n) {
      std::string name = names[n];
      chop_whitespace(name);
      if (name[0] == '!') {
        continue;
      }

      int idx = find_string(var_names1, name, interFace.nocase_var_names);
      if (idx >= 0) {
        if (interFace.summary_flag ||
            find_string(var_names2, name, interFace.nocase_var_names) >= 0) {
          tols[tmp_list.size()] = tols[n];
          tmp_list.push_back(var_names1[idx]);
        }
        else {
          *diff_found = true;
          if (!interFace.quiet_flag) {
            std::ostringstream diff;
            fmt::print(diff,
                       "exodiff: DIFFERENCE .. The {} variable \"{}\" is not in the second file.\n",
                       type, name);
            DIFF_OUT(diff);
          }
        }
      }
      else {
        *diff_found = true;
        if (!interFace.quiet_flag) {
          std::ostringstream diff;
          fmt::print(
              diff,
              "exodiff: DIFFERENCE .. Specified {} variable \"{}\" is not in the first file.\n",
              type, name);
          DIFF_OUT(diff);
        }
      }
    }
    names = tmp_list;
  }

  template <typename INT>
  void build_truth_table(EXOTYPE type, const char *label, std::vector<std::string> &names,
                         size_t num_entity, ExoII_Read<INT> &file1, ExoII_Read<INT> &file2,
                         const std::vector<std::string> &var_names1,
                         const std::vector<std::string> &var_names2, std::vector<int> &truth_tab,
                         bool quiet_flag, bool *diff_found)
  {
    if (!names.empty()) {
      int num_vars = names.size();

      truth_tab.resize(num_vars * num_entity);
      for (int i = num_vars * num_entity - 1; i >= 0; --i) {
        truth_tab[i] = 0;
      }

      for (size_t b = 0; b < num_entity; ++b) {
        Exo_Entity *set1 = file1.Get_Entity_by_Index(type, b);
        Exo_Entity *set2 = nullptr;
        if (interFace.by_name) {
          set2 = file2.Get_Entity_by_Name(type, set1->Name());
        }
        else {
          set2 = file2.Get_Entity_by_Id(type, set1->Id());
        }

        if (set2 == nullptr) {
          if (interFace.map_flag != PARTIAL) {
            *diff_found = true;
            std::ostringstream diff;
            fmt::print(diff,
                       "exodiff: DIFFERENCE {} id {} exists in first file but not the second...\n",
                       label, set1->Id());
            DIFF_OUT(diff);
          }
          continue;
        }

        for (int out_idx = 0; out_idx < num_vars; ++out_idx) {
          const std::string &name = names[out_idx];
          int                idx1 = find_string(var_names1, name, interFace.nocase_var_names);
          int                idx2 = find_string(var_names2, name, interFace.nocase_var_names);
          if (idx1 < 0 || idx2 < 0) {
            Error(fmt::format("Unable to find variable named '{}' on database.\n", name));
            exit(1);
          }

          if (set1->is_valid_var(idx1)) {
            if (set2->is_valid_var(idx2)) {
              truth_tab[b * num_vars + out_idx] = 1;
            }
            else if (!quiet_flag) {
              std::ostringstream diff;
              fmt::print(diff,
                         "exodiff: INFO {0} variable \"{1}\" is not saved for {0}"
                         " id {2} in the second file but is in the first (by virtue of the truth "
                         "tables).  "
                         "This variable won't be considered for this {0}.\n",
                         label, name, set1->Id());
              DIFF_OUT(diff, fmt::color::yellow);
            }
          }
          else if (set2->is_valid_var(idx2) && !quiet_flag) {
            std::ostringstream diff;
            fmt::print(
                diff,
                "exodiff: INFO {0} variable \"{1}\" is not saved for {0}"
                " id {2} in the first file but is in the second (by virtue of the truth tables).  "
                "This variable won't be considered for this {0}.\n",
                label, name, set1->Id());
            DIFF_OUT(diff, fmt::color::yellow);
          }
        }
      }
    }
  }
} // End of namespace

template int  Create_File(ExoII_Read<int> &file1, ExoII_Read<int> &file2,
                          const std::string &diffile_name, bool *diff_found);
template void Build_Variable_Names(ExoII_Read<int> &file1, ExoII_Read<int> &file2,
                                   bool *diff_found);

template int  Create_File(ExoII_Read<int64_t> &file1, ExoII_Read<int64_t> &file2,
                          const std::string &diffile_name, bool *diff_found);
template void Build_Variable_Names(ExoII_Read<int64_t> &file1, ExoII_Read<int64_t> &file2,
                                   bool *diff_found);
