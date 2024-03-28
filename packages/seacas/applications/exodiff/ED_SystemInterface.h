// Copyright(C) 1999-2020, 2022, 2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once

#include "GetLongOpt.h" // for GetLongOption
#include "Tolerance.h"  // for Tolerance, etc
#include "util.h"
#include <cmath>
#include <string>  // for string
#include <utility> // for pair
#include <vector>  // for vector

class SystemInterface
{
public:
  SystemInterface();

  bool parse_options(int argc, char **argv);

  static void show_version();

  void Parse_Command_File();
  void Set_Max_Names(int size);

  // Program parameters.
  Tolerance coord_tol{ToleranceMode::ABSOLUTE_, 1.0e-6, 0.0};
  Tolerance time_tol{ToleranceMode::RELATIVE_, 1.0e-6, 1.0e-15};
  Tolerance final_time_tol{ToleranceMode::RELATIVE_, 0.0, 0.0};
  Tolerance default_tol{ToleranceMode::RELATIVE_, 1.0e-6, 0.0};
  Tolerance ss_df_tol{ToleranceMode::RELATIVE_, 1.0e-6, 0.0};

  // These should correspond to the values specified during parsing of
  // coordinate tolerance.
  // Offset of timesteps between first and second databases.
  int time_step_offset{0};
  int time_step_start{1};     // First step to compare (1-based)
  int time_step_stop{-1};     // Last step to compare
  int time_step_increment{1}; // Step increment

  // Offset/Scale time values on input database -- Time_mod = Scale * Time_db + Offset
  double time_value_offset{0.0}; // Add offset to time values on first database
  double time_value_scale{1.0};  // Scale time values on first database

  std::pair<int, int> explicit_steps{}; // Only compare these two steps (db1:db2) if nonzero.

  int max_warnings{100};

  std::vector<std::string> glob_var_names{};
  Tolerance                glob_var_default{ToleranceMode::RELATIVE_, 1.0e-6, 0.0};
  std::vector<Tolerance>   glob_var{};

  std::vector<std::string> node_var_names{};
  Tolerance                node_var_default{ToleranceMode::RELATIVE_, 1.0e-6, 0.0};
  std::vector<Tolerance>   node_var{};

  std::vector<std::string> elmt_var_names{};
  Tolerance                elmt_var_default{ToleranceMode::RELATIVE_, 1.0e-6, 0.0};
  std::vector<Tolerance>   elmt_var{};

  std::vector<std::string> elmt_att_names{};
  Tolerance                elmt_att_default{ToleranceMode::RELATIVE_, 1.0e-6, 0.0};
  std::vector<Tolerance>   elmt_att{};

  std::vector<std::string> ns_var_names{};
  Tolerance                ns_var_default{ToleranceMode::RELATIVE_, 1.0e-6, 0.0};
  std::vector<Tolerance>   ns_var{};

  std::vector<std::string> ss_var_names{};
  Tolerance                ss_var_default{ToleranceMode::RELATIVE_, 1.0e-6, 0.0};
  std::vector<Tolerance>   ss_var{};

  std::vector<std::string> eb_var_names{};
  Tolerance                eb_var_default{ToleranceMode::RELATIVE_, 1.0e-6, 0.0};
  std::vector<Tolerance>   eb_var{};

  std::vector<std::string> fb_var_names{};
  Tolerance                fb_var_default{ToleranceMode::RELATIVE_, 1.0e-6, 0.0};
  std::vector<Tolerance>   fb_var{};

  // time step exclusion data
  std::vector<int> exclude_steps{};

  std::string file1{};
  std::string file2{};
  std::string diff_file{};
  std::string command_file{};

  bool quiet_flag{false};     // By default, warnings and other info is produced
  bool show_all_diffs{false}; // Be default, show only maximum diff for each variable;
                              // if set, show all diff that exceeds tolerance.
  ToleranceMode output_type{
      ToleranceMode::ABSOLUTE_};           // By default, output file diffs are absolute.
  MapType map_flag{MapType::USE_FILE_IDS}; // By default, no searching is done to match
                                           // nodes & elements.
  bool nsmap_flag{true};                   // By default, nodeset nodelist match is off
  bool ssmap_flag{true};                   // By default, sideset elem/side match is off
  bool short_block_check{true};            // By default, element block compares are
                                           // case in-sensitive and full.  This switch
                                           // checks only up to the shortest string length.
  bool nocase_var_names{true};             // By default, variable name compares are
                                           // case sensitive and full.  This switch
                                           // ignores case when comparing.
  bool summary_flag{false};                // By default, summary mode is not in effect.
  bool ignore_maps{false};                 // By default, use the node and element number
                                           // maps to report node and element ids.
  bool ignore_nans{false};                 // Don't check for NaNs
  bool ignore_dups{false}; // If two elements/nodes in same location in map or partial map
                           // case, just return first match instead of aborting.
  bool ignore_steps{false};
  bool ignore_attributes{false}; // Don't compare attributes...
  bool ignore_sideset_df{false}; // Don't compare sideset df

  bool ints_64_bits{false};

  bool coord_sep{false};
  bool exit_status_switch{true};
  bool dump_mapping{false};         // By default, mappings are not printed.
  bool show_unmatched{false};       // Show elements not matched in partial mode
  bool noSymmetricNameCheck{false}; // By default, the second file's variable
  bool allowNameMismatch{false};    // By default, name in 1st db must be in second also.
  bool doL1Norm{false};
  bool doL2Norm{false};
  bool pedantic{false}; // Be most picky on what is different (not fully picky yet)

  bool interpolating{false}; // Interpolate times on file2 to match times on file1;
  bool by_name{false};       // Match entities by name instead of by id.

  bool glob_var_do_all_flag{false};
  bool node_var_do_all_flag{false};
  bool elmt_var_do_all_flag{false};
  bool elmt_att_do_all_flag{false};
  bool ns_var_do_all_flag{false};
  bool ss_var_do_all_flag{false};
  bool eb_var_do_all_flag{false};
  bool fb_var_do_all_flag{false};

private:
  void          enroll_options();
  GetLongOption options_{}; //!< Options parsing
};

extern SystemInterface interFace;
