// Copyright(C) 2008 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#ifndef Sierra_SystemInterface_h
#define Sierra_SystemInterface_h

#include "GetLongOpt.h" // for GetLongOption
#include "Tolerance.h"  // for Tolerance, etc
#include "terminal_color.h"
#include <cmath>
#include <string>  // for string
#include <utility> // for pair
#include <vector>  // for vector

#define DEFAULT_MAX_NUMBER_OF_NAMES 1000

#if 0
#define ERROR(x) std::cerr << trmclr::red << "exodiff: ERROR: " << x << trmclr::normal
#else
#define ERROR(x) std::cerr << "exodiff: ERROR: " << x
#endif

class SystemInterface
{
public:
  SystemInterface();
  ~SystemInterface();

  bool parse_options(int argc, char **argv);

  static void show_version();

  void Parse_Command_File();
  void Set_Max_Names(int size);

  // Program parameters.
  bool quiet_flag;                 // By default, warnings and other info is produced
  bool show_all_diffs;             // Be default, show only maximum diff for each variable;
                                   // if set, show all diff that exceeds tolerance.
  TOLERANCE_TYPE_enum output_type; // By default, output file diffs are absolute.
  MAP_TYPE_enum       map_flag;    // By default, no searching is done to match
                                   // nodes & elements.
  bool nsmap_flag;                 // By default, nodeset nodelist match is off
  bool ssmap_flag;                 // By default, sideset elem/side match is off
  bool short_block_check;          // By default, element block compares are
                                   // case in-sensitive and full.  This switch
                                   // checks only up to the shortest string length.
  bool nocase_var_names;           // By default, variable name compares are
                                   // case sensitive and full.  This switch
                                   // ignores case when comparing.
  bool summary_flag;               // By default, summary mode is not in effect.
  bool ignore_maps;                // By default, use the node and element number
                                   // maps to report node and element ids.
  bool ignore_nans;                // Don't check for NaNs
  bool ignore_dups;                // If two elements/nodes in same location in map or partial map
                                   // case, just return first match instead of aborting.

  bool ignore_attributes; // Don't compare attributes...
  bool ignore_sideset_df; // Don't compare sideset df

  bool ints_64_bits;

  bool coord_sep;
  bool exit_status_switch;
  bool dump_mapping;         // By default, mappings are not printed.
  bool show_unmatched;       // Show elements not matched in partial mode
  bool noSymmetricNameCheck; // By default, the second file's variable
  bool allowNameMismatch;    // By default, name in 1st db must be in second also.
  bool doL1Norm;
  bool doL2Norm;
  bool pedantic; // Be most picky on what is different (not fully picky yet)

  bool interpolating; // Interpolate times on file2 to match times on file1;
  bool by_name;       // Match entities by name instead of by id.

  // These should correspond to the values specified during parsing of
  // coordinate tolerance.
  Tolerance coord_tol;
  Tolerance time_tol;
  Tolerance final_time_tol;

  // Offset of timesteps between first and second databases.
  int time_step_offset;
  int time_step_start;     // First step to compare (1-based)
  int time_step_stop;      // Last step to compare
  int time_step_increment; // Step increment

  std::pair<int, int> explicit_steps; // Only compare these two steps (db1:db2) if nonzero.

  int max_number_of_names;

  Tolerance default_tol;

  std::vector<std::string> glob_var_names;
  bool                     glob_var_do_all_flag;
  Tolerance                glob_var_default;
  std::vector<Tolerance>   glob_var;

  std::vector<std::string> node_var_names;
  bool                     node_var_do_all_flag;
  Tolerance                node_var_default;
  std::vector<Tolerance>   node_var;

  std::vector<std::string> elmt_var_names;
  bool                     elmt_var_do_all_flag;
  Tolerance                elmt_var_default;
  std::vector<Tolerance>   elmt_var;

  std::vector<std::string> elmt_att_names;
  bool                     elmt_att_do_all_flag;
  Tolerance                elmt_att_default;
  std::vector<Tolerance>   elmt_att;

  std::vector<std::string> ns_var_names;
  bool                     ns_var_do_all_flag;
  Tolerance                ns_var_default;
  std::vector<Tolerance>   ns_var;

  std::vector<std::string> ss_var_names;
  bool                     ss_var_do_all_flag;
  Tolerance                ss_var_default;
  std::vector<Tolerance>   ss_var;

  Tolerance ss_df_tol;

  // time step exclusion data
  std::vector<int> exclude_steps;

  std::string file1;
  std::string file2;
  std::string diff_file;
  std::string command_file;

private:
  void          enroll_options();
  GetLongOption options_; //!< Options parsing
};

extern SystemInterface interface;
#endif
