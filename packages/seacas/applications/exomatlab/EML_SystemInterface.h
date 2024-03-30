// Copyright(C) 1999-2022, 2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once

#include <EML_CodeTypes.h> // for StringIdVector
#include <GetLongOpt.h>    // for GetLongOption
#include <string>          // for string

class SystemInterface
{
public:
  SystemInterface();

  bool parse_options(int argc, char **argv);

  char field_suffix() const { return fieldSuffix_; }

  StringIdVector global_var_names() const { return globalVarNames_; }
  StringIdVector vars_to_list() const { return varsToList_; }
  bool           list_vars() const { return listVars_; }

  std::string input_file() const { return inputFile_; }
  std::string output_file() const { return outputFile_; }

  double minimum_time() const { return minimumTime_; }
  double maximum_time() const { return maximumTime_; }

  //! Dumps representation of data in this class to cerr

  static void show_version();

private:
  void enroll_options();

  double minimumTime_{0.0};
  double maximumTime_{-1.0};

  GetLongOption options_; //!< Options parsing

  std::string inputFile_{};
  std::string outputFile_{};

  StringIdVector globalVarNames_;
  StringIdVector nodeVarNames_;
  StringIdVector elemVarNames_;
  StringIdVector nsetVarNames_;
  StringIdVector ssetVarNames_;
  StringIdVector varsToList_;

  bool listVars_{false};
  char fieldSuffix_{0};
};
