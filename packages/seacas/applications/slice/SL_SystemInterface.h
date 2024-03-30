// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once

#include "GetLongOpt.h"

#include <iosfwd>
#include <string>
#include <vector>
using StringVector   = std::vector<std::string>;
using Omissions      = std::vector<StringVector>;
using StringIdVector = std::vector<std::pair<std::string, int>>;

class SystemInterface
{
public:
  SystemInterface();

  bool parse_options(int argc, char **argv);

  size_t             processor_count() const { return processorCount_; }
  const std::string &decomposition_method() const { return decompMethod_; }
  const std::string &decomposition_file() const { return decompFile_; }
  const std::string &decomposition_variable() const { return decompVariable_; }
  const std::string &output_path() const { return outputPath_; }

  int debug() const { return debugLevel_; }
  int screen_width() const { return screenWidth_; }

  int step_min() const { return stepMin_; }
  int step_max() const { return stepMax_; }
  int step_interval() const { return stepInterval_; }

  size_t max_files() const { return maxFiles_; }
  size_t partial() const { return partialReadCount_; }
  bool   contiguous_decomposition() const { return contig_; }

  const StringIdVector &global_var_names() const { return globalVarNames_; }
  const StringIdVector &node_var_names() const { return nodeVarNames_; }
  const StringIdVector &elem_var_names() const { return elemVarNames_; }
  const StringIdVector &nset_var_names() const { return nsetVarNames_; }
  const StringIdVector &sset_var_names() const { return ssetVarNames_; }

  const Omissions &block_omissions() const { return blockOmissions_; }

  //! Dumps representation of data in this class to cerr
  void dump(std::ostream &str) const;

  static void show_version();

  // Make this private eventually...
  std::string inputFile_;
  std::string inputFormat_;
  std::string nemesisFile_;

private:
  void enroll_options();
  void parse_exclude(const char *list);

  /*! The defined formats for the count attribute are:<br>
    - <missing> -- default -- 1 <= count <= oo  (all steps)
    - "X"                  -- X <= count <= X  (just step X) (if X == -1, do last step only)
    - "X:Y"                -- X to Y by 1
    - "X:"                 -- X to oo by 1
    - ":Y"                 -- 1 to Y by 1
    - "::Z"                -- 1 to oo by Z

    The count and step must always be >= 0
  */
  void parse_step_option(const char *tokens);

  GetLongOption options_; //!< Options parsing

  std::string decompMethod_{"linear"};
  std::string decompFile_;
  std::string decompVariable_{"processor_id"};
  std::string outputPath_;

  size_t partialReadCount_{1'000'000'000};
  size_t maxFiles_{1'020};
  int    processorCount_{1};
  int    debugLevel_{0};
  int    screenWidth_{0};
  int    stepMin_{1};
  int    stepMax_{1 << 30};
  int    stepInterval_{1};

public:
  int         compressionLevel_{0};
  bool        shuffle_{false};
  bool        ints64Bit_{false};
  bool        netcdf4_{false};
  bool        netcdf5_{false};
  bool        disableFieldRecognition_{false};
  bool        szip_{false};
  bool        zlib_{true};
  bool        outputDecompMap_{false};
  bool        outputDecompField_{false};
  bool        ignore_x_{false};
  bool        ignore_y_{false};
  bool        ignore_z_{false};
  bool        lineDecomp_{false};
  std::string lineSurfaceList_{};

private:
  bool contig_{false};

  Omissions blockOmissions_;

  StringIdVector globalVarNames_;
  StringIdVector nodeVarNames_;
  StringIdVector elemVarNames_;
  StringIdVector nsetVarNames_;
  StringIdVector ssetVarNames_;
};
