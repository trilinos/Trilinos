// Copyright(C) 2016-2017 National Technology & Engineering Solutions of
// Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above
//   copyright notice, this list of conditions and the following
//   disclaimer in the documentation and/or other materials provided
//   with the distribution.
//
// * Neither the name of NTESS nor the names of its
//   contributors may be used to endorse or promote products derived
//   from this software without specific prior written permission.
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

#include "GetLongOpt.h"

#include <iosfwd>
#include <string>
#include <vector>
using StringVector = std::vector<std::string>;
using Omissions    = std::vector<StringVector>;
typedef std::vector<std::pair<std::string, int>> StringIdVector;

class SystemInterface
{
public:
  SystemInterface();
  ~SystemInterface();

  bool parse_options(int argc, char **argv);

  size_t             processor_count() const { return processorCount_; }
  const std::string &decomposition_method() const { return decompMethod_; }
  const std::string &decomposition_file() const { return decompFile_; }
  const std::string &output_path() const { return outputPath_; }

  int debug() const { return debugLevel_; }
  int screen_width() const { return screenWidth_; }

  int step_min() const { return stepMin_; }
  int step_max() const { return stepMax_; }
  int step_interval() const { return stepInterval_; }

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
  std::string outputPath_;

  size_t partialReadCount_{1000000000};
  int    processorCount_{1};
  int    debugLevel_{0};
  int    screenWidth_{0};
  int    stepMin_{1};
  int    stepMax_{INT_MAX};
  int    stepInterval_{1};
  bool   omitNodesets_{false};
  bool   omitSidesets_{false};
  bool   disableFieldRecognition_{false};
  bool   contig_{false};

  Omissions blockOmissions_;

  StringIdVector globalVarNames_;
  StringIdVector nodeVarNames_;
  StringIdVector elemVarNames_;
  StringIdVector nsetVarNames_;
  StringIdVector ssetVarNames_;
};
#endif
