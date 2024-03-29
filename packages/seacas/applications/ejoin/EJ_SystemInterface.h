// Copyright(C) 1999-, 20232023,  National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once

#include "EJ_CodeTypes.h" // for StringIdVector, Omissions, etc
#include "EJ_vector3d.h"  // for vector3d
#include "GetLongOpt.h"   // for GetLongOption
#include <climits>        // for INT_MAX
#include <string>         // for string
#include <vector>         // for vector

class SystemInterface
{
public:
  SystemInterface();

  bool parse_options(int argc, char **argv);

  int debug() const { return debugLevel_; }

  double tolerance() const { return tolerance_; }
  bool   match_node_ids() const { return matchNodeIds_; }
  bool   match_node_xyz() const { return matchNodeXYZ_; }
  bool   match_elem_ids() const { return matchElemIds_; }
  bool   omit_nodesets() const { return omitNodesets_; }
  bool   omit_sidesets() const { return omitSidesets_; }
  bool   omit_assemblies() const { return omitAssemblies_; }
  bool   convert_nodes_to_nodesets(int part_number) const;
  bool   create_assemblies() const { return createAssemblies_; }

  bool disable_field_recognition() const { return disableFieldRecognition_; }
  bool ints64bit() const { return ints64bit_; }
  bool use_netcdf4() const { return useNetcdf4_; }
  bool ignore_element_ids() const { return ignoreElementIds_; }

  int  compression_level() const { return compressionLevel_; }
  bool zlib() const { return zlib_; }
  bool szip() const { return szip_; }
  int  step_min() const { return stepMin_; }
  int  step_max() const { return stepMax_; }
  int  step_interval() const { return stepInterval_; }

  vector3d                offset() const { return offset_; }
  const std::vector<int> &information_record_parts() const { return infoRecordParts_; }
  const StringIdVector   &global_var_names() const { return globalVarNames_; }
  const StringIdVector   &node_var_names() const { return nodeVarNames_; }
  const StringIdVector   &elem_var_names() const { return elemVarNames_; }
  const StringIdVector   &nset_var_names() const { return nsetVarNames_; }
  const StringIdVector   &sset_var_names() const { return ssetVarNames_; }

  const Omissions &block_inclusions() const { return blockInclusions_; }
  const Omissions &block_omissions() const { return blockOmissions_; }
  const Omissions &nset_omissions() const { return nsetOmissions_; }
  const Omissions &sset_omissions() const { return ssetOmissions_; }
  const Omissions &assembly_omissions() const { return assemblyOmissions_; }

  const std::string &block_prefix() const { return blockPrefix_; }

  //! Dumps representation of data in this class to cerr
  static void show_version();

  // Make this private eventually...
  StringVector inputFiles_;
  std::string  outputName_;

private:
  void enroll_options();

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

  int  debugLevel_{0};
  int  stepMin_{1};
  int  stepMax_{INT_MAX};
  int  stepInterval_{1};
  int  compressionLevel_{0};
  bool omitNodesets_{false};
  bool omitSidesets_{false};
  bool omitAssemblies_{false};
  bool matchNodeIds_{false};
  bool matchNodeXYZ_{false};
  bool matchElemIds_{false};
  bool disableFieldRecognition_{false};
  bool ints64bit_{false};
  bool useNetcdf4_{false};
  bool ignoreElementIds_{false};
  bool zlib_{true};
  bool szip_{false};
  bool createAssemblies_{true};

  std::string blockPrefix_{"p"};

  vector3d offset_;
  double   tolerance_{0.0};

  Omissions blockInclusions_;
  Omissions blockOmissions_;
  Omissions assemblyOmissions_;
  Omissions nsetOmissions_;
  Omissions ssetOmissions_;

  std::vector<int> nodesetConvertParts_;
  std::vector<int> infoRecordParts_;

  StringIdVector globalVarNames_;
  StringIdVector nodeVarNames_;
  StringIdVector elemVarNames_;
  StringIdVector nsetVarNames_;
  StringIdVector ssetVarNames_;
};
