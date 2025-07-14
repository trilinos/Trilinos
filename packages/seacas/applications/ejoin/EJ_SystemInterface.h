// Copyright(C) 1999-2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once

#include "EJ_CodeTypes.h"
#include "EJ_vector3d.h"
#include "GetLongOpt.h"
#include <climits>
#include <string>
#include <vector>

class SystemInterface
{
public:
  SystemInterface();

  bool parse_options(int argc, char **argv);

  int debug() const { return debugLevel_; }

  double tolerance() const { return tolerance_; }
  bool   match_node_ids() const { return matchNodeIds_; }
  bool   match_node_xyz() const { return matchNodeXYZ_; }
  bool   match_nodeset_nodes() const { return !nodesetMatch_.empty(); }
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

  int  quantize_nsd() const { return quantizeNSD_; }
  bool quantize() const { return quantizeNSD_ > 0; }
  int  compression_level() const { return compressionLevel_; }
  bool zlib() const { return zlib_; }
  bool szip() const { return szip_; }
  bool zstd() const { return zstd_; }
  bool bz2() const { return bz2_; }

  int step_min() const { return stepMin_; }
  int step_max() const { return stepMax_; }
  int step_interval() const { return stepInterval_; }

  bool               combine_nodesets() const { return combineNodesets_; }
  bool               combine_sidesets() const { return combineSidesets_; }
  bool               combine_element_blocks() const { return combineElementBlocks_; }
  const std::string &elementblock_combines() const { return elementBlockCombines_; }
  const std::string &nodeset_combines() const { return nodesetCombines_; }
  const std::string &sideset_combines() const { return sidesetCombines_; }

  const std::vector<vector3d> &offset() const { return offset_; }
  const vector3d              &offset(size_t part) const { return offset_[part]; }
  const std::vector<vector3d> &scale() const { return scale_; }
  const vector3d              &scale(size_t part) const { return scale_[part]; }
  const std::vector<int>      &information_record_parts() const { return infoRecordParts_; }
  const StringIdVector        &global_var_names() const { return globalVarNames_; }
  const StringIdVector        &node_var_names() const { return nodeVarNames_; }
  const StringIdVector        &elem_var_names() const { return elemVarNames_; }
  const StringIdVector        &nodeset_var_names() const { return nodesetVarNames_; }
  const StringIdVector        &sideset_var_names() const { return sidesetVarNames_; }

  const Omissions &block_inclusions() const { return blockInclusions_; }
  const Omissions &block_omissions() const { return blockOmissions_; }
  const Omissions &nodeset_omissions() const { return nodesetOmissions_; }
  const Omissions &sideset_omissions() const { return sidesetOmissions_; }
  const Omissions &assembly_omissions() const { return assemblyOmissions_; }
  const Omissions &nodeset_match() const { return nodesetMatch_; }

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
  int  quantizeNSD_{0};
  bool zlib_{false};
  bool szip_{false};
  bool zstd_{false};
  bool bz2_{false};

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

  bool combineElementBlocks_{false};
  bool combineNodesets_{false};
  bool combineSidesets_{false};

  bool createAssemblies_{true};

  std::string blockPrefix_{"p"};

  std::vector<vector3d> offset_;
  std::vector<vector3d> scale_;
  double                tolerance_{0.0};

  Omissions blockInclusions_;
  Omissions blockOmissions_;
  Omissions assemblyOmissions_;
  Omissions nodesetOmissions_;
  Omissions sidesetOmissions_;

  Omissions nodesetMatch_;

  std::string elementBlockCombines_{};
  std::string nodesetCombines_{};
  std::string sidesetCombines_{};

  std::vector<int> nodesetConvertParts_;
  std::vector<int> infoRecordParts_;

  StringIdVector globalVarNames_;
  StringIdVector nodeVarNames_;
  StringIdVector elemVarNames_;
  StringIdVector nodesetVarNames_;
  StringIdVector sidesetVarNames_;
};
