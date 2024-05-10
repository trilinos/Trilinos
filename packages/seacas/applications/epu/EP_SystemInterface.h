/*
 * Copyright(C) 1999-2022, 2024 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once

#include "GetLongOpt.h" // for GetLongOption
#include <climits>      // for INT_MAX
#include <iosfwd>       // for ostream
#include <string>       // for string
#include <utility>      // for pair
#include <vector>       // for vector

namespace Excn {
  using StringIdVector = std::vector<std::pair<std::string, int>>;

  class SystemInterface
  {
  public:
    explicit SystemInterface(int rank = 0);

    bool parse_options(int argc, char **argv);

    int debug() const { return debugLevel_; }
    int processor_count() const { return processorCount_; }
    int start_part() const { return startPart_; }
    int part_count() const;
    int screen_width() const { return screenWidth_; }

    int step_min() const { return stepMin_; }
    int step_max() const { return stepMax_; }
    int step_interval() const { return stepInterval_; }
    int subcycle() const { return subcycle_; }
    int cycle() const { return cycle_; }

    void subcycle_join(bool tf) { subcycleJoin_ = tf; }
    void subcycle(int cycles) { subcycle_ = cycles; }
    void processor_count(int count) { processorCount_ = count; }
    void step_min(int my_step_min) { stepMin_ = my_step_min; }
    void step_max(int my_step_max) { stepMax_ = my_step_max; }
    void step_interval(int interval) { stepInterval_ = interval; }
    void set_output_filename(const std::string &filename) const { outputFilename_ = filename; }
    std::string output_filename() const { return outputFilename_; }

    std::string cwd() const { return cwd_; }
    std::string basename() const { return basename_; }
    std::string exodus_suffix() const { return inExtension_; }
    std::string output_suffix() const;

    std::string root_dir() const { return rootDirectory_; }
    std::string sub_dir() const { return subDirectory_; }

    bool add_nodal_communication_map() const { return addNodalCommunicationMap_; }
    bool add_processor_id_field() const { return addProcessorIdField_; }
    bool add_processor_id_map() const { return addProcessorIdMap_; }
    bool sum_shared_nodes() const { return sumSharedNodes_; }
    bool use_netcdf4() const { return useNetcdf4_; }
    bool use_netcdf5() const { return useNetcdf5_; }
    void set_use_netcdf4() const { useNetcdf4_ = true; }
    bool append() const { return append_; }
    bool map_element_ids() const { return mapIds_; }
    bool map_edge_ids() const { return mapEdgeIds_; }
    bool map_face_ids() const { return mapFaceIds_; }
    bool omit_nodesets() const { return omitNodesets_; }
    bool omit_sidesets() const { return omitSidesets_; }
    bool omit_edgeblocks() const { return omitEdgeBlocks_; }
    bool omit_faceblocks() const { return omitFaceBlocks_; }
    bool int64() const { return intIs64Bit_; }
    void set_int64() const { intIs64Bit_ = true; }
    int  compress_data() const { return compressData_; }
    bool zlib() const { return zlib_; }
    bool szip() const { return szip_; }
    bool subcycle_join() const { return subcycleJoin_; }
    bool output_shared_nodes() const { return outputSharedNodes_; }
    bool is_auto() const { return auto_; }
    bool keep_temporary() const { return keepTemporary_; }
    bool remove_file_per_rank_files() const;
    bool verify_valid_file() const { return verifyValidFile_; }
    int  max_open_files() const
    {
      return maxOpenFiles_;
    } // Used to test auto-subcyling without thousands of files...

    StringIdVector global_var_names() const { return globalVarNames_; }
    StringIdVector node_var_names() const { return nodeVarNames_; }
    StringIdVector elem_var_names() const { return elemVarNames_; }
    StringIdVector nset_var_names() const { return nsetVarNames_; }
    StringIdVector sset_var_names() const { return ssetVarNames_; }
    StringIdVector edblk_var_names() const { return edblkVarNames_; }
    StringIdVector fablk_var_names() const { return fablkVarNames_; }

    //! Dumps representation of data in this class to cerr
    void dump(std::ostream &str) const;

    static void show_version(int rank = 0);

  private:
    void enroll_options();
    bool decompose_filename(const std::string &cs);

    /*! The defined formats for the count attribute are:<br>
        - <missing> -- default -- 1 <= count <= oo  (all steps)
        - "X"                  -- X <= count <= X  (just step X). -1 for last step.
        - "X:Y"                -- X to Y by 1
        - "X:"                 -- X to oo by 1
        - ":Y"                 -- 1 to Y by 1
        - "::Z"                -- 1 to oo by Z

     The count and step must always be >= 0
    */
    void parse_step_option(const char *tokens);

    GetLongOption options_; //!< Options parsing

    std::string inExtension_{};
    std::string outExtension_{};
    std::string cwd_{};
    std::string rootDirectory_{};
    std::string subDirectory_{};
    std::string basename_{};

    // Used for a storage area only.  Needed for subcyle and auto-join option
    // Not directly settable through the user-interFace (maybe should be?)
    mutable std::string outputFilename_{};

    int          myRank_{0};
    int          processorCount_{1};
    int          startPart_{};
    int          partCount_{-1};
    int          debugLevel_{};
    int          screenWidth_{};
    int          stepMin_{1};
    int          stepMax_{INT_MAX};
    int          stepInterval_{1};
    int          subcycle_{-1};
    int          cycle_{-1};
    int          compressData_{0};
    int          maxOpenFiles_{0};
    bool         zlib_{true};
    bool         szip_{false};
    bool         sumSharedNodes_{false};
    bool         addProcessorIdField_{false};
    bool         addProcessorIdMap_{false};
    bool         mapIds_{true};
    bool         mapEdgeIds_{true};
    bool         mapFaceIds_{true};
    bool         omitNodesets_{false};
    bool         omitSidesets_{false};
    bool         omitEdgeBlocks_{false};
    bool         omitFaceBlocks_{false};
    mutable bool useNetcdf4_{false};
    bool         useNetcdf5_{false};
    bool         append_{false};
    mutable bool intIs64Bit_{false};
    bool         subcycleJoin_{false};
    bool         outputSharedNodes_{false};
    bool         auto_{false};
    bool         keepTemporary_{false};
    bool         removeFilePerRankFiles_{false};
    bool         verifyValidFile_{false};
    bool         addNodalCommunicationMap_{false};

    StringIdVector globalVarNames_{};
    StringIdVector nodeVarNames_{};
    StringIdVector elemVarNames_{};
    StringIdVector nsetVarNames_{};
    StringIdVector ssetVarNames_{};
    StringIdVector edblkVarNames_{};
    StringIdVector fablkVarNames_{};
  };

  inline int SystemInterface::part_count() const
  {
    return partCount_ > 0 ? partCount_ : processorCount_;
  }
} // namespace Excn
