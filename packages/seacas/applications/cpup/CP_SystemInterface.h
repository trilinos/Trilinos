/*
 * Copyright(C) 1999-2024 National Technology & Engineering Solutions
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
#include <vector>       // for vector

namespace Cpup {
  using StringVector = std::vector<std::string>;

  class SystemInterface
  {
  public:
    explicit SystemInterface(int rank = 0);

    bool parse_options(int argc, char **argv);

    int processor_count() const { return processorCount_; }
    int screen_width() const { return screenWidth_; }

    int step_min() const { return stepMin_; }
    int step_max() const { return stepMax_; }
    int step_interval() const { return stepInterval_; }

    void step_min(int my_step_min) { stepMin_ = my_step_min; }
    void step_max(int my_step_max) { stepMax_ = my_step_max; }
    void step_interval(int interval) { stepInterval_ = interval; }
    void set_output_filename(const std::string &filename) const { outputFilename_ = filename; }
    std::string output_filename() const;

    std::string cwd() const { return cwd_; }
    std::string basename() const { return basename_; }
    std::string cgns_suffix() const { return inExtension_; }
    std::string output_suffix() const;

    std::string root_dir() const { return rootDirectory_; }
    std::string sub_dir() const { return subDirectory_; }

    bool minimize_open_files() const { return minimizeOpenFiles_; }
    bool is_auto() const { return auto_; }
    int  debug() const { return debugLevel_; }
#if 0
    int start_part() const { return startPart_; }
    int part_count() const;
    int subcycle() const { return subcycle_; }
    int cycle() const { return cycle_; }

    void subcycle_join(bool tf) { subcycleJoin_ = tf; }
    void subcycle(int cycles) { subcycle_ = cycles; }
    void processor_count(int count) { processorCount_ = count; }
    bool add_nodal_communication_map() const { return addNodalCommunicationMap_; }
    bool add_processor_id_field() const { return addProcessorIdField_; }
    bool add_processor_id_map() const { return addProcessorIdMap_; }
    bool disable_field_recognition() const { return false; }
    bool sum_shared_nodes() const { return sumSharedNodes_; }
    bool append() const { return append_; }
    bool map_element_ids() const { return mapIds_; }
    bool omit_nodesets() const { return omitNodesets_; }
    bool omit_sidesets() const { return omitSidesets_; }
    int  compress_data() const { return compressData_; }
    bool zlib() const { return zlib_; }
    bool szip() const { return szip_; }
    bool subcycle_join() const { return subcycleJoin_; }
    bool keep_temporary() const { return keepTemporary_; }
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
#endif
    const StringVector &var_names() const { return varNames_; }
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

    StringVector varNames_;
    std::string  inExtension_{};
    std::string  outExtension_{};
    std::string  cwd_{};
    std::string  rootDirectory_{};
    std::string  subDirectory_{};
    std::string  basename_{};

    // Used for a storage area only.  Needed for subcyle and auto-join option
    // Not directly settable through the user-interFace (maybe should be?)
    mutable std::string outputFilename_{};

    int  myRank_{0};
    int  processorCount_{1};
    int  screenWidth_{};
    int  stepMin_{1};
    int  stepMax_{INT_MAX};
    int  stepInterval_{1};
    int  debugLevel_{};
    bool auto_{false};
    bool minimizeOpenFiles_{false};
#if 0
    int  startPart_{};
    int  partCount_{-1};
    int  subcycle_{-1};
    int  cycle_{-1};
    int  compressData_{0};
    int  maxOpenFiles_{0};
    bool zlib_{true};
    bool szip_{false};
    bool sumSharedNodes_{false};
    bool addProcessorIdField_{false};
    bool addProcessorIdMap_{false};
    bool mapIds_{true};
    bool omitNodesets_{false};
    bool omitSidesets_{false};
    bool append_{false};
    bool subcycleJoin_{false};
    bool keepTemporary_{false};
    bool verifyValidFile_{false};
    bool addNodalCommunicationMap_{false};

    StringIdVector nodeVarNames_;
    StringIdVector elemVarNames_;
    StringIdVector nsetVarNames_;
    StringIdVector ssetVarNames_;
#endif
  };

#if 0
  inline int SystemInterface::part_count() const
  {
    return partCount_ > 0 ? partCount_ : processorCount_;
  }
#endif
} // namespace Cpup
