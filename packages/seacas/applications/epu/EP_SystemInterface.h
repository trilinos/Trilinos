/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef Sierra_SystemInterface_h
#define Sierra_SystemInterface_h

#include "GetLongOpt.h" // for GetLongOption
#include <climits>      // for INT_MAX
#include <iosfwd>       // for ostream
#include <string>       // for string
#include <utility>      // for pair
#include <vector>       // for vector

namespace Excn {
  typedef std::vector<std::pair<std::string, int>> StringIdVector;

  class SystemInterface
  {
  public:
    explicit SystemInterface(int rank = 0);
    ~SystemInterface();

    bool parse_options(int argc, char **argv);

    int debug() const { return debugLevel_; }
    int raid_offset() const { return raidOffset_; }
    int raid_count() const { return raidCount_; }
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

    bool add_processor_id_field() const { return addProcessorId_; }
    bool sum_shared_nodes() const { return sumSharedNodes_; }
    bool use_netcdf4() const { return useNetcdf4_; }
    void set_use_netcdf4() const { useNetcdf4_ = true; }
    bool append() const { return append_; }
    bool map_element_ids() const { return mapIds_; }
    bool omit_nodesets() const { return omitNodesets_; }
    bool omit_sidesets() const { return omitSidesets_; }
    bool int64() const { return intIs64Bit_; }
    void set_int64() const { intIs64Bit_ = true; }
    int  compress_data() const { return compressData_; }
    bool subcycle_join() const { return subcycleJoin_; }
    bool output_shared_nodes() const { return outputSharedNodes_; }
    bool is_auto() const { return auto_; }
    bool keep_temporary() const { return keepTemporary_; }
    int  max_open_files() const
    {
      return maxOpenFiles_;
    } // Used to test auto-subcyling without thousands of files...

    StringIdVector global_var_names() const { return globalVarNames_; }
    StringIdVector node_var_names() const { return nodeVarNames_; }
    StringIdVector elem_var_names() const { return elemVarNames_; }
    StringIdVector nset_var_names() const { return nsetVarNames_; }
    StringIdVector sset_var_names() const { return ssetVarNames_; }

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
    int          raidOffset_{};
    int          raidCount_{};
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
    bool         sumSharedNodes_{false};
    bool         addProcessorId_{false};
    bool         mapIds_{true};
    bool         omitNodesets_{false};
    bool         omitSidesets_{false};
    mutable bool useNetcdf4_{false};
    bool         append_{false};
    mutable bool intIs64Bit_{false};
    bool         subcycleJoin_{false};
    bool         outputSharedNodes_{false};
    bool         auto_{false};
    bool         keepTemporary_{false};

    StringIdVector globalVarNames_;
    StringIdVector nodeVarNames_;
    StringIdVector elemVarNames_;
    StringIdVector nsetVarNames_;
    StringIdVector ssetVarNames_;
  };

  inline int SystemInterface::part_count() const
  {
    return partCount_ > 0 ? partCount_ : processorCount_;
  }
} // namespace Excn
#endif
