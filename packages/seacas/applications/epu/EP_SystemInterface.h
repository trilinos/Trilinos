/*
 * Copyright(C) 2010 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
#ifndef Sierra_SystemInterface_h
#define Sierra_SystemInterface_h

#include "GetLongOpt.h" // for GetLongOption
#include <iosfwd>       // for ostream
#include <string>       // for string
#include <utility>      // for pair
#include <vector>       // for vector

namespace Excn {
  typedef std::vector<std::pair<std::string, int>> StringIdVector;

  class SystemInterface
  {
  public:
    SystemInterface();
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

    void subcycle(int cycles) { subcycle_ = cycles; }
    void processor_count(int count) { processorCount_ = count; }
    void step_min(int my_step_min) { stepMin_ = my_step_min; }
    void step_max(int my_step_max) { stepMax_ = my_step_max; }
    void step_interval(int interval) { stepInterval_ = interval; }

    std::string cwd() const { return cwd_; }
    std::string basename() const { return basename_; }
    std::string exodus_suffix() const { return inExtension_; }
    std::string output_suffix() const;

    std::string root_dir() const { return rootDirectory_; }
    std::string sub_dir() const { return subDirectory_; }

    bool           add_processor_id_field() const { return addProcessorId_; }
    bool           sum_shared_nodes() const { return sumSharedNodes_; }
    bool           use_netcdf4() const { return useNetcdf4_; }
    bool           append() const { return append_; }
    bool           map_element_ids() const { return mapIds_; }
    bool           omit_nodesets() const { return omitNodesets_; }
    bool           omit_sidesets() const { return omitSidesets_; }
    bool           int64() const { return intIs64Bit_; }
    void           set_int64() const { intIs64Bit_ = true; }
    int            compress_data() const { return compressData_; }
    bool           subcycle_join() const { return subcycleJoin_; }
    bool           output_shared_nodes() const { return outputSharedNodes_; }
    StringIdVector global_var_names() const { return globalVarNames_; }
    StringIdVector node_var_names() const { return nodeVarNames_; }
    StringIdVector elem_var_names() const { return elemVarNames_; }
    StringIdVector nset_var_names() const { return nsetVarNames_; }
    StringIdVector sset_var_names() const { return ssetVarNames_; }

    //! Dumps representation of data in this class to cerr
    void dump(std::ostream &str) const;

    static void show_version();

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

    std::string inExtension_;
    std::string outExtension_;
    std::string cwd_;
    std::string rootDirectory_;
    std::string subDirectory_;
    std::string basename_;

    int          raidOffset_;
    int          raidCount_;
    int          processorCount_;
    int          startPart_;
    int          partCount_;
    int          debugLevel_;
    int          screenWidth_;
    int          stepMin_;
    int          stepMax_;
    int          stepInterval_;
    int          subcycle_;
    int          cycle_;
    int          compressData_;
    bool         sumSharedNodes_;
    bool         addProcessorId_;
    bool         mapIds_;
    bool         omitNodesets_;
    bool         omitSidesets_;
    bool         useNetcdf4_;
    bool         append_;
    mutable bool intIs64Bit_;
    bool         subcycleJoin_;
    bool         outputSharedNodes_;

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
