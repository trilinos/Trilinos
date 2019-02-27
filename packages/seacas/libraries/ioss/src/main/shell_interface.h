/*
 * Copyright(C) 1999-2017 National Technology & Engineering Solutions
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
 */
#ifndef shell_SystemInterface_h
#define shell_SystemInterface_h

#include "Ioss_GetLongOpt.h"

#include <iosfwd>
#include <string>
#include <vector>

/** \brief A special namespace for the io_shell demonstration program interface.
 */
namespace IOShell {
  class Interface
  {
  public:
    Interface();
    ~Interface();

    bool parse_options(int argc, char **argv);

    //! Dumps representation of data in this class to cerr

    void enroll_options();

    Ioss::GetLongOption options_;

    std::vector<std::string> inputFile;
    std::string              outputFile;
    std::string              inFiletype{"unknown"};
    std::string              outFiletype{"unknown"};
    std::string              groupName;
    std::string              decomp_method;
    std::string              compose_output{"default"};
    double                   maximum_time{std::numeric_limits<double>::max()};
    double                   minimum_time{-std::numeric_limits<double>::max()};
    double                   append_time{std::numeric_limits<double>::max()};
    double                   timestep_delay{0.0};
    int                      append_step{std::numeric_limits<int>::max()};
    int                      surface_split_type{1};
    int                      data_storage_type{0};
    int                      compression_level{0};
    int                      serialize_io_size{0};
    int                      flush_interval{0};

    //! If non-zero, then put `split_times` timesteps in each file. Then close file and start new
    //! file.
    // If `split_cyclic == 0`, then filenames will be
    // filename.e-s000X; if `split_cyclic > 0`, then filenames will be
    // filename.e.{A|B|C...}
    int split_times{0};
    //! If non-zero, then the `split_times` timesteps will be put into
    // `split_cyclic` files and then recycle filenames.  For example:
    // `split_times=1` and `split_cyclic=3`, t=1 ->file.A, t=2
    // ->file.B, t=3 -> file.C, t=4 -> file.A If `split_times=2` and
    // `split_cyclic=2`, then t=1,2 -> file.A, t=3,4 -> file.B, t=5,6
    // -> file.A, t=7,8 -> file.B
    int  split_cyclic{0};
    bool shuffle{false};
    bool debug{false};
    bool statistics{false};
    bool memory_statistics{false};
    bool do_transform_fields{false};
    bool ints_64_bit{false};
    bool ints_32_bit{false};
    bool reals_32_bit{false};
    bool netcdf4{false};
    bool netcdf5{false};
    bool quiet{false};
    bool in_memory_read{false};
    bool in_memory_write{false};
    bool lower_case_variable_names{true};
    bool delete_timesteps{false};
    bool minimize_open_files{false};
    bool disable_field_recognition{false};
    bool file_per_state{
        false}; // Put transient data for each timestep in separate file (EXPERMENTAL)
    char fieldSuffixSeparator{'_'};
  };
} // namespace IOShell
#endif
