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
#ifndef skinner_SystemInterface_h
#define skinner_SystemInterface_h

#include "Ioss_GetLongOpt.h" // for GetLongOption
#include <iosfwd>            // for ostream
#include <string>            // for string

/** \brief A special namespace for the skinner demonstration program interface.
 */
namespace Skinner {
  class Interface
  {
  public:
    Interface();
    ~Interface();

    bool parse_options(int argc, char **argv);

    bool ints_64_bit() const { return ints64Bit_; }

    bool no_output() const { return noOutput_; }

    std::string input_filename() const { return inputFile_; }
    std::string output_filename() const { return outputFile_; }
    std::string input_type() const { return inFiletype_; }
    std::string output_type() const { return outFiletype_; }

  private:
    void enroll_options();

    Ioss::GetLongOption options_;

    std::string inputFile_;
    std::string outputFile_;
    std::string inFiletype_;
    std::string outFiletype_;

  public:
    std::string decomp_method;
    std::string compose_output;
    int         compression_level;
    bool        shuffle;
    bool        debug;
    bool        statistics;
    bool        ints64Bit_;
    bool        netcdf4;
    bool        ignoreFaceIds_;
    bool        noOutput_;
  };
} // namespace Skinner
#endif
