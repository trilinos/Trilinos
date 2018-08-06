// Copyright(C) 2009-2010-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of NTESS nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
#ifndef Sierra_SystemInterface_h
#define Sierra_SystemInterface_h

#include "CJ_CodeTypes.h"
#include "GetLongOpt.h" // for GetLongOption
#include <iosfwd>       // for ostream
#include <string>       // for string

namespace Excn {
  class SystemInterface
  {
  public:
    SystemInterface();
    ~SystemInterface();

    bool parse_options(int argc, char **argv);

    int debug() const { return debugLevel_; }
    int screen_width() const { return screenWidth_; }
    int compress_data() const { return compressionLevel_; }

    bool omit_nodesets() const { return omitNodesets_; }
    bool omit_sidesets() const { return omitSidesets_; }
    bool ints_64_bit() const { return ints64Bit_; }
    bool ignore_coordinates() const { return ignoreCoordinates_; }
    bool use_netcdf4() const { return useNetcdf4_; }

    double alive_value() const { return aliveValue_; }

    double         interpart_minimum_time_delta() const { return interpartMinimumTimeDelta_; }
    StringIdVector global_var_names() const { return globalVarNames_; }
    StringIdVector node_var_names() const { return nodeVarNames_; }
    StringIdVector elem_var_names() const { return elemVarNames_; }
    StringIdVector nset_var_names() const { return nsetVarNames_; }
    StringIdVector sset_var_names() const { return ssetVarNames_; }

    const std::string &element_status_variable() const { return elementStatusVariable_; }
    const std::string &nodal_status_variable() const { return nodalStatusVariable_; }
    const std::string &combined_mesh_status_variable() const { return meshCombineStatusVariable_; }

    //! Dumps representation of data in this class to cerr
    void dump(std::ostream &str) const;

    static void show_version();

    // Make this private eventually...
    StringVector inputFiles_;
    std::string  outputName_;

  private:
    void enroll_options();

    GetLongOption options_; //!< Options parsing

    int  debugLevel_{0};
    int  screenWidth_{0};
    int  compressionLevel_{0};
    bool omitNodesets_{false};
    bool omitSidesets_{false};
    bool ints64Bit_{false};
    bool ignoreCoordinates_{false};
    bool useNetcdf4_{false};

    double aliveValue_{-1.0};
    double interpartMinimumTimeDelta_{0.0};

    std::string elementStatusVariable_;
    std::string nodalStatusVariable_;

    // Mesh status variable to combine with elementStatusVariable_
    std::string meshCombineStatusVariable_;

    StringIdVector globalVarNames_;
    StringIdVector nodeVarNames_;
    StringIdVector elemVarNames_;
    StringIdVector nsetVarNames_;
    StringIdVector ssetVarNames_;
  };
} // namespace Excn
#endif
