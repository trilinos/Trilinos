// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
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
    bool sort_times() const { return sortTimes_; }

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
    bool sortTimes_{false};

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
