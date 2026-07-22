// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef __CATALYST_CGNS_MESH_BASE_H
#define __CATALYST_CGNS_MESH_BASE_H

#ifndef __CATALYST_PLUGIN_BUILD
#include "iovs_export.h"
#else
#define IOVS_EXPORT
#endif

#include <string>
#include <vector>

namespace Iovs_cgns {

  class IOVS_EXPORT CatalystCGNSMeshBase
  {

  public:
    CatalystCGNSMeshBase() {};
    virtual ~CatalystCGNSMeshBase() {};

    // Description:
    // Calls the ParaView Catalyst pipeline to run co-processing for this time iteration.
    virtual void PerformCoProcessing(std::vector<int>         &error_and_warning_codes,
                                     std::vector<std::string> &error_and_warning_messages) = 0;

    // Description:
    // Sets time data for this ParaView Catalyst co-processing iteration.
    // currentTime is the current Ioss simulation time and timeStep is
    // the current time iteration count.
    virtual void SetTimeData(double currentTime, int timeStep) = 0;

    // Description:
    // Clears all nodal and element variables
    // Clears the global vtkPoints.
    virtual void ReleaseMemory() = 0;

    // Description:
    // Collects memory usage information from all processors and
    // writes the min, max, and mean to the log file.  Also writes the
    // min, max, and mean of the elapsed time since this method was
    // last called.
    virtual void logMemoryUsageAndTakeTimerReading() = 0;

    virtual void Delete() = 0;

    struct ZoneData
    {
      int         zone_id;
      std::string zone_name;
      std::string data_name;
      int         ni;
      int         nj;
      int         nk;
      int         comp_count;
      bool        is_cell_field;
      enum { T_DOUBLE, T_INT, T_INT64 } data_type;
      union {
        int     *p_int;
        int64_t *p_int64;
        double  *p_double;
      } data;
      int size;
    };

    virtual void AddStructuredZoneData(const ZoneData &zoneData) = 0;
  };

} // namespace Iovs_cgns

#endif // __CATALYST_CGNS_MESH_BASE_H
