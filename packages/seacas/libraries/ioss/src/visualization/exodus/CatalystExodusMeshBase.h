// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef __CATALYST_EXODUS_MESH_BASE_H
#define __CATALYST_EXODUS_MESH_BASE_H

#ifndef __CATALYST_PLUGIN_BUILD
#include "iovs_export.h"
#else
#define IOVS_EXPORT
#endif

#include <string>
#include <vector>

namespace Iovs_exodus {

  class IOVS_EXPORT CatalystExodusMeshBase
  {

  public:
    CatalystExodusMeshBase() {};
    virtual ~CatalystExodusMeshBase() {};

    // Description:
    // Calls the ParaView Catalyst pipeline to run co-processing
    // for this time iteration.
    virtual void PerformCoProcessing(std::vector<int>         &error_and_warning_codes,
                                     std::vector<std::string> &error_and_warning_messages) = 0;

    // Description:
    // Sets time data for this ParaView Catalyst co-processing iteration.
    // currentTime is the current Ioss simulation time and timeStep is
    // the current time iteration count.
    virtual void SetTimeData(double currentTime, int timeStep) = 0;

    // Description:
    // Collects memory usage information from all processors and
    // writes the min, max, and mean to the log file.  Also writes the
    // min, max, and mean of the elapsed time since this method was
    // last called.
    virtual void logMemoryUsageAndTakeTimerReading() = 0;

    // Description:
    // Clears memory buffers used in mesh construction and field data
    // for the current time-step. Clears the global vtkPoints.
    virtual void ReleaseMemory() = 0;

    virtual void Delete() = 0;

    // Description:
    // Creates a global variable
    // Creates the global variable on all element blocks.
    virtual void CreateGlobalVariable(const std::string &variable_name, int num_comps,
                                      const double *data) = 0;

    // Description:
    // Creates a global variable
    // Creates the global variable on all element blocks.
    virtual void CreateGlobalVariable(const std::string &variable_name, int num_comps,
                                      const int *data) = 0;

    // Description:
    // Creates a global variable
    // Creates the global variable on all element blocks.
    virtual void CreateGlobalVariable(const std::string &variable_name, int num_comps,
                                      const int64_t *data) = 0;

    // Description:
    // Initializes the global array of points
    // defined by num_points, dimension (2,3), and data.  Clears any existing data.
    virtual void InitializeGlobalPoints(int num_points, int dimension, const double *data) = 0;

    // Description:
    // Initializes the element blocks to NULL data sets with ids and names in elemBlkIdNameList.
    // This method must be called first.
    using ElementBlockIdNameList = std::vector<std::pair<int, std::string>>;
    virtual void InitializeElementBlocks(const ElementBlockIdNameList &elemBlkIdNameList) = 0;

    // Description:
    // Creates a vtkUnstructuredGrid
    // that represents and element block in the Exodus II data.  The global_points
    // array contains all of the points in the Exodus II file.
    virtual void CreateElementBlock(const char *elem_block_name, int elem_block_id,
                                    const std::string &elem_type, int nodes_per_elem, int num_elem,
                                    const int64_t *global_elem_ids, int *connectivity) = 0;

    // Description:
    // Creates a vtkUnstructuredGrid
    // that represents and element block in the Exodus II data.  The global_points
    // array contains all of the points in the Exodus II file.
    virtual void CreateElementBlock(const char *elem_block_name, int elem_block_id,
                                    const std::string &elem_type, int nodes_per_elem, int num_elem,
                                    const int64_t *global_elem_ids, int64_t *connectivity) = 0;

    // Description:
    // Creates an element variable
    virtual void CreateElementVariable(const std::string &variable_name, int num_comps,
                                       int elem_block_id, const double *data) = 0;

    // Description:
    // Creates an element variable
    virtual void CreateElementVariable(const std::string &variable_name, int num_comps,
                                       int elem_block_id, const int *data) = 0;

    // Description:
    // Creates an element variable
    virtual void CreateElementVariable(const std::string &variable_name, int num_comps,
                                       int elem_block_id, const int64_t *data) = 0;

    // Description:
    // Creates a nodal variable
    virtual void CreateNodalVariable(const std::string &variable_name, int num_comps,
                                     const double *data) = 0;

    // Description:
    // Creates a nodal variable
    virtual void CreateNodalVariable(const std::string &variable_name, int num_comps,
                                     const int *data) = 0;

    // Description:
    // Creates a nodal variable
    virtual void CreateNodalVariable(const std::string &variable_name, int num_comps,
                                     const int64_t *data) = 0;
  };

} // namespace Iovs_exodus

#endif // __CATALYST_EXODUS_MESH_BASE_H
