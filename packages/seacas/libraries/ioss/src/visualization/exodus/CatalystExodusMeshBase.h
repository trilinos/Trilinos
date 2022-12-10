// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef __CATALYST_EXODUS_MESH_BASE_H
#define __CATALYST_EXODUS_MESH_BASE_H

#include "iovs_export.h"

#include <string>
#include <vector>

namespace Iovs_exodus {

  class IOVS_EXPORT CatalystExodusMeshBase
  {

  public:
    CatalystExodusMeshBase(){};
    virtual ~CatalystExodusMeshBase(){};

    // Description:
    // Calls the ParaView Catalyst pipeline to run co-processing
    // for this time iteration.
    virtual void PerformCoProcessing(std::vector<int> &        error_and_warning_codes,
                                     std::vector<std::string> &error_and_warning_messages) = 0;

    // Description:
    // Sets time data for this ParaView Catalyst co-processing iteration.
    // currentTime is the current Ioss simulation time and timeStep is
    // the current time iteration count.
    virtual void SetTimeData(double currentTime, int timeStep) = 0;

    // Description:
    // Clears all nodal and element variables from the vtkMultiBlockDataSet.
    // Clears the global vtkPoints.
    virtual void ReleaseMemory() = 0;

    // Description:
    // Collects memory usage information from all processors and
    // writes the min, max, and mean to the log file.  Also writes the
    // min, max, and mean of the elapsed time since this method was
    // last called.
    virtual void logMemoryUsageAndTakeTimerReading() = 0;

    virtual void Delete() = 0;

    // Description:
    // Creates a global variable on the vtkExodusIIMultiBlockDataSet.
    // Creates the global variable on all element blocks.
    virtual void CreateGlobalVariable(std::vector<std::string> &component_names,
                                      const double *            data) = 0;

    // Description:
    // Creates a global variable on the vtkExodusIIMultiBlockDataSet.
    // Creates the global variable on all element blocks.
    virtual void CreateGlobalVariable(std::vector<std::string> &component_names,
                                      const int *               data) = 0;

    // Description:
    // Initializes the vtkMultiBlockDataSet with a global array of points
    // defined by num_points, dimension (2,3), and data.  Clears any existing data.
    virtual void InitializeGlobalPoints(int num_points, int dimension, const double *data) = 0;

    // Description:
    // Initializes the element blocks to NULL data sets with ids in element_block_id_list.
    // This method must be called first.
    virtual void InitializeElementBlocks(const std::vector<int> &element_block_id_list) = 0;

    // Description:
    // Creates a vtkUnstructuredGrid on the vtkExodusIIMultiBlockDataSet
    // that represents and element block in the Exodus II data.  The global_points
    // array contains all of the points in the Exodus II file.
    virtual void CreateElementBlock(const char *elem_block_name, int elem_block_id,
                                    const std::string &elem_type, int nodes_per_elem, int num_elem,
                                    const int64_t *global_elem_ids, int *connectivity) = 0;

    // Description:
    // Creates a vtkUnstructuredGrid on the vtkExodusIIMultiBlockDataSet
    // that represents and element block in the Exodus II data.  The global_points
    // array contains all of the points in the Exodus II file.
    virtual void CreateElementBlock(const char *elem_block_name, int elem_block_id,
                                    const std::string &elem_type, int nodes_per_elem, int num_elem,
                                    const int64_t *global_elem_ids, int64_t *connectivity) = 0;

    // Description:
    // Creates a vtkUnstructuredGrid representing the node set in the Exodus II
    // data. Node sets are arbitrary lists of mesh point ids.
    virtual void CreateNodeSet(const char *node_set_name, int node_set_id, int num_ids,
                               const int *data) = 0;

    // Description:
    // Creates a vtkUnstructuredGrid representing the node set in the Exodus II
    // data. Node sets are arbitrary lists of mesh point ids.
    virtual void CreateNodeSet(const char *node_set_name, int node_set_id, int num_ids,
                               const int64_t *data) = 0;

    // Description:
    // Creates a vtkUnstructuredGrid representing the side set (also Side Block) in
    // the Exodus II data. Side sets are collections of element faces and edges.
    virtual void CreateSideSet(const char *ss_owner_name, int side_set_id, int num_ids,
                               const int *element_ids, const int *face_ids) = 0;

    // Description:
    // Creates a vtkUnstructuredGrid representing the side set (also Side Block) in
    // the Exodus II data. Side sets are collections of element faces and edges.
    virtual void CreateSideSet(const char *ss_owner_name, int side_set_id, int num_ids,
                               const int64_t *element_ids, const int64_t *face_ids) = 0;

    // Description:
    // Creates an element variable the vtkExodusIIMultiBlockDataSet.
    virtual void CreateElementVariable(std::vector<std::string> &component_names, int elem_block_id,
                                       const double *data) = 0;

    // Description:
    // Creates an element variable the vtkExodusIIMultiBlockDataSet.
    virtual void CreateElementVariable(std::vector<std::string> &component_names, int elem_block_id,
                                       const int *data) = 0;

    // Description:
    // Creates an element variable the vtkExodusIIMultiBlockDataSet.
    virtual void CreateElementVariable(std::vector<std::string> &component_names, int elem_block_id,
                                       const int64_t *data) = 0;

    // Description:
    // Creates a nodal variable the vtkExodusIIMultiBlockDataSet.
    virtual void CreateNodalVariable(std::vector<std::string> &component_names,
                                     const double *            data) = 0;

    // Description:
    // Creates a nodal variable the vtkExodusIIMultiBlockDataSet.
    virtual void CreateNodalVariable(std::vector<std::string> &component_names,
                                     const int *               data) = 0;

    // Description:
    // Creates a nodal variable the vtkExodusIIMultiBlockDataSet.
    virtual void CreateNodalVariable(std::vector<std::string> &component_names,
                                     const int64_t *           data) = 0;
  };

} // namespace Iovs_exodus

#endif // __CATALYST_EXODUS_MESH_BASE_H
