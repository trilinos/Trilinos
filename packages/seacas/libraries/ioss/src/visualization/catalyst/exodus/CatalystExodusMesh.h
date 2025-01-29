// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef __CATALYST_EXODUS_MESH_H
#define __CATALYST_EXODUS_MESH_H

#include "vtkNew.h"
#include "vtkPartitionedDataSetCollection.h"
#include <map>
#include <vector>
#include <visualization/catalyst/manager/CatalystManager.h>
#include <visualization/exodus/CatalystExodusMeshBase.h>

class vtkUnstructuredGrid;
class vtkVariant;
class vtkPoints;

namespace Iovs_exodus {

  class CatalystExodusMesh : public CatalystExodusMeshBase
  {

    using CatalystPipelineInfo = Iovs::CatalystManager::CatalystPipelineInfo;

  public:
    CatalystExodusMesh(Iovs::CatalystManager *cm, CatalystPipelineInfo &catalystPipelineInfo);

    ~CatalystExodusMesh();

    void PerformCoProcessing(std::vector<int>         &error_and_warning_codes,
                             std::vector<std::string> &error_and_warning_messages);

    void SetTimeData(double currentTime, int timeStep);

    void logMemoryUsageAndTakeTimerReading();

    void ReleaseMemory();

    void Delete();

    void CreateGlobalVariable(const std::string &variable_name, int num_comps, const double *data);

    void CreateGlobalVariable(const std::string &variable_name, int num_comps, const int *data);

    void CreateGlobalVariable(const std::string &variable_name, int num_comps, const int64_t *data);

    void InitializeGlobalPoints(int num_points, int dimension, const double *data);

    void InitializeElementBlocks(const ElementBlockIdNameList &elemBlockNameIdList);

    void CreateElementBlock(const char *elem_block_name, int elem_block_id,
                            const std::string &elem_type, int nodes_per_elem, int num_elem,
                            const int64_t *global_elem_ids, int *connectivity);

    void CreateElementBlock(const char *elem_block_name, int elem_block_id,
                            const std::string &elem_type, int nodes_per_elem, int num_elem,
                            const int64_t *global_elem_ids, int64_t *connectivity);

    void CreateElementVariable(const std::string &variable_name, int num_comps, int elem_block_id,
                               const double *data);

    void CreateElementVariable(const std::string &variable_name, int num_comps, int elem_block_id,
                               const int *data);

    void CreateElementVariable(const std::string &variable_name, int num_comps, int elem_block_id,
                               const int64_t *data);

    void CreateNodalVariable(const std::string &variable_name, int num_comps, const double *data);

    void CreateNodalVariable(const std::string &variable_name, int num_comps, const int *data);

    void CreateNodalVariable(const std::string &variable_name, int num_comps, const int64_t *data);

    // Description:
    // If true (the default), vector variables will contain a
    // trailing underscore in their name.  The default behavior
    // is consistent with the ParaView Exodus II file reader.
    bool UnderscoreVectorsON();
    void SetUnderscoreVectors(bool status);

    // Description:
    // If true (the default), displacements will be applied to the
    // mesh nodes before being sent to the in-situ pipeline.  The node
    // displacement variable is called either DISPL or displ.  The
    // default behavior is consistent with the ParaView Exodus II
    // file reader.
    bool ApplyDisplacementsON();
    void SetApplyDisplacements(bool status);

    vtkPartitionedDataSetCollection *getPartitionedDataSetCollection();

  private:
    std::map<int, std::map<int, int>> ebmap;
    std::map<int, std::map<int, int>> ebmap_reverse;
    std::map<int, std::map<int, int>> global_elem_id_map;
    std::vector<int>                  global_point_id_to_global_elem_id;
    std::map<int, unsigned int>       ebidmap;
    double                            GetArrayValue(vtkVariant &v, const void *data, int index);
    void                              ReleaseGlobalPoints();
    vtkPoints                        *global_points;
    int                               num_global_points;
    bool                              writeCatalystMesh;
    std::string                       catalystMeshFilePrefix;

    void CreateElementBlockInternal(const char *elem_block_name, int elem_block_id,
                                    const std::string &elem_type, int nodes_per_elem, int num_elem,
                                    vtkVariant &v, const int64_t *global_elem_ids,
                                    void *connectivity);

    void CreateGlobalVariableVariant(const std::string &variable_name, int num_comps, vtkVariant &v,
                                     const void *data);
    void CreateGlobalVariableInternal(const std::string &variable_name, int num_comps,
                                      vtkUnstructuredGrid *ug, vtkVariant &v, const void *data);

    void CreateNodalVariableVariant(const std::string &variable_name, int num_comps, vtkVariant &v,
                                    const void *data);
    void CreateNodalVariableInternal(const std::string &variable_name, int num_comps,
                                     vtkUnstructuredGrid *ug, int element_block_id,
                                     std::map<int, std::map<int, int>> &point_map, vtkVariant &v,
                                     const void *data);

    void CreateElementVariableVariant(const std::string &variable_name, int num_comps,
                                      int elem_block_id, vtkVariant &v, const void *data);
    void CreateElementVariableInternal(const std::string &variable_name, int num_comps,
                                       vtkUnstructuredGrid *ug, vtkVariant &v, const void *data);

    void ReleaseMemoryInternal(vtkUnstructuredGrid *ug);

    CatalystExodusMesh();
    CatalystExodusMesh(const CatalystExodusMesh &)            = delete;
    CatalystExodusMesh &operator=(const CatalystExodusMesh &) = delete;

    vtkNew<vtkPartitionedDataSetCollection> vpdc;
    Iovs::CatalystManager                  *catManager = nullptr;
    bool                                    UnderscoreVectors;
    bool                                    ApplyDisplacements;
    CatalystPipelineInfo                    catalystPipelineInfo;
    const std::string                       ASSEMBLY_LABEL              = "label";
    const std::string                       ASSEMBLY_ROOT_NAME          = "IOSS";
    const std::string                       ASSEMBLY_ELEMENT_BLOCKS     = "element_blocks";
    const int                               PDS_UNSTRUCTURED_GRID_INDEX = 0;
    int                                     getElementBlocksAssemblyNode();
    vtkUnstructuredGrid                    *getUnstructuredGrid(int blockId);
  };

} // namespace Iovs_exodus

#endif // __CATALYST_EXODUS_MESH_H
