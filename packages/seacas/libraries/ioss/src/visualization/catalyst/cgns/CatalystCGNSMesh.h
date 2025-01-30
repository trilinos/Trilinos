// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef __CATALYST_CGNS_MESH_H
#define __CATALYST_CGNS_MESH_H

#include "vtkNew.h"
#include "vtkPartitionedDataSetCollection.h"
#include <map>
#include <vector>
#include <visualization/catalyst/manager/CatalystManager.h>
#include <visualization/cgns/CatalystCGNSMeshBase.h>

class vtkStructuredGrid;

namespace Iovs_cgns {

  class CatalystCGNSMesh : public CatalystCGNSMeshBase
  {

    using CatalystPipelineInfo = Iovs::CatalystManager::CatalystPipelineInfo;

  public:
    CatalystCGNSMesh(Iovs::CatalystManager *cm, CatalystPipelineInfo &catalystPipelineInfo);

    ~CatalystCGNSMesh();

    void PerformCoProcessing(std::vector<int>         &error_and_warning_codes,
                             std::vector<std::string> &error_and_warning_messages);

    void SetTimeData(double currentTime, int timeStep);

    void ReleaseMemory();

    void logMemoryUsageAndTakeTimerReading();

    void Delete();

    void AddStructuredZoneData(const ZoneData &zoneData);

    vtkPartitionedDataSetCollection *getPartitionedDataSetCollection();

  private:
    CatalystCGNSMesh();
    CatalystCGNSMesh(const CatalystCGNSMesh &)            = delete;
    CatalystCGNSMesh &operator=(const CatalystCGNSMesh &) = delete;

    std::map<int, int> zone_id_to_zone_location_map;

    vtkNew<vtkPartitionedDataSetCollection> vpdc;
    Iovs::CatalystManager                  *catManager = nullptr;
    bool                                    writeCatalystMesh;
    std::string                             catalystMeshFilePrefix;
    CatalystPipelineInfo                    catalystPipelineInfo;
    const std::string                       ASSEMBLY_LABEL             = "label";
    const std::string                       ASSEMBLY_ROOT_NAME         = "IOSS";
    const std::string                       ASSEMBLY_STRUCTURED_BLOCKS = "structured_blocks";
    const int                               PDS_STRUCTURED_GRID_INDEX  = 0;
    int                                     getStructuredBlocksAssemblyNode();
    void               createPartitionedDataSet(const ZoneData &zoneData, vtkStructuredGrid *sg);
    vtkStructuredGrid *getStucturedGrid(const ZoneData &zoneData);
  };

} // namespace Iovs_cgns

#endif // __CATALYST_CGNS_MESH_H
