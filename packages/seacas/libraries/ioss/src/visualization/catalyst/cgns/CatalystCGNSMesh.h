// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef __CATALYST_CGNS_MESH_H
#define __CATALYST_CGNS_MESH_H

#include "CatalystCGNSMeshBase.h"
#include "CatalystManager.h"
#include <vector>
#include <map>

class vtkMultiBlockDataSet;

namespace Iovs_cgns {

class CatalystCGNSMesh : public CatalystCGNSMeshBase {

    using CatalystPipelineInfo = Iovs::CatalystManager::CatalystPipelineInfo;

public:

    CatalystCGNSMesh(Iovs::CatalystManager *cm,
        CatalystPipelineInfo& catalystPipelineInfo);

    ~CatalystCGNSMesh();

    void PerformCoProcessing(std::vector<int> &error_and_warning_codes,
                             std::vector<std::string> &error_and_warning_messages);

    void SetTimeData(double currentTime, int timeStep);

    void ReleaseMemory();

    void logMemoryUsageAndTakeTimerReading();

    void Delete();

    void CreateBase(int base_id,
                    const std::string& base_name);

    void AddStructuredZoneData(int base_id,
                               int zone_id,
                               const std::string& zone_name,
                               const std::string& data_name,
                               int ni,
                               int nj,
                               int nk,
                               int comp_count,
                               bool is_cell_field,
                               char field_suffix_separator,
                               double* data,
                               int size);

    vtkMultiBlockDataSet* getMultiBlockDataSet();

private:

    const unsigned int BASES_BLOCK_ID   = 0;
    const char *       BASES_BLOCK_NAME = "Bases";

    const unsigned int ZONES_BLOCK_ID   = 0;
    const char *       ZONES_BLOCK_NAME = "Zones";

    CatalystCGNSMesh();
    CatalystCGNSMesh(const CatalystCGNSMesh &) = delete;
    CatalystCGNSMesh &operator=(const CatalystCGNSMesh &) = delete;

    std::string createFieldVariableName(std::string fieldNamePrefix,
        char fieldSuffixSeparator, int componentIndex, int componentCount);  

    struct base {
        int base_location;
        std::map<int, int> zone_id_to_zone_location_map;
    };

    std::map<int, base> base_id_to_base_map;

    vtkMultiBlockDataSet* multiBlock = nullptr;
    Iovs::CatalystManager* catManager = nullptr;
    bool writeCatalystMesh;
    std::string catalystMeshFilePrefix;
    CatalystPipelineInfo catalystPipelineInfo;
};

} // namespace Iovs_cgns

#endif // __CATALYST_CGNS_MESH_H
