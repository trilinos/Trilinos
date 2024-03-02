#ifndef KrinoMeshAdapt_hpp
#define KrinoMeshAdapt_hpp

#include <stk_mesh/base/Types.hpp>
#include <stk_util/parallel/Parallel.hpp>

namespace krino { struct MeshAdaptAlgorithmParameters; }

namespace krino
{

void refine_mesh_with_params(stk::mesh::BulkData & mesh,
    const MeshAdaptAlgorithmParameters &algParams,
    const stk::ParallelMachine comm);

}

#endif
