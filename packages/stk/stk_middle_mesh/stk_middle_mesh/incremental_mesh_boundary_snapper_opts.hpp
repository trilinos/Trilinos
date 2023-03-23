#ifndef INCREMENTAL_MESH_BOUNDARY_SNAPPER_OPTS_H
#define INCREMENTAL_MESH_BOUNDARY_SNAPPER_OPTS_H

#include "mesh_quality_improver_opts.hpp"

namespace stk {
namespace middle_mesh {

struct IncrementalBoundarySnapperOpts
{
    int boundarySnapNsteps = 5; // The IncrementalMeshBoundarySnapper does an alternating
                                // pattern of applying a portion of the boundary snap and
                                // quality improvement.  This parameter determines the number
                                // steps.

    MeshQualityImproverOpts qualityImproverOpts;
};

} // namespace middle_mesh
} // namespace stk

#endif