#ifndef NONCONFORMAL_STANDARD_H
#define NONCONFORMAL_STANDARD_H

#include "nonconformal_abstract.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

struct NonconformalOpts
{
    bool enableVolumeSnap = false; // run MeshSnapper (snaps all entities,
                                   // not just boundary)
    double volumeSnapTol = 1e-8;

    // MeshQualityImprover options
    bool enableSnapAndQuality = true; // enable both boundary snapping
                                      // and subsequent mesh quality
                                      // improver
    int nlayers      = -1;
    double maxDeltaX = 1e-2;
    int itermax      = 500;

    // quality metric options
    double delta = 1e-3;

    // PointClassification options
    double eps = 1e-12;
};

std::shared_ptr<nonconformal::impl::NonconformalAbstract>
create_nonconformal_standard(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2,
                             const NonconformalOpts& opts);

} // namespace impl
} // namespace middle_mesh
} // namespace stk
#endif
