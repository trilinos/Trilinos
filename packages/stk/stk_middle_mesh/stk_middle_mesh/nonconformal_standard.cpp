
#include "nonconformal_standard.hpp"
#include "boundary_fixture.hpp"
#include "incremental_mesh_boundary_snapper.hpp"
#include "mesh_snapper.hpp"
#include <stdexcept>
#include "nonconformal4.hpp"
#include "stk_util/parallel/Parallel.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

std::shared_ptr<nonconformal::impl::NonconformalAbstract>
create_nonconformal_standard(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2,
                             const NonconformalOpts& opts)
{
  if (opts.enableVolumeSnap)
  {
    mesh::impl::MeshSnapper snapper({opts.volumeSnapTol});
    snapper.snap(mesh1, mesh2);
  }

  if (opts.enableSnapAndQuality)
  {
    IncrementalBoundarySnapperOpts snapperOpts;
    snapperOpts.qualityImproverOpts = {opts.nlayers, opts.maxDeltaX, opts.itermax, opts.delta};

    auto fixer = mesh::impl::make_incremental_boundary_snapper(mesh1, mesh2, parallel_machine_world(), snapperOpts);

    fixer->snap();
    if (fixer->get_mesh1_quality_improver()->count_invalid_points() > 0 ||
        fixer->get_mesh2_quality_improver()->count_invalid_points() > 0)
    {
      throw std::runtime_error("unable to make interface boundaries conform and valid");
    }
  }

  NormalProjectionOpts normalOpts;
  normalOpts.classifierTolerances = PointClassifierNormalWrapperTolerances(opts.eps);
  normalOpts.edgeTracerTolerances = middle_mesh::impl::EdgeTracerTolerances(opts.eps);

  auto maker = std::make_shared<nonconformal4::impl::Nonconformal4>(mesh1, mesh2, normalOpts);
  maker->create();

  return maker;
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
