#include "mesh_quality_improver.h"
#include "mesh_io.h"
#include "optimization_step.h"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

void MeshQualityImprover::run()
{
  // improve overall quality
  std::cout << "\nrunning final quality improvement" << std::endl;
  std::vector<opt::impl::ActiveVertData*> verts;
  get_all_verts(verts);
  RunOpts opts{m_distTol, m_itermax, false, true};
  run_single(m_quality, verts, opts);
}

void MeshQualityImprover::run_single(std::shared_ptr<opt::impl::OptimizationStep> step,
                                     std::vector<opt::impl::ActiveVertData*>& verts, const RunOpts& opts)
{
  for (int i = 0; i < opts.itermax; ++i)
  {
    double deltaX1 = step->improve_quality(verts.begin(), verts.end(), opts.trimValid);
    double deltaX2 = step->improve_quality(verts.rbegin(), verts.rend(), opts.trimValid);
    double deltaX  = std::max(deltaX1, deltaX2);
    std::cout << "iteration " << i << " max delta_x = " << deltaX << std::endl;
    if (deltaX < opts.maxDeltaX && !(opts.requireValid && count_invalid_points() != 0))
    {
      std::cout << "dist_tol satisfied" << std::endl;
      break;
    }
  }
}

void MeshQualityImprover::get_invalid_verts(std::shared_ptr<opt::impl::OptimizationStep> opt,
                                            std::vector<opt::impl::ActiveVertData*>& invalidVerts)
{
  invalidVerts.clear();
  for (auto& active : m_activeVerts)
  {
    if (opt->has_invalid(active))
    {
      invalidVerts.push_back(&active);
    }
  }
}

void MeshQualityImprover::get_all_verts(std::vector<opt::impl::ActiveVertData*>& verts)
{
  verts.clear();
  for (auto& active : m_activeVerts)
    verts.push_back(&active);
}

int MeshQualityImprover::count_invalid_points()
{
  int count = 0;
  for (auto& active : m_activeVerts)
  {
    if (m_quality->has_invalid(active))
    {
      count += 1;
      continue;
    }
  }

  return count;
}

} // namespace impl
} // namespace mesh
} // namespace middle_mesh
} // namespace stk
