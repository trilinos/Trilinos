#include "mesh_quality_improver.hpp"
#include "mesh_io.hpp"
#include "optimization_step.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

void MeshQualityImprover::run()
{
  // improve overall quality
  std::vector<opt::impl::ActiveVertData*> verts;
  get_all_verts(verts);
  RunOpts opts{m_distTol, m_itermax, false, true};
  run_single(m_quality, verts, opts);
}

void MeshQualityImprover::run_single(std::shared_ptr<opt::impl::OptimizationStep> step,
                                     std::vector<opt::impl::ActiveVertData*>& verts, const RunOpts& opts)
{
  m_activeVertContainer.update_remote_coords();

  bool amIRoot = utils::impl::comm_rank(m_mesh->get_comm()) == 0;
  for (int i = 0; i < opts.itermax; ++i)
  {
    //std::cout << "\ndoing forward step" << std::endl;
    double deltaX1 = step->improve_quality(verts.begin(), verts.end(), opts.trimValid);
    m_activeVertContainer.update_remote_coords();

    //std::cout << "\ndoing backwards step" << std::endl;
    double deltaX2 = step->improve_quality(verts.rbegin(), verts.rend(), opts.trimValid);
    m_activeVertContainer.update_remote_coords();

    double deltaXLocal  = std::max(deltaX1, deltaX2), deltaXGlobal=0;
    MPI_Allreduce(&deltaXLocal, &deltaXGlobal, 1, MPI_DOUBLE, MPI_MAX, m_mesh->get_comm());

    if (amIRoot && verbose_output())
      std::cout << "iteration " << i << " max delta_x = " << deltaXGlobal << std::endl;

    if (deltaXGlobal < opts.maxDeltaX && !(opts.requireValid && count_invalid_points() != 0))
    {
      if (amIRoot)
        std::cout << "dist_tol satisfied" << std::endl;
      break;
    }
  }
}

void MeshQualityImprover::get_invalid_verts(std::shared_ptr<opt::impl::OptimizationStep> opt,
                                            std::vector<opt::impl::ActiveVertData*>& invalidVerts)
{
  invalidVerts.clear();
  for (auto& active : m_activeVertContainer.get_active_verts())
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
  for (auto& active : m_activeVertContainer.get_active_verts())
    verts.push_back(&active);
}

int MeshQualityImprover::count_invalid_points()
{
  int count = 0;
  for (auto& active : m_activeVertContainer.get_active_verts())
  {
    if (m_quality->has_invalid(active))
    {
      count += 1;
      continue;
    }
  }

  int countGlobal = 0;
  MPI_Allreduce(&count, &countGlobal, 1, MPI_INT, MPI_SUM, m_mesh->get_comm());

  return countGlobal;
}

} // namespace impl
} // namespace mesh
} // namespace middle_mesh
} // namespace stk
