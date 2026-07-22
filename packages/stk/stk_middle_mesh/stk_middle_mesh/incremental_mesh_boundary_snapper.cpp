#include "incremental_mesh_boundary_snapper.hpp"
#include "create_mesh_quality_improver.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

void IncrementalMeshBoundarySnapper::snap()
{

  m_snapper.snap(m_mesh1, m_mesh2, m_unionComm);

  if (m_mesh1)
    record_final_position_and_restore_initial_position(m_mesh1, m_mesh1Data);

  if (m_mesh2)
    record_final_position_and_restore_initial_position(m_mesh2, m_mesh2Data);

  if (m_mesh1)
  {
    bool amIRoot = utils::impl::comm_rank(m_mesh1->get_comm()) == 0;

    if (m_improver1->verbose_output() && amIRoot) {
      std::cout << "applying displacement to mesh1" << std::endl;
    }
    apply_displacement(m_mesh1, m_mesh1Data, m_improver1, "mesh2");
  }

  if (m_mesh2)
  {
    bool amIRoot = utils::impl::comm_rank(m_mesh2->get_comm()) == 0;

    if (m_improver2->verbose_output() && amIRoot) {
      std::cout << "\napplying displacement to mesh2" << std::endl;
    }
    apply_displacement(m_mesh2, m_mesh2Data, m_improver2, "mesh2");
  }
}

void IncrementalMeshBoundarySnapper::record_initial_position(MeshPtr mesh, FieldPtr field)
{
  auto& data = *field;
  for (auto& vert : mesh->get_vertices())
    if (vert)
      data(vert, 0, 0).ptOrig = vert->get_point_orig(0);
}

void IncrementalMeshBoundarySnapper::record_final_position_and_restore_initial_position(MeshPtr mesh, FieldPtr field)
{
  auto& data = *field;
  for (auto& vert : mesh->get_vertices())
    if (vert)
    {
      data(vert, 0, 0).ptDest = vert->get_point_orig(0);
      vert->set_point_orig(0, data(vert, 0, 0).ptOrig);
    }
}

void IncrementalMeshBoundarySnapper::apply_displacement(MeshPtr mesh, FieldPtr field,
                                                        std::shared_ptr<MeshQualityImprover> improver,
                                                        const std::string& /*prefix*/)
{
  bool amIRoot = utils::impl::comm_rank(mesh->get_comm()) == 0;
  auto& data = *field;
  for (int step = 0; step < m_nsteps; ++step)
  {
    if (improver->verbose_output() && amIRoot) {
      std::cout << "\nApplying mesh snap step " << step << " / " << m_nsteps << std::endl;
    }
    for (auto& vert : mesh->get_vertices())
      if (vert)
      {
        utils::Point displacement = (data(vert, 0, 0).ptDest - data(vert, 0, 0).ptOrig) / m_nsteps;
        auto pt                   = vert->get_point_orig(0);
        vert->set_point_orig(0, pt + displacement);
      }

    //if (improver->verbose_output()) {
    //  std::cout << "before smoothing, number of invalid points = " << improver->countInvalidPoints() << std::endl;
    //}

    improver->run();
    //if (improver->verbose_output()) {
    //  std::cout << "after smoothing, number of invalid points = " << improver->countInvalidPoints() << std::endl;
    //}
  }
}

std::shared_ptr<IncrementalMeshBoundarySnapper>
make_incremental_boundary_snapper(MeshPtr mesh1, MeshPtr mesh2, MPI_Comm unionComm, const IncrementalBoundarySnapperOpts& opts)
{
  std::shared_ptr<MeshQualityImprover> fixer1, fixer2;
  if (mesh1)
  {
    BoundaryFixture filter1(mesh1);
    fixer1 = make_standard_improver(mesh1, filter1, opts.qualityImproverOpts);
  }

  if (mesh2)
  {
    BoundaryFixture filter2(mesh2);
    fixer2 = make_standard_improver(mesh2, filter2, opts.qualityImproverOpts);
  }

  return std::make_shared<IncrementalMeshBoundarySnapper>(mesh1, fixer1, mesh2, fixer2, opts.boundarySnapNsteps, unionComm);
}

} // namespace impl
} // namespace mesh
} // namespace middle_mesh
} // namespace stk
