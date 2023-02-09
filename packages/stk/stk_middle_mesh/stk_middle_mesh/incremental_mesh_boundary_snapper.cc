#include "incremental_mesh_boundary_snapper.h"
#include "create_mesh_quality_improver.h"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

void IncrementalMeshBoundarySnapper::snap()
{
  m_snapper.snap(m_mesh1, m_mesh2);

  record_final_position_and_restore_initial_position(m_mesh1, m_mesh1Data);
  record_final_position_and_restore_initial_position(m_mesh2, m_mesh2Data);

  std::cout << "applying displacement to mesh1" << std::endl;
  apply_displacement(m_mesh1, m_mesh1Data, m_improver1, "mesh2");

  std::cout << "\napplying displacement to mesh2" << std::endl;
  apply_displacement(m_mesh2, m_mesh2Data, m_improver2, "mesh2");
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
                                                        const std::string& prefix)
{
  auto& data = *field;
  for (int step = 0; step < m_nsteps; ++step)
  {
    std::cout << "\nApplying mesh snap step " << step << " / " << m_nsteps << std::endl;
    for (auto& vert : mesh->get_vertices())
      if (vert)
      {
        utils::Point displacement = (data(vert, 0, 0).ptDest - data(vert, 0, 0).ptOrig) / m_nsteps;
        auto pt                   = vert->get_point_orig(0);
        vert->set_point_orig(0, pt + displacement);
      }

    // std::cout << "before smoothing, number of invalid points = " << improver->countInvalidPoints() << std::endl;

    improver->run();
    // std::cout << "after smoothing, number of invalid points = " << improver->countInvalidPoints() << std::endl;
  }
}

std::shared_ptr<IncrementalMeshBoundarySnapper>
make_incremental_boundary_snapper(MeshPtr mesh1, MeshPtr mesh2, const IncrementalBoundarySnapperOpts& opts)
{
  BoundaryFixture filter1(mesh1);
  auto fixer1 = make_standard_improver(mesh1, filter1, opts.qualityImproverOpts);

  BoundaryFixture filter2(mesh2);
  auto fixer2 = make_standard_improver(mesh2, filter2, opts.qualityImproverOpts);

  return std::make_shared<IncrementalMeshBoundarySnapper>(mesh1, fixer1, mesh2, fixer2, opts.boundarySnapNsteps);
}

} // namespace impl
} // namespace mesh
} // namespace middle_mesh
} // namespace stk
