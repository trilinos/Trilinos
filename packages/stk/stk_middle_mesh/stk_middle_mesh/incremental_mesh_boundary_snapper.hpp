#ifndef INCREMENTAL_MESH_BOUNDARY_SNAPPER_H
#define INCREMENTAL_MESH_BOUNDARY_SNAPPER_H

#include "boundary_fixture.hpp"
#include "create_mesh_quality_improver.hpp"
#include "incremental_mesh_boundary_snapper_opts.hpp"
#include "mesh.hpp"
#include "mesh_boundary_snapper.hpp"
#include "mesh_quality_improver.hpp"

#include "mesh_quality_statistics.hpp" //TODO: DEBUGGING
#include "patch_distortion_objective.hpp"
#include "regularized_distortion_metric.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

struct SnapData
{
    utils::Point ptOrig;
    utils::Point ptDest;
};

using MeshPtr = std::shared_ptr<Mesh>;

// TODO: it would be better if the size of each step was determined by the minimum edge length connected to the vert
//       Would need to void stalling when edges get very small
class IncrementalMeshBoundarySnapper
{
  public:
    IncrementalMeshBoundarySnapper(MeshPtr mesh1, std::shared_ptr<MeshQualityImprover> improver1, MeshPtr mesh2,
                                   std::shared_ptr<MeshQualityImprover> improver2, int nsteps)
      : m_mesh1(mesh1)
      , m_mesh2(mesh2)
      , m_mesh1Data(create_field<impl::SnapData>(mesh1, FieldShape(1, 0, 0), 1))
      , m_mesh2Data(create_field<impl::SnapData>(mesh2, FieldShape(1, 0, 0), 1))
      , m_improver1(improver1)
      , m_improver2(improver2)
      , m_nsteps(nsteps)
    {
      record_initial_position(mesh1, m_mesh1Data);
      record_initial_position(mesh2, m_mesh2Data);
    }

    void snap();

    std::shared_ptr<MeshQualityImprover> get_mesh1_quality_improver() { return m_improver1; }

    std::shared_ptr<MeshQualityImprover> get_mesh2_quality_improver() { return m_improver2; }

  private:
    using FieldPtr = std::shared_ptr<Field<impl::SnapData>>;

    void record_initial_position(MeshPtr mesh, FieldPtr field);

    void record_final_position_and_restore_initial_position(MeshPtr mesh, FieldPtr field);

    void apply_displacement(MeshPtr mesh, FieldPtr field, std::shared_ptr<MeshQualityImprover> improver,
                            const std::string& prefix);

    MeshPtr m_mesh1;
    MeshPtr m_mesh2;
    MeshBoundarySnapper m_snapper;
    FieldPtr m_mesh1Data;
    FieldPtr m_mesh2Data;
    std::shared_ptr<MeshQualityImprover> m_improver1;
    std::shared_ptr<MeshQualityImprover> m_improver2;
    int m_nsteps;
};

std::shared_ptr<IncrementalMeshBoundarySnapper>
make_incremental_boundary_snapper(MeshPtr mesh1, MeshPtr mesh2,
                                  const IncrementalBoundarySnapperOpts& opts = IncrementalBoundarySnapperOpts());

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
