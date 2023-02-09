#ifndef APPLICATION_INTERFACE_FAKE_PARALLEL_H
#define APPLICATION_INTERFACE_FAKE_PARALLEL_H

#include "application_interface.h"
#include "field.h"
#include "mesh.h"
#include "variable_size_field.h"

namespace stk {
namespace middle_mesh {
namespace impl {

class ApplicationInterfaceFakeParallelImpl
{
  public:
    // union_comm is a MPI_Comm that contains at least the union of all the processes on
    // the mesh1 and mesh2 comms.  Having other processes is allowed.
    ApplicationInterfaceFakeParallelImpl(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2,
                                         MPI_Comm unionComm, const ParallelSearchOpts& parallelSearchOpts,
                                         const VolumeSnapOpts& volumeSnapOpts,
                                         const BoundarySnapAndQualityImprovementOpts& boundarySnapOpts,
                                         const MiddleGridOpts& middleGridOpts)
      : m_mesh1Parallel(mesh1)
      , m_mesh2Parallel(mesh2)
      , m_unionComm(unionComm)
      , m_parallelSearchOpts(parallelSearchOpts)
      , m_volumeSnapOpts(volumeSnapOpts)
      , m_boundarySnapOpts(boundarySnapOpts)
      , m_middleGridOpts(middleGridOpts)
    {
      if (unionComm == MPI_COMM_NULL)
        throw std::runtime_error("union communicator cannot be null");

      check_union_comm_size();
      check_both_meshes_exist_somewhere();

#ifndef NDEBUG
      if (mesh1)
        mesh::check_topology(mesh1);

      if (mesh2)
        mesh::check_topology(mesh2);
#endif
    }

    void create_middle_grid();

    std::shared_ptr<mesh::Mesh> get_middle_grid_for_mesh1();

    std::shared_ptr<mesh::Mesh> get_middle_grid_for_mesh2();

    mesh::FieldPtr<mesh::MeshEntityPtr> get_mesh1_classification();

    mesh::FieldPtr<mesh::MeshEntityPtr> get_mesh2_classification();

    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> compute_mesh1_inverse_classification();

    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> compute_mesh2_inverse_classification();

  private:
    void check_union_comm_size();

    void check_both_meshes_exist_somewhere();

    void gather_meshes_to_root();

    void do_volume_snap();

    void do_boundary_snap();

    void create_serial_middle_grid();

    void scatter_meshes();

    void scatter_mesh(mesh::FieldPtr<mesh::MeshEntityPtr> meshClassificationSerial,
                      mesh::FieldPtr<mesh::RemoteSharedEntity> elementOrigins,
                      std::shared_ptr<mesh::Mesh> inputMeshParallel, std::shared_ptr<mesh::Mesh>& middleGridParallel,
                      mesh::FieldPtr<mesh::MeshEntityPtr>& meshClassificationParallel);

    mesh::FieldPtr<int> get_element_destinations(std::shared_ptr<mesh::Mesh> middleGridSerial,
                                                 mesh::FieldPtr<mesh::MeshEntityPtr> meshClassification,
                                                 mesh::FieldPtr<mesh::RemoteSharedEntity> elementOrigins);

    std::shared_ptr<mesh::Mesh> m_mesh1Parallel;
    std::shared_ptr<mesh::Mesh> m_mesh2Parallel;
    MPI_Comm m_unionComm;
    ParallelSearchOpts m_parallelSearchOpts;
    VolumeSnapOpts m_volumeSnapOpts;
    BoundarySnapAndQualityImprovementOpts m_boundarySnapOpts;
    MiddleGridOpts m_middleGridOpts;

    std::shared_ptr<mesh::Mesh> m_mesh1Serial;
    std::shared_ptr<mesh::Mesh> m_mesh2Serial;
    mesh::FieldPtr<mesh::RemoteSharedEntity> m_mesh1ElementOrigins;
    mesh::FieldPtr<mesh::RemoteSharedEntity> m_mesh2ElementOrigins;

    std::shared_ptr<mesh::Mesh> m_middleGridSerial;
    mesh::FieldPtr<mesh::MeshEntityPtr> m_mesh1ClassificationSerial;
    mesh::FieldPtr<mesh::MeshEntityPtr> m_mesh2ClassificationSerial;

    std::shared_ptr<mesh::Mesh> m_middleGridParallel1;
    std::shared_ptr<mesh::Mesh> m_middleGridParallel2;
    mesh::FieldPtr<mesh::MeshEntityPtr> m_mesh1ClassificationParallel;
    mesh::FieldPtr<mesh::MeshEntityPtr> m_mesh2ClassificationParallel;
};

std::pair<int, int> get_comm_sizes_on_root(MPI_Comm comm1, MPI_Comm comm2, MPI_Comm unionComm, int rootRank);

} // namespace impl

} // namespace middle_mesh
} // namespace stk

#endif