#ifndef APPLICATION_INTERFACE_FAKE_PARALLEL_H
#define APPLICATION_INTERFACE_FAKE_PARALLEL_H

#include "application_interface.hpp"
#include "field.hpp"
#include "mesh.hpp"
#include "variable_size_field.hpp"

#include "utils.hpp"  //TODO: move to src file

namespace stk {
namespace middle_mesh {
namespace impl {

class ApplicationInterfaceFakeParallelImpl : public ApplicationInterface
{
  public:
    // union_comm is a MPI_Comm that contains at least the union of all the processes on
    // the mesh1 and mesh2 comms.  Having other processes is allowed.
    ApplicationInterfaceFakeParallelImpl(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2,
                                         MPI_Comm unionComm, const ParallelSearchOpts& parallelSearchOpts,
                                         const VolumeSnapOpts& volumeSnapOpts,
                                         const BoundarySnapAndQualityImprovementOpts& boundarySnapOpts,
                                         const MiddleGridOpts& middleGridOpts,
                                         std::shared_ptr<XiCoordinates> xiPts)
      : m_mesh1Parallel(mesh1)
      , m_mesh2Parallel(mesh2)
      , m_unionComm(unionComm)
      , m_parallelSearchOpts(parallelSearchOpts)
      , m_volumeSnapOpts(volumeSnapOpts)
      , m_boundarySnapOpts(boundarySnapOpts)
      , m_middleGridOpts(middleGridOpts)
      , m_xiPts(xiPts)
      , m_rootRankOnUnionComm(decide_root_rank(unionComm, mesh1, mesh2))
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

    void create_middle_grid() override;

    std::shared_ptr<mesh::Mesh> get_middle_grid_for_mesh1() override;

    std::shared_ptr<mesh::Mesh> get_middle_grid_for_mesh2() override;

    mesh::FieldPtr<mesh::MeshEntityPtr> get_mesh1_classification() override;

    mesh::FieldPtr<mesh::MeshEntityPtr> get_mesh2_classification() override;

    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> compute_mesh1_inverse_classification() override;

    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> compute_mesh2_inverse_classification() override;

    mesh::FieldPtr<mesh::RemoteSharedEntity> get_remote_info_mesh_one_to_two() override;

    mesh::FieldPtr<mesh::RemoteSharedEntity> get_remote_info_mesh_two_to_one() override;

    mesh::FieldPtr<utils::Point> get_xi_points_on_mesh1() override;

    mesh::FieldPtr<utils::Point> get_xi_points_on_mesh2() override;

  private:
    void check_union_comm_size();

    void check_both_meshes_exist_somewhere();

    void gather_meshes_to_root();

    void do_volume_snap();

    void do_boundary_snap();

    void create_serial_middle_grid();

    void scatter_meshes();

    void scatter_remote_info();

    void scatter_xi_points();

    void scatter_xi_points(std::shared_ptr<mesh::Mesh> middleMeshParallel,
                           mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> entityDestinations,
                           mesh::FieldPtr<utils::Point> xiPointsSerial,
                           mesh::FieldPtr<utils::Point>& xiPointsParallel);

    void scatter_mesh(mesh::FieldPtr<mesh::MeshEntityPtr> meshClassificationSerial,
                      mesh::FieldPtr<mesh::RemoteSharedEntity> elementOrigins,                      
                      std::shared_ptr<mesh::Mesh> inputMeshParallel, 
                      mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity>& middleGridSerialEntityDestinations,
                      std::shared_ptr<mesh::Mesh>& middleGridParallel,
                      mesh::FieldPtr<mesh::MeshEntityPtr>& meshClassificationParallel);

    mesh::FieldPtr<int> get_element_destinations(std::shared_ptr<mesh::Mesh> middleGridSerial,
                                                 mesh::FieldPtr<mesh::MeshEntityPtr> meshClassification,
                                                 mesh::FieldPtr<mesh::RemoteSharedEntity> elementOrigins);

    int decide_root_rank(MPI_Comm unionComm, std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2);

    std::shared_ptr<mesh::Mesh> m_mesh1Parallel;
    std::shared_ptr<mesh::Mesh> m_mesh2Parallel;
    MPI_Comm m_unionComm;
    ParallelSearchOpts m_parallelSearchOpts;
    VolumeSnapOpts m_volumeSnapOpts;
    BoundarySnapAndQualityImprovementOpts m_boundarySnapOpts;
    MiddleGridOpts m_middleGridOpts;
    std::shared_ptr<XiCoordinates> m_xiPts; 

    std::shared_ptr<mesh::Mesh> m_mesh1Serial;
    std::shared_ptr<mesh::Mesh> m_mesh2Serial;
    mesh::FieldPtr<mesh::RemoteSharedEntity> m_mesh1ElementOrigins;
    mesh::FieldPtr<mesh::RemoteSharedEntity> m_mesh2ElementOrigins;

    std::shared_ptr<mesh::Mesh> m_middleGridSerial;
    mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> m_middleGridSerialEntityDestinations1;
    mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> m_middleGridSerialEntityDestinations2;
    mesh::FieldPtr<mesh::MeshEntityPtr> m_mesh1ClassificationSerial;
    mesh::FieldPtr<mesh::MeshEntityPtr> m_mesh2ClassificationSerial;
    mesh::FieldPtr<utils::Point> m_mesh1XiPointsSerial;
    mesh::FieldPtr<utils::Point> m_mesh2XiPointsSerial;

    std::shared_ptr<mesh::Mesh> m_middleGridParallel1;
    std::shared_ptr<mesh::Mesh> m_middleGridParallel2;
    mesh::FieldPtr<mesh::MeshEntityPtr> m_mesh1ClassificationParallel;
    mesh::FieldPtr<mesh::MeshEntityPtr> m_mesh2ClassificationParallel;
    mesh::FieldPtr<mesh::RemoteSharedEntity> m_remoteInfoMeshOneToTwo;
    mesh::FieldPtr<mesh::RemoteSharedEntity> m_remoteInfoMeshTwoToOne;
    mesh::FieldPtr<utils::Point> m_mesh1XiPointsParallel;
    mesh::FieldPtr<utils::Point> m_mesh2XiPointsParallel;
    bool m_middleGridCreated = false;

    const int m_rootRankOnUnionComm = 0;
};

std::pair<int, int> get_comm_sizes_on_root(MPI_Comm comm1, MPI_Comm comm2, MPI_Comm unionComm, int rootRank);

} // namespace impl

} // namespace middle_mesh
} // namespace stk

#endif