#ifndef STK_MIDDLE_MESH_APPLICATION_INTERFACE_PARALLEL_H
#define STK_MIDDLE_MESH_APPLICATION_INTERFACE_PARALLEL_H

#include "application_interface.hpp"
#include "mesh.hpp"
#include "mesh_entity.hpp"
#include "mesh_scatter.hpp"
#include "mesh_scatter_spec.hpp"
#include "stk_middle_mesh/variable_size_field.hpp"
#include "stk_middle_mesh/mesh_relational_data.hpp"
#include "stk_middle_mesh/bounding_box_search.hpp"


namespace stk {
namespace middle_mesh {
namespace impl {

class ApplicationInterfaceParallel : public ApplicationInterface
{
  public:
    ApplicationInterfaceParallel(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2,
                                 MPI_Comm unionComm, const ParallelSearchOpts& parallelSearchOpts,
                                 const VolumeSnapOpts& volumeSnapOpts,
                                 const BoundarySnapAndQualityImprovementOpts& boundarySnapOpts,
                                 const MiddleGridOpts& middleGridOpts,
                                 std::shared_ptr<XiCoordinates> xiPts);

    // creates the middle grid
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

    using BoundingBoxSearch = search::ElementToElementBoundingBoxSearch;

    void check_inputs();

    void check_meshes_exist();

    void check_mesh_comm_size(std::shared_ptr<mesh::Mesh> mesh, const std::string& name);

    void check_mesh_element_count(std::shared_ptr<mesh::Mesh> mesh, const std::string& name);

    void do_volume_snap();

    void do_boundary_snap();

    void create_scatter_spec(std::shared_ptr<mesh::Mesh> sendMesh, std::shared_ptr<mesh::Mesh> recvMesh,
                             std::shared_ptr<mesh::impl::MeshScatterSpec> scatterSpecSendToRecv,
                             std::shared_ptr<mesh::impl::MeshScatterSpec> scatterSpecRecvToSend);

    void invert_coarse_search_result(std::shared_ptr<mesh::Mesh> sendMesh,
                                     std::shared_ptr<mesh::Mesh> recvMesh,
                                     const BoundingBoxSearch::EntityProcRelationVec& meshRecvToSendRelations,
                                     std::shared_ptr<mesh::impl::MeshScatterSpec> scatterSpecRecvToSend);

    int translate_comm_rank(MPI_Comm inputComm, MPI_Comm outputComm, int rankOnInputComm);  

    void scatter_mesh_1to2(std::shared_ptr<mesh::impl::MeshScatterSpec> scatterSpec);

    void scatter_mesh_2to1(std::shared_ptr<mesh::impl::MeshScatterSpec> scatterSpec);

    void do_mesh_projections();

    void scatter_mesh_relational_data_2to1();

    void create_middle_mesh_verts_and_edges();

    void create_middle_mesh_triangles();

    void project_xi_points_onto_input_meshes();

    void apply_geometry_improvers();

    void scatter_meshin_1to2();

    void create_meshin_remote_info(const mesh::impl::MeshScatter& scatterer);

    mesh::FieldPtr<mesh::RemoteSharedEntity> extract_element_remote_info(std::shared_ptr<mesh::Mesh> meshIn,
                                                                         mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> entityRemoteInfoPtr);

    void scatter_meshin_fields_1to2();


    bool m_middleMeshCreated = false;
    std::shared_ptr<mesh::Mesh> m_mesh1;
    std::shared_ptr<mesh::Mesh> m_mesh2;
    MPI_Comm m_unionComm;
    ParallelSearchOpts m_parallelSearchOpts;
    VolumeSnapOpts m_volumeSnapOpts;
    BoundarySnapAndQualityImprovementOpts m_boundarySnapOpts;
    MiddleGridOpts m_middleGridOpts;
    std::shared_ptr<XiCoordinates> m_xiPts;   

    std::shared_ptr<mesh::Mesh> m_mesh1ScatteredToMesh2;
    std::shared_ptr<mesh::Mesh> m_mesh2ScatteredToMesh1;

    //TODO: see if we need to keep these around
    mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> m_mesh1EntityOrigins;
    mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> m_mesh1EntityDestinations;
    mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> m_mesh2EntityOrigins;
    mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> m_mesh2EntityDestinations;

    std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> m_pointClassifierOnMesh2;
    std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> m_pointClassifierOnMesh1;
    std::shared_ptr<nonconformal4::impl::MeshRelationalData> m_meshRelationalDataOnMesh2;
    std::shared_ptr<nonconformal4::impl::MeshRelationalData> m_meshRelationalDataOnMesh1;
    std::shared_ptr<mesh::Mesh> m_meshInOnMesh1Procs;
    std::shared_ptr<mesh::Mesh> m_meshInOnMesh2Procs;


    mesh::FieldPtr<mesh::MeshEntityPtr> m_mesh1Classification;
    mesh::FieldPtr<mesh::MeshEntityPtr> m_mesh2Classification;
    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> m_mesh1InverseClassification;
    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> m_mesh2InverseClassification;
    mesh::FieldPtr<mesh::RemoteSharedEntity> m_meshInRemoteInfo1to2;
    mesh::FieldPtr<mesh::RemoteSharedEntity> m_meshInRemoteInfo2to1;
    mesh::FieldPtr<utils::Point> m_xiPtsProjectedOntoMesh1;
    mesh::FieldPtr<utils::Point> m_xiPtsProjectedOntoMesh2ScatteredToMesh1; 
    mesh::FieldPtr<utils::Point> m_xiPtsProjectedOntoMesh2;
};

}
}
}

#endif