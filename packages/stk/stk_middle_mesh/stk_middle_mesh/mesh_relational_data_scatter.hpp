#ifndef STK_MIDDLE_MESH_MESH_RELATIONAL_DATA_SCATTER_H
#define STK_MIDDLE_MESH_MESH_RELATIONAL_DATA_SCATTER_H

#include "mesh_entity.hpp"
#include "mesh_relational_data.hpp"
#include "predicates/intersection_common.hpp"
#include "predicates/quad_to_triangles.hpp"
#include "stk_middle_mesh/mesh.hpp"
#include "stk_middle_mesh/variable_size_field.hpp"
#include "stk_middle_mesh/entity_sorted_by_owner.hpp"
#include "stk_util/parallel/CommBuffer.hpp"
#include "stk_util/parallel/DataExchangeKnownPatternNonBlockingBuffer.hpp"
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlockingBuffer.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_middle_mesh/utils.hpp"


namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

struct MeshRelationalDataScatterInput
{
  std::shared_ptr<mesh::Mesh> mesh1;
  std::shared_ptr<mesh::Mesh> mesh2;
  std::shared_ptr<mesh::Mesh> mesh1ScatteredToMesh2;
  std::shared_ptr<mesh::Mesh> mesh2ScatteredToMesh1;
  mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> mesh1EntityOrigins;
  mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> mesh1EntityDestinations;
  mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> mesh2EntityOrigins;
  mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> mesh2EntityDestinations;
};


class MeshRelationalDataScatter
{
  public:
    // given a MeshRelationalData, whose fields are defined on mesh2 and mesh1ScatteredToMesh2, send the fields to
    // mesh1 and mesh2ScatteredToMesh1
    // mesh1EntityOrigins maps mesh1ScatteredToMesh2 to mesh1
    // mesh2EntityDestinations maps mesh2 to mesh2Scattered to mesh1
    MeshRelationalDataScatter(const MeshRelationalDataScatterInput& scatterInput,     
                              std::shared_ptr<MeshRelationalData> meshRelationalData,
                              std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> pointClassifierInput,
                              std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> pointClassifierOutput,           
                              MPI_Comm unionComm);

    std::shared_ptr<mesh::Mesh> get_middle_mesh() { return m_meshIn; }

    std::shared_ptr<MeshRelationalData> scatter();

  private:
    using Exchanger        = stk::DataExchangeKnownPatternNonBlockingCommBuffer;
    using ExchangerUnknown = stk::DataExchangeUnknownPatternNonBlockingCommBuffer;


    void pack_mesh1_vert_fields(Exchanger& exchanger);

    void size_recv_buffer_mesh1_vert_fields(Exchanger& exchanger);

    void unpack_mesh1_vert_fields(Exchanger& exchanger);

    void pack_mesh2_vert_fields(ExchangerUnknown& exchanger);

    void pack(stk::CommBuffer& buf, int destRank, const predicates::impl::PointRecord& record);

    void pack(stk::CommBuffer& buf, const predicates::impl::PointRecordForTriangle& record);

    int get_id_on_dest(mesh::MeshEntityPtr el1, int destRank);

    predicates::impl::PointRecord unpack_point_record(stk::CommBuffer& buf);

    predicates::impl::PointRecordForTriangle unpack_point_record_for_triangle(stk::CommBuffer& buf, mesh::MeshEntityPtr parentEl);

    void unpack_mesh2_vert_fields(ExchangerUnknown& exchanger);

    void pack_mesh1_edge_fields(ExchangerUnknown& exchanger);

    void get_owned_splits(mesh::MeshEntityPtr edge, std::vector<int>& ownedSplits);

    void pack_edge_splits(stk::CommBuffer& buf, const mesh::RemoteSharedEntity& dest,
                          mesh::MeshEntityPtr edge, const std::vector<int> ownedSplits);

    void unpack_mesh1_edge_fields(ExchangerUnknown& exchanger);

    void unpack_edge_split(stk::CommBuffer& buf, int senderRank);

    void pack_mesh2_edge_fields(ExchangerUnknown& exchanger);

    void unpack_mesh2_edge_fields(ExchangerUnknown& exchanger);

    void unpack_vert_on_edge(stk::CommBuffer& buf, int senderRank);

    void check_fields_are_on_correct_meshes();
    
    std::shared_ptr<mesh::Mesh> m_mesh1;
    std::shared_ptr<mesh::Mesh> m_mesh2;
    std::shared_ptr<mesh::Mesh> m_mesh1ScatteredToMesh2;
    std::shared_ptr<mesh::Mesh> m_mesh2ScatteredToMesh1;
    mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> m_mesh1EntityOrigins;
    mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> m_mesh1EntityDestinations;
    mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> m_mesh2EntityOrigins;
    mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> m_mesh2EntityDestinations;
    std::shared_ptr<mesh::Mesh> m_meshIn;
    std::shared_ptr<MeshRelationalData> m_meshRelationDataInput;
    std::shared_ptr<MeshRelationalData> m_meshRelationDataOutput;
    std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> m_pointClassifierInput;
    std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> m_pointClassifierOutput;
    predicates::impl::QuadToTriangles* m_quadToTrianglesOrigin;
    predicates::impl::QuadToTriangles* m_quadToTrianglesDest;

    FakeVertGenerator m_fakeVertGenerator;
    std::vector<std::map<int, FakeVert>> m_senderFakeVertsToLocal;
    MPI_Comm m_unionComm;

};

}
}
}
}

#endif