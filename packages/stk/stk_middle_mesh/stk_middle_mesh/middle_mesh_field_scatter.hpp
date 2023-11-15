#ifndef STK_MIDDLE_MESH_MIDDLE_MESH_FIELD_SCATTER
#define STK_MIDDLE_MESH_MIDDLE_MESH_FIELD_SCATTER

#include "application_interface.hpp"
#include "variable_size_field.hpp"
#include "stk_util/parallel/DataExchangeKnownPatternNonBlockingBuffer.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

class MiddleMeshFieldScatter
{
  public:
    MiddleMeshFieldScatter(std::shared_ptr<mesh::Mesh> meshInOnMesh1Procs,
                           std::shared_ptr<mesh::Mesh> meshInOnMesh2Procs,
                           std::shared_ptr<mesh::Mesh> mesh2,
                           std::shared_ptr<XiCoordinates> xiPts,
                           mesh::FieldPtr<mesh::MeshEntityPtr> meshInClassificationOnMesh2ScatteredToMesh1,
                           mesh::FieldPtr<utils::Point> xiPtsProjectedOntoMesh2ScatteredToMesh1,
                           mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> mesh2EntityOrigins,
                           mesh::FieldPtr<mesh::RemoteSharedEntity> meshInRemoteInfo1to2,
                           mesh::FieldPtr<mesh::RemoteSharedEntity> meshInRemoteInfo2to1,
                           MPI_Comm unionComm);

    void scatter();

    mesh::FieldPtr<mesh::MeshEntityPtr> get_meshin_classification_on_mesh2()
    {
      return m_meshInClassificationOnMesh2;
    }

    mesh::FieldPtr<utils::Point> get_xi_pts_projected_onto_mesh2() const
    {
      return m_xiPtsProjectedOntoMesh2;
    }

  private:

    using Exchanger = stk::DataExchangeKnownPatternNonBlockingCommBuffer;

    void pack_send_bufs(Exchanger& exchanger);

    void size_recv_bufs(Exchanger& exchanger);

    void unpack_recv_bufs(Exchanger& exchanger);

    std::shared_ptr<mesh::Mesh> m_meshInOnMesh1Procs;
    std::shared_ptr<mesh::Mesh> m_meshInOnMesh2Procs;
    std::shared_ptr<mesh::Mesh> m_mesh2;
    std::shared_ptr<XiCoordinates> m_xiPts;   

    mesh::FieldPtr<mesh::MeshEntityPtr> m_meshInClassificationOnMesh2ScatteredToMesh1;
    mesh::FieldPtr<utils::Point> m_xiPtsProjectedOntoMesh2ScatteredToMesh1;

    mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> m_mesh2EntityOrigins;
    mesh::FieldPtr<mesh::RemoteSharedEntity> m_meshInRemoteInfo1to2;
    mesh::FieldPtr<mesh::RemoteSharedEntity> m_meshInRemoteInfo2to1;
    int m_numXiPts = 0;

    mesh::FieldPtr<mesh::MeshEntityPtr> m_meshInClassificationOnMesh2;
    mesh::FieldPtr<utils::Point> m_xiPtsProjectedOntoMesh2;

    MPI_Comm m_unionComm;
};

}
}
}
}

#endif