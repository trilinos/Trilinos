#include "middle_mesh_field_scatter.hpp"
#include "stk_middle_mesh/utils.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

MiddleMeshFieldScatter::MiddleMeshFieldScatter(std::shared_ptr<mesh::Mesh> meshInOnMesh1Procs,
                        std::shared_ptr<mesh::Mesh> meshInOnMesh2Procs,
                        std::shared_ptr<mesh::Mesh> mesh2,
                        std::shared_ptr<XiCoordinates> xiPts,
                        mesh::FieldPtr<mesh::MeshEntityPtr> meshInClassificationOnMesh2ScatteredToMesh1,
                        mesh::FieldPtr<utils::Point> xiPtsProjectedOntoMesh2ScatteredToMesh1,
                        mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> mesh2EntityOrigins,
                        mesh::FieldPtr<mesh::RemoteSharedEntity> meshInRemoteInfo1to2,
                        mesh::FieldPtr<mesh::RemoteSharedEntity> meshInRemoteInfo2to1,
                        MPI_Comm unionComm) :
  m_meshInOnMesh1Procs(meshInOnMesh1Procs),
  m_meshInOnMesh2Procs(meshInOnMesh2Procs),
  m_mesh2(mesh2),
  m_xiPts(xiPts),
  m_meshInClassificationOnMesh2ScatteredToMesh1(meshInClassificationOnMesh2ScatteredToMesh1),
  m_xiPtsProjectedOntoMesh2ScatteredToMesh1(xiPtsProjectedOntoMesh2ScatteredToMesh1),
  m_mesh2EntityOrigins(mesh2EntityOrigins),
  m_meshInRemoteInfo1to2(meshInRemoteInfo1to2),
  m_meshInRemoteInfo2to1(meshInRemoteInfo2to1),
  m_unionComm(unionComm)
{
  if (m_xiPts)
  {
    m_numXiPts = std::max(m_xiPts->get_xi_coords(mesh::MeshEntityType::Triangle).size(),
                    m_xiPts->get_xi_coords(mesh::MeshEntityType::Quad).size());
  }
}

void MiddleMeshFieldScatter::scatter()
{
  Exchanger exchanger(m_unionComm);

  pack_send_bufs(exchanger);
  size_recv_bufs(exchanger);
  exchanger.allocate_send_buffers();
  exchanger.allocate_recv_buffers();
  pack_send_bufs(exchanger);


  unpack_recv_bufs(exchanger);
}


void MiddleMeshFieldScatter::pack_send_bufs(Exchanger& exchanger)
{
  if (m_meshInOnMesh1Procs)
  {
    auto& meshInToMesh2Els         = *m_meshInClassificationOnMesh2ScatteredToMesh1;
    auto& mesh2EntityOrigins       = *m_mesh2EntityOrigins;  
    auto& meshInEntityDestinations = *m_meshInRemoteInfo1to2;
    for (auto& elIn : m_meshInOnMesh1Procs->get_elements())
    {
      if (elIn)
      {
        mesh::RemoteSharedEntity elInDestination = meshInEntityDestinations(elIn, 0, 0);
        mesh::MeshEntityPtr el2ScatteredToMesh1  = meshInToMesh2Els(elIn, 0, 0);
        mesh::RemoteSharedEntity el2Origin       = mesh2EntityOrigins(el2ScatteredToMesh1, 0, 0); 

        assert(elInDestination.remoteRank == el2Origin.remoteRank);
        int remoteRank = elInDestination.remoteRank;
        exchanger.get_send_buf(remoteRank).pack(elInDestination.remoteId);
        exchanger.get_send_buf(remoteRank).pack(el2Origin.remoteId);
        if (m_numXiPts > 0)
        {
          auto& xiPtsProjectedOntoMesh2ScatteredToMesh1 = *m_xiPtsProjectedOntoMesh2ScatteredToMesh1;
          for (int i=0; i < m_numXiPts; ++i)
          {
            exchanger.get_send_buf(remoteRank).pack(xiPtsProjectedOntoMesh2ScatteredToMesh1(elIn, i, 0));
          }
        }
      }
    }
  }
}

void MiddleMeshFieldScatter::size_recv_bufs(Exchanger& exchanger)
{
  if (m_meshInOnMesh2Procs)
  {
    auto& meshInRemoteInfo2to1 = *m_meshInRemoteInfo2to1;
    for (auto& elIn : m_meshInOnMesh2Procs->get_elements())
    {
      if (elIn)
      {
        mesh::RemoteSharedEntity elInOnMesh1 = meshInRemoteInfo2to1(elIn, 0, 0);
        int idVal = -1;
        utils::Point xiCoords;
        exchanger.get_recv_buf(elInOnMesh1.remoteRank).pack(idVal);
        exchanger.get_recv_buf(elInOnMesh1.remoteRank).pack(idVal);
        for (int i=0; i < m_numXiPts; ++i)
        {
          exchanger.get_recv_buf(elInOnMesh1.remoteRank).pack(xiCoords);
        }
      }
    }
  }

  for (int rank=0; rank < utils::impl::comm_size(m_unionComm); ++rank)
  {
    auto& buf = exchanger.get_recv_buf(rank);
    exchanger.set_recv_buffer_size(rank, buf.size());
  }
}

void MiddleMeshFieldScatter::unpack_recv_bufs(Exchanger& exchanger)
{
  if (m_meshInOnMesh2Procs)
  {
    m_meshInClassificationOnMesh2 = mesh::create_field<mesh::MeshEntityPtr>(m_meshInOnMesh2Procs, mesh::FieldShape(0, 0, 1), 1, nullptr);
    if (m_xiPts)
    {
      m_xiPtsProjectedOntoMesh2 = mesh::create_field<utils::Point>(m_meshInOnMesh2Procs, mesh::FieldShape(0, 0, m_numXiPts), 1);
    }
  }

  exchanger.start_nonblocking();

  auto unpacker = [&](int /*rank*/, stk::CommBuffer& buf)
  {
    auto& mesh2Classification = *m_meshInClassificationOnMesh2;

    while (buf.remaining() > 0)
    {
      int elInId, el2Id;
      buf.unpack(elInId);
      buf.unpack(el2Id);
      mesh::MeshEntityPtr elIn = m_meshInOnMesh2Procs->get_elements()[elInId];
      mesh::MeshEntityPtr el2  = m_mesh2->get_elements()[el2Id];
      mesh2Classification(elIn, 0, 0) = el2;

      for (int i=0; i < m_numXiPts; ++i)
      {
        auto& xiPtsProjectedOntoMesh2 = *m_xiPtsProjectedOntoMesh2;
        utils::Point xiPt;
        buf.unpack(xiPt);
        xiPtsProjectedOntoMesh2(elIn, i, 0) = xiPt;
      }
    }
  };

  exchanger.complete_receives(unpacker);   
}


}
}
}
}