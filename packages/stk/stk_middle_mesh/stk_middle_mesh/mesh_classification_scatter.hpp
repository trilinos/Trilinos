#ifndef MESH_CLASSIFICATION_SCATTER_H
#define MESH_CLASSIFICATION_SCATTER_H

#include "field.hpp"
#include "variable_size_field.hpp"
#include "mesh.hpp"
#include "mesh_entity.hpp"
#include "parallel_exchange.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

class MeshClassificationScatter
{
  public:
    MeshClassificationScatter(MPI_Comm inputComm, 
                              std::shared_ptr<Mesh> middleGridSerial,
                              VariableSizeFieldPtr<RemoteSharedEntity> middleGridSerialEntityDestinations,
                              std::shared_ptr<Mesh> middleGridParallel,
                              std::shared_ptr<Mesh> inputMeshParallel,
                              FieldPtr<RemoteSharedEntity> inputMeshElementOrigins,
                              FieldPtr<MeshEntityPtr> middleGridClassificationSerial)
      : m_inputComm(inputComm)
      , m_middleGridSerial(middleGridSerial)
      , m_middleGridSerialEntityDestinations(middleGridSerialEntityDestinations)
      , m_middleGridParallel(middleGridParallel)
      , m_inputMeshParallel(inputMeshParallel)
      , m_inputMeshElementOrigins(inputMeshElementOrigins)
      , m_middleGridClassificationSerial(middleGridClassificationSerial)
    {
      m_amIRootRank   = middleGridSerial != nullptr;
      m_amIOnMeshComm = middleGridParallel != nullptr;

      m_rootRankOnInputComm = get_root_rank_on_input_comm();
      if (m_amIOnMeshComm)
      {
        m_middleGridClassificationParallel =
            create_field<MeshEntityPtr>(m_middleGridParallel, FieldShape(0, 0, 1), 1, nullptr);
      }
    }

    FieldPtr<MeshEntityPtr> scatter()
    {
      send_classification_info();
      return m_middleGridClassificationParallel;
    }

  private:
    struct LocalIdPair
    {
      int middleMeshElementLocalId;
      int inputMeshElementLocalId;
    };

    int get_root_rank_on_input_comm()
    {
      int dataInput  = m_amIRootRank ? utils::impl::comm_rank(m_inputComm) : 0;
      int dataOutput = -1;
      MPI_Allreduce(&dataInput, &dataOutput, 1, MPI_INT, MPI_SUM, m_inputComm);

      return dataOutput;
    }

    void send_classification_info()
    {
      utils::impl::ParallelExchange<LocalIdPair> exchanger(m_inputComm, 81);
      if (m_amIRootRank)
        pack_data(exchanger);

      if (m_amIOnMeshComm)
        exchanger.set_recv_buffer_size(m_rootRankOnInputComm, count_valid(m_middleGridParallel->get_elements()));

      exchanger.start_sends();
      exchanger.start_recvs();
      exchanger.complete_recvs();

      if (m_amIOnMeshComm)
        unpack_data(exchanger);

      exchanger.complete_sends();
    }

    void pack_data(utils::impl::ParallelExchange<LocalIdPair>& exchanger)
    {
      assert(m_amIRootRank);

      auto& middleGridClassificationSerial = *m_middleGridClassificationSerial;
      auto& middleGridEntityDestinations = *m_middleGridSerialEntityDestinations;
      auto& inputMeshElementOrigins        = *m_inputMeshElementOrigins;
      for (auto& el : m_middleGridSerial->get_elements())
        if (el)
        {
          assert(middleGridEntityDestinations.get_num_comp(el, 0) == 1);
          int elLocalIdOnDestination         = middleGridEntityDestinations(el, 0, 0).remoteId;
          MeshEntityPtr inputMeshSerialEl    = middleGridClassificationSerial(el, 0, 0);
          RemoteSharedEntity inputMeshRemote = inputMeshElementOrigins(inputMeshSerialEl, 0, 0);
          exchanger.get_send_buffer(inputMeshRemote.remoteRank).push_back(LocalIdPair{elLocalIdOnDestination, inputMeshRemote.remoteId});
        }
    }

    void unpack_data(utils::impl::ParallelExchange<LocalIdPair>& exchanger)
    {
      auto& middleGridClassificationParallel = *m_middleGridClassificationParallel;
      auto& recvBuf                          = exchanger.get_recv_buffer(m_rootRankOnInputComm);
      for (const LocalIdPair& idPair : recvBuf)
      {
        MeshEntityPtr elOnMiddleGridParallel = m_middleGridParallel->get_elements()[idPair.middleMeshElementLocalId];
        MeshEntityPtr elOnInputGridParallel  = m_inputMeshParallel->get_elements()[idPair.inputMeshElementLocalId];

        middleGridClassificationParallel(elOnMiddleGridParallel, 0, 0) = elOnInputGridParallel;
      }
    }

    bool m_amIRootRank;
    bool m_amIOnMeshComm;

    MPI_Comm m_inputComm;
    int m_rootRankOnInputComm;
    // MPI_Comm m_comm;
    std::shared_ptr<Mesh> m_middleGridSerial;
    VariableSizeFieldPtr<RemoteSharedEntity> m_middleGridSerialEntityDestinations;
    std::shared_ptr<Mesh> m_middleGridParallel;
    std::shared_ptr<Mesh> m_inputMeshParallel;
    FieldPtr<RemoteSharedEntity> m_inputMeshElementOrigins;
    FieldPtr<MeshEntityPtr> m_middleGridClassificationSerial;
    FieldPtr<MeshEntityPtr> m_middleGridClassificationParallel;
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif