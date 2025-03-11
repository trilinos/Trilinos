#ifndef STK_MIDDLE_MESH_FIELD_SCATTER
#define STK_MIDDLE_MESH_FIELD_SCATTER

#include "mesh.hpp"
#include "field.hpp"
#include "variable_size_field.hpp"
#include "utils.hpp"
#include "point.hpp"
#include "stk_util/parallel/DataExchangeKnownPatternNonBlockingBuffer.hpp"
#include "stk_util/util/ReportHandler.hpp"

struct Test
{
  int foo;
  double bar;
};

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

template <typename T>
class FieldScatterFromRoot
{
  public:
    FieldScatterFromRoot(MPI_Comm unionComm, int rootRankOnUnionComm,
                         VariableSizeFieldPtr<RemoteSharedEntity> entityDestinations,
                         FieldPtr<T> fieldSrc, FieldPtr<T> fieldDest) :
      m_unionComm(unionComm),
      m_rootRankOnUnionComm(rootRankOnUnionComm),
      m_entityDestinations(entityDestinations),
      m_fieldSrc(fieldSrc),
      m_fieldDest(fieldDest)
    {
      if (utils::impl::comm_rank(unionComm) == rootRankOnUnionComm)
        STK_ThrowRequireMsg(fieldSrc != nullptr, "fieldSrc must be provided on rootRankOnUnionComm");


      if (fieldSrc)
      {
        STK_ThrowRequireMsg(entityDestinations != nullptr, "the entity destination field must be defined on meshSerial");
        STK_ThrowRequireMsg(entityDestinations->get_mesh() == entityDestinations->get_mesh(), 
                        "the entity destination field must be defined on the same mesh as the fieldSrc");
      }

      check_fieldshapes_same_debug_only();
    }

    void scatter()
    {
      do_scatter();
    }
      
  private:
    void check_fieldshapes_same_debug_only()
    {
#ifndef NDEBUG      
      std::array<int, 4> fieldShapePlusNumComp;
      if (m_fieldSrc)
      {
        for (int i=0; i < 3; ++i)
          fieldShapePlusNumComp[i] = m_fieldSrc->get_field_shape().count[i];
        fieldShapePlusNumComp[3] = m_fieldSrc->get_num_comp();
      }

      MPI_Bcast(fieldShapePlusNumComp.data(), 4, MPI_INT, m_rootRankOnUnionComm, m_unionComm);

      if (m_fieldDest)
      {
        for (int i=0; i < 3; ++i)
          STK_ThrowRequireMsg(fieldShapePlusNumComp[i] == m_fieldDest->get_field_shape().count[i], 
                          "FieldShape must be same on source and destination fields");

        STK_ThrowRequireMsg(fieldShapePlusNumComp[3] == m_fieldDest->get_num_comp(), "Number of Field components must be same on source and destination fields");
      }
#endif
    }

    using Exchanger = stk::DataExchangeKnownPatternNonBlockingCommBuffer;
    void do_scatter()
    {
      Exchanger exchanger(m_unionComm);

      if (m_fieldDest)
      {
        set_recv_buffer_sizes(exchanger);
      }
      exchanger.allocate_recv_buffers();

      if (m_fieldSrc)
      {
        pack_data(exchanger);
      }

      exchanger.allocate_send_buffers();

      if (m_fieldSrc)
      {
        pack_data(exchanger);
      }

      exchanger.start_nonblocking();

      auto f = [&](int /*rank*/, stk::CommBuffer& buf)
      { 
        unpack_data(buf);
      };

      exchanger.complete_receives(f);
      exchanger.complete_sends();
    }

    void set_recv_buffer_sizes(Exchanger& exchanger)
    {
      STK_ThrowRequire(m_fieldDest != nullptr);

      auto meshDest     = m_fieldDest->get_mesh();
      FieldShape fshape = m_fieldDest->get_field_shape();
      auto& buf = exchanger.get_recv_buf(m_rootRankOnUnionComm);
      for (int dim=0; dim < 3; ++dim)
      {
        int numEntities = mesh::count_valid(meshDest->get_mesh_entities(dim));
        for (int i=0; i < numEntities; ++i)
        {
          buf.template pack<int>(0);
          buf.template pack<int>(0);
          for (int j=0; j < fshape.count[dim]; ++j)
            for (int k=0; k < m_fieldDest->get_num_comp(); ++k)
              buf.template pack<T>(T());
        }
      }

      exchanger.set_recv_buffer_size(m_rootRankOnUnionComm, buf.size());
    }

    void pack_data(Exchanger& exchanger)
    {
      STK_ThrowRequire(m_fieldSrc != nullptr);

      auto& entityDestinations = *m_entityDestinations;
      auto& fieldSrc           = *m_fieldSrc;
      FieldShape fshape        = fieldSrc.get_field_shape();
      auto meshSerial          = fieldSrc.get_mesh();
      for (int dim=0; dim < 3; ++dim)
        for (auto& entity : meshSerial->get_mesh_entities(dim))
          if (entity)
          {
            for (int i=0; i < entityDestinations.get_num_comp(entity, 0); ++i)
            {
              int destRank = entityDestinations(entity, 0, i).remoteRank;
              int destLocalId = entityDestinations(entity, 0, i).remoteId;

              auto& buf = exchanger.get_send_buf(destRank);
              buf.pack(dim);
              buf.pack(destLocalId);

              for (int j=0; j < fshape.count[dim]; ++j)
                for (int k=0; k < fieldSrc.get_num_comp(); ++k)
                {
                  buf.pack(fieldSrc(entity, j, k));
                }
            }   
          }
    }

    void unpack_data(stk::CommBuffer& buf)
    {
      STK_ThrowRequire(m_fieldDest);

      auto& fieldDest   = *m_fieldDest;
      FieldShape fshape = fieldDest.get_field_shape();
      while (buf.remaining() > 0)
      {
        int dim, localId;
        buf.unpack<int>(dim);
        buf.unpack<int>(localId);
        MeshEntityPtr entity = fieldDest.get_mesh()->get_mesh_entities(dim)[localId];

        for (int i=0; i < fshape.count[dim]; ++i)
          for (int j=0; j < fieldDest.get_num_comp(); ++j)
            buf.unpack<T>(fieldDest(entity, i, j));
      }
    }
    
    MPI_Comm m_unionComm;
    int m_rootRankOnUnionComm;
    VariableSizeFieldPtr<RemoteSharedEntity> m_entityDestinations;
    FieldPtr<T> m_fieldSrc;
    FieldPtr<T> m_fieldDest;

};


}
}
}
}

#endif
