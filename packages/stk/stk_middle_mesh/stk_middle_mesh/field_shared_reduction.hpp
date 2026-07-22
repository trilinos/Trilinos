#ifndef STK_MIDDLE_MESH_FIELD_SHARED_REDUCTION_H
#define STK_MIDDLE_MESH_FIELD_SHARED_REDUCTION_H

#include "field.hpp"
#include "stk_util/parallel/DataExchangeKnownPatternNonBlockingBuffer.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {

template <typename T>
class ReductionOp
{
  public:
    // compute v = op(a, b)
    virtual T operator()(const T& a, const T& b) const = 0;

    virtual ~ReductionOp() = default;
};

template <typename T>
class ReductionOpSum : public ReductionOp<T>
{
  public:
    // compute v = op(a, b)
    T operator()(const T& a, const T& b) const override { return a + b; }
};

template <typename T>
class ReductionOpLOR : public ReductionOp<T>
{
  public:
    // compute v = op(a, b)
    T operator()(const T& a, const T& b) const override { return a || b; }
};


// applies a reduction operation to shared entities, resulting in a field
// that is parallel consistent.  No guarantees are made about what order
// the reduction operation is applied to the values
template <typename T>
class FieldSharedReduction 
{
  template <typename T2>
  using Exchanger = stk::DataExchangeKnownPatternNonBlockingBuffer<T2>;

  public:
    FieldSharedReduction(FieldPtr<T> field, const ReductionOp<T>& op) :
      m_field(field),
      m_op(op)
    {}

    void reduce()
    {
      FieldShape fshape = m_field->get_field_shape();
      for (int dim=0; dim < 2; ++dim)
      {
        if (fshape.count[dim] > 0)
        {
          Exchanger<int> indexExchanger(m_field->get_mesh()->get_comm());
          pack_indices(dim, indexExchanger);

          indexExchanger.start_nonblocking();

          Exchanger<T> dataExchanger(m_field->get_mesh()->get_comm());
          pack_buffers(dim, dataExchanger);
          dataExchanger.start_nonblocking();
          unpack_buffers(dim, indexExchanger, dataExchanger);
        }
      }
    }


  private:
    void pack_indices(int dim, Exchanger<int>& indexExchanger)
    {
      for (auto& entity : m_field->get_mesh()->get_mesh_entities(dim))
        if (entity)
        {
          for (int i=0; i < entity->count_remote_shared_entities(); ++i)
          {
            RemoteSharedEntity remote = entity->get_remote_shared_entity(i);
            indexExchanger.get_send_buf(remote.remoteRank).push_back(remote.remoteId);
            indexExchanger.get_recv_buf(remote.remoteRank).push_back(0);
          }
        }
    }

    void pack_buffers(int dim, Exchanger<T>& exchanger)
    {
      FieldShape fshape = m_field->get_field_shape();
      auto& field = *m_field;
      for (auto& entity : m_field->get_mesh()->get_mesh_entities(dim))
        if (entity)
        {
          for (int i=0; i < entity->count_remote_shared_entities(); ++i)
          {
            RemoteSharedEntity remote = entity->get_remote_shared_entity(i);
            for (int j=0; j < fshape.count[dim]; ++j)
              for (int k=0; k < field.get_num_comp(); ++k)
              {
                exchanger.get_send_buf(remote.remoteRank).push_back(field(entity, j, k));
                exchanger.get_recv_buf(remote.remoteRank).push_back(T());
              }
          }
        }
    }

    void unpack_buffers(int dim, Exchanger<int>& indexExchanger, Exchanger<T>& exchanger)
    {

      auto f = [&](int /*rank*/, const std::vector<int>& /*buf*/) {};
      indexExchanger.complete_receives(f);

      auto unpackData = [&](int rank, const std::vector<T>& buf)
      {
        unpack_buffer(dim, rank, buf, indexExchanger.get_recv_buf(rank));
      };

      exchanger.complete_receives(unpackData);
    }

    void unpack_buffer(int dim, int /*rank*/, const std::vector<T>& buf, const std::vector<int>& indices)
    {
      FieldShape fshape = m_field->get_field_shape();
      assert(buf.size() == fshape.count[dim] * m_field->get_num_comp() * indices.size());
      
      auto& entities = m_field->get_mesh()->get_mesh_entities(dim);
      int idx = 0;
      auto& field = *m_field;
      for (auto& localId : indices)
      {
        MeshEntityPtr entity = entities[localId];
        for (int i=0; i < fshape.count[dim]; ++i)
          for (int j=0; j < field.get_num_comp(); ++j)
          {
            field(entity, i, j) = m_op(field(entity, i, j), buf[idx++]);
          }
      }
    }


    FieldPtr<T> m_field;
    const ReductionOp<T>& m_op;
};

}
}
}

#endif