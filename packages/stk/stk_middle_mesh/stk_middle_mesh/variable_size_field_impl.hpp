#ifndef VARIABLE_SIZE_FIELD_IMPL_H
#define VARIABLE_SIZE_FIELD_IMPL_H

#include "mesh_entity.hpp"
#include <vector>

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

template <typename T>
class VariableSizeFieldForDimension
{
  public:
    VariableSizeFieldForDimension(int entityDimension, int numNodesPerEntity, int numEntities)
      : m_entityDimension(entityDimension)
      , m_numNodesPerEntity(numNodesPerEntity)
      , m_indices(numEntities * numNodesPerEntity, Indices{0, 0})
    {}

    T& operator()(MeshEntityPtr entity, int node, int component)
    {
      assert(get_type_dimension(entity->get_type()) == m_entityDimension);
      int idx = get_value_idx(entity, node, component);
      return m_values[idx];
    }

    const T& operator()(MeshEntityPtr entity, int node, int component) const
    {
      assert(get_type_dimension(entity->get_type()) == m_entityDimension);
      int idx = get_value_idx(entity, node, component);
      return m_values[idx];
    }

    void insert(MeshEntityPtr entity, int node, const T& val = T())
    {
      assert(m_values.size() == m_freeIndices.size());
      Indices& idxs = m_indices[get_indices_idx(entity, node)];
      int ncomp     = get_num_comp(idxs);
      if (get_num_comp(idxs) == 0)
      {
        m_values.push_back(val);
        m_freeIndices.push_back(false);
        idxs.start         = m_values.size() - 1;
        idxs.onePastTheEnd = m_values.size();

      } else
      {
        int candidateIdx = get_value_idx(entity, node, ncomp - 1) + 1;
        if (size_t(candidateIdx) < m_values.size() && m_freeIndices[candidateIdx])
        {
          m_values[candidateIdx]      = val;
          m_freeIndices[candidateIdx] = false;
          idxs.onePastTheEnd++;
        } else
        {
          move_values_to_end(idxs);
          m_values.push_back(val);
          m_freeIndices.push_back(false);
          idxs.onePastTheEnd++;
        }
      }
    }

    void set(const T& init)
    {
      for (auto& val : m_values)
        val = init;
    }

    void clear()
    {
      m_values.clear();
      m_freeIndices.clear();
      std::fill(m_indices.begin(), m_indices.end(), Indices{0, 0});
    }

    int get_num_nodes() const { return m_numNodesPerEntity; }

    int get_num_comp(MeshEntityPtr entity, int node) const
    {
      auto& indices = m_indices[get_indices_idx(entity, node)];
      return get_num_comp(indices);
    }

    void add_entity() { m_indices.push_back(Indices{0, 0}); }

    void condense_array(const std::vector<MeshEntityPtr>& entities)
    {
      if (entities.size() == 0)
        return;

      unsigned int offset = 0;
      for (unsigned int i = 0; i < entities.size(); ++i)
      {
        MeshEntityPtr e = entities[i];
        if (e)
        {
          for (int node = 0; node < m_numNodesPerEntity; ++node)
          {
            int idx                 = get_indices_idx(e, node);
            m_indices[idx - offset] = m_indices[idx];
          }
        } else
          offset += m_numNodesPerEntity;
      }

      m_indices.resize(m_indices.size() - offset);
      // TODO: maybe squeeze out empty space in m_values here
    }

  private:
    struct Indices
    {
        int start;
        int onePastTheEnd;
    };

    int get_num_comp(const Indices& indices) const { return indices.onePastTheEnd - indices.start; }

    int get_value_idx(MeshEntityPtr entity, int node, int component) const
    {
      const Indices& idxs = m_indices[get_indices_idx(entity, node)];
      assert(component >= 0 && component < get_num_comp(idxs));
      // int idx = idxs.start + component
      int start = idxs.start;
      int idx   = start + component;
      assert(size_t(idx) >= 0 && size_t(idx) < m_values.size());

      return idx;
    }

    int get_indices_idx(MeshEntityPtr entity, int node) const
    {
      assert(node >= 0 && node < m_numNodesPerEntity);
      return entity->get_id() * m_numNodesPerEntity + node;
    }

    void move_values_to_end(Indices& idxs)
    {
      int newStart = m_values.size();
      for (int i = idxs.start; i < idxs.onePastTheEnd; ++i)
      {
        m_values.push_back(m_values[i]);
        m_freeIndices[i] = false;
        m_freeIndices.push_back(false);
      }

      idxs.start         = newStart;
      idxs.onePastTheEnd = m_values.size();
    }

    int m_entityDimension;
    int m_numNodesPerEntity;
    std::vector<T> m_values;
    std::vector<Indices> m_indices;
    std::vector<bool> m_freeIndices;
};

} // namespace impl
} // namespace mesh

} // namespace middle_mesh
} // namespace stk
#endif