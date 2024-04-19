#ifndef VARIABLE_SIZE_FIELD_IMPL_H
#define VARIABLE_SIZE_FIELD_IMPL_H

#include "mesh_entity.hpp"
#include <vector>

namespace stk {
namespace middle_mesh {
namespace mesh {

struct StorageUsage
{
  size_t numUsed = 0;
  size_t numAllocated = 0;
};

inline StorageUsage operator+(const StorageUsage& lhs, const StorageUsage& rhs)
{
  return {lhs.numUsed + rhs.numUsed, lhs.numAllocated + rhs.numAllocated};
}

inline std::ostream& operator<<(std::ostream& os, const StorageUsage& usage)
{
  os << "used: " << usage.numUsed << ", allocated: " << usage.numAllocated;
  return os;
}

namespace impl {

template <typename T>
class IteratorRange
{
  public:
    IteratorRange(T* begin, T* end) :
      m_begin(begin),
      m_end(end)
    {}

    T* begin() { return m_begin;}

    const T* cbegin() const { return m_begin;}

    T* end() { return m_end; }

    const T* cend() {return m_end; }

  private:
    T* m_begin;
    T* m_end;
};

template <typename T>
class ConstIteratorRange
{
  public:
    ConstIteratorRange(const T* begin, const T* end) :
      m_begin(begin),
      m_end(end)
    {}
    
    const T* begin() const { return m_begin;}

    const T* cbegin() const { return m_begin;}

    const T* end() const {return m_end; }

    const T* cend() const {return m_end; }


  private:
    const T* m_begin;
    const T* m_end;
}; 



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

    IteratorRange<T> operator()(MeshEntityPtr entity, int node)
    {
      assert(get_type_dimension(entity->get_type()) == m_entityDimension);
      const Indices& idxs = m_indices[get_indices_idx(entity, node)];
      return {&(m_values[0]) + idxs.start, &(m_values[0]) + idxs.onePastTheEnd};
    }    

    const T& operator()(MeshEntityPtr entity, int node, int component) const
    {
      assert(get_type_dimension(entity->get_type()) == m_entityDimension);
      int idx = get_value_idx(entity, node, component);
      return m_values[idx];
    }

    ConstIteratorRange<T> operator()(MeshEntityPtr entity, int node) const
    {
      assert(get_type_dimension(entity->get_type()) == m_entityDimension);
      const Indices& idxs = m_indices[get_indices_idx(entity, node)];
      return {&(m_values[0]) + idxs.start, &(m_values[0]) + idxs.onePastTheEnd};
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
        size_t candidateIdx = get_value_idx(entity, node, ncomp - 1) + 1;
        bool segmentIsAtTheEnd = candidateIdx == m_values.size();
        bool nextIndexFree     = candidateIdx < m_values.size() && m_freeIndices[candidateIdx];
        if (segmentIsAtTheEnd)
        {
          m_values.push_back(val);
          m_freeIndices.push_back(false);
          idxs.onePastTheEnd++;
        } else if (nextIndexFree)
        {
          m_values[candidateIdx]      = val;
          m_freeIndices[candidateIdx] = false;
          idxs.onePastTheEnd++;
        } else
        {
          int moveTowardsFrontDistance = compute_max_move_towards_front(idxs, std::max(ncomp, 4));
          if (moveTowardsFrontDistance > 0)
          {
            move_values_towards_front(idxs, moveTowardsFrontDistance);
            assert(m_freeIndices[idxs.onePastTheEnd]);
            m_values[idxs.onePastTheEnd]      = val;
            m_freeIndices[idxs.onePastTheEnd] = false;
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
    }

    void resize(MeshEntityPtr entity, int node, int newSize, const T& val=T())
    {
      assert(get_type_dimension(entity->get_type()) == m_entityDimension);
      assert(newSize >= 0);

      Indices& idxs = m_indices[get_indices_idx(entity, node)];
      int currSize = get_num_comp(idxs);
      if (newSize > currSize)
      {
        int numNewEntries = newSize - currSize;
        for (int i=0; i < numNewEntries; ++i)
          insert(entity, node, val);
      } else
      {
        int newEnd = idxs.start + newSize;
        for (int i=newEnd; i < idxs.onePastTheEnd; ++i)
          m_freeIndices[i] = true;

        idxs.onePastTheEnd = newEnd;
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

    StorageUsage get_storage_usage() const
    {
      size_t numUsed = 0;
      for (const Indices& idxs : m_indices)
      {
        numUsed += idxs.onePastTheEnd - idxs.start;
      }

      return {numUsed, m_values.size()};
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
      assert(idx >= 0 && size_t(idx) < m_values.size());

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
        m_freeIndices[i] = true;
        m_freeIndices.push_back(false);
      }

      idxs.start         = newStart;
      idxs.onePastTheEnd = m_values.size();
    }

    int compute_max_move_towards_front(const Indices& idxs, int maxDist)
    {
      maxDist = std::min(maxDist, idxs.start);
      int endIdx = idxs.start - maxDist;
      for (int i=idxs.start-1; i >= endIdx; --i)
      {
        if (!m_freeIndices[i])
        {
          return idxs.start - i - 1;
        }
      }

      return maxDist;
    }

    void move_values_towards_front(Indices& idxs, int n)
    {
      assert(idxs.start >= n);
      for (int srcIdx = idxs.start; srcIdx < idxs.onePastTheEnd; ++srcIdx)
      {
        int destIdx = srcIdx - n;
        assert(m_freeIndices[destIdx]);

        m_values[destIdx]      = m_values[srcIdx];
        m_freeIndices[destIdx] = false;
        m_freeIndices[srcIdx]  = true;
      }

      idxs.start         -= n;
      idxs.onePastTheEnd -= n;
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
