#ifndef FIELD_MANAGER_H
#define FIELD_MANAGER_H

#include "field_base.hpp"
#include <algorithm>
#include <memory>
#include <stdexcept>

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

class FieldManager
{
  public:
    void add_field(std::shared_ptr<FieldBase> field)
    {
      auto isEqual = [&](const std::weak_ptr<FieldBase>& valWp) -> bool {
        auto val = valWp.lock();
        return val ? val == field : false;
      };

      if (std::find_if(m_fields.begin(), m_fields.end(), isEqual) != m_fields.end())
        throw std::invalid_argument("fields must be unique");

      m_fields.emplace_back(field);
    }

    void add_entity(const int dim)
    {
      for (auto& fieldWp : m_fields)
        if (auto field = fieldWp.lock())
          field->add_entity(dim);
    }

    void condense_arrays(const std::vector<MeshEntityPtr>& verts, const std::vector<MeshEntityPtr>& edges,
                         const std::vector<MeshEntityPtr>& elements)
    {
      prune();

      for (auto& fieldWp : m_fields)
        if (auto field = fieldWp.lock())
          field->condense_arrays(verts, edges, elements);
    }

    // remove fields that have been deleted
    void prune()
    {
      if (m_fields.size() > 0)
        for (int i = m_fields.size() - 1; i >= 0; --i)
          if (m_fields[i].expired())
            m_fields.erase(m_fields.begin() + i);
    }

  private:
    std::vector<std::weak_ptr<FieldBase>> m_fields;
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
