#ifndef FIELD_H
#define FIELD_H

#include "field_base.hpp"
#include "mesh.hpp"
#include "mesh_entity.hpp"
#include <array>
#include <vector>

namespace stk {
namespace middle_mesh {
namespace mesh {

template <typename T>
class Field : public impl::FieldBase
{
  public:
    Field(const FieldShape& fshape, int ncomp, const impl::EntityCount count, std::shared_ptr<Mesh> mesh,
          const T& init = T(), bool isFieldHeldByMesh=false)
      : m_fshape(fshape)
      , m_ncomp(ncomp)
      , m_meshSharedPtr(mesh)
      , m_mesh(mesh.get())
    {
      for (int dim = 0; dim < 3; ++dim)
        m_data[dim].resize(fshape.count[dim] * count.count[dim] * ncomp);

      set(init);

      if (isFieldHeldByMesh)
        m_meshSharedPtr = nullptr; // dont create reference cycle
    }

    Field(const Field<T>&) = delete;

    Field<T>& operator=(const Field<T>&) = delete;

    using value_type = T;

    virtual ~Field() {}

    T& operator()(MeshEntityPtr e, int node, int component)
    {
      assert(is_entity_on_mesh(e));
      int dim = get_type_dimension(e->get_type());
      int id  = e->get_id();
      return m_data[dim][get_idx(id, node, component, dim)];
    }

    int get_num_nodes(const int dim) const { return m_fshape.count[dim]; }

    int get_num_comp() const { return m_ncomp; }

    void set(const T& init)
    {
      /*
            for (int dim =0; dim < 3; ++dim)
              for ( auto& val : m_data[dim])
                val = init;
      */
      for (int dim = 0; dim < 3; ++dim)
        for (size_t i = 0; i < m_data[dim].size(); ++i)
          m_data[dim][i] = init;
    }

    const FieldShape& get_field_shape() const { return m_fshape; }

    std::shared_ptr<Mesh> get_mesh() const
    { 
      assert(m_meshSharedPtr);
      return m_meshSharedPtr;
    }

  protected:
    void add_entity(const int dim) override
    {
      if (get_num_nodes(dim) > 0)
      {
        int startidx = m_data[dim].size();
        int nentries = get_num_nodes(dim) * get_num_comp();
        m_data[dim].resize(startidx + nentries);

        for (int i = startidx; i < startidx + nentries; ++i)
          m_data[dim][i] = T();
      }
    }

    void condense_arrays(const std::vector<MeshEntityPtr>& verts, const std::vector<MeshEntityPtr>& edges,
                         const std::vector<MeshEntityPtr>& elements) override
    {
      if (m_fshape.count[0] > 0)
        condense_array(verts, m_data[0]);

      if (m_fshape.count[1] > 0)
        condense_array(edges, m_data[1]);

      if (m_fshape.count[2] > 0)
        condense_array(elements, m_data[2]);
    }

    void condense_array(const std::vector<MeshEntityPtr>& entities, std::vector<T>& data)
    {
      if (entities.size() == 0)
        return;

      int dim = 0;
      for (auto& e : entities)
        if (e)
        {
          dim = get_type_dimension(e->get_type());
          break;
        }

      unsigned int offset = 0;
      for (unsigned int i = 0; i < entities.size(); ++i)
      {
        MeshEntityPtr e = entities[i];
        if (e)
        {
          int id = e->get_id();
          for (int node = 0; node < get_num_nodes(dim); ++node)
            for (int n = 0; n < get_num_comp(); ++n)
              data[get_idx(id - offset, node, n, dim)] = data[get_idx(id, node, n, dim)];
        } else
          offset += 1;
      }

      data.resize((data.size() - offset) * get_num_nodes(dim) * get_num_comp());
    }

  private:
    int get_idx(const int id, const int node, const int component, const int dim)
    {
      // boundscheck
      assert(m_fshape.count[dim] > 0);
      assert(id < static_cast<int>(m_data[dim].size()));
      assert(node >= 0 && node < m_fshape.count[dim]);
      assert(component >= 0 && component < m_ncomp);

      int offset = get_num_nodes(dim) * m_ncomp;

      return id * offset + node * m_ncomp + component;
    }

    bool is_entity_on_mesh(MeshEntityPtr entity)
    {
      int dim               = get_type_dimension(entity->get_type());
      MeshEntityPtr entity2 = m_mesh->get_mesh_entities(dim)[entity->get_id()];
      return entity == entity2;
    }

    std::array<std::vector<T>, 3> m_data;
    FieldShape m_fshape;
    int m_ncomp;
    std::shared_ptr<Mesh> m_meshSharedPtr;
    Mesh* m_mesh;
};

template <typename T>
using FieldPtr = std::shared_ptr<Field<T>>;

template <typename T>
FieldPtr<T> create_field(std::shared_ptr<Mesh> mesh, const FieldShape& fshape, const int ncomp,
                         const T& init = T(), bool isFieldHeldByMesh=false)
{
  impl::EntityCount count(mesh->get_vertices().size(), mesh->get_edges().size(), mesh->get_elements().size());
  auto field = std::make_shared<Field<T>>(fshape, ncomp, count, mesh, init, isFieldHeldByMesh);
  mesh->attach_field(field);

  return field;
}

} // namespace mesh

} // namespace middle_mesh
} // namespace stk
#endif
