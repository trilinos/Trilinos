#ifndef MESH_AGGLOMERATOR_H
#define MESH_AGGLOMERATOR_H

#include "field.hpp"
#include "mesh.hpp"
#include <ostream>
#include <queue>
#include <set>
#include <vector>

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

class MeshAgglomerator
{
  public:
    template <typename T, typename Tremap>
    MeshAgglomerator(std::shared_ptr<Mesh> meshIn, T func, Tremap entityRemap)
    {
      setup_groups(meshIn, func);
      convert_entities(entityRemap);
    }

    template <typename T>
    MeshAgglomerator(std::shared_ptr<Mesh> meshIn, T func)
      : MeshAgglomerator(meshIn, func, [](MeshEntityPtr entity) { return entity; })
    {}

    template <typename T>
    using SetType = std::set<T, MeshEntityCompare>;

    template <typename T>
    using VectorType = std::vector<T>;

    // gets the indices of all groups that have at least nthres verts
    // from verts_in
    std::vector<int> get_group_idxs(SetType<MeshEntityPtr>& vertsIn, const int nthres);

    int get_num_groups() const { return m_elements.size(); }

    const VectorType<MeshEntityPtr>& get_group_elements(const int i) const { return m_elements[i]; }
    const VectorType<MeshEntityPtr>& get_group_verts(const int i) const { return m_verts[i]; }

  private:
    // func is a function f(MeshEntityPtr) -> bool. True if an edge
    // is not crossable
    // entity_remap is a function f(MeshEntityPtr) -> MeshEntityPtr that maps vertices and
    // element from mesh_in to whatever mesh entities will be passed into getGroupEntities()
    template <typename T>
    void setup_groups(std::shared_ptr<Mesh> meshIn, T func);

    template <typename T>
    void create_group(std::shared_ptr<Mesh> meshIn, MeshEntityPtr elStart, T func);

    void get_verts(std::shared_ptr<Mesh> meshIn, const VectorType<MeshEntityPtr>& els, VectorType<MeshEntityPtr>& verts);

    template <typename T>
    void convert_entities(T entityRemap);

    bool contains_entity_sorted(const VectorType<MeshEntityPtr>& entities, MeshEntityPtr entity);

    std::vector<VectorType<MeshEntityPtr>> m_elements;
    std::vector<VectorType<MeshEntityPtr>> m_verts;
};

template <typename T>
void MeshAgglomerator::setup_groups(std::shared_ptr<Mesh> meshIn, T func)
{
  using Bool = int_least8_t;
  int nelementsSeen = 0;
  int nelem         = count_valid(meshIn->get_elements());
  size_t elStartIdx = 0;
  FieldPtr<Bool> seenElsFieldPtr = create_field<Bool>(meshIn, FieldShape(0, 0, 1), 1, false);
  auto& seenElsField = *seenElsFieldPtr;

  while (nelementsSeen != nelem)
  {
    // find element not yet seen
    MeshEntityPtr elStart = nullptr;
    for (size_t i=elStartIdx; i < meshIn->get_elements().size(); ++i)
    {
      MeshEntityPtr el1 = meshIn->get_elements()[i];
      if (el1 && !(seenElsField(el1, 0, 0)))
      {
        elStart    = el1;
        elStartIdx = i;
        break;
      }
    }

    assert(elStart);

    create_group(meshIn, elStart, func);
    for (auto& el : m_elements.back())
    {
      seenElsField(el, 0, 0) = true;
    }
    nelementsSeen += m_elements.back().size();
    elStartIdx++;
  }
}

template <typename T>
void MeshAgglomerator::create_group(std::shared_ptr<Mesh> meshIn, MeshEntityPtr elStart, T func)
{
  using Bool = int_least8_t;
  FieldPtr<Bool> seenEntitiesFieldPtr = create_field<Bool>(meshIn, FieldShape(0, 0, 1), 1, false);
  auto& seenEntitiesField = *seenEntitiesFieldPtr;

  m_elements.emplace_back();
  auto& els = m_elements.back();
  els.push_back(elStart);
  seenEntitiesField(elStart, 0, 0) = true;

  std::queue<MeshEntityPtr> que;
  que.push(elStart);

  while (que.size() > 0)
  {
    auto el = que.front();
    que.pop();
    for (int i = 0; i < el->count_down(); ++i)
    {
      MeshEntityPtr edge = el->get_down(i);
      if (edge->count_up() == 1)
        continue;

      if (!func(edge))
      {
        auto otherEl = edge->get_up(0) == el ? edge->get_up(1) : edge->get_up(0);
        if (!(seenEntitiesField(otherEl, 0, 0)))
        {
          els.push_back(otherEl);
          que.push(otherEl);
          seenEntitiesField(otherEl, 0, 0) = true;
        }
      }
    }
  }

  m_verts.emplace_back();
  get_verts(meshIn, els, m_verts.back());
}

template <typename T>
void MeshAgglomerator::convert_entities(T entityRemap)
{
  for (int i = 0; i < get_num_groups(); ++i)
  {
    SetType<MeshEntityPtr> verts;
    for (auto& vert : m_verts[i])
      vert = entityRemap(vert);

    SetType<MeshEntityPtr> els;
    for (auto& el : m_elements[i])
      el = entityRemap(el);

    std::sort(m_verts[i].begin(), m_verts[i].end(),       MeshEntityCompare());
    std::sort(m_elements[i].begin(), m_elements[i].end(), MeshEntityCompare());
  }
}

std::ostream& operator<<(std::ostream& os, const MeshAgglomerator& agg);

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
