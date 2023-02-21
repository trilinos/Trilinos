#ifndef MESH_AGGLOMERATOR_H
#define MESH_AGGLOMERATOR_H

#include "field.h"
#include "mesh.h"
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

    // gets the indices of all groups that have at least nthres verts
    // from verts_in
    std::vector<int> get_group_idxs(SetType<MeshEntityPtr>& vertsIn, const int nthres);

    int get_num_groups() const { return m_elements.size(); }

    const SetType<MeshEntityPtr>& get_group_elements(const int i) const { return m_elements[i]; }
    const SetType<MeshEntityPtr>& get_group_verts(const int i) const { return m_verts[i]; }

  private:
    // func is a function f(MeshEntityPtr) -> bool. True if an edge
    // is not crossable
    // entity_remap is a function f(MeshEntityPtr) -> MeshEntityPtr that maps vertices and
    // element from mesh_in to whatever mesh entities will be passed into getGroupEntities()
    template <typename T>
    void setup_groups(std::shared_ptr<Mesh> meshIn, T func);

    template <typename T>
    void create_group(MeshEntityPtr elStart, T func);

    template <typename T>
    void convert_entities(T entityRemap);

    std::vector<SetType<MeshEntityPtr>> m_elements;
    std::vector<SetType<MeshEntityPtr>> m_verts;
};

template <typename T>
void MeshAgglomerator::setup_groups(std::shared_ptr<Mesh> meshIn, T func)
{
  int nelementsSeen = 0;
  int nelem         = count_valid(meshIn->get_elements());

  while (nelementsSeen != nelem)
  {
    // find element not yet seen
    MeshEntityPtr elStart = nullptr;
    for (auto& el1 : meshIn->get_elements())
    {
      if (el1)
      {
        bool found = false;
        for (int i = 0; i < get_num_groups(); ++i)
          if (m_elements[i].count(el1) > 0)
          {
            found = true;
            break;
          }

        if (found)
          continue;
        else
        {
          elStart = el1;
          break;
        }
      }
    }

    assert(elStart);

    create_group(elStart, func);
    nelementsSeen += m_elements.back().size();
  }
}

template <typename T>
void MeshAgglomerator::create_group(MeshEntityPtr elStart, T func)
{
  m_elements.emplace_back();
  auto& els = m_elements.back();
  els.insert(elStart);

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
        if (els.count(otherEl) == 0)
        {
          els.insert(otherEl);
          que.push(otherEl);
        }
      }
    }
  }

  // get the verts
  m_verts.emplace_back();
  auto& vertSet = m_verts.back();
  MeshEntityPtr verts[MAX_DOWN];
  for (auto& el : els)
  {
    int nverts = get_downward(el, 0, verts);
    for (int i = 0; i < nverts; ++i)
      vertSet.insert(verts[i]);
  }
}

template <typename T>
void MeshAgglomerator::convert_entities(T entityRemap)
{
  for (int i = 0; i < get_num_groups(); ++i)
  {
    // TODO: could replace the entities one by one in the same set
    SetType<MeshEntityPtr> verts;
    for (auto& vert : m_verts[i])
      verts.insert(entityRemap(vert));
    m_verts[i] = verts;

    SetType<MeshEntityPtr> els;
    for (auto& el : m_elements[i])
      els.insert(entityRemap(el));
    m_elements[i] = els;
  }
}

std::ostream& operator<<(std::ostream& os, const MeshAgglomerator& agg);

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
