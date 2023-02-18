#ifndef MESH_LAYERS_H
#define MESH_LAYERS_H

#include <limits>
#include <queue>
#include <vector>

#include "field.h"
#include "mesh.h"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

class MeshLayers
{
    using Bool = int_least8_t;

  public:
    // gets all the mesh vertices that are reachable from the given
    // roots.  This typically gets all vertices, unless there are vertices
    // enclosed by vertices excluded by func
    template <typename Tfunc>
    void get_all_layers(std::shared_ptr<Mesh> mesh, Tfunc func, std::vector<MeshEntityPtr>& roots,
                        std::vector<MeshEntityPtr>& output)
    {
      get_layers(mesh, func, roots, std::numeric_limits<int>::max() - 1, output);
    }

    // func is a function f(MeshEntityPtr) -> bool. returning false
    // if a mesh entity should be excluded
    // roots: the set of starting points.  Not all of roots needs to be
    //        f(root) -> true
    // nlayers: all verts that can be connected to roots by nlayers
    //          edges or less will be returned
    template <typename Tfunc>
    void get_layers(std::shared_ptr<Mesh> mesh, Tfunc func, std::vector<MeshEntityPtr>& roots, const int nlayers,
                    std::vector<MeshEntityPtr>& output)
    {
      assert(roots.size() > 0);
      m_seenEntities = create_field<Bool>(mesh, FieldShape(1, 0, 0), 1, false);
      std::queue<MeshEntityPtr> que1, que2;
      initialize_que(roots, que1);

      output.clear();
      std::queue<MeshEntityPtr>*queCurr = &que1, *queNext = &que2;
      for (int i = 0; i < nlayers + 1; ++i)
      {
        get_layer(queCurr, queNext, func, output, i == 0);

        if (queNext->size() == 0)
          break;

        // shuffle queues
        auto queTmp = queCurr;
        queCurr     = queNext;
        queNext     = queTmp;
      }
    }

    template <typename Tfunc>
    void get_layer(std::queue<MeshEntityPtr>* queCurr, std::queue<MeshEntityPtr>* queNext, Tfunc func,
                   std::vector<MeshEntityPtr>& output, const bool queExcludedAdjacent = false)
    {
      while (queCurr->size() != 0)
      {
        MeshEntityPtr v = queCurr->front();
        queCurr->pop();
        if (func(v))
        {
          output.push_back(v);
          que_adjacent_verts(v, queNext);
        } else if (queExcludedAdjacent)
        {
          que_adjacent_verts(v, queNext);
        }
      }
    }

    void initialize_que(std::vector<MeshEntityPtr>& roots, std::queue<MeshEntityPtr>& que);

    void que_adjacent_verts(MeshEntityPtr v, std::queue<MeshEntityPtr>* que);

    void mark_entity_seen(MeshEntityPtr e);

    bool is_entity_seen(MeshEntityPtr e);

  private:
    std::shared_ptr<Field<Bool>> m_seenEntities;
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
