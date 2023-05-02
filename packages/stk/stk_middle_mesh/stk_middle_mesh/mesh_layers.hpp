#ifndef MESH_LAYERS_H
#define MESH_LAYERS_H

#include <limits>
#include <queue>
#include <vector>

#include "field.hpp"
#include "mesh.hpp"

#include "stk_util/parallel/DataExchangeUnknownPatternBlockingBuffer.hpp"
#include "stk_util/parallel/ParallelReduceBool.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

class MeshLayers
{
    using Bool = int_least8_t;

  public:
    explicit MeshLayers(std::shared_ptr<Mesh> mesh) :
      m_mesh(mesh),
      m_seenEntities(create_field<Bool>(mesh, FieldShape(1, 0, 0), 1, false)),
      m_exchanger(mesh->get_comm())
    {}
    // gets all the mesh vertices that are reachable from the given
    // roots.  This typically gets all vertices, unless there are vertices
    // enclosed by vertices excluded by func
    template <typename Tfunc>
    void get_all_layers(Tfunc func, std::vector<MeshEntityPtr>& roots,
                        std::vector<MeshEntityPtr>& output)
    {
      get_layers(func, roots, std::numeric_limits<int>::max() - 1, output);
    }

    // func is a function f(MeshEntityPtr) -> bool. returning false
    // if a mesh entity should be excluded
    // roots: the set of starting points.  Not all of roots needs to be
    //        f(root) -> true
    // nlayers: all verts that can be connected to roots by nlayers
    //          edges or less will be returned
    template <typename Tfunc>
    void get_layers(Tfunc func, std::vector<MeshEntityPtr>& roots, const int nlayers,
                    std::vector<MeshEntityPtr>& output)
    {
      std::vector<MeshEntityPtr> que1, que2;
      initialize_que(roots, que1);

      output.clear();
      std::vector<MeshEntityPtr>*queCurr = &que1, *queNext = &que2;
      for (int i = 0; i < nlayers + 1; ++i)
      {
        get_layer(*queCurr, *queNext, func, output, i == 0);

        if (stk::is_true_on_all_procs(m_mesh->get_comm(), queNext->size() == 0))
          break;

        // shuffle queues
        auto queTmp = queCurr;
        queCurr     = queNext;
        queNext     = queTmp;
        queNext->clear();
      }
    }

  private:
    enum class AddToQue : int
    {
      Output,
      Next
    };

    struct RemoteQueueUpdate
    {
      int localId;
      AddToQue queue;
    };

    template <typename Tfunc>
    void get_layer(const std::vector<MeshEntityPtr>& queCurr, std::vector<MeshEntityPtr>& queNext, Tfunc func,
                   std::vector<MeshEntityPtr>& output, const bool queExcludedAdjacent = false)
    {
      for (auto& v : queCurr)
      {
        if (func(v))
        {          
          add_vert_to_output(v, output);
          que_adjacent_verts(v, queNext);
        } else if (queExcludedAdjacent)
        {
          que_adjacent_verts(v, queNext);
        }   
      }

      process_remote_updates(queNext, output);
    }

    void initialize_que(std::vector<MeshEntityPtr>& roots, std::vector<MeshEntityPtr>& que);

    void add_vert_to_output(MeshEntityPtr v, std::vector<MeshEntityPtr>& que);

    void que_adjacent_verts(MeshEntityPtr v, std::vector<MeshEntityPtr>& que, std::set<int>* queIds = nullptr);

    void process_remote_updates(std::vector<MeshEntityPtr>& queNext, std::vector<MeshEntityPtr>& output);

    void mark_entity_seen(MeshEntityPtr e);

    bool is_entity_seen(MeshEntityPtr e);

  private:
    std::shared_ptr<Mesh> m_mesh;
    std::shared_ptr<Field<Bool>> m_seenEntities;
    stk::DataExchangeUnknownPatternBlockingBuffer<RemoteQueueUpdate> m_exchanger;
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
