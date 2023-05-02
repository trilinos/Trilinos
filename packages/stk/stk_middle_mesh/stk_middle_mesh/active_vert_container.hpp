#ifndef STK_MIDDLE_MESH_ACTIVE_VERT_CONTAINER_H
#define STK_MIDDLE_MESH_ACTIVE_VERT_CONTAINER_H

#include "mesh.hpp"
#include "mesh_layers.hpp"
#include "active_vert_data.hpp"
#include "utils.hpp"
#include "remote_coordinate_updator.hpp"

#include "stk_util/parallel/DataExchangeKnownPatternNonBlockingBuffer.hpp"


namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {


class ActiveVertContainer
{
  public:

    // Creates an ActiveVertData for each selected vertex.  The selection process
    // is combination of 2 factors
    //   1. The vertex must have filter(MeshEntityPtr) return true
    //   2. The vertex must be reachable from within nlayers edges from a vertex that
    //      is either
    //     a. has filter(MeshEntityPtr) returning false or 
    //     b. has isValid(MeshEntityPtr) return false
    // Setting nlayers = -1 means there is no limit on the distance between the
    // root vertices and the selected vertices
    template <typename Tfilter, typename Tvalid>
    ActiveVertContainer(std::shared_ptr<Mesh> mesh, Tfilter filter, Tvalid isValid, const int nlayers=-1) :
      m_mesh(mesh),
      m_recvPatchMapping(utils::impl::comm_size(mesh->get_comm())),
      m_coordExchanger(mesh->get_comm()),
      m_coordUpdator(mesh)
    {
      get_active_verts(filter, isValid, nlayers);
    }

    std::vector<opt::impl::ActiveVertData>& get_active_verts() { return m_activeVerts; }

    void update_remote_coords()
    {
      start_coord_update();
      finish_coord_update();
    }

  private:
    using Exchanger = stk::DataExchangeUnknownPatternBlockingCommBuffer;

    struct RemoteActiveVertData
    {
      std::vector<RemoteSharedEntity> vertIds;
      std::vector<int> triVertIndices;
    };

    struct PatchRemoteInfo
    {
      explicit PatchRemoteInfo(int patchIdx_) :
        patchIdx(patchIdx_)
      {}

      int patchIdx;
      std::vector<int> uniqueVertIdxsOnDest;
    };

    template <typename Tfilter, typename Tvalid>
    void get_active_verts(Tfilter filter, Tvalid isValid, const int nlayers)
    {
      int myRank = utils::impl::comm_rank(m_mesh->get_comm());
      
      std::vector<MeshEntityPtr> roots, verts;
      get_roots(filter, isValid, roots);

      MeshLayers layers(m_mesh);
      if (nlayers == -1)
        layers.get_all_layers(filter, roots, verts);
      else
        layers.get_layers(filter, roots, nlayers, verts);

      Exchanger exchanger(m_mesh->get_comm());
      for (auto& v : verts)
        if (v && filter(v))
        {
          if (mesh::get_owner(m_mesh, v) == myRank)
          {
            m_activeVerts.emplace_back(m_mesh, v);
            assert(m_activeVerts.back().get_num_elements() > 0);
          }
          else
          {
            auto& patch = m_nonOwnedActiveVerts.emplace_back(m_mesh, v);
            assert(patch.get_num_elements() != 0);
            set_patch_destinations(patch, exchanger);
          }
        }

      exchanger.allocate_send_buffers();
      for (auto& patch : m_nonOwnedActiveVerts)
        set_patch_destinations(patch, exchanger);

      exchanger.execute();
      get_remote_patch_contributions(exchanger);
      update_remote_coords();
      set_remote_coords_orig();
    }

    template <typename Tfilter, typename Tvalid>
    void get_roots(Tfilter filter, Tvalid isValid, std::vector<MeshEntityPtr>& roots)
    {
      roots.clear();
      for (auto v : m_mesh->get_vertices())
      {
        if (v && (!filter(v) || !isValid(v)))
        {
            roots.push_back(v);
        }
      }

      if (stk::is_true_on_all_procs(m_mesh->get_comm(), roots.size() == 0))
        throw std::runtime_error("some roots must be provided for mesh quality optimization");
    }

    void set_patch_destinations(opt::impl::ActiveVertData& patch, Exchanger& exchanger);

    void get_remote_patch_contributions(Exchanger& exchanger);

    void compute_vertex_to_patch_map(std::map<int, int>& vertexToPatch);


    void pack(Exchanger& exchanger, int destRank, const RemoteActiveVertData& remoteData);

    RemoteActiveVertData unpack(Exchanger& exchanger, int sendRank); 

    void merge_patches(int patchIdx, int senderRank, RemoteActiveVertData& remoteData);

    void start_coord_update();

    void finish_coord_update();

    void unpack_buffer(int rank, const std::vector<utils::Point>& vertCoords);

    void set_remote_coords_orig();

    std::shared_ptr<Mesh> m_mesh;
    std::vector<opt::impl::ActiveVertData> m_activeVerts;
    std::vector<opt::impl::ActiveVertData> m_nonOwnedActiveVerts;

    std::vector<std::vector<PatchRemoteInfo>> m_recvPatchMapping;
    stk::DataExchangeKnownPatternNonBlockingBuffer<utils::Point> m_coordExchanger;
    stk::middle_mesh::mesh::impl::RemoteCoordinateUpdator m_coordUpdator;
};


}
}
}
}


#endif