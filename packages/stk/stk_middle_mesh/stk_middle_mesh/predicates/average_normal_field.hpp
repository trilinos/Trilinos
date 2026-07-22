#ifndef STK_MIDDLE_MESH_AVERAGE_NORMAL_FIELD_H
#define STK_MIDDLE_MESH_AVERAGE_NORMAL_FIELD_H

#include "stk_middle_mesh/mesh_entity.hpp"
#include "stk_middle_mesh/field.hpp"
#include "stk_middle_mesh/utils.hpp"
#include "stk_util/parallel/DataExchangeKnownPatternNonBlockingBuffer.hpp"

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

class AveragedNormalField
{
  public:
    AveragedNormalField(std::shared_ptr<mesh::Mesh> mesh) :
      m_mesh(mesh),
      m_normalField(mesh::create_field<utils::Point>(mesh, mesh::FieldShape(1, 0, 0), 1, {0, 0, 0}))
    {
      auto edgeInfoFieldPtr = mesh::create_field<EdgeInfo>(m_mesh, mesh::FieldShape(1, 0, 0), 1);

      compute_averaged_normal_vectors(edgeInfoFieldPtr);
      compute_edge_lengths(edgeInfoFieldPtr);
      parallel_sum_edge_info(edgeInfoFieldPtr);
      apply_scale_factor(edgeInfoFieldPtr);
    }

    mesh::FieldPtr<utils::Point> get_field() const { return m_normalField; }

  private:
    struct EdgeInfo
    {
      int count = 0;
      utils::Point normal = {0, 0, 0};
      double edgeLen = 0;

      EdgeInfo& operator+=(const EdgeInfo& other)
      {
        count += other.count;
        normal += other.normal;
        edgeLen   += other.edgeLen;

        return *this;
      }
    };

    struct EntityAndLength
    {
      int entityId;
      EdgeInfo lenAndCount;
    };


    void compute_averaged_normal_vectors(mesh::FieldPtr<EdgeInfo> edgeInfoFieldPtr)
    {
      auto& edgeInfoField = *edgeInfoFieldPtr;

      std::vector<mesh::MeshEntityPtr> elements;
      std::array<mesh::MeshEntityPtr, 4> verts;
      int myrank = utils::impl::comm_rank(m_mesh->get_comm());
      for (auto& vert : m_mesh->get_vertices())
      {
        if (vert)
        {
          elements.clear();
          get_upward(vert, 2, elements);

          for (auto el : elements)
          {
            if (mesh::get_owner(m_mesh, el) == myrank)
            {
              // get index of vert
              int nverts  = get_downward(el, 0, verts.data());
              int vertIdx = get_vert_idx(verts, nverts, vert);
              assert(vertIdx >= 0 && vertIdx < nverts);

              // get two adjacent vertices
              mesh::MeshEntityPtr vertNext = verts[(vertIdx + 1) % nverts];
              mesh::MeshEntityPtr vertPrev = verts[(vertIdx - 1 + nverts) % nverts];

              // form form vector v_i+1 - v_i and v_i-1 - v_i
              utils::Point b1 = vertNext->get_point_orig(0) - vert->get_point_orig(0);
              utils::Point b2 = vertPrev->get_point_orig(0) - vert->get_point_orig(0);

              // compute normal
              utils::Point normal = cross(b1, b2);
              edgeInfoField(vert, 0, 0).normal += normal;
            }
          }
        }
      }
    }


    int get_vert_idx(const std::array<mesh::MeshEntityPtr, 4>& verts, int nverts, mesh::MeshEntityPtr vert)
    {
      for (int i = 0; i < nverts; ++i)
        if (verts[i] == vert)
          return i;

      return -1;
    }


    void compute_edge_lengths(mesh::FieldPtr<EdgeInfo> edgeInfoFieldPtr)
    {
      int myrank = utils::impl::comm_rank(m_mesh->get_comm());
      auto& edgeInfoField = *edgeInfoFieldPtr;

      for (auto& vert : m_mesh->get_vertices())
      {
        if (vert)
        {
          double totalEdgeLength = 0;
          int numOwnedEdges = 0;
          for (int i = 0; i < vert->count_up(); ++i)
          {
            mesh::MeshEntityPtr edge = vert->get_up(i);
            if (mesh::get_owner(m_mesh, edge) == myrank)
            {
              utils::Point disp        = edge->get_down(1)->get_point_orig(0) - edge->get_down(0)->get_point_orig(0);
              totalEdgeLength += std::sqrt(dot(disp, disp));
              numOwnedEdges++;
            }
          }

          edgeInfoField(vert, 0, 0).count = numOwnedEdges;
          edgeInfoField(vert, 0, 0).edgeLen = totalEdgeLength;          
        }
      }   
    }


    void parallel_sum_edge_info(mesh::FieldPtr<EdgeInfo> edgeInfoFieldPtr)
    {
      auto& normalField = *m_normalField;
      stk::DataExchangeKnownPatternNonBlockingBuffer<EntityAndLength> exchanger(m_mesh->get_comm());
      auto& edgeInfoField = *edgeInfoFieldPtr;

      for (auto& vert : m_mesh->get_vertices())
      {
        if (vert)
        {
          for (int i=0; i < vert->count_remote_shared_entities(); ++i)
          {
            mesh::RemoteSharedEntity remote = vert->get_remote_shared_entity(i);
            exchanger.get_send_buf(remote.remoteRank).push_back({remote.remoteId, edgeInfoField(vert, 0, 0)});
            exchanger.get_recv_buf(remote.remoteRank).push_back(EntityAndLength{});
          }
        }
      }

      exchanger.start_nonblocking();

      auto unpacker = [&](int /*rank*/, const std::vector<EntityAndLength>& buf)
      {
        for (EntityAndLength entityAndLength : buf)
        {
          EdgeInfo edgeInfo = entityAndLength.lenAndCount;

          mesh::MeshEntityPtr vert = m_mesh->get_vertices()[entityAndLength.entityId];
          edgeInfoField(vert, 0, 0) += edgeInfo;
          normalField(vert, 0, 0) = edgeInfoField(vert, 0, 0).normal;
        }
      };

      exchanger.complete_receives(unpacker);
    }


    void apply_scale_factor(mesh::FieldPtr<EdgeInfo> edgeInfoFieldPtr)
    {
      auto& normalField   = *m_normalField;
      auto& edgeInfoField = *edgeInfoFieldPtr;
      for (auto& vert : m_mesh->get_vertices())
      {
        if (vert)
        {
          EdgeInfo edgeInfo = edgeInfoField(vert, 0, 0);
          double avgEdgeLength = edgeInfo.edgeLen / edgeInfo.count;

          normalField(vert, 0, 0) = scale_normal(edgeInfo.normal, avgEdgeLength);
        }  
      } 
    }     

    utils::Point scale_normal(const utils::Point& normal, double length)
    {
      double normalLength = std::sqrt(dot(normal, normal));
      if (normalLength < 1e-13)
        throw std::runtime_error("normal vector has zero length");

      return length * normal / normalLength;
    }


    std::shared_ptr<mesh::Mesh> m_mesh;
    mesh::FieldPtr<utils::Point> m_normalField;
};

}
}
}
}

#endif