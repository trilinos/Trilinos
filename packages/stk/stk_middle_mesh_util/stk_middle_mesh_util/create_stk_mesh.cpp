#include "stk_io/WriteMesh.hpp"
#include "create_stk_mesh.hpp"
#include "stk_middle_mesh/field.hpp"
#include "constants.hpp"
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_io/DatabasePurpose.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include "Ioss_Property.h"


namespace stk {
namespace middle_mesh {
namespace stk_interface {

std::string StkMeshCreator::vertex_field_name()
{
  return impl::NAME_PREFIX + "nonconformal_interface_vert_field";
}

void StkMeshCreator::declare_stk_vert_field()
{
  m_metaDataPtr->enable_late_fields();
  m_stkNodeField = &(m_metaDataPtr->declare_field<VertIdType>(stk::topology::NODE_RANK, vertex_field_name()));
}

void StkMeshCreator::load_mesh(const std::string& fname)
{
  stk::io::StkMeshIoBroker reader(m_bulkDataPtr->parallel());
  if (m_autodecompMethod != "NONE")
  {
    reader.property_add(Ioss::Property("DECOMPOSITION_METHOD", m_autodecompMethod));
  }
  reader.set_bulk_data(*m_bulkDataPtr);
  reader.add_mesh_database(fname, stk::io::READ_MESH);
  reader.create_input_mesh();
  reader.populate_bulk_data();
}

MeshPart StkMeshCreator::create_mesh_from_part(const std::string& name)
{
  std::shared_ptr<mesh::Mesh> mesh = mesh::make_empty_mesh(m_bulkDataPtr->parallel());
  MeshFieldPtr stkEls =
      mesh::create_field<MeshFieldPtr::element_type::value_type>(mesh, mesh::impl::FieldShape(0, 0, 1), 1);
  m_part = m_metaDataPtr->get_part(name);
  stk::mesh::put_field_on_mesh(*m_stkNodeField, *m_part, 1, 0);

  create_nodes(mesh);
  setup_vert_sharing(mesh);

  if (m_bulkDataPtr->does_sideset_exist(*m_part)) {
    create_edges(mesh, stk::topology::FACE_RANK);
    create_faces_from_sideset(mesh, stkEls);
  }
  else {
    create_edges(mesh, stk::topology::ELEM_RANK);
    create_faces_from_shells(mesh, stkEls);
  }

  setup_edge_sharing(mesh, stkEls);

  return {mesh, stkEls, m_stkNodeField, m_part};
}

// copies the coordinates from the Mesh to the STK mesh
void StkMeshCreator::write_back_coords(std::shared_ptr<mesh::Mesh> mesh, const std::string& name)
{
  m_part                = m_metaDataPtr->get_part(name);
  const auto& meshVerts = mesh->get_vertices();

  // TODO: distinguish between NODE_RANK and VERT_RANK
  const stk::mesh::FieldBase& coordField = *(m_metaDataPtr->coordinate_field());
  stk::mesh::Selector select             = *m_part;
  const stk::mesh::BucketVector& buckets = m_bulkDataPtr->get_buckets(stk::topology::NODE_RANK, select);

  for (stk::mesh::Bucket* bucket : buckets)
    for (auto& vert : *bucket)
    {
      VertIdType* vertIdxs = (stk::mesh::field_data(*m_stkNodeField, vert));
      mesh::MeshEntityPtr vert2 = meshVerts[vertIdxs[0]];
      utils::Point pt           = vert2->get_point_orig(0);

      double* coordsV = reinterpret_cast<double*>(stk::mesh::field_data(coordField, vert));
      coordsV[0]      = pt.x;
      coordsV[1]      = pt.y;
      coordsV[2]      = pt.z;
    }
}

void StkMeshCreator::create_nodes(std::shared_ptr<mesh::Mesh> mesh)
{
  // TODO: distinguish between NODE_RANK and VERT_RANK
  const stk::mesh::FieldBase& coordField = *(m_metaDataPtr->coordinate_field());
  const stk::mesh::Part& owned = m_metaDataPtr->locally_owned_part();
  const stk::mesh::Part& shared = m_metaDataPtr->globally_shared_part();
  stk::mesh::Selector select             = *m_part & (owned | shared);
  const stk::mesh::BucketVector& buckets = m_bulkDataPtr->get_buckets(stk::topology::NODE_RANK, select);

  for (stk::mesh::Bucket* bucket : buckets)
    for (auto& vert : *bucket)
    {
      const double* coordsV = reinterpret_cast<const double*>(stk::mesh::field_data(coordField, vert));

      auto vert2 = mesh->create_vertex(coordsV[0], coordsV[1], coordsV[2]);

      // record relation between meshes
      VertIdType* nodeFieldV = stk::mesh::field_data(*m_stkNodeField, vert);
      nodeFieldV[0]               = vert2->get_id();
    }
}

void StkMeshCreator::create_faces_from_sideset(std::shared_ptr<mesh::Mesh> mesh, MeshFieldPtr stkEls)
{
  stk::mesh::Entity invalidEntity;
  assert(m_bulkDataPtr->does_sideset_exist(*m_part));
  const stk::mesh::SideSet& sideset = m_bulkDataPtr->get_sideset(*m_part);
  const auto& meshVerts             = mesh->get_vertices();

  std::array<mesh::MeshEntityPtr, 4> verts;
  for (auto sIt = sideset.begin(); sIt != sideset.end(); ++sIt)
  {
    stk::mesh::SideSetEntry entry = *sIt;
    stk::mesh::Entity face        = get_side_entity_for_elem_side_pair(*m_bulkDataPtr, entry.element, entry.side);
    assert(face != invalidEntity); // check entity is valid
    int i = 0;
    for (auto it = m_bulkDataPtr->begin_nodes(face); it != m_bulkDataPtr->end_nodes(face); ++it)
    {
      VertIdType* vertIdxs = (stk::mesh::field_data(*m_stkNodeField, *it));
      verts.at(i)               = meshVerts[vertIdxs[0]];
      ++i;
    }

    assert(i == 4);
    mesh::MeshEntityPtr el = mesh->create_quad_from_verts(verts[0], verts[1], verts[2], verts[3]);
    (*stkEls)(el, 0, 0)    = entry;
  }
}

void StkMeshCreator::create_faces_from_shells(std::shared_ptr<mesh::Mesh> mesh, MeshFieldPtr stkEls)
{
  const auto& meshVerts = mesh->get_vertices();

  const stk::mesh::Part& owned = m_metaDataPtr->locally_owned_part();
  stk::mesh::Selector select   = *m_part & owned;
  const stk::mesh::BucketVector& buckets = m_bulkDataPtr->get_buckets(stk::topology::ELEMENT_RANK, select);

  stk::topology::topology_t topo = m_part->topology();
  if (!(topo == stk::topology::SHELL_TRI_3 || topo == stk::topology::SHELL_QUAD_4))
    throw std::runtime_error(
        "part topology is not SHELL_TRI_3 or SHELL_QUAD_4.  These are the only supported topologies");

  std::array<mesh::MeshEntityPtr, 4> verts;
  for (stk::mesh::Bucket* bucket : buckets)
    for (auto& stkEl : *bucket)
    {
      size_t idx = 0;
      for (auto it = m_bulkDataPtr->begin_nodes(stkEl); it != m_bulkDataPtr->end_nodes(stkEl); ++it)
      {
        assert(idx < verts.size());
        VertIdType* vertIdxs = (stk::mesh::field_data(*m_stkNodeField, *it));
        verts[idx]                = meshVerts[vertIdxs[0]];
        ++idx;
      }

      mesh::MeshEntityPtr meshEl = nullptr;
      if (topo == stk::topology::SHELL_TRI_3)
        meshEl = mesh->create_triangle_from_verts(verts[0], verts[1], verts[2]);
      else if (topo == stk::topology::SHELL_QUAD_4)
        meshEl = mesh->create_quad_from_verts(verts[0], verts[1], verts[2], verts[3]);
      else
        throw std::runtime_error("found element that is not triangle or quad");

      (*stkEls)(meshEl, 0, 0) = stk::mesh::SideSetEntry(stkEl);
    }
}

void StkMeshCreator::setup_vert_sharing(std::shared_ptr<mesh::Mesh> mesh)
{
  const stk::mesh::Part& shared = m_metaDataPtr->globally_shared_part();
  stk::mesh::Selector select             = *m_part & shared;
  const stk::mesh::BulkData& bulk = *m_bulkDataPtr;
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::NODE_RANK, select);
  FieldType& stkNodeField = *m_stkNodeField;

  stk::CommSparse commSparse(m_bulkDataPtr->parallel());
  stk::pack_and_communicate(commSparse, [&commSparse, &buckets, &bulk, &stkNodeField]() {
    std::vector<int> sharingProcs;
    for(const stk::mesh::Bucket* bptr : buckets) {
      for(stk::mesh::Entity node : *bptr) {
        bulk.comm_shared_procs(node, sharingProcs);
        VertIdType* vertIdx = (stk::mesh::field_data(stkNodeField, node));
        for(int shProc : sharingProcs) {
          commSparse.send_buffer(shProc).pack<stk::mesh::EntityKey>(bulk.entity_key(node));
          int vertId = *vertIdx;
          commSparse.send_buffer(shProc).pack<int>(vertId);
        }
      }
    }
  });

  const std::vector<mesh::MeshEntityPtr>& verts = mesh->get_vertices();

  for(int p=0; p<commSparse.parallel_size(); ++p) {
    stk::CommBuffer& buf = commSparse.recv_buffer(p);
    while(buf.remaining()) {
      stk::mesh::EntityKey key;
      int remoteVertId;
      buf.unpack<stk::mesh::EntityKey>(key);
      buf.unpack<int>(remoteVertId);
      stk::mesh::Entity node = bulk.get_entity(key);
      STK_ThrowRequireMsg(bulk.is_valid(node), "StkMeshCreator::setup_vert_sharing failed to find local node for "<<key<<" recvd from P"<<p);

      VertIdType* localVertId = stk::mesh::field_data(stkNodeField, node);
      mesh::MeshEntityPtr localVert = verts[*localVertId];
      STK_ThrowRequireMsg(localVert, "StkMeshCreator::setup_vert_sharing null vert for localVertId="<<*localVertId);

      localVert->add_remote_shared_entity({p, remoteVertId});
    }
  }
}

void StkMeshCreator::create_edges(std::shared_ptr<mesh::Mesh> mesh, stk::mesh::EntityRank rank)
{
  const stk::mesh::Part& owned = m_metaDataPtr->locally_owned_part();
  stk::mesh::Selector select             = *m_part & owned;
  const stk::mesh::BulkData& bulk = *m_bulkDataPtr;
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(rank, select);
  FieldType& stkNodeField = *m_stkNodeField;
  constexpr unsigned maxNumNodesPerEdge = 2;

  std::vector<stk::mesh::Entity> edgeNodes(maxNumNodesPerEdge);
  const std::vector<mesh::MeshEntityPtr>& verts = mesh->get_vertices();

  for(const stk::mesh::Bucket* bptr : buckets) {
    stk::topology stkTopo = bptr->topology();
    edgeNodes.resize(stkTopo.num_edges());

    for(stk::mesh::Entity stkEl : *bptr) {
      const stk::mesh::Entity* nodes = bulk.begin_nodes(stkEl);

      for(unsigned edgeOrdinal=0; edgeOrdinal<stkTopo.num_edges(); ++edgeOrdinal) {
        stk::topology stkEdgeTopo = stkTopo.edge_topology(edgeOrdinal);
        STK_ThrowRequireMsg(stkEdgeTopo.num_nodes() == maxNumNodesPerEdge,
                            "StkMeshCreator::create_edges ERROR, edges must have 2 vertices");
        stkTopo.edge_nodes(nodes, edgeOrdinal, edgeNodes.data());

        int vert1Idx = *(stk::mesh::field_data(stkNodeField, edgeNodes[0]));
        int vert2Idx = *(stk::mesh::field_data(stkNodeField, edgeNodes[1]));

        if (bulk.identifier(edgeNodes[1]) < bulk.identifier(edgeNodes[0])) {
          std::swap(vert1Idx, vert2Idx);
        }

        if (!mesh::get_common_edge(verts[vert1Idx], verts[vert2Idx])) {
          mesh->create_edge(verts[vert1Idx], verts[vert2Idx]);
        }
      }
    }
  }
}

int find_remote_vert(mesh::MeshEntityPtr vert, int remoteProc)
{
  int numRemote = vert->count_remote_shared_entities();
  int remoteVertId = -1;
  for(int i=0; i<numRemote; ++i) {
    if (vert->get_remote_shared_entity(i).remoteRank == remoteProc) {
      remoteVertId = vert->get_remote_shared_entity(i).remoteId;
    }
  }

  STK_ThrowRequireMsg(remoteVertId != -1, "StkMeshCreator failed to find remote vert for shared-proc "<<remoteProc);
  return remoteVertId;
}

void StkMeshCreator::setup_edge_sharing(std::shared_ptr<mesh::Mesh> mesh, MeshFieldPtr stkElsField)
{
  const stk::mesh::Part& shared = m_metaDataPtr->globally_shared_part();
  stk::mesh::Selector select             = *m_part & shared;
  const stk::mesh::BulkData& bulk = *m_bulkDataPtr;

  stk::CommSparse commSparse(m_bulkDataPtr->parallel());
  stk::pack_and_communicate(commSparse, [&commSparse, &bulk, &mesh, &stkElsField]() {
    std::vector<int> sharingProcs;
    constexpr unsigned maxNumEdgeNodes = 3;
    std::vector<stk::mesh::Entity> edgeNodes(maxNumEdgeNodes);
    std::vector<mesh::MeshEntityPtr> edgeVerts(maxNumEdgeNodes);
    
    const std::vector<mesh::MeshEntityPtr>& surfaceElems = mesh->get_elements();
    for(const mesh::MeshEntityPtr& elem : surfaceElems) {
      if (elem) {
        const stk::mesh::SideSetEntry& ssetEntry = (*stkElsField)(elem, 0, 0);
        stk::mesh::Entity stkEl = ssetEntry.element;
        
        const bool stkElemIsFace = ssetEntry.side != stk::mesh::INVALID_CONNECTIVITY_ORDINAL;
        if (stkElemIsFace) {
          stkEl = stk::mesh::get_side_entity_for_elem_side_pair(bulk, stkEl, ssetEntry.side);
        }
        
        stk::topology stkTopo = bulk.bucket(stkEl).topology();
        
        const stk::mesh::Entity* nodes = bulk.begin_nodes(stkEl);

        for(int dn=0; dn<elem->count_down(); ++dn) {
          mesh::MeshEntityPtr edgeEnt = elem->get_down(dn);
          STK_ThrowRequire((edgeEnt && edgeEnt->get_type() == mesh::MeshEntityType::Edge));
          edgeNodes.resize(edgeEnt->count_down());
          stkTopo.edge_nodes(nodes, dn, edgeNodes.data());
          
          edgeVerts.resize(edgeEnt->count_down());

          for(int n=0; n<edgeEnt->count_down(); ++n) {
            mesh::MeshEntityPtr vert = edgeEnt->get_down(n);
            STK_ThrowRequire((vert && vert->get_type() == mesh::MeshEntityType::Vertex));
            edgeVerts[n] = vert;
          }

          bulk.shared_procs_intersection(edgeNodes, sharingProcs);
          if (sharingProcs.empty()) {
            continue;
          }

          for(int remoteProc : sharingProcs) {
            commSparse.send_buffer(remoteProc).pack<int>(find_remote_vert(edgeVerts[0], remoteProc));
            commSparse.send_buffer(remoteProc).pack<int>(find_remote_vert(edgeVerts[1], remoteProc));
            commSparse.send_buffer(remoteProc).pack<int>(edgeEnt->get_id());
          }
        }
      }
    }
  });

  const std::vector<mesh::MeshEntityPtr>& verts = mesh->get_vertices();

  for(int p=0; p<commSparse.parallel_size(); ++p) {
    stk::CommBuffer& buf = commSparse.recv_buffer(p);
    while(buf.remaining()) {
      int vert1;
      int vert2;
      int remoteEdge;
      buf.unpack<int>(vert1);
      buf.unpack<int>(vert2);
      buf.unpack<int>(remoteEdge);

      mesh::MeshEntityPtr edge = mesh::get_common_edge(verts[vert1], verts[vert2]);
      if (edge)
      {
        edge->add_remote_shared_entity({p, remoteEdge});
      }
    }
  }
}

} // namespace stk_interface
} // namespace middle_mesh
} // namespace stk

