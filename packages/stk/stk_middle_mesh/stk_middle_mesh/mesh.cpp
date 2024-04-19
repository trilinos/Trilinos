#include "mesh.hpp"
#include "field.hpp"
#include <iostream>
#include <set>
#include <unordered_set>

#include "mesh_entity.hpp"
#include "newton2.hpp"
#include "parallel_exchange.hpp"
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlockingBuffer.hpp"


namespace stk {
namespace middle_mesh {
namespace mesh {

const std::vector<utils::Point> NO_NODES;
const std::vector<MeshEntityPtr> NO_DOWN;

Mesh::~Mesh()
{
  for (auto& v : m_verts)
    if (v)
      delete v;

  for (auto& edge : m_edges)
    if (edge)
      delete edge;

  for (auto& el : m_elements)
    if (el)
      delete el;
}

MeshEntityPtr Mesh::create_vertex(const double x, const double y, const double z)
{
  auto id = m_verts.size();
  std::vector<utils::Point> nodes{utils::Point(x, y, z)};
  auto vert = new MeshEntity(MeshEntityType::Vertex, id, nodes);
  m_verts.push_back(vert);
  m_fieldManager.add_entity(0);
  return vert;
}

MeshEntityPtr Mesh::create_edge(MeshEntityPtr vert1, MeshEntityPtr vert2)
{
  if (does_edge_exist(vert1, vert2))
    throw std::invalid_argument("creating duplicate edges not allowed");

  auto id = m_edges.size();
  std::vector<MeshEntityPtr> down{vert1, vert2};

  if (!is_unique(down))
    throw std::invalid_argument("entities must be unique");

  if (is_null(down))
    throw std::invalid_argument("entities must not be null");

  auto edge = new MeshEntity(MeshEntityType::Edge, id, NO_NODES, down);
  for (auto e : down)
    e->set_up(edge);

  m_edges.push_back(edge);
  m_fieldManager.add_entity(1);

  return edge;
}

MeshEntityPtr Mesh::create_triangle(MeshEntityPtr edge1, MeshEntityPtr edge2, MeshEntityPtr edge3,
                                    EntityOrientation edge1Orient)
{
  assert(edge1->get_type() == MeshEntityType::Edge);
  assert(edge2->get_type() == MeshEntityType::Edge);
  assert(edge3->get_type() == MeshEntityType::Edge);

  if (does_face_exist(edge1, edge2, edge3))
    throw std::invalid_argument("creating duplicate faces not allowed");

  auto id = m_elements.size();
  std::vector<MeshEntityPtr> down{edge1, edge2, edge3};
  auto orient = get_edges_orientation(down, edge1Orient);

  if (!is_unique(down))
    throw std::invalid_argument("entities must be unique");

  if (is_null(down))
    throw std::invalid_argument("entities must not be null");

  auto el = new MeshEntity(MeshEntityType::Triangle, id, NO_NODES, down, orient);

  for (auto e : down)
    e->set_up(el);

  m_elements.push_back(el);
  m_fieldManager.add_entity(2);

  return el;
}

// create element from vertices, creating edges as neeeded
MeshEntityPtr Mesh::create_triangle_from_verts(MeshEntityPtr vert1, MeshEntityPtr vert2, MeshEntityPtr vert3)
{
  assert(vert1->get_type() == MeshEntityType::Vertex);
  assert(vert2->get_type() == MeshEntityType::Vertex);
  assert(vert3->get_type() == MeshEntityType::Vertex);

  std::vector<MeshEntityPtr> down{vert1, vert2, vert3};
  if (!is_unique(down))
    throw std::invalid_argument("entities must be unique");

  if (is_null(down))
    throw std::invalid_argument("entities must not be null");

  // orientation: edge1 is v1 - v2, edge2 is v2 to v3, edge3 is v3 to v1
  MeshEntityPtr verts[4] = {vert1, vert2, vert3, vert1};
  MeshEntityPtr edges[3];

  for (int i = 0; i < 3; ++i)
  {
    // find/create edge
    auto edge = get_common_edge(verts[i], verts[i + 1]);
    if (!edge)
      edge = create_edge(verts[i], verts[i + 1]);

    edges[i] = edge;
  }

  EntityOrientation orient;
  if (edges[0]->get_down(0) == vert1)
    orient = EntityOrientation::Standard;
  else
    orient = EntityOrientation::Reversed;

  return create_triangle(edges[0], edges[1], edges[2], orient);
}

MeshEntityPtr Mesh::create_quad(MeshEntityPtr edge1, MeshEntityPtr edge2, MeshEntityPtr edge3, MeshEntityPtr edge4,
                                EntityOrientation edge1Orient)
{
  if (does_face_exist(edge1, edge2, edge3, edge4))
    throw std::invalid_argument("creating duplicate faces not allowed");

  auto id = m_elements.size();
  std::vector<MeshEntityPtr> down{edge1, edge2, edge3, edge4};
  auto orient = get_edges_orientation(down, edge1Orient);

  if (!is_unique(down))
    throw std::invalid_argument("entities must be unique");

  if (is_null(down))
    throw std::invalid_argument("entities must not be null");

  auto el = new MeshEntity(MeshEntityType::Quad, id, NO_NODES, down, orient);

  for (auto& e : down)
    e->set_up(el);

  m_elements.push_back(el);
  m_fieldManager.add_entity(2);

  return el;
}

MeshEntityPtr Mesh::create_quad_from_verts(MeshEntityPtr vert1, MeshEntityPtr vert2, MeshEntityPtr vert3,
                                           MeshEntityPtr vert4)
{
  assert(vert1->get_type() == MeshEntityType::Vertex);
  assert(vert2->get_type() == MeshEntityType::Vertex);
  assert(vert3->get_type() == MeshEntityType::Vertex);
  assert(vert4->get_type() == MeshEntityType::Vertex);

  std::vector<MeshEntityPtr> down{vert1, vert2, vert3, vert4};
  if (!is_unique(down))
    throw std::invalid_argument("entities must be unique");

  if (is_null(down))
    throw std::invalid_argument("entities must not be null");

  // orientation: edge1 is v1 - v2, edge2 is v2 to v3, edge3 is v3 to v4, edge4
  // is v4 to v1
  MeshEntityPtr verts[5] = {vert1, vert2, vert3, vert4, vert1};
  MeshEntityPtr edges[4];

  for (int i = 0; i < 4; ++i)
  {
    // find/create edge

    auto edge = get_common_edge(verts[i], verts[i + 1]);
    if (!edge)
      edge = create_edge(verts[i], verts[i + 1]);

    edges[i] = edge;
  }

  EntityOrientation orient;
  if (edges[0]->get_down(0) == vert1)
    orient = EntityOrientation::Standard;
  else
    orient = EntityOrientation::Reversed;

  return create_quad(edges[0], edges[1], edges[2], edges[3], orient);
}

const impl::GeoClassification& Mesh::get_geo_class(MeshEntityPtr e)
{
  return (*m_geoClass)(e, 0, 0);
}

void Mesh::set_geo_class(MeshEntityPtr e, const impl::GeoClassification& classI)
{
  (*m_geoClass)(e, 0, 0) = classI;
}

MeshEntityPtr Mesh::split_edge_broken(MeshEntityPtr e, double xi, MeshEntityPtr* edgesNew)
{
  double eps = 1e-12;
  assert(e->get_type() == MeshEntityType::Edge);
  // assert(xi >= -eps && xi <= 1 + eps);

  // dont create zero length edges
  if (xi < eps || xi > 1 - eps)
    throw std::invalid_argument("splitting at xi = " + std::to_string(xi) + " would create a zero length edge");

  utils::Point pt = compute_edge_coords_orig(e, xi);
  auto v1         = e->get_down(0);
  auto v3         = e->get_down(1);
  auto v2         = create_vertex(pt);

  auto e1 = create_edge(v1, v2);
  auto e2 = create_edge(v2, v3);

  v1->delete_up(e);
  v3->delete_up(e);

  if (edgesNew)
  {
    edgesNew[0] = e1;
    edgesNew[1] = e2;
  }

  // this is where fixing the downward adjacencies of faces that contained
  // edge e would go
  for (int i = 0; i < e->count_up(); ++i)
    e->get_up(i)->delete_down(e);

  delete_edge(e);

  return v2;
}

void Mesh::split_edge(MeshEntityPtr edge, double xi, MeshEntityPtr* faces)
{
  assert(edge->get_type() == MeshEntityType::Edge);
  bool haveTwo = edge->count_up() == 2;

  MeshEntityPtr elL, elR = nullptr;

  elL = edge->get_up(0);
  assert(elL->get_type() == MeshEntityType::Triangle);
  if (haveTwo)
  {
    elR = edge->get_up(1);
    assert(elR->get_type() == MeshEntityType::Triangle);
  }

  MeshEntityPtr vertsL[MAX_DOWN], vertsR[MAX_DOWN];
  get_downward(elL, 0, vertsL);
  if (haveTwo)
    get_downward(elR, 0, vertsR);

  // get local index of edge on each element
  int idxL = -1, idxR = -1;
  for (int i = 0; i < elL->count_down(); ++i)
    if (elL->get_down(i) == edge)
    {
      idxL = i;
      break;
    }
  assert(idxL >= 0);

  if (haveTwo)
  {
    for (int i = 0; i < elR->count_down(); ++i)
      if (elR->get_down(i) == edge)
      {
        idxR = i;
        break;
      }
    assert(idxR >= 0);
  }

  auto newVert = split_edge_broken(edge, xi);

  // remove elL and elR from adjacency data structure
  for (int i = 0; i < elL->count_down(); ++i)
    if (elL->get_down(i))
      elL->get_down(i)->delete_up(elL);

  if (haveTwo)
    for (int i = 0; i < elR->count_down(); ++i)
      if (elR->get_down(i))
        elR->get_down(i)->delete_up(elR);

  p_delete_face(elL);
  if (haveTwo)
    p_delete_face(elR);

  // create new left elements
  MeshEntityPtr elL1, elL2, elR1, elR2;
  elL1 = create_triangle_from_verts(vertsL[idxL], newVert, vertsL[(idxL + 2) % 3]);
  elL2 = create_triangle_from_verts(newVert, vertsL[(idxL + 1) % 3], vertsL[(idxL + 2) % 3]);

  // create new right elements
  if (haveTwo)
  {
    elR1 = create_triangle_from_verts(vertsR[idxR], newVert, vertsR[(idxR + 2) % 3]);
    elR2 = create_triangle_from_verts(newVert, vertsR[(idxR + 1) % 3], vertsR[(idxR + 2) % 3]);
  }

  if (faces)
  {
    faces[0] = elL1;
    faces[1] = elL2;
    if (haveTwo)
    {
      faces[2] = elR1;
      faces[3] = elR2;
    }
  }
}

void Mesh::delete_face(MeshEntityPtr face)
{
  // collect edges to be deleted
  int nedges = 0;
  MeshEntityPtr edges[MAX_DOWN];
  for (int i = 0; i < face->count_down(); ++i)
  {
    MeshEntityPtr edge = face->get_down(i);
    if (edge->count_up() == 1) // if the topology data structure is correct,
    {                          // this upward adjacency must be face
      edges[nedges++] = edge;
    }
  }

  // delete verts
  MeshEntityPtr verts[MAX_DOWN];
  int nverts = get_downward(face, 0, verts);
  for (int i = 0; i < nverts; ++i)
  {
    MeshEntityPtr vert = verts[i];

    // delete verts that only have edges to be deleted as upward adjacencies
    bool shouldDelete = true;
    for (int j = 0; j < vert->count_up(); ++j)
    {
      bool found = false;
      for (int k = 0; k < nedges; ++k)
        if (vert->get_up(j) == edges[k])
        {
          found = true;
          break;
        }

      if (!found)
      {
        shouldDelete = false;
        break;
      }
    }

    if (shouldDelete)
    {
      for (int j = 0; j < vert->count_up(); ++j)
        vert->get_up(j)->delete_down(vert);
      // no need to fix upward adjacencies because we are going to delete them
      delete_vert(vert);
    }
  }

  // delete edges
  for (int i = 0; i < nedges; ++i)
  {
    MeshEntityPtr edgeI = edges[i];
    for (int j = 0; j < edgeI->count_down(); ++j)
    {
      MeshEntityPtr vJ = edgeI->get_down(j);
      if (vJ)
        vJ->delete_up(edgeI);
    }
    face->delete_down(edgeI);
    delete_edge(edgeI);
  }

  // delete face
  for (int j = 0; j < face->count_down(); ++j)
  {
    MeshEntityPtr edgeJ = face->get_down(j);
    if (edgeJ)
      edgeJ->delete_up(face);
  }
  p_delete_face(face);
}

void Mesh::set_geo_classification(std::shared_ptr<Field<impl::GeoClassification>> field)
{
  m_geoClass = field;
}

void Mesh::condense_arrays()
{
  m_fieldManager.condense_arrays(m_verts, m_edges, m_elements);
  condense_array(m_verts);
  condense_array(m_edges);
  condense_array(m_elements);
}

const Mesh::TStorage& Mesh::get_mesh_entities(const int dim)
{
  if (dim == 0)
    return m_verts;
  else if (dim == 1)
    return m_edges;
  else if (dim == 2)
    return m_elements;
  else
    throw std::invalid_argument("invalid dimension");
}

void Mesh::condense_array(std::vector<MeshEntityPtr>& entities)
{
  unsigned int offset = 0;
  for (unsigned int i = 0; i < entities.size(); ++i)
  {
    MeshEntityPtr e = entities[i];
    if (e)
    {
      entities[i - offset] = e;
      e->set_id(i - offset);
    } else
      offset += 1;
  }

  entities.resize(entities.size() - offset);
}

std::vector<EntityOrientation> Mesh::get_edges_orientation(std::vector<MeshEntityPtr> edges,
                                                           EntityOrientation edge1Orient)
{
  // figure out edge orientations
  std::vector<EntityOrientation> orient{edge1Orient};
  MeshEntityPtr down2[MAX_DOWN], down3[MAX_DOWN];
  for (unsigned int i = 1; i < edges.size(); ++i)
  {
    get_downward(edges[i - 1], 0, down2);
    apply_orientation(orient[i - 1], down2, 2);
    get_downward(edges[i], 0, down3);

    if (down2[1] == down3[0])
      orient.push_back(EntityOrientation::Standard);
    else
      orient.push_back(EntityOrientation::Reversed);
  }

  return orient;
}

bool Mesh::does_edge_exist(MeshEntityPtr vert1, MeshEntityPtr vert2)
{
  for (int i = 0; i < vert1->count_up(); ++i)
    for (int j = 0; j < vert2->count_up(); ++j)
      if (vert1->get_up(i) == vert2->get_up(j))
        return true;

  return false;
}

bool Mesh::does_face_exist(MeshEntityPtr edge1, MeshEntityPtr edge2, MeshEntityPtr edge3)
{
  for (int i = 0; i < edge1->count_up(); ++i)
    for (int j = 0; j < edge2->count_up(); ++j)
      for (int k = 0; k < edge3->count_up(); ++k)
        if (edge1->get_up(i) == edge2->get_up(j) && edge1->get_up(i) == edge3->get_up(k))
          return true;

  return false;
}

bool Mesh::does_face_exist(MeshEntityPtr edge1, MeshEntityPtr edge2, MeshEntityPtr edge3, MeshEntityPtr edge4)
{
  for (int i = 0; i < edge1->count_up(); ++i)
    for (int j = 0; j < edge2->count_up(); ++j)
      for (int k = 0; k < edge3->count_up(); ++k)
        for (int p = 0; p < edge4->count_up(); ++p)
          if (edge1->get_up(i) == edge2->get_up(j) && edge1->get_up(i) == edge3->get_up(k) &&
              edge1->get_up(i) == edge4->get_up(p))

            return true;

  return false;
}

int count_valid(const std::vector<MeshEntityPtr>& entities)
{
  int count = 0;
  for (auto e : entities)
    if (e)
      count += 1;

  return count;
}

int check_angles(std::shared_ptr<Mesh> mesh, const double thetaMin, const double thetaMax)
{
  double pi        = std::atan(1) * 4;
  double thetaMinR = thetaMin * pi / 180.0;
  double thetaMaxR = thetaMax * pi / 180.0;
  int count        = 0;
  double angles[MAX_DOWN];
  for (auto& el : mesh->get_elements())
  {
    if (!el)
      continue;

    int nangles = compute_angles(el, angles);
    for (int i = 0; i < nangles; ++i)
      if (angles[i] < thetaMinR || angles[i] > thetaMaxR)
        ++count;
  }

  return count;
}

void check_topology(std::shared_ptr<Mesh> mesh, int maxDim)
{
  if (maxDim >= 1)
  {
    check_topology_down(mesh->get_edges(), mesh->get_vertices());
  }

  if (maxDim >= 2)
  {
    check_topology_down(mesh->get_elements(), mesh->get_edges());
  }

  // check every vertex has at least two upward adjacencies, and edges
  // have one
  if (maxDim >= 1)
  {
    check_topology_up(mesh->get_vertices(), mesh->get_edges(), 2);
  }

  if (maxDim >= 2)
  {
    check_topology_up(mesh->get_edges(), mesh->get_elements(), 1);
  }

  check_remotes_unique(mesh);
  check_remotes_symmetric(mesh);

  if (maxDim >= 2)
  {
    check_edge_orientation_parallel(mesh);
  }
}

void check_topology_down(const std::vector<MeshEntityPtr>& entities, const std::vector<MeshEntityPtr>& entitiesDown)
{
  // check all downward adjacencies are not null
  std::array<int, 4> expectedNumDown = {0, 2, 3, 4};
  for (auto& entity : entities)
    if (entity)
    { 
      int expectedDown = expectedNumDown.at(static_cast<int>(entity->get_type()));
      if (entity->count_down() != expectedDown)
      {
        std::stringstream ss;
        ss << "number of downward adjacencies is incorrect for " << entity << " expected " 
           << expectedDown << ", found " << entity->count_down();
        throw std::runtime_error(ss.str());
      }
      
      
      for (int i = 0; i < entity->count_down(); ++i)
      {
        auto entityDown = entity->get_down(i);
        if (!entityDown)
          throw std::runtime_error("entity has deleted downward adjacency");

        // attempt to verify that the downward adjacency entity still exists
        // in the mesh
        if (entityDown != entitiesDown[entityDown->get_id()])
          throw std::runtime_error("downward entity ID is inconsistent");
      }

      std::array<MeshEntityPtr, MAX_DOWN> down;
      int entityDim = get_type_dimension(entity->get_type());
      for (int dim=0; dim < entityDim; ++dim)
      {
        int ndown = get_downward(entity, dim, down.data());
        std::sort(down.begin(), down.begin() + ndown);
        for (int i=1; i < ndown; ++i)
          if (down[i] == down[i-1])
            throw std::runtime_error("duplicate downward entity");
      }
    }
}

void check_topology_up(const std::vector<MeshEntityPtr>& entities, const std::vector<MeshEntityPtr>& entitiesUp,
                       const int minUp)
{
  std::vector<MeshEntityPtr> retrievedUpEntities;
  for (auto& entity : entities)
  {
    if (!entity)
      continue;

    if (entity->count_up() < minUp)
      throw std::runtime_error("entity has insufficient upward adjacencies");

    for (int i = 0; i < entity->count_up(); ++i)
    {
      auto entityUp = entity->get_up(i);

      if (!entityUp)
        throw std::runtime_error("vertex has null upward edge");

      if (entityUp != entitiesUp[entityUp->get_id()])
        throw std::runtime_error("upward entity ID is inconsistent");
    }


    int entityDim = get_type_dimension(entity->get_type());
    for (int dim=entityDim+1; dim <= 2; ++dim)
    {
      retrievedUpEntities.clear();
      int nup = get_upward(entity, dim, retrievedUpEntities);
      std::sort(retrievedUpEntities.begin(), retrievedUpEntities.end());
      for (int i=1; i < nup; ++i)
        if (retrievedUpEntities[i-1] == retrievedUpEntities[i])
          throw std::runtime_error("duplicate upward entity");
    }
  }
}

namespace {
struct EdgeAndVertIds
{
    int edgeId;
    int vert1Id;
    int vert2Id;
};
} // namespace

void check_edge_orientation_parallel(std::shared_ptr<Mesh> mesh)
{
  utils::impl::ParallelExchange<EdgeAndVertIds> exchanger(mesh->get_comm(), 77);
  std::vector<int> recvSizes(utils::impl::comm_size(mesh->get_comm()), 0);

  int myrank = utils::impl::comm_rank(mesh->get_comm());
  for (auto& edge : mesh->get_edges())
    if (edge && edge->count_remote_shared_entities() > 0)
    {
      int ownerRank = get_owner(mesh, edge);
      if (ownerRank == myrank)
      {
        MeshEntityPtr v1 = edge->get_down(0);
        MeshEntityPtr v2 = edge->get_down(1);

        for (int i = 0; i < edge->count_remote_shared_entities(); ++i)
        {
          const RemoteSharedEntity& edgeRemote = edge->get_remote_shared_entity(i);
          RemoteSharedEntity v1Remote          = get_remote_shared_entity(v1, edgeRemote.remoteRank);
          RemoteSharedEntity v2Remote          = get_remote_shared_entity(v2, edgeRemote.remoteRank);

          EdgeAndVertIds ids{edgeRemote.remoteId, v1Remote.remoteId, v2Remote.remoteId};
          exchanger.get_send_buffer(edgeRemote.remoteRank).push_back(ids);
        }
      } else
      {
        recvSizes[ownerRank]++;
      }
    }

  for (size_t i = 0; i < recvSizes.size(); ++i)
    if (recvSizes[i] > 0)
      exchanger.set_recv_buffer_size(i, recvSizes[i]);

  exchanger.start_sends();
  exchanger.start_recvs();
  exchanger.complete_recvs();

  for (size_t sendRank = 0; sendRank < recvSizes.size(); ++sendRank)
    if (recvSizes[sendRank] > 0)
    {
      const auto& recvBuf = exchanger.get_recv_buffer(sendRank);
      for (const EdgeAndVertIds& ids : recvBuf)
      {
        MeshEntityPtr edge = mesh->get_edges()[ids.edgeId];
        MeshEntityPtr v1   = edge->get_down(0);
        MeshEntityPtr v2   = edge->get_down(1);
        if (v1->get_id() != ids.vert1Id || v2->get_id() != ids.vert2Id)
          throw std::runtime_error("edge has incorrect orientation or vertex RemoteSharedEntities are inconsistent with the edge RemoteSharedEntity");
      }
    }

  exchanger.complete_sends();
}

struct RemoteData
{
  int remoteId;
  int localId;
  int localRank;
};

void check_remotes_symmetric(std::shared_ptr<Mesh> mesh)
{
  int myRank = utils::impl::comm_rank(mesh->get_comm());
  for (int dim=0; dim <= 1; ++dim)
  {
    stk::DataExchangeUnknownPatternNonBlockingBuffer<RemoteData> exchanger(mesh->get_comm());

    for (auto& entity : mesh->get_mesh_entities(dim))
      if (entity)
      {
        for (int i=0; i < entity->count_remote_shared_entities(); ++i)
        {
          const RemoteSharedEntity& remote = entity->get_remote_shared_entity(i);
          RemoteData data{remote.remoteId, entity->get_id(), myRank};
          exchanger.get_send_buf(remote.remoteRank).push_back(data);
        }
      }

    exchanger.start_nonblocking();
    exchanger.post_nonblocking_receives();

    auto f = [](int rank, const std::vector<RemoteData>& buf) {};
    exchanger.complete_receives(f);
    exchanger.complete_sends();

    for (int rank=0; rank < utils::impl::comm_size(mesh->get_comm()); ++rank)
    {
      for (const RemoteData& data : exchanger.get_recv_buf(rank))
      {
        MeshEntityPtr entity = mesh->get_mesh_entities(dim)[data.remoteId];
        RemoteSharedEntity remote = get_remote_shared_entity(entity, data.localRank);

        if (remote.remoteId != data.localId)
        {
          std::stringstream ss;
          ss << "entity " << data.localId << " on rank " << data.localRank << " has a remote entity " 
             << data.remoteId << " on rank " << myRank << ", but the remote does not have a symmetric RemoteSharedEntity";
          throw std::runtime_error(ss.str());
        }
      }
    }
  }
}

void check_remotes_unique(std::shared_ptr<Mesh> mesh)
{
  std::vector<RemoteSharedEntity> remotes;
  auto cmp = [](const RemoteSharedEntity& lhs, const RemoteSharedEntity& rhs)
  {
    if (lhs.remoteRank != rhs.remoteRank)
      return lhs.remoteRank < rhs.remoteRank;
    else
      return lhs.remoteId < rhs.remoteId;
  };

  for (int dim=0; dim < 2; ++dim)
    for (MeshEntityPtr entity : mesh->get_mesh_entities(dim))
      if (entity)
      {
        int nremotes = entity->count_remote_shared_entities();
        remotes.clear();
        for (int i=0; i < nremotes; ++i)
          remotes.push_back(entity->get_remote_shared_entity(i));

        std::sort(remotes.begin(), remotes.end(), cmp);
        for (size_t i=1; i < remotes.size(); ++i)
          if (!cmp(remotes[i-1], remotes[i]) && !cmp(remotes[i], remotes[i-1]))
            throw std::runtime_error(std::string("detected duplicate remote for entity ") +
                                     std::to_string(entity->get_id()) + " of dimensions " + std::to_string(dim));
      }
}


void check_coordinate_field(std::shared_ptr<Mesh> mesh)
{
  for (int dim = 0; dim < 3; ++dim)
    for (auto& entity : mesh->get_mesh_entities(dim))
      for (int i = 0; i < entity->count_points(); ++i)
      {
        auto pt = entity->get_point_orig(i);
        if (!std::isfinite(pt.get_x()) || !std::isfinite(pt.get_y()) || !std::isfinite(pt.get_z()))
          throw std::runtime_error("original coord has non-finite value");
      }
}

void apply_orientation(EntityOrientation flag, MeshEntityPtr* down, const int n)
{
  if (flag == EntityOrientation::Reversed)
    for (int i = 0; i < n / 2; ++i)
    {
      int i2   = n - i - 1;
      auto tmp = down[i2];
      down[i2] = down[i];
      down[i]  = tmp;
    }
}

void reverse_edge(MeshEntityPtr edge)
{
  assert(edge->get_type() == MeshEntityType::Edge);

  std::vector<MeshEntityPtr> els;
  get_upward(edge, 2, els);

  std::vector<int> edgeIdxs(els.size());
  for (size_t i=0; i < els.size(); ++i)
  {
    for (int j=0; j < els[i]->count_down(); ++j)
    {
      if (els[i]->get_down(j) == edge)
      {
        edgeIdxs[i] = j;
      }
    }
  }

  mesh::MeshEntityPtr v0 = edge->get_down(0);
  mesh::MeshEntityPtr v1 = edge->get_down(1);

  edge->replace_down(0, v1);
  edge->replace_down(1, v0);

  for (size_t i=0; i < els.size(); ++i)
  {
    EntityOrientation orientOld = els[i]->get_down_orientation(edgeIdxs[i]);
    EntityOrientation orientNew = reverse(orientOld);
    els[i]->set_down_orientation(edgeIdxs[i], orientNew);
  }
}


int get_vertices(MeshEntityPtr e, MeshEntityPtr* verts)
{
  int nnodes;
  switch (e->get_type())
  {
    case MeshEntityType::Vertex: {
      verts[0] = e;
      nnodes   = 1;
      break;
    }
    case MeshEntityType::Edge: {
      verts[0] = e->get_down(0);
      verts[1] = e->get_down(1);
      nnodes   = 2;
      break;
    }
    case MeshEntityType::Triangle: {
      // TODO: is there a way to build apply_orientation into the get_downward
      //        framework?
      MeshEntityPtr down[MAX_DOWN];
      int ndown = get_downward(e->get_down(0), 0, down);
      apply_orientation(e->get_down_orientation(0), down, ndown);

      verts[0] = down[0];
      verts[1] = down[1];

      ndown = get_downward(e->get_down(1), 0, down);
      apply_orientation(e->get_down_orientation(1), down, ndown);
      verts[2] = down[1];
      nnodes   = 3;
      break;
    }
    case MeshEntityType::Quad: {
      MeshEntityPtr down[MAX_DOWN];
      int ndown = get_downward(e->get_down(0), 0, down);
      apply_orientation(e->get_down_orientation(0), down, ndown);

      verts[0] = down[0];
      verts[1] = down[1];

      ndown = get_downward(e->get_down(1), 0, down);
        
      apply_orientation(e->get_down_orientation(1), down, ndown);
      verts[2] = down[1];

      ndown = get_downward(e->get_down(2), 0, down);
      apply_orientation(e->get_down_orientation(2), down, ndown);
      verts[3] = down[1];
      nnodes   = 4;
      break;
    }
    default:
      throw std::invalid_argument("unrecognized MeshEntityType");
  }

  return nnodes;
}

int get_edges(MeshEntityPtr e, MeshEntityPtr* edges)
{
  int nnodes;
  switch (e->get_type())
  {
    case MeshEntityType::Edge: {
      edges[0] = e;
      nnodes   = 1;
      break;
    }
    case MeshEntityType::Triangle: {
      edges[0] = e->get_down(0);
      edges[1] = e->get_down(1);
      edges[2] = e->get_down(2);
      nnodes   = 3;
      break;
    }
    case MeshEntityType::Quad: {
      edges[0] = e->get_down(0);
      edges[1] = e->get_down(1);
      edges[2] = e->get_down(2);
      edges[3] = e->get_down(3);
      nnodes   = 4;
      break;
    }
    default:
      throw std::invalid_argument("unsupported MeshEntityType");
  }

  return nnodes;
}

int get_faces(MeshEntityPtr e, MeshEntityPtr* faces)
{
  if (get_type_dimension(e->get_type()) == 2)
  {
    faces[0] = e;
    return 1;
  } else
    throw std::invalid_argument("unsupported MeshEntityType");
}

int get_downward(MeshEntityPtr e, int dim, MeshEntityPtr* down)
{
  assert(dim < get_type_dimension(e->get_type()));

  switch (dim)
  {
    case 0: {
      return get_vertices(e, down);
    }
    case 1: {
      return get_edges(e, down);
    }
    case 2: {
      return get_faces(e, down);
    }
    default:
      throw std::invalid_argument("unsupported dimension");
  }
}

// gets one level upward adjacencies
int get_up_one(MeshEntityPtr e, std::vector<MeshEntityPtr>& up)
{
  up.resize(e->count_up());
  for (int i = 0; i < e->count_up(); ++i)
    up[i] = e->get_up(i);
  return e->count_up();
}

int get_vertex_up(MeshEntityPtr e, int dim, std::vector<MeshEntityPtr>& up)
{
  assert(e->get_type() == MeshEntityType::Vertex);
  static std::set<MeshEntityPtr, MeshEntityCompare> tmp;

  switch (dim)
  {
    case 1: {
      return get_up_one(e, up);
    }

    case 2: {
      tmp.clear();
      for (int i = 0; i < e->count_up(); ++i)
      {
        MeshEntityPtr edge = e->get_up(i);
        for (int j = 0; j < edge->count_up(); ++j)
          // according to the internet, using an unordered_set is the
          // fastest method of making the elements unique
          tmp.insert(edge->get_up(j));
      }
      up.assign(tmp.begin(), tmp.end());
      return up.size();
    }

    default:
      throw std::invalid_argument("unsupported dimension");
  }
}

int get_edge_up(MeshEntityPtr e, int dim, std::vector<MeshEntityPtr>& up)
{
  assert(e->get_type() == MeshEntityType::Edge);

  if (dim == 2)
    return get_up_one(e, up);
  else
    throw std::invalid_argument("unsupported dimension");
}

int get_upward(MeshEntityPtr e, int dim, std::vector<MeshEntityPtr>& up)
{
  int edim = get_type_dimension(e->get_type());
  assert(dim > edim);

  switch (edim)
  {
    case 0: {
      return get_vertex_up(e, dim, up);
    }
    case 1: {
      return get_edge_up(e, dim, up);
    }
    default:
      throw std::invalid_argument("unsupported dimension");
  }
}

// gets second order adjacencies, ie.
// entities of dimension target_dim that are connected to e via
// entities of dimension via_dim
int get_bridge_adjacent(MeshEntityPtr e, const int viaDim, const int targetDim, std::vector<MeshEntityPtr>& entities)
{
  // this is slightly inefficient because it may require two de-duplication
  // steps: one inside get_upward and one here

  int edim = get_type_dimension(e->get_type());
  assert(edim != viaDim);
  assert(viaDim != targetDim);

  static std::vector<MeshEntityPtr> viaEntities;
  static std::vector<MeshEntityPtr> targetEntities;
  static std::set<MeshEntityPtr> tmp; // for de-duplication

  // get via entities
  int nVia = 0;
  if (viaDim < edim)
  {
    viaEntities.resize(MAX_DOWN);
    nVia = get_downward(e, viaDim, viaEntities.data());
  } else // via_dim > edim
    nVia = get_upward(e, viaDim, viaEntities);

  // get target entities
  tmp.clear();
  if (targetDim < viaDim)
  {
    targetEntities.resize(MAX_DOWN);
    for (int i = 0; i < nVia; ++i)
    {
      int nentities = get_downward(viaEntities[i], targetDim, targetEntities.data());
      tmp.insert(targetEntities.begin(), targetEntities.begin() + nentities);
    }
  } else
  {
    for (int i = 0; i < nVia; ++i)
    {
      get_upward(viaEntities[i], targetDim, targetEntities);
      tmp.insert(targetEntities.begin(), targetEntities.end());
    }
  }

  entities.assign(tmp.begin(), tmp.end());

  return entities.size();
}

int get_owner(std::shared_ptr<Mesh> mesh, MeshEntityPtr entity)
{
  int owner = utils::impl::comm_rank(mesh->get_comm());
  for (int i = 0; i < entity->count_remote_shared_entities(); ++i)
    owner = std::min(owner, entity->get_remote_shared_entity(i).remoteRank);

  return owner;
}

int get_local_id(MeshEntityPtr higherDimensionEntity, MeshEntityPtr lowerDimensionEntity)
{
  int outputDim = get_type_dimension(lowerDimensionEntity->get_type());
  assert(get_type_dimension(higherDimensionEntity->get_type()) > outputDim);

  std::array<MeshEntityPtr, MAX_DOWN> entities;
  int nentities = get_downward(higherDimensionEntity, outputDim, entities.data());
  for (int i=0; i < nentities; ++i)
  {
    if (entities[i] == lowerDimensionEntity)
    {
      return i;
    }
  }

  throw std::runtime_error("unable to find entity");
}


RemoteSharedEntity get_owner_remote(std::shared_ptr<Mesh> mesh, MeshEntityPtr entity)
{
  int myrank = utils::impl::comm_rank(mesh->get_comm());
  RemoteSharedEntity remote{myrank, entity->get_id()};
  for (int i=0; i < entity->count_remote_shared_entities(); ++i)
  {
    RemoteSharedEntity remote_i = entity->get_remote_shared_entity(i);
    if (remote_i.remoteRank < remote.remoteRank)
      remote = remote_i;
  }

  return remote;
}


bool check_is_entity_owner(std::shared_ptr<Mesh> mesh, MeshEntityPtr entity)
{
  return get_owner(mesh, entity) == utils::impl::comm_rank(mesh->get_comm());
}

RemoteSharedEntity get_remote_shared_entity(MeshEntityPtr entity, int rank)
{
  for (int i = 0; i < entity->count_remote_shared_entities(); ++i)
  {
    if (entity->get_remote_shared_entity(i).remoteRank == rank)
    {
      return entity->get_remote_shared_entity(i);
    }
  }

  throw std::runtime_error("unable to find remote on given rank");
}

MeshEntityPtr get_other_vert(MeshEntityPtr v, MeshEntityPtr edge)
{
  assert(get_type_dimension(edge->get_type()) == 1);
  assert(get_type_dimension(v->get_type()) == 0);
  assert(edge->get_down(0) == v || edge->get_down(1) == v);

  return edge->get_down(0) == v ? edge->get_down(1) : edge->get_down(0);
}

// returns the edge connecting vert1 and vert2, or nullptr if no
// such edge exists
MeshEntityPtr get_common_edge(MeshEntityPtr vert1, MeshEntityPtr vert2)
{
  for (int i = 0; i < vert1->count_up(); ++i)
    for (int j = 0; j < vert2->count_up(); ++j)
    {
      if (vert1->get_up(i) == vert2->get_up(j))
      {
        return vert1->get_up(i);
      }
    }

  return nullptr;
}

utils::Point compute_edge_coords_orig(MeshEntityPtr edge, const double xi)
{
  assert(edge->get_type() == MeshEntityType::Edge);
  auto p1 = edge->get_down(0)->get_point_orig(0);
  auto p2 = edge->get_down(1)->get_point_orig(0);
  return compute_edge_coords(p1, p2, xi);
}

utils::Point compute_edge_coords(const utils::Point& p1, const utils::Point& p2, const double xi)
{
  double dx = p2.get_x() - p1.get_x();
  double dy = p2.get_y() - p1.get_y();
  double dz = p2.get_z() - p1.get_z();
  double x  = xi * dx + p1.get_x();
  double y  = xi * dy + p1.get_y();
  double z  = xi * dz + p1.get_z();

  return utils::Point(x, y, z);
}

void compute_lagrange_vals(const double xi, double vals[2])
{
  vals[0] = -(xi - 1);
  vals[1] = xi;
}

void compute_lagrange_derivs(const double xi, double derivs[2])
{
  derivs[0] = -1;
  derivs[1] = 1;
}

utils::Point compute_quad_coords_from_xi(const double xi[2], const utils::Point& pt1, const utils::Point& pt2,
                                         const utils::Point& pt3, const utils::Point& pt4)
{
  std::array<utils::Point, 4> pts = {pt1, pt2, pt3, pt4};
  auto pt = interpolate_quad(pts.data(), utils::Point(xi[0], xi[1]));
  pt.z = 0;
  return pt;
}

utils::Point compute_quad_coords_from_xi_3d(MeshEntityPtr quad, const utils::Point& ptXi)
{
  assert(quad->get_type() == MeshEntityType::Quad);

  std::array<MeshEntityPtr, MAX_DOWN> verts;
  get_downward(quad, 0, verts.data());

  std::array<double, 2> ptXiArray = {ptXi.x, ptXi.y};
  return compute_quad_coords_from_xi_3d(ptXiArray.data(), verts[0]->get_point_orig(0), verts[1]->get_point_orig(0),
                                                          verts[2]->get_point_orig(0), verts[3]->get_point_orig(0));
}

utils::Point compute_quad_coords_from_xi_3d(const double xi[2], const utils::Point& pt1, const utils::Point& pt2,
                                            const utils::Point& pt3, const utils::Point& pt4)
{
  std::array<utils::Point, 4> pts = {pt1, pt2, pt3, pt4};
  auto pt = interpolate_quad(pts.data(), utils::Point(xi[0], xi[1]));
  return pt;
}

utils::Point compute_quad_centroid(MeshEntityPtr quad)
{
  assert(quad->get_type() == MeshEntityType::Quad);

  return compute_quad_centroid_3d(quad);
}

utils::Point compute_quad_centroid_3d(MeshEntityPtr quad)
{
  assert(quad->get_type() == MeshEntityType::Quad);

  MeshEntityPtr verts[MAX_DOWN];
  get_downward(quad, 0, verts);
  auto n1 = verts[0]->get_point_orig(0);
  auto n2 = verts[1]->get_point_orig(0);
  auto n3 = verts[2]->get_point_orig(0);
  auto n4 = verts[3]->get_point_orig(0);

  double xi[2] = {0.5, 0.5};
  return compute_quad_coords_from_xi_3d(xi, n1, n2, n3, n4);
}

utils::Point compute_quad_coords_from_xi(MeshEntityPtr quad, const utils::Point& pt)
{
  assert(quad->get_type() == MeshEntityType::Quad);

  MeshEntityPtr verts[MAX_DOWN];
  get_downward(quad, 0, verts);

  double xi[2] = {pt.x, pt.y};
  return compute_quad_coords_from_xi(xi, verts[0]->get_point_orig(0), verts[1]->get_point_orig(0),
                                     verts[2]->get_point_orig(0), verts[3]->get_point_orig(0));
}

// pt is the xi coordinates
// they are the coordinates associated with verts 1 and 2, with verts 0 being the
// dependent coordinate
utils::Point compute_tri_coords_from_xi(MeshEntityPtr tri, const utils::Point& pt)
{
  assert(tri->get_type() == MeshEntityType::Triangle);
  return compute_tri_coords_from_xi_3d(tri, pt);
}

utils::Point compute_tri_coords_from_xi_3d(MeshEntityPtr tri, const utils::Point& ptXi)
{
  assert(tri->get_type() == MeshEntityType::Triangle);

  MeshEntityPtr verts[MAX_DOWN];
  get_downward(tri, 0, verts);
  auto n1 = verts[0]->get_point_orig(0);
  auto n2 = verts[1]->get_point_orig(0);
  auto n3 = verts[2]->get_point_orig(0);

  std::array<utils::Point, 3> pts = {n1, n2, n3};
  return interpolate_triangle(pts.data(), ptXi);
}

utils::Point compute_tri_centroid(MeshEntityPtr tri)
{
  utils::Point pt = compute_tri_coords_from_xi(tri, utils::Point(1.0 / 3.0, 1.0 / 3.0));

  return pt;
}

utils::Point compute_tri_centroid_3d(MeshEntityPtr tri)
{
  utils::Point pt = compute_tri_coords_from_xi_3d(tri, utils::Point(1.0 / 3.0, 1.0 / 3.0));

  return pt;
}

utils::Point compute_edge_centroid3_d(MeshEntityPtr edge)
{
  assert(edge->get_type() == MeshEntityType::Edge);

  return (edge->get_down(0)->get_point_orig(0) + edge->get_down(1)->get_point_orig(0)) / 2;
}

/*
//TODO: Copied
utils::Point computeXiCoords(MeshEntityPtr el, const utils::Point& pt)
{
  switch(el->get_type())
  {
    case MeshEntityType::Triangle: { return computeTriXiCoords(el, pt);}
    case MeshEntityType::Quad:     { return computeQuadXiCoords(el, pt);}
    default:
      throw std::invalid_argument("unsupported entity type");
  }
}
*/

utils::Point compute_centroid(MeshEntityPtr el)
{
  switch (el->get_type())
  {
    case MeshEntityType::Vertex: {
      return el->get_point_orig(0);
    }
    case MeshEntityType::Edge: {
      return compute_edge_centroid3_d(el);
    }
    case MeshEntityType::Triangle: {
      return compute_tri_centroid(el);
    }
    case MeshEntityType::Quad: {
      return compute_quad_centroid(el);
    }
    default:
      throw std::invalid_argument("unsupported entity type");
  }
}

utils::Point compute_centroid_3d(MeshEntityPtr el)
{
  switch (el->get_type())
  {
    case MeshEntityType::Vertex: {
      return el->get_point_orig(0);
    }
    case MeshEntityType::Edge: {
      return compute_edge_centroid3_d(el);
    }
    case MeshEntityType::Triangle: {
      return compute_tri_centroid_3d(el);
    }
    case MeshEntityType::Quad: {
      return compute_quad_centroid_3d(el);
    }
    default:
      throw std::invalid_argument("unsupported entity type");
  }
}

utils::Point compute_coords_from_xi(MeshEntityPtr el, const utils::Point& pt)
{
  switch (el->get_type())
  {
    case MeshEntityType::Triangle: {
      return compute_tri_coords_from_xi(el, pt);
    }
    case MeshEntityType::Quad: {
      return compute_quad_coords_from_xi(el, pt);
    }
    default:
      throw std::invalid_argument("unsupported entity type");
  }
}

utils::Point compute_coords_from_xi_3d(MeshEntityPtr el, const utils::Point& pt)
{
  switch (el->get_type())
  {
    case MeshEntityType::Triangle: {
      return compute_tri_coords_from_xi_3d(el, pt);
    }
    case MeshEntityType::Quad: {
      return compute_quad_coords_from_xi_3d(el, pt);
    }
    default:
      throw std::invalid_argument("unsupported entity type");
  }
}

double convert_xi_coords_from_range(double xiStart, double xiEnd, double xi)
{
  assert(xiEnd > xiStart);
  return (xi - xiStart)/(xiEnd - xiStart);
}

utils::Point convert_xi_coords_from_range(double xiStart, double xiEnd, const utils::Point& ptXi)
{
  return {convert_xi_coords_from_range(xiStart, xiEnd, ptXi.x),
          convert_xi_coords_from_range(xiStart, xiEnd, ptXi.y)};
}

double convert_xi_coords_to_range(double xiStart, double xiEnd, double xi)
{
  assert(xiEnd > xiStart);
  return xi * (xiEnd - xiStart) + xiStart;
}

utils::Point convert_xi_coords_to_range(double xiStart, double xiEnd, const utils::Point& ptXi)
{
  return {convert_xi_coords_to_range(xiStart, xiEnd, ptXi.x),
          convert_xi_coords_to_range(xiStart, xiEnd, ptXi.y)};
}

int compute_angles(MeshEntityPtr el, double angles[MAX_DOWN])
{
  MeshEntityPtr verts[MAX_DOWN];
  int nverts = get_downward(el, 0, verts);

  for (int i = 0; i < nverts; ++i)
  {
    utils::Point v0 = verts[(i - 1 + nverts) % nverts]->get_point_orig(0);
    utils::Point v1 = verts[i]->get_point_orig(0);
    utils::Point v2 = verts[(i + 1) % nverts]->get_point_orig(0);

    utils::Point r1 = v0 - v1;
    utils::Point r2 = v2 - v1;
    double arg      = dot(r2, r1) / (std::sqrt(dot(r1, r1)) * std::sqrt(dot(r2, r2)));
    angles[i]       = std::acos(arg);
  }

  return nverts;
}

bool is_unique(std::vector<MeshEntityPtr> entities)
{
  // this is only used for small vectors, so use naive algorithm
  for (unsigned int i = 0; i < entities.size(); ++i)
    for (unsigned int j = 0; j < entities.size(); ++j)
      if (i != j && entities[i] == entities[j])
        return false;

  return true;
}

bool is_null(std::vector<MeshEntityPtr> entities)
{
  for (auto& el : entities)
    if (!el)
      return true;

  return false;
}

void check_vertices_null(std::vector<MeshEntityPtr> entities)
{
  MeshEntityPtr verts[MAX_DOWN];
  for (auto& e : entities)
    if (e)
    {
      int nverts = get_downward(e, 0, verts);
      for (int i = 0; i < nverts; ++i)
        if (!verts[i])
          throw std::invalid_argument("found null vertex");
    }
}

std::shared_ptr<Mesh> make_empty_mesh(MPI_Comm comm)
{
  Mesh* m = new Mesh(comm);
  std::shared_ptr<Mesh> mp(m);
  auto f = create_field<impl::GeoClassification>(mp, FieldShape(1, 1, 1), 1, impl::GeoClassification(), true);
  mp->set_geo_classification(f);
  // auto mesh = std::make_shared<Mesh>();
  return mp;
}

int count_entities_of_type(std::shared_ptr<mesh::Mesh> mesh, MeshEntityType type)
{
  int count = 0;
  int dim = get_type_dimension(type);
  for (auto& entity : mesh->get_mesh_entities(dim))
    if (entity && entity->get_type() == type)
      count++;

  return count;
}

} // namespace mesh

} // namespace middle_mesh
} // namespace stk
