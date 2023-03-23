#include "element_mesh_triangulator.hpp"

// TODO: can go in src file
#include "CDTInterface.hpp"
#include "mesh_io.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

using predicates::impl::PointClassification;

void ElementMeshTriangulator::triangulate(ElementMeshData elementMeshData)
{
  m_elementMeshData      = elementMeshData;
  int numConstraintEdges = count_valid(m_elementMeshData.elementMeshIn->get_edges());
  perturb_boundary_nodes(1);
  triangulate_element_mesh();
  perturb_boundary_nodes(-1);
  fix_triangulation(numConstraintEdges);
  fix_triangulation_with_boundary_ordering();
#ifndef NDEBUG
  check_for_cap_elements();
#endif
}

void ElementMeshTriangulator::perturb_boundary_nodes(double fac)
{
  // when fac > 0, offsets all points classified on edges such that the
  // CDT predicate thinks the points is on the interior.  When this is true,
  // the fixTriangulationD will be able to delete the spurious triangles
  // Call again with fac < 0 to undo the perturbation
  if (M_OUTPUT)
  {
    std::cout << "\nEntered offsetEdgePoints" << std::endl;
    std::cout << "fac = " << fac << std::endl;
  }

  double eps = fac * 1e-16;
  std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
  int nverts            = get_downward(m_elementMeshData.el1, 0, verts.data());
  utils::Point centroid = compute_centroid_3d(m_elementMeshData.el1);
  utils::Point offsets[mesh::MAX_DOWN];
  for (int i = 0; i < nverts; ++i)
  {
    utils::Point pt1       = verts[i]->get_point_orig(0);
    utils::Point pt2       = verts[(i + 1) % nverts]->get_point_orig(0);
    utils::Point midpoint  = (pt1 + pt2) / 2;
    utils::Point direction = centroid - midpoint;
    offsets[i]             = eps * (direction / std::sqrt(dot(direction, direction)));
  }

  for (auto& vertIn : m_elementMeshData.elementMeshIn->get_vertices())
  {
    predicates::impl::PointRecord& record = (*m_elementMeshData.elementMeshInVertClass)(vertIn, 0, 0);
    if (record.type == PointClassification::Edge)
      vertIn->set_point_orig(0, vertIn->get_point_orig(0) + offsets[record.id]);
  }
}

void ElementMeshTriangulator::triangulate_element_mesh()
{
  assert(m_elementMeshData.elementMeshIn->get_elements().size() == 0);

  std::array<mesh::MeshEntityPtr, 4> el1Verts;
  int nverts = get_downward(m_elementMeshData.el1, 0, el1Verts.data());
  utils::impl::Projection proj(el1Verts[0]->get_point_orig(0), el1Verts[1]->get_point_orig(0),
                               el1Verts[nverts - 1]->get_point_orig(0));

  CDTInterface cdt(m_elementMeshData.elementMeshIn);
  cdt.triangulate(proj);

  if (M_OUTPUT)
  {
    mesh::impl::print_vert_edges("mesh_in_triangulated", m_elementMeshData.elementMeshIn);
  }
}

void ElementMeshTriangulator::fix_triangulation(int numConstraintEdges)
{
  // find and delete "cap" elements: elements on the boundary of the
  // triangulation where the edge on the boundary is not a constrained edge.
  // These can occur when three points on the boundary are *almost* colinear,
  // but the middle point is slightly inside the line formed by the other two.
  // The Delauney triangulation has to form an edge between the other two
  // vertices to satisfy the convex hull property, but we don't want that
  // element to exist for our purposes
  int ndeleted = -1;
  while (ndeleted != 0)
  {
    // TODO: do adjacency search to find all elements to delete in one go
    ndeleted = 0;
    for (auto& edge : m_elementMeshData.elementMeshIn->get_edges())
    {
      if (edge && edge->get_id() >= numConstraintEdges && edge->count_up() == 1)
      {
        m_elementMeshData.elementMeshIn->delete_face(edge->get_up(0));
        ndeleted += 1;
      }
    }
  }
}

bool ElementMeshTriangulator::should_delete_edges(const predicates::impl::PointRecord& class1,
                                                  const predicates::impl::PointRecord& class2,
                                                  mesh::MeshEntityPtr v1ElementMeshIn,
                                                  mesh::MeshEntityPtr v2ElementMeshIn)
{
  assert(v1ElementMeshIn != v2ElementMeshIn);

  auto& mesh1EdgesToSplit             = *(m_relationalData->mesh1EdgesToSplit);
  auto& fakeVertsToRealVertsIn        = m_relationalData->fakeVertsToVertsIn;
  auto& elementMeshVertsToMeshInVerts = *m_elementMeshData.elementMeshEntitiesToMeshInEntities;

  auto edge1               = get_entity(class1);
  auto edge2               = get_entity(class2);
  mesh::MeshEntityPtr vIn1 = elementMeshVertsToMeshInVerts(v1ElementMeshIn, 0, 0);
  mesh::MeshEntityPtr vIn2 = elementMeshVertsToMeshInVerts(v2ElementMeshIn, 0, 0);

  if (edge1 != edge2)
    return false;

  int idx1 = -1, idx2 = -1;
  for (int i = 0; i < mesh1EdgesToSplit.get_num_comp(edge1, 0); ++i)
  {
    EdgeSplitRecord edgeSplit  = mesh1EdgesToSplit(edge1, 0, i);
    mesh::MeshEntityPtr vertIn = fakeVertsToRealVertsIn[edgeSplit.vert.id];

    if (vertIn == vIn1)
    {
      assert(idx1 == -1);
      idx1 = i;
    }

    if (vertIn == vIn2)
    {
      assert(idx2 == -1);
      idx2 = i;
    }

    if (idx1 != -1 && idx2 != -1)
      break;
  }

  assert(idx1 != -1);
  assert(idx2 != -1);

  return std::abs(idx1 - idx2) != 1;
}

bool ElementMeshTriangulator::should_delete_edge_between_vert_and_edge(const predicates::impl::PointRecord& classV,
                                                                       const predicates::impl::PointRecord& classE,
                                                                       mesh::MeshEntityPtr elementMeshVertOnMesh1Edge)
{
  assert(classV.type == PointClassification::Vert);
  assert(classE.type == PointClassification::Edge);

  auto& elementMeshVertsToMeshInVerts = *m_elementMeshData.elementMeshEntitiesToMeshInEntities;
  auto& fakeVertsToRealVertsIn        = m_relationalData->fakeVertsToVertsIn;
  auto& mesh1EdgesToSplit             = *(m_relationalData->mesh1EdgesToSplit);

  mesh::MeshEntityPtr meshInVertOnMesh1Edge = elementMeshVertsToMeshInVerts(elementMeshVertOnMesh1Edge, 0, 0);

  int id1     = classV.id;
  auto edge2  = get_entity(classE);
  int id2     = classE.id;
  int nsplits = mesh1EdgesToSplit.get_num_comp(edge2, 0);

  // figure out which edge split caused mesh_in_vert_on_mesh1_edge to be created
  int idx2 = -1;
  for (int i = 0; i < nsplits; ++i)
  {
    EdgeSplitRecord edgeSplit = mesh1EdgesToSplit(edge2, 0, i);
    if (fakeVertsToRealVertsIn[edgeSplit.vert.id] == meshInVertOnMesh1Edge)
    {
      idx2 = i;
      break;
    }
  }

  assert(idx2 != -1);

  // figure out if v1 is first or last vertex on edge2, in the
  // coordinate system of the edge, not el1
  // first -> index of -1, last -> index splits.size()
  int idx1;
  if (m_elementMeshData.el1->get_down_orientation(id2) == mesh::EntityOrientation::Standard)
    idx1 = id1 == id2 ? -1 : nsplits;
  else
    idx1 = id1 == id2 ? nsplits : -1;

  return std::abs(idx1 - idx2) != 1;
}

void ElementMeshTriangulator::fix_triangulation_with_boundary_ordering()
{
  auto& mesh1EdgesToSplit = *(m_relationalData->mesh1EdgesToSplit);

  if (M_OUTPUT)
    std::cout << "\nEntered fixTriangulationD" << std::endl;
  // TODO: is this necessary?  With the edge offsetting the simpler test used
  //       in fixTriangulation might be sufficient
  int ndeleted = -1;
  while (ndeleted != 0)
  {
    ndeleted = 0;
    for (auto& edge : m_elementMeshData.elementMeshIn->get_edges())
    {
      if (edge && edge->count_up() == 1)
      {
        auto v1 = edge->get_down(0);
        auto v2 = edge->get_down(1);

        auto class1 = (*m_elementMeshData.elementMeshInVertClass)(v1, 0, 0);
        auto class2 = (*m_elementMeshData.elementMeshInVertClass)(v2, 0, 0);

        std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> elVerts;
        get_downward(edge->get_up(0), 0, elVerts.data());

        mesh::MeshEntityPtr vertsIn[mesh::MAX_DOWN];
        get_downward(edge->get_up(0), 0, vertsIn);

        if (class1.type == PointClassification::Vert && class2.type == PointClassification::Vert)
        {
          if (std::abs(class1.id - class2.id) != 1)
            continue;

          assert(std::abs(class1.id - class2.id) == 1);
          int minId = std::min(class1.id, class2.id);
          // id of left vertex == edge id for triangles and quads
          auto edge1 = m_elementMeshData.el1->get_down(minId);

          if (mesh1EdgesToSplit.get_num_comp(edge1, 0) != 0)
          {
            // if we get here, then: there is an edge between two vertices of the mesh1 element even
            // though that edge was split.  Thus the face connected to that edge should be deleted
            m_elementMeshData.elementMeshIn->delete_face(edge->get_up(0));
            ndeleted += 1;
          }

          continue;
        }

        if (class1.type == PointClassification::Edge && class2.type == PointClassification::Edge)
        {
          if (should_delete_edges(class1, class2, v1, v2))
          {
            m_elementMeshData.elementMeshIn->delete_face(edge->get_up(0));
            ndeleted += 1;
          }

          continue;
        }

        // only one is on a vertex
        bool shouldDelete = false;
        if (class1.type == PointClassification::Vert)
          shouldDelete = should_delete_edge_between_vert_and_edge(class1, class2, v2);
        else // class2_type == Edge
          shouldDelete = should_delete_edge_between_vert_and_edge(class2, class1, v1);

        if (shouldDelete)
        {
          m_elementMeshData.elementMeshIn->delete_face(edge->get_up(0));
          ndeleted += 1;
        }
      }
    }
  }
}

void ElementMeshTriangulator::check_for_cap_elements()
{
  auto& mesh1EdgesToSplit = *(m_relationalData->mesh1EdgesToSplit);

  for (auto& edgeIn : m_elementMeshData.elementMeshIn->get_edges())
    if (edgeIn)
    {
      auto v1     = edgeIn->get_down(0);
      auto v2     = edgeIn->get_down(1);
      auto class1 = (*m_elementMeshData.elementMeshInVertClass)(v1, 0, 0);
      auto class2 = (*m_elementMeshData.elementMeshInVertClass)(v2, 0, 0);

      if (class1.type == PointClassification::Interior || class2.type == PointClassification::Interior)
      {
        continue;
      }

      if (class1.type == PointClassification::Vert && class2.type == PointClassification::Vert)
      {
        mesh::MeshEntityPtr edge1 = get_common_edge(get_entity(class1), get_entity(class2));
        if (edge1)
        {
          if (mesh1EdgesToSplit.get_num_comp(edge1, 0) != 0)
            throw std::runtime_error("found cap element");
        }
      } else if (class1.type == PointClassification::Edge && class2.type == PointClassification::Edge)
      {
        mesh::MeshEntityPtr vert1Edge = get_entity(class1);
        mesh::MeshEntityPtr vert2Edge = get_entity(class2);
        if (vert1Edge == vert2Edge)
        {
          int idx1 = get_edge_split_index(vert1Edge, v1);
          int idx2 = get_edge_split_index(vert2Edge, v2);
          if (std::abs(idx1 - idx2) != 1)
            throw std::runtime_error("found cap element");
        }
      } else
      {
        // if we get here, 1 vert is classified on an edge and the other is on a vert
        // mesh::MeshEntityPtr vert_on_vert = class1.type == PointClassification::Vert ? v1 : v2;
        mesh::MeshEntityPtr vertOnEdge                = class1.type == PointClassification::Vert ? v2 : v1;
        predicates::impl::PointRecord vertOnVertClass = class1.type == PointClassification::Vert ? class1 : class2;
        predicates::impl::PointRecord vertOnEdgeClass = class1.type == PointClassification::Vert ? class2 : class1;
        assert(vertOnEdgeClass.type == PointClassification::Edge);

        mesh::MeshEntityPtr vert1 = get_entity(vertOnVertClass);
        mesh::MeshEntityPtr edge1 = get_entity(vertOnEdgeClass);
        int idxEdge               = get_edge_split_index(edge1, vertOnEdge);
        int idxVert;
        if (edge1->get_down(0) == vert1)
          idxVert = -1;
        else if (edge1->get_down(1) == vert1)
          idxVert = mesh1EdgesToSplit.get_num_comp(edge1, 0);
        else
        {
          continue;
        }

        if (std::abs(idxVert - idxEdge) != 1)
          throw std::runtime_error("found cap element");
      }
    }
}

int ElementMeshTriangulator::get_edge_split_index(mesh::MeshEntityPtr edge1, mesh::MeshEntityPtr elementMeshVert)
{
  auto& mesh1EdgesToSplit                  = *(m_relationalData->mesh1EdgesToSplit);
  auto& elementMeshEntitesToMeshInEntities = *(m_elementMeshData.elementMeshEntitiesToMeshInEntities);
  mesh::MeshEntityPtr vertIn               = elementMeshEntitesToMeshInEntities(elementMeshVert, 0, 0);

  for (int i = 0; i < mesh1EdgesToSplit.get_num_comp(edge1, 0); ++i)
  {
    FakeVert fv                 = mesh1EdgesToSplit(edge1, 0, i).vert;
    mesh::MeshEntityPtr vertInI = m_relationalData->fakeVertsToVertsIn[fv.id];
    if (vertIn == vertInI)
      return i;
  }

  throw std::runtime_error("unable to find edge split for given vertex");
}

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
