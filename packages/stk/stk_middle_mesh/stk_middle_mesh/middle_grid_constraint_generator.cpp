#include "middle_grid_constraint_generator.hpp"
#include "predicates/point_classifier_normal_wrapper.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

void MiddleGridConstraintGenerator::generate()
{
  create_mesh1_vertices();
  create_mesh2_interior_vertices();
  create_mesh1_edges();
  split_edges();
  create_internal_edges();
}

void MiddleGridConstraintGenerator::create_mesh1_vertices()
{
  if (m_output)
    std::cout << "creating mesh1 vertices" << std::endl;

  auto& verts1ToFakeVerts   = *(m_relationalData->verts1ToFakeVerts);
  auto& vertsInClassOnMesh1 = *(m_relationalData->vertsInClassOnMesh1);
  auto& fakeVertsToVertsIn  = m_relationalData->fakeVertsToVertsIn;
  auto& mesh1ElsToVertsIn   = *(m_relationalData->mesh1ElsToVertsIn);

  std::vector<mesh::MeshEntityPtr> els1;

  for (auto& vert1 : m_mesh1->get_vertices())
    if (vert1)
    {
      FakeVert fv = verts1ToFakeVerts(vert1, 0, 0);
      auto vertIn = m_meshIn->create_vertex(fv.pt);

      if (m_output && m_vertIds.count(vertIn->get_id()) > 0)
        std::cout << "created vert_in id " << vertIn->get_id() << " from fakevert id " << fv.id << std::endl;

      mesh::MeshEntityPtr el1 = vert1->get_up(0)->get_up(0);
      int localId             = predicates::impl::get_entity_id(el1, vert1);
      vertsInClassOnMesh1(vertIn, 0, 0) = m_pointClassifier->create_vert_record(el1, localId);
          //predicates::impl::PointRecord(predicates::impl::PointClassification::Vert, localId, el1);
      fakeVertsToVertsIn[fv.id] = vertIn;

      int nels = get_upward(vert1, 2, els1);
      for (int i = 0; i < nels; ++i)
      {
        mesh1ElsToVertsIn.insert(els1[i], 0, vertIn);
        if (m_output && m_vertIds.count(vertIn->get_id()) > 0)
          std::cout << "adding vert_in id " << vertIn->get_id() << " to el1 with id " << els1[i]->get_id() << std::endl;
      }
    }
}

void MiddleGridConstraintGenerator::create_mesh2_interior_vertices()
{
  if (m_output)
    std::cout << "creating mesh2 interior vertices" << std::endl;
  auto& verts2ClassOnMesh1  = *(m_relationalData->verts2ClassOnMesh1);
  auto& verts2ToFakeVerts   = *(m_relationalData->verts2ToFakeVerts);
  auto& vertsInClassOnMesh1 = *(m_relationalData->vertsInClassOnMesh1);
  auto& fakeVertsToVertsIn  = m_relationalData->fakeVertsToVertsIn;
  auto& mesh1ElsToVertsIn   = *(m_relationalData->mesh1ElsToVertsIn);

  for (auto& vert2 : m_mesh2->get_vertices())
    if (vert2)
    {
      predicates::impl::PointRecord& record = verts2ClassOnMesh1(vert2, 0, 0);
      if (record.type == predicates::impl::PointClassification::Interior)
      {
        FakeVert fv = verts2ToFakeVerts(vert2, 0, 0);
        auto vertIn = m_meshIn->create_vertex(fv.pt);

        if (m_output && m_vertIds.count(vertIn->get_id()) > 0)
          std::cout << "created vert_in id " << vertIn->get_id() << " from fakevert id " << fv.id << std::endl;

        vertsInClassOnMesh1(vertIn, 0, 0) = record;
        fakeVertsToVertsIn[fv.id]         = vertIn;
        mesh1ElsToVertsIn.insert(record.el, 0, vertIn);
        if (m_output && m_vertIds.count(vertIn->get_id()) > 0)
          std::cout << "adding vert_in id " << vertIn->get_id() << " to el1 with id " << record.el->get_id()
                    << std::endl;
      }
    }
}

void MiddleGridConstraintGenerator::create_mesh1_edges()
{
  auto& verts1ToFakeVerts = *(m_relationalData->verts1ToFakeVerts);
  for (auto& edge : m_mesh1->get_edges())
    if (edge)
    {
      FakeVert fv1             = verts1ToFakeVerts(edge->get_down(0), 0, 0);
      FakeVert fv2             = verts1ToFakeVerts(edge->get_down(1), 0, 0);
      mesh::MeshEntityPtr v1In = m_relationalData->fakeVertsToVertsIn[fv1.id];
      mesh::MeshEntityPtr v2In = m_relationalData->fakeVertsToVertsIn[fv2.id];
      assert(!stk::middle_mesh::mesh::get_common_edge(v1In, v2In));

      m_meshIn->create_edge(v1In, v2In);
    }
}

void MiddleGridConstraintGenerator::split_edges()
{
  if (m_output)
    std::cout << "splitting edges" << std::endl;
  auto& mesh1EdgesToSplit   = *(m_relationalData->mesh1EdgesToSplit);
  auto& fakeVertsToVertsIn  = m_relationalData->fakeVertsToVertsIn;
  auto& vertsInClassOnMesh1 = *(m_relationalData->vertsInClassOnMesh1);
  auto& mesh1ElsToVertsIn   = *(m_relationalData->mesh1ElsToVertsIn);

  for (auto& edge1 : m_mesh1->get_edges())
    if (edge1 && mesh1EdgesToSplit.get_num_comp(edge1, 0) > 0)
    {
      sort_edge_splits(edge1);

      mesh::MeshEntityPtr el1 = edge1->get_up(0);
      int localId             = predicates::impl::get_entity_id(el1, edge1);

      auto currEdgeIn = get_mesh_in_edge_from_mesh1_edge(edge1);
      double xiStart  = 0;
      mesh::MeshEntityPtr newEdges[2];

      // std::cout << "about to split edge at xi: ";
      // for (int i=0; i < mesh1_edges_to_split.get_num_comp(edge1, 0); ++i)
      //{
      //   std::cout << mesh1_edges_to_split(edge1, 0, i).xi << ", ";
      // }
      // std::cout << std::endl;

      for (int i = 0; i < mesh1EdgesToSplit.get_num_comp(edge1, 0); ++i)
      {
        double xi   = mesh1EdgesToSplit(edge1, 0, i).xi;
        FakeVert fv = mesh1EdgesToSplit(edge1, 0, i).vert;

        // map xi from edge1 to curr_edge_in
        double xiNew = (xi - xiStart) / (1 - xiStart);
        auto newVert = m_meshIn->split_edge_broken(currEdgeIn, xiNew, newEdges);

        if (m_output && (m_vertIds.count(newVert->get_id()) > 0))
          std::cout << "created vert_in id " << newVert->get_id() << " from fakevert id " << fv.id << std::endl;

        fakeVertsToVertsIn[fv.id] = newVert;

        double edgeXiOnReferenceEl = el1->get_down_orientation(localId) == mesh::EntityOrientation::Standard ? xi : 1 - xi;
        vertsInClassOnMesh1(newVert, 0, 0) = m_pointClassifier->create_edge_record(el1, localId, edgeXiOnReferenceEl);
        for (int j = 0; j < edge1->count_up(); ++j)
        {
          mesh1ElsToVertsIn.insert(edge1->get_up(j), 0, newVert);
          if (m_output && m_vertIds.count(newVert->get_id()) > 0)
            std::cout << "adding vert_in id " << newVert->get_id() << " to el1 with id " << edge1->get_up(j)->get_id()
                      << std::endl;
        }

        // update for next iteration
        xiStart    = xi;
        currEdgeIn = newEdges[1];
      }
    }
}

void MiddleGridConstraintGenerator::sort_edge_splits(mesh::MeshEntityPtr edge1)
{
  auto& mesh1EdgesToSplit = *(m_relationalData->mesh1EdgesToSplit);
  auto cmp = [](const EdgeSplitRecord& split1, const EdgeSplitRecord& split2) { return split1.xi < split2.xi; };

  int nsplits                  = mesh1EdgesToSplit.get_num_comp(edge1, 0);
  EdgeSplitRecord* splitsStart = &(mesh1EdgesToSplit(edge1, 0, 0));
  EdgeSplitRecord* splitsEnd   = &(mesh1EdgesToSplit(edge1, 0, nsplits - 1)) + 1;
  std::sort(splitsStart, splitsEnd, cmp);
}

mesh::MeshEntityPtr MiddleGridConstraintGenerator::get_mesh_in_edge_from_mesh1_edge(mesh::MeshEntityPtr edge1)
{
  auto& verts1ToFakeVerts = *(m_relationalData->verts1ToFakeVerts);

  FakeVert fv1             = verts1ToFakeVerts(edge1->get_down(0), 0, 0);
  FakeVert fv2             = verts1ToFakeVerts(edge1->get_down(1), 0, 0);
  mesh::MeshEntityPtr v1In = m_relationalData->fakeVertsToVertsIn[fv1.id];
  mesh::MeshEntityPtr v2In = m_relationalData->fakeVertsToVertsIn[fv2.id];

  return get_common_edge(v1In, v2In);
}

void MiddleGridConstraintGenerator::create_internal_edges()
{
  auto& edges2ToFakeVertsIn = *(m_relationalData->edges2ToFakeVertsIn);

  for (auto& edge2 : m_mesh2->get_edges())
    if (edge2)
    {
      assert(edges2ToFakeVertsIn.get_num_comp(edge2, 0) >= 2);
      // sortVertsOnEdge2(edge2);  //TODO: make this an assert?

      for (int i = 0; i < edges2ToFakeVertsIn.get_num_comp(edge2, 0) - 1; ++i)
      {
        FakeVert fv1             = edges2ToFakeVertsIn(edge2, 0, i).vert;
        FakeVert fv2             = edges2ToFakeVertsIn(edge2, 0, i + 1).vert;
        mesh::MeshEntityPtr v1In = m_relationalData->fakeVertsToVertsIn[fv1.id];
        mesh::MeshEntityPtr v2In = m_relationalData->fakeVertsToVertsIn[fv2.id];

        if (m_output && (v1In->get_id() == 1916 || v2In->get_id() == 1916))
        {
          std::cout << "found possible edge involving vert 1916" << std::endl;
          std::cout << "vert ids = " << v1In->get_id() << ", " << v2In->get_id() << std::endl;

          for (int j = 0; j < edges2ToFakeVertsIn.get_num_comp(edge2, 0); ++j)
            std::cout << "fake vert " << j << " on specified edge has id " << edges2ToFakeVertsIn(edge2, 0, j).vert.id
                      << ", xi = " << edges2ToFakeVertsIn(edge2, 0, j).xi << std::endl;
        }

        if (!get_common_edge(v1In, v2In))
          m_meshIn->create_edge(v1In, v2In);
      }
    }
}

void MiddleGridConstraintGenerator::sort_verts_on_edge2(mesh::MeshEntityPtr edge2)
{
  auto& edges2ToFakeVertsIn = *(m_relationalData->edges2ToFakeVertsIn);
  auto cmp                  = [](const VertOnEdge& v1, const VertOnEdge& v2) { return v1.xi < v2.xi; };

  int nsplits            = edges2ToFakeVertsIn.get_num_comp(edge2, 0);
  VertOnEdge* vertsStart = &(edges2ToFakeVertsIn(edge2, 0, 0));
  VertOnEdge* vertsEnd   = &(edges2ToFakeVertsIn(edge2, 0, nsplits - 1)) + 1;
  std::sort(vertsStart, vertsEnd, cmp);
}

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
