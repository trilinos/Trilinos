#include "mesh_projection_calculator.h"

#include "adjacency_search.h"
#include "mesh_entity.h"
#include <limits>

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

using predicates::impl::PointClassification;

std::shared_ptr<mesh::Mesh> MeshProjectionCalculator::project()
{
  create_mesh1_fake_verts();
  project_mesh2_vertices_onto_mesh1();
  intersect_mesh2_edges_with_mesh1();
  // sortEdge2FakeVerts();

  m_relationalData->fakeVertsToVertsIn.resize(m_fakeVertGen.get_num_verts());

  return m_meshIn;
}

void MeshProjectionCalculator::create_mesh1_fake_verts()
{
  if (m_output)
    std::cout << "\ncreating mesh1 fake verts" << std::endl;
  // TODO: this is not great: this takes a bunch of memory that
  //       is a copy of existing data.  Maybe create the mesh_in
  //       verts here and find a way to use them when a mesh2
  //       edge overlaps a mesh1 vertex?
  auto& verts1ToFakeVerts = *(m_relationalData->verts1ToFakeVerts);
  for (auto& vert1 : m_mesh1->get_vertices())
    if (vert1)
    {
      FakeVert fv                    = m_fakeVertGen.get_vert(vert1->get_point_orig(0));
      verts1ToFakeVerts(vert1, 0, 0) = fv;

      if (m_output && m_vertIds.count(fv.id) > 0)
        std::cout << "created fakevert id " << fv.id << " from mesh1 vert id " << vert1->get_id() << std::endl;
    }
}

void MeshProjectionCalculator::project_mesh2_vertices_onto_mesh1()
{
  if (m_output)
    std::cout << "\nprojecting mesh1 vertices onto mesh1 " << std::endl;
  // TODO: for memory savings, could
  //   1. put first result into the regular field
  //   2. when second result is found, copy them both to variable size field
  //   3. Need to initialize regular field to type == Exterior, so we can detect
  //      vertices that don't project to anywhere
  mesh::VariableSizeFieldPtr<predicates::impl::PointRecord> verts2AllClassificationsPtr =
      mesh::create_variable_size_field<predicates::impl::PointRecord>(m_mesh2, mesh::impl::FieldShape(1, 0, 0));

  mesh::MeshEntityPtr el1;
  std::vector<mesh::MeshEntityPtr> mesh2Els;
  SetType<mesh::MeshEntityPtr> seenVerts; // TODO: this isn't quite what we want: there can
                                          //       be multiple equivalent result if a vertex
                                          //       projects onto a vertex or edge
  auto& verts2AllClassifications = *verts2AllClassificationsPtr;
  auto& mesh1ElToMesh2Els        = *(m_relationalData->mesh1ElementsToMesh2Elements);
  mesh::impl::AdjacencySearch search(m_mesh1, m_mesh2);
  while ((el1 = search.get_next(mesh2Els)))
  {
    seenVerts.clear();
    std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
    for (auto& el2 : mesh2Els)
    {
      int nverts = get_downward(el2, 0, verts.data());
      for (int i = 0; i < nverts; ++i)
        if (seenVerts.count(verts[i]) == 0)
        {
          predicates::impl::PointRecord record = m_pointClassifier->classify(el1, el2, verts[i]->get_point_orig(0));

          if (record.type != PointClassification::Exterior)
          {
            verts2AllClassifications.insert(verts[i], 0, record);
            if (m_output && verts[i]->get_id() == 1001)
            {
              std::cout << "classified vertex with id " << verts[i]->get_id() << " on element " << el1
                        << ", record = " << record << std::endl;
              utils::Point ptOrig = verts[i]->get_point_orig(0);
              utils::Point ptProj = m_pointClassifier->compute_xyz_coords(record);
              std::cout << "vertex coords = " << ptOrig << ", projected coordinates = " << ptProj << std::endl;
              std::cout << "distance = " << std::sqrt(dot(ptOrig - ptProj, ptOrig - ptProj)) << std::endl;
            }
          }

          seenVerts.insert(verts[i]);
        }
    }

    for (auto& el2 : mesh2Els)
      mesh1ElToMesh2Els.insert(el1, 0, el2);
  }

  choose_unique_vert_projection(verts2AllClassificationsPtr);
}

void MeshProjectionCalculator::choose_unique_vert_projection(
    mesh::VariableSizeFieldPtr<predicates::impl::PointRecord> verts2AllClassificationsPtr)
{
  auto& verts2AllClassifications = *verts2AllClassificationsPtr;
  for (auto& vert2 : m_mesh2->get_vertices())
    if (vert2)
    {
      int nprojections = verts2AllClassifications.get_num_comp(vert2, 0);
      if (nprojections > 0)
      {
        utils::Point closestPt;
        int minIdx = get_closest_projection(verts2AllClassificationsPtr, vert2, closestPt);
        record_vert2_classification(vert2, verts2AllClassifications(vert2, 0, minIdx), closestPt);
      } else
        throw std::runtime_error("vertex does not projection onto any mesh1 element");
    }
}

void MeshProjectionCalculator::record_vert2_classification(mesh::MeshEntityPtr vert2,
                                                           const predicates::impl::PointRecord& record,
                                                           const utils::Point& pt)
{
  auto& verts2ClassOnMesh1  = *(m_relationalData->verts2ClassOnMesh1);
  auto& verts1ToFakeVerts   = *(m_relationalData->verts1ToFakeVerts);
  auto& verts2ToFakeVerts   = *(m_relationalData->verts2ToFakeVerts);
  auto& edges2ToFakeVertsIn = *(m_relationalData->edges2ToFakeVertsIn);
  auto& mesh1EdgesToSplit   = *(m_relationalData->mesh1EdgesToSplit);

  FakeVert fv;
  if (record.type == PointClassification::Vert)
    fv = verts1ToFakeVerts(get_entity(record), 0, 0);
  else
    fv = m_fakeVertGen.get_vert(pt);

  if (m_output && m_vertIds.count(fv.id) > 0)
  {
    std::cout << "created fakevert id " << fv.id << " from mesh2 vert id " << vert2->get_id() << std::endl;
    std::cout << "record = " << record << ", el = " << record.el << std::endl;
  }

  verts2ClassOnMesh1(vert2, 0, 0) = record;
  verts2ToFakeVerts(vert2, 0, 0)  = fv;
  record_edge2_vertex_association(edges2ToFakeVertsIn, vert2, fv);
  /*
    for (int i=0; i < vert2->count_up(); ++i)
    {
      mesh::MeshEntityPtr edge2 = vert2->get_up(i);
      double xi = edge2->get_down(0) == vert2 ? 0 : 1;
      edges2_to_fake_verts_in.insert(edge2, 0, {fv, xi});

      if (edges2_to_fake_verts_in.get_num_comp(edge2, 0) == 0)
      {
        edges2_to_fake_verts_in.insert(edge2, 0);
        edges2_to_fake_verts_in.insert(edge2, 0);
      }

      int idx = edge2->get_down(0) == vert2 ? 0 : 1;
      assert(edge2->get_down(idx) == vert2);
      edges2_to_fake_verts_in(edge2, 0, idx) = VertOnEdge{fv, xi};

    }
  */

  if (record.type == PointClassification::Edge)
  {
    double xi                = m_pointClassifier->get_edge_xi(record);
    mesh::MeshEntityPtr edge = get_entity(record);
    mesh1EdgesToSplit.insert(edge, 0, EdgeSplitRecord{fv, xi});

    if (m_output && m_vertIds.count(fv.id) > 0)
    {
      std::cout << "fakevert from vert-edge intersection" << std::endl;
      std::cout << "created fakevert id " << fv.id << std::endl;
      std::cout << "vert2 = " << vert2 << std::endl;
      std::cout << "point record = " << record << ", el id " << record.el->get_id() << std::endl;
    }
  }
}

void MeshProjectionCalculator::record_edge2_vertex_association(mesh::VariableSizeField<VertOnEdge>& edges2ToFakeVertsIn,
                                                               mesh::MeshEntityPtr vert2, FakeVert fv)
{
  for (int i = 0; i < vert2->count_up(); ++i)
  {
    mesh::MeshEntityPtr edge2 = vert2->get_up(i);
    double xi                 = edge2->get_down(0) == vert2 ? 0 : 1;

    if (edges2ToFakeVertsIn.get_num_comp(edge2, 0) == 0)
    {
      edges2ToFakeVertsIn.insert(edge2, 0, {});
      edges2ToFakeVertsIn.insert(edge2, 0, {});
    }

    int idx = edge2->get_down(0) == vert2 ? 0 : 1;
    assert(edge2->get_down(idx) == vert2);
    edges2ToFakeVertsIn(edge2, 0, idx) = VertOnEdge{fv, xi};
  }
}

int MeshProjectionCalculator::get_closest_projection(
    mesh::VariableSizeFieldPtr<predicates::impl::PointRecord> verts2AllClassificationsPtr, mesh::MeshEntityPtr vert2,
    utils::Point& closestPt)
{
  auto& verts2AllClassifications = *verts2AllClassificationsPtr;
  int nprojections               = verts2AllClassifications.get_num_comp(vert2, 0);
  assert(nprojections >= 1);

  if (nprojections == 1)
  {
    closestPt = m_pointClassifier->compute_xyz_coords(verts2AllClassifications(vert2, 0, 0));
    return 0;
  }

  // choose closest projection in xyz space
  utils::Point srcPt = vert2->get_point_orig(0);
  double minDist     = std::numeric_limits<double>::max();
  int minIdx         = -1;
  for (int i = 0; i < nprojections; ++i)
  {
    utils::Point destPt = m_pointClassifier->compute_xyz_coords(verts2AllClassifications(vert2, 0, i));
    utils::Point disp   = destPt - srcPt;
    double distSquared  = dot(disp, disp);
    if (distSquared < minDist)
    {
      minDist   = distSquared;
      minIdx    = i;
      closestPt = destPt;
    }
  }

  assert(minIdx != -1);
  return minIdx;
}

void MeshProjectionCalculator::intersect_mesh2_edges_with_mesh1()
{
  // FieldPtr<Bool> seen_verts = create_field<Bool>(m_mesh1, FieldShape(1, 0, 0), 1, false);
  auto& verts2ClassOnMesh1 = *(m_relationalData->verts2ClassOnMesh1);
  auto& verts2ToFakeVerts  = *(m_relationalData->verts2ToFakeVerts); // TODO: DEBUGGING
  std::vector<EdgeIntersection> intersections;
  for (auto& edge : m_mesh2->get_edges())
    if (edge)
    {
      predicates::impl::PointRecord& recordStart = verts2ClassOnMesh1(edge->get_down(0), 0, 0);
      predicates::impl::PointRecord& recordEnd   = verts2ClassOnMesh1(edge->get_down(1), 0, 0);

      if (m_output)
      {
        std::cout << "\ntracing edge " << edge << std::endl;
        std::cout << "edge vertex ids: " << edge->get_down(0)->get_id() << ", " << edge->get_down(1)->get_id()
                  << std::endl;
        std::cout << "edge = " << edge << std::endl;
        std::cout << "record_start = " << recordStart << ", el = " << recordStart.el << std::endl;
        std::cout << "record_end = " << recordEnd << ", el = " << recordEnd.el << std::endl;
      }

      m_edgeTracer.trace_edge(edge, recordStart, recordEnd, intersections);

      if (m_output && (m_vertIds.count(verts2ToFakeVerts(edge->get_down(0), 0, 0).id) > 0 ||
                       m_vertIds.count(verts2ToFakeVerts(edge->get_down(1), 0, 0).id) > 0))
      {
        std::cout << "edge " << edge << " has vertex of interest" << std::endl;
        std::cout << "edge vert ids = " << (verts2ToFakeVerts(edge->get_down(0), 0, 0).id) << ", "
                  << (verts2ToFakeVerts(edge->get_down(1), 0, 0).id) << std::endl;
        std::cout << "specified edge has " << intersections.size() << " intersections" << std::endl;
        std::cout << "first vert record = " << recordStart << ", el id = " << recordStart.el->get_id() << std::endl;
        std::cout << "second vert record = " << recordEnd << ", el id = " << recordEnd.el->get_id() << std::endl;
        std::cout << "second vert entity id = " << get_entity(recordEnd)->get_id() << std::endl;
        for (size_t i = 0; i < intersections.size(); ++i)
        {
          std::cout << "intersection " << i << " with mesh1 entity " << get_entity(intersections[i].record)
                    << std::endl;
          std::cout << "  alpha = " << intersections[i].alpha << ", " << intersections[i].beta << std::endl;
        }
      }

      for (auto& intersection : intersections)
      {
        assert(intersection.record.type == PointClassification::Vert ||
               intersection.record.type == PointClassification::Edge);

        if (intersection.record.type == PointClassification::Vert)
          record_mesh2_edge_on_mesh1_vert(edge, intersection);
        else
        {
          assert(intersection.record.type == PointClassification::Edge);
          record_mesh2_edge_intersects_mesh1_edge(edge, intersection);
        }
      }

      sort_edge2_fake_verts(edge);
    }
}

void MeshProjectionCalculator::record_mesh2_edge_on_mesh1_vert(mesh::MeshEntityPtr edge2,
                                                               const EdgeIntersection& intersection)
{
  assert(edge2->get_type() == stk::middle_mesh::mesh::MeshEntityType::Edge);
  assert(intersection.record.type == PointClassification::Vert);

  auto& edges2ToFakeVertsIn = *(m_relationalData->edges2ToFakeVertsIn);
  auto& verts1ToFakeVerts   = *(m_relationalData->verts1ToFakeVerts);

  mesh::MeshEntityPtr vert1 = get_entity(intersection.record);
  assert(vert1->get_type() == stk::middle_mesh::mesh::MeshEntityType::Vertex);

  auto& fakeVert1 = verts1ToFakeVerts(vert1, 0, 0);
  edges2ToFakeVertsIn.insert(edge2, 0, {fakeVert1, intersection.alpha});
}

void MeshProjectionCalculator::record_mesh2_edge_intersects_mesh1_edge(mesh::MeshEntityPtr edge2,
                                                                       const EdgeIntersection& intersection)
{
  // what happens if two successive intersections are on mesh1 vertices that share an edge?
  assert(intersection.record.type == PointClassification::Edge);
  auto& mesh1EdgesToSplit   = *(m_relationalData->mesh1EdgesToSplit);
  auto& edges2ToFakeVertsIn = *(m_relationalData->edges2ToFakeVertsIn);
  // auto& verts2_to_fake_verts = *(m_relationalData->verts2_to_fake_verts); // TODO: DEBUGGING

  // auto& verts2_class_on_mesh1 = *(m_relationalData->verts2_class_on_mesh1); //TODO: DEBUGGING

  mesh::MeshEntityPtr edge1 = get_entity(intersection.record);
  double xi                 = intersection.beta;
  auto pt                   = compute_edge_coords_orig(edge1, xi);

  FakeVert fv = m_fakeVertGen.get_vert(pt);

  // if (m_output && (m_vertIds.count(verts2_to_fake_verts(edge2->get_down(0), 0, 0).id) > 0 ||
  //     m_vertIds.count(verts2_to_fake_verts(edge2->get_down(1), 0, 0).id) > 0))
  if (m_output && m_vertIds.count(fv.id) > 0)
  {
    std::cout << "edge between specified vertices intersects a mesh1 edge" << std::endl;
    std::cout << "mesh1 edge = " << get_entity(intersection.record) << std::endl;
    std::cout << "mesh2 edge = " << edge2 << std::endl;
    std::cout << "alpha = " << intersection.alpha << ", beta = " << intersection.beta << std::endl;
    std::cout << "created fakevert id = " << fv.id << std::endl;
  }
  /*
    if (m_vertIds.count(fv.id) > 0)
    {
      //std::cout << "\nfakevert from edge-edge intersection" << std::endl;
      //std::cout << "created fakevert id " << fv.id << std::endl;
      //std::cout << "pointrrecord = " << intersection.record
      //          << ", alpha = " << intersection.alpha << ", beta = " << intersection.beta << std::endl;
      //std::cout << "edge1 = " << edge1 << std::endl;
      //std::cout << "edge2 = " << edge2 << std::endl;

      //mesh::MeshEntityPtr v1 = edge2->get_down(0);
      //mesh::MeshEntityPtr v2 = edge2->get_down(1);

      //std::cout << "edge vertex ids = " << v1->get_id() << ", " << v2->get_id() << std::endl;
      //std::cout << "edge2 first vert record = " << verts2_class_on_mesh1(v1, 0, 0) << ", el = " <<
    verts2_class_on_mesh1(v1, 0, 0).el << std::endl;
      //std::cout << "edge2 second vert record = " << verts2_class_on_mesh1(v2, 0, 0) << ", el = " <<
    verts2_class_on_mesh1(v2, 0, 0).el << std::endl;

      //std::cout << "edge2 element ids: ";
      //for (int i=0; i < edge2->count_up(); ++i)
      //  std::cout << edge2->get_up(i)->get_id() << ", ";
      //std::cout << std::endl;
    }
  */
  mesh1EdgesToSplit.insert(edge1, 0, EdgeSplitRecord{fv, xi});
  edges2ToFakeVertsIn.insert(edge2, 0, {fv, intersection.alpha});
}

void MeshProjectionCalculator::sort_edge2_fake_verts(mesh::MeshEntityPtr edge2)
{
  // the other functions in this class store the fake verts on each edge
  // in a particular order: (first vert, last vert, verts from edge intersections in sorted order)
  // So all we have to do here is move the last vert to the end and shift the data from edge
  // intersections left.
  auto& edges2ToFakeVertsIn = *(m_relationalData->edges2ToFakeVertsIn);

  assert(edges2ToFakeVertsIn.get_num_comp(edge2, 0) >= 2);
  int ncomp = edges2ToFakeVertsIn.get_num_comp(edge2, 0);
  if (ncomp > 2)
  {
    auto lastVert = edges2ToFakeVertsIn(edge2, 0, 1);
    for (int i = 2; i < ncomp; ++i)
      edges2ToFakeVertsIn(edge2, 0, i - 1) = edges2ToFakeVertsIn(edge2, 0, i);
    edges2ToFakeVertsIn(edge2, 0, ncomp - 1) = lastVert;
  }
}

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
