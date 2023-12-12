#include "stk_middle_mesh/predicates/point_classifier_normal_wrapper.hpp"
#include "stk_middle_mesh/field.hpp"
#include "stk_middle_mesh/plane_projection.hpp"

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

PointRecord PointClassifierNormalWrapper::classify(mesh::MeshEntityPtr destFace, mesh::MeshEntityPtr srcFace,
                                                   const utils::Point& pt)
{
  assert(is_entity_on_mesh(srcFace));

  utils::Point normal = compute_normal_vector(srcFace, pt);
  if (destFace->get_type() == mesh::MeshEntityType::Triangle)
  {
    // compute xi coordinates of pt in src_face
    // compute normal vector at that point
    auto r = m_normalClassifier.classify(destFace, pt, normal);
    return PointRecord(r.type, r.id, r.el, r);
  } else // quad
  {
    m_quadToTriangles.set_triangles(destFace);

    auto r1 = m_normalClassifier.classify(m_quadToTriangles.el1, pt, normal);
    auto r2 = m_normalClassifier.classify(m_quadToTriangles.el2, pt, normal);

    auto quadRecord = m_quadToTriangles.get_quad_record(destFace, r1, r2);

    return quadRecord;
  }
}

PointRecord PointClassifierNormalWrapper::classify_reverse(mesh::MeshEntityPtr destFace, const utils::Point& pt,
                                                           bool logicalResultOnly)
{
  assert(is_entity_on_mesh(destFace));

  // utils::Point normal = computeNormalVector(src_face, pt);
  auto normals = get_normal_vectors(destFace);
  if (destFace->get_type() == mesh::MeshEntityType::Triangle)
  {
    // compute xi coordinates of pt in src_face
    // compute normal vector at that point
    auto r = m_normalClassifier.classify_reverse(destFace, pt, normals, logicalResultOnly);
    return PointRecord(r.type, r.id, r.el, r);
  } else // quad
  {
    m_quadToTriangles.set_triangles(destFace);
    std::array<utils::Point, 4> triNormals = {normals[m_quadToTriangles.VERTMAP_TRI1_TO_QUAD[0]],
                                              normals[m_quadToTriangles.VERTMAP_TRI1_TO_QUAD[1]],
                                              normals[m_quadToTriangles.VERTMAP_TRI1_TO_QUAD[2]]};
    auto r1 = m_normalClassifier.classify_reverse(m_quadToTriangles.el1, pt, triNormals, logicalResultOnly);

    for (int i = 0; i < 3; ++i)
      triNormals[i] = normals[m_quadToTriangles.VERTMAP_TRI2_TO_QUAD[i]];
    auto r2 = m_normalClassifier.classify_reverse(m_quadToTriangles.el2, pt, triNormals, logicalResultOnly);

    auto quadRecord = m_quadToTriangles.get_quad_record(destFace, r1, r2);

    return quadRecord;
  }
}

utils::Point PointClassifierNormalWrapper::compute_xyz_coords(const PointRecord& record, bool allowExterior)
{
  if (!allowExterior)
  {
    assert(record.type != PointClassification::Exterior);
  }

  if (record.el->get_type() == mesh::MeshEntityType::Triangle)
  {
    return m_triangleCoordUtils.compute_xyz_coords(record.m_r1, allowExterior);
  } else if (record.el->get_type() == mesh::MeshEntityType::Quad)
  {
    m_quadToTriangles.set_triangles(record.el);
    return m_quadToTriangles.compute_xyz_coords(record, allowExterior);
    //const PointRecordForTriangle& r = record.m_r1.type != PointClassification::Exterior ? record.m_r1 : record.m_r2;
    //return m_triangleCoordUtils.compute_xyz_coords(r);
  } else
    throw std::runtime_error("element must be triangle or quad");
}

utils::Point PointClassifierNormalWrapper::compute_xi_coords(const PointRecord& record, bool allowExterior)
{
  if (!allowExterior)
  {
    assert(record.type != PointClassification::Exterior);
  }

  if (record.el->get_type() == mesh::MeshEntityType::Triangle)
  {
    return m_triangleCoordUtils.compute_xi_coords(record.m_r1, allowExterior);
  } else if (record.el->get_type() == mesh::MeshEntityType::Quad)
  {
    m_quadToTriangles.set_triangles(record.el);
    return m_quadToTriangles.get_quad_xi_coords(record, allowExterior);
  } else
    throw std::runtime_error("element must be triangle or quad");
}


double PointClassifierNormalWrapper::get_edge_xi(const PointRecord& record)
{
  assert(record.type == PointClassification::Edge);

  if (record.el->get_type() == mesh::MeshEntityType::Triangle)
    return m_triangleCoordUtils.get_edge_xi(record.m_r1);
  else // quad
  {
    m_quadToTriangles.set_triangles(record.el);
    const PointRecordForTriangle& r = record.m_r1.type == PointClassification::Edge ? record.m_r1 : record.m_r2;
    assert(r.type == PointClassification::Edge);

    double xi = m_triangleCoordUtils.get_edge_xi(r);
    if (record.el->get_down_orientation(record.id) == mesh::EntityOrientation::Reversed)
      xi = 1 - xi;

    return xi;
  }
}

double PointClassifierNormalWrapper::compute_orthogonal_dist(const PointRecord& record1, const int id)
{
  if (record1.el->get_type() == mesh::MeshEntityType::Triangle)
    return m_triangleCoordUtils.compute_orthogonal_dist(record1.m_r1, id);
  else // quad
  {
    int triId                 = m_quadToTriangles.EDGEMAP_QUAD_TO_TRI[id].first;
    mesh::MeshEntityPtr triEl = triId == 0 ? m_quadToTriangles.el1 : m_quadToTriangles.el2;
    int triEdgeId             = m_quadToTriangles.EDGEMAP_QUAD_TO_TRI[id].second;

    const PointRecordForTriangle& r = record1.m_r1.el == triEl ? record1.m_r1 : record1.m_r2;
    return m_triangleCoordUtils.compute_orthogonal_dist(r, triEdgeId);
  }
}

PointRecord PointClassifierNormalWrapper::create_vert_record(mesh::MeshEntityPtr el, int vertId)
{
  assert(el->get_type() == mesh::MeshEntityType::Triangle || el->get_type() == mesh::MeshEntityType::Quad);

  if (el->get_type() == mesh::MeshEntityType::Triangle) {
    PointRecordForTriangle record = m_triangleCoordUtils.create_record(el, vertId);
    return PointRecord(PointClassification::Vert, vertId, el, record);
  } else // quad
    return m_quadToTriangles.create_record(el, vertId);
}

PointRecord PointClassifierNormalWrapper::create_edge_record(mesh::MeshEntityPtr el, int edgeId, double edgeXiOnReferenceEl)
{
  assert(el->get_type() == mesh::MeshEntityType::Triangle || el->get_type() == mesh::MeshEntityType::Quad);

  if (el->get_type() == mesh::MeshEntityType::Triangle) {
    PointRecordForTriangle record = m_triangleCoordUtils.create_record(el, edgeId, edgeXiOnReferenceEl);
    return PointRecord(PointClassification::Edge, edgeId, el, record);
  } else // quad
    return m_quadToTriangles.create_record(el, edgeId, edgeXiOnReferenceEl);
}

PointRecord PointClassifierNormalWrapper::classify_onto(const PointRecord& record, mesh::MeshEntityPtr el)
{
  assert(el->get_type() == mesh::MeshEntityType::Triangle || el->get_type() == mesh::MeshEntityType::Quad);

  if (el->get_type() == mesh::MeshEntityType::Triangle) {
    PointRecordForTriangle recordTri = m_triangleCoordUtils.classify_onto(record.m_r1, el);
    return PointRecord(record.type, record.id, el, recordTri);
  } else // quad
    return m_quadToTriangles.classify_onto(record, el) ;
}



std::array<utils::Point, 4> PointClassifierNormalWrapper::get_normal_vectors(mesh::MeshEntityPtr face)
{
  std::array<mesh::MeshEntityPtr, 4> verts;
  int nverts = get_downward(face, 0, verts.data());

  auto& field = *m_normalField;
  std::array<utils::Point, 4> normals;
  for (int i = 0; i < nverts; ++i)
    normals[i] = field(verts[i], 0, 0);

  return normals;
}

utils::Point PointClassifierNormalWrapper::compute_normal_vector(mesh::MeshEntityPtr face, const utils::Point& pt)
{
  if (face->get_type() == mesh::MeshEntityType::Triangle)
    return compute_normal_vector_triangle(face, pt);
  else if (face->get_type() == mesh::MeshEntityType::Quad)
    return compute_normal_vector_quad(face, pt);
  else
    throw std::runtime_error("unrecognized type");
}

utils::Point PointClassifierNormalWrapper::compute_normal_vector_quad(mesh::MeshEntityPtr face, const utils::Point& pt)
{
  assert(face->get_type() == mesh::MeshEntityType::Quad);

  std::array<mesh::MeshEntityPtr, 4> verts;
  get_downward(face, 0, verts.data());

  std::array<utils::Point, 3> tri1Pts = {verts[m_quadToTriangles.VERTMAP_TRI1_TO_QUAD[0]]->get_point_orig(0),
                                         verts[m_quadToTriangles.VERTMAP_TRI1_TO_QUAD[1]]->get_point_orig(0),
                                         verts[m_quadToTriangles.VERTMAP_TRI1_TO_QUAD[2]]->get_point_orig(0)};

  std::array<utils::Point, 3> tri2Pts = {verts[m_quadToTriangles.VERTMAP_TRI2_TO_QUAD[0]]->get_point_orig(0),
                                         verts[m_quadToTriangles.VERTMAP_TRI2_TO_QUAD[1]]->get_point_orig(0),
                                         verts[m_quadToTriangles.VERTMAP_TRI2_TO_QUAD[2]]->get_point_orig(0)};
  auto tri1Xi                         = compute_triangle_xi_coords(tri1Pts, pt);
  auto tri2Xi                         = compute_triangle_xi_coords(tri2Pts, pt);

  int triId = -1;
  if (in_range(tri1Xi.x, 0, 1, m_tolerances.quadBreakupTolerance) &&
      in_range(tri1Xi.y, 0, 1, m_tolerances.quadBreakupTolerance) &&
      in_range(1 - tri1Xi.x - tri1Xi.y, 0, 1, m_tolerances.quadBreakupTolerance))
  {
    triId = 1;
  } else if (in_range(tri2Xi.x, 0, 1, m_tolerances.quadBreakupTolerance) &&
             in_range(tri2Xi.y, 0, 1, m_tolerances.quadBreakupTolerance) &&
             in_range(1 - tri2Xi.x - tri2Xi.y, 0, 1, m_tolerances.quadBreakupTolerance))
  {
    triId = 2;
  } else
    throw std::runtime_error("could not find point in triangle");
  // TODO: figure out which element the point is "most" inside of

  // std::array<utils::Point, 3>& pts = tri_id == 1 ? tri1_pts : tri2_pts;
  utils::Point& ptXi = triId == 1 ? tri1Xi : tri2Xi;
  auto& vertmap      = triId == 1 ? m_quadToTriangles.VERTMAP_TRI1_TO_QUAD : m_quadToTriangles.VERTMAP_TRI2_TO_QUAD;

  auto& field = *m_normalField;
  utils::Point normal(0, 0, 0);
  std::array<double, 3> xi = {1 - ptXi.x - ptXi.y, ptXi.x, ptXi.y};
  for (int i = 0; i < 3; ++i)
    normal += xi[i] * field(verts[vertmap[i]], 0, 0);

  return normal;
}

utils::Point PointClassifierNormalWrapper::compute_normal_vector_triangle(mesh::MeshEntityPtr face,
                                                                          const utils::Point& pt)
{
  assert(face->get_type() == mesh::MeshEntityType::Triangle);

  std::array<mesh::MeshEntityPtr, 3> verts;
  [[maybe_unused]] int nverts = get_downward(face, 0, verts.data());
  assert(nverts == 3);

  auto ptXi   = compute_triangle_xi_coords(verts, pt);
  auto& field = *m_normalField;
  utils::Point normal(0, 0, 0);
  std::array<double, 3> xi = {1 - ptXi.x - ptXi.y, ptXi.x, ptXi.y};
  for (int i = 0; i < 3; ++i)
    normal += xi[i] * field(verts[i], 0, 0);

  return normal;
}

utils::Point PointClassifierNormalWrapper::compute_triangle_xi_coords(std::array<mesh::MeshEntityPtr, 3>& verts,
                                                                      const utils::Point& pt)
{
  std::array<utils::Point, 3> vertPts;
  for (int i = 0; i < 3; ++i)
  {
    assert(verts[i]->get_type() == mesh::MeshEntityType::Vertex);
    vertPts[i] = verts[i]->get_point_orig(0);
  }

  return compute_triangle_xi_coords(vertPts, pt);
}

// takes the verts of a triangle (in 3d) and a point that lies on the triangle,
// and returns the xi coordinates of the point
utils::Point PointClassifierNormalWrapper::compute_triangle_xi_coords(std::array<utils::Point, 3>& verts,
                                                                      const utils::Point& pt)
{
  // figure out the best 2 coordinates to use
  utils::Point delta1 = verts[1] - verts[0];
  utils::Point delta2 = verts[2] - verts[0];

  double detXy       = delta1.x * delta2.y - delta1.y * delta2.x;
  double detYz       = delta1.y * delta2.z - delta1.z * delta2.y;
  double detZx       = delta1.z * delta2.x - delta1.x * delta2.z;
  bool shouldReverse = false;

  int coord1 = -1, coord2 = -1;
  if (std::abs(detXy) > std::abs(detYz) && std::abs(detXy) > std::abs(detZx))
  {
    coord1        = 0;
    coord2        = 1;
    shouldReverse = detXy < 0;
  } else if (std::abs(detYz) > std::abs(detXy) && std::abs(detYz) > std::abs(detZx))
  {
    coord1        = 1;
    coord2        = 2;
    shouldReverse = detYz < 0;
  } else // zx is best
  {
    coord1        = 2;
    coord2        = 0;
    shouldReverse = detZx < 0;
  }

  if (shouldReverse)
  {
    int tmpCoord = coord1;
    coord1       = coord2;
    coord2       = tmpCoord;
  }

  std::array<utils::Point, 3> vertsProjected;
  for (int i = 0; i < 3; ++i)
  {
    vertsProjected[i] = utils::Point(verts[i][coord1], verts[i][coord2]);
  }
  utils::Point ptProjected(pt[coord1], pt[coord2]);

  mesh::impl::ElementOperations2D elemOps;
  return elemOps.compute_tri_xi_coords(vertsProjected, ptProjected);
}


bool PointClassifierNormalWrapper::is_entity_on_mesh(mesh::MeshEntityPtr entity)
{
  int id          = entity->get_id();
  auto& entityVec = m_mesh->get_mesh_entities(get_type_dimension(entity->get_type()));

  return id < int(entityVec.size()) && entityVec[id] == entity;
}

} // namespace impl

} // namespace predicates
} // namespace middle_mesh
} // namespace stk
