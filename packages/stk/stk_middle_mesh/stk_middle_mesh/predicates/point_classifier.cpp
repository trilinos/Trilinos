#include "stk_middle_mesh/predicates/point_classifier.hpp"
#include <limits>

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

PointClassifier::PointClassifier(const double eps, const double smallSlopeTol)
  : m_classifier(eps, smallSlopeTol)
  , m_defaultProjection(utils::Point(0, 0, 0), utils::Point(1, 0, 0), utils::Point(0, 1, 0))
{
  set_projection(&m_defaultProjection);
}

void PointClassifier::set_projection(utils::impl::Projection* proj)
{
  m_elemOps.set_projection(proj);
  m_classifier.set_projection(proj);
}

PointRecord PointClassifier::classify(mesh::MeshEntityPtr face, const utils::Point& pt)
{
  if (face->get_type() == mesh::MeshEntityType::Triangle)
  {
    auto r = m_classifier.classify(face, pt);
    return PointRecord(r.type, r.id, r.el, r);
  } else // quad
  {
    m_quadToTriangles.set_triangles(face);
    auto r1 = m_classifier.classify(m_quadToTriangles.el1, pt);
    auto r2 = m_classifier.classify(m_quadToTriangles.el2, pt);

    auto quadRecord = m_quadToTriangles.get_quad_record(face, r1, r2);

    return quadRecord;
  }
}

double PointClassifier::get_edge_xi(const PointRecord& record)
{
  if (record.el->get_type() == mesh::MeshEntityType::Triangle)
    return m_classifier.get_edge_xi(record.m_r1);
  else // quad
  {
    m_quadToTriangles.set_triangles(record.el);
    double xi = m_classifier.get_edge_xi(record.m_r1);
    return convert_edge_xi(xi, record.el, record.id);
  }
}

// TODO: replace std::pair with struct
PossibleEdgeIntersection PointClassifier::get_edge_xi(const PointRecord& record1, const PointRecord& record2)
{
  if (record1.el->get_type() == mesh::MeshEntityType::Triangle)
  {
    auto p = m_classifier.get_edge_xi(record1.m_r1, record2.m_r1);
    return {convert_point_record(p.record), p.edgeXi};
  } else // quad
  {
    m_quadToTriangles.set_triangles(record1.el);
    const PointRecord& rI = record1.type == PointClassification::Interior ? record1 : record2;
    const PointRecord& rE = record1.type == PointClassification::Interior ? record2 : record1;

    // get PointRecordForTriangle for record2 corresponding to record1.m_r1 and m_r2
    const PointRecordForTriangle& rSame  = rI.m_r1.el == rE.m_r1.el ? rE.m_r1 : rE.m_r2;
    const PointRecordForTriangle& rOther = rI.m_r2.el == rE.m_r2.el ? rE.m_r2 : rE.m_r1;

    // check if the point is on the interior edge
    if (rI.m_r1.type == PointClassification::Edge)
    {
      PossibleEdgeIntersectionForTriangle p = m_classifier.get_edge_intersection_xi(rI.m_r1, rSame);

      if (p.record.id == -1)
        p = m_classifier.get_edge_intersection_xi(rI.m_r2, rOther);

      // try again without excluding verts on the end of the edge
      if (p.record.id == -1)
      {
        p = m_classifier.get_edge_intersection_xi(rI.m_r1, rSame, false);

        if (p.record.id == -1)
          p = m_classifier.get_edge_intersection_xi(rI.m_r2, rOther, false);
      }

      if (p.record.id != -1)
      {
        return get_exterior_quad_point_record(record1.el, p.record, p.edgeXi);
      } else
      {
        throw std::runtime_error("unable to find intersection");
      }
    }

    // else point is on interior of tri1_id
    auto p = m_classifier.get_edge_xi(rI.m_r1, rSame);
    if (!(p.record.type == PointClassification::Edge && p.record.id == m_quadToTriangles.INTERIOR_EDGE))
    {
      PointRecordForTriangle r2Unused{PointClassification::Exterior};
      PointRecord recordQuad = m_quadToTriangles.get_quad_record(record1.el, p.record, r2Unused);
      double quadEdgeXi      = convert_edge_xi(p.edgeXi, record1.el, recordQuad.id);
      return {recordQuad, quadEdgeXi};
    }

    // else the intersection is across the internal edge, need to
    // consider corner intersection with the other triangle
    CornerRecordForTriangle corner = m_classifier.get_edge_xi_corner(rI.m_r2, rOther);
    return get_exterior_quad_point_record(record1.el, corner);
  }
}

CornerRecord PointClassifier::get_edge_xi_corner(const PointRecord& record1, const PointRecord& record2)
{
  if (record1.el->get_type() == mesh::MeshEntityType::Triangle)
  {
    CornerRecordForTriangle c = m_classifier.get_edge_xi_corner(record1.m_r1, record2.m_r1);
    return CornerRecord(c.hasIntersection, convert_point_record(c.record1), convert_point_record(c.record2), c.xi1,
                        c.xi2);
  } else // quad
  {
    m_quadToTriangles.set_triangles(record1.el);

    // get PointRecordForTriangle for record2 corresponding to record1.m_r1 and m_r2
    const PointRecordForTriangle& rSame  = record1.m_r1.el == record2.m_r1.el ? record2.m_r1 : record2.m_r2;
    const PointRecordForTriangle& rOther = record1.m_r2.el == record2.m_r2.el ? record2.m_r2 : record2.m_r1;
    CornerRecordForTriangle c1           = m_classifier.get_edge_xi_corner(record1.m_r1, rSame);
    CornerRecordForTriangle c2           = m_classifier.get_edge_xi_corner(record1.m_r2, rOther);

    // there can be an inconsistent case due to floating point math:
    // the edge can intersect the interior edge of one triangle but not the
    // interior edge of the other triangle (even though they are the same edge).
    // This can only happen if the edge
    // is very close to being outside the element.  In this case, pretend there
    // is no intersection
    // TODO: consider using bitwise operations to avoid short circuiting causing additional
    //      branches
    bool doesC1IntersectMiddleEdge =
        (c1.record1.type == PointClassification::Edge && c1.record1.id == m_quadToTriangles.INTERIOR_EDGE) ||
        (c1.record2.type == PointClassification::Edge && c1.record2.id == m_quadToTriangles.INTERIOR_EDGE);
    bool doesC2IntersectMiddleEdge =
        (c2.record1.type == PointClassification::Edge && c2.record1.id == m_quadToTriangles.INTERIOR_EDGE) ||
        (c2.record2.type == PointClassification::Edge && c2.record2.id == m_quadToTriangles.INTERIOR_EDGE);
    bool isInconsistentCase = (c1.hasIntersection && doesC1IntersectMiddleEdge && !c2.hasIntersection) ||
                              (c2.hasIntersection && doesC2IntersectMiddleEdge && !c1.hasIntersection);
    if ((!c1.hasIntersection && !c2.hasIntersection) || isInconsistentCase)
      return CornerRecord(false);

    return convert_tri_corners_to_quad_corner(record1.el, c1, c2);
  }
}

PossibleEdgeIntersection PointClassifier::get_edge_intersection_xi(const PointRecord& record1,
                                                                   const PointRecord& record2)
{
  assert(record1.type == PointClassification::Exterior || record2.type == PointClassification::Exterior);
  assert(record1.type != PointClassification::Interior);
  assert(record2.type != PointClassification::Interior);

  if (record1.el->get_type() == mesh::MeshEntityType::Triangle)
  {
    auto p = m_classifier.get_edge_intersection_xi(record1.m_r1, record2.m_r1);
    return {PointRecord(p.record.type, p.record.id, p.record.el, p.record), p.edgeXi};
  } else // quad
  {
    m_quadToTriangles.set_triangles(record1.el);
    const PointRecord& recordI = record1.type == PointClassification::Exterior ? record2 : record1;
    const PointRecord& recordE = record1.type == PointClassification::Exterior ? record1 : record2;

    assert(recordI.type != PointClassification::Exterior);
    bool isPointOnSharedVert = recordI.m_r1.type == PointClassification::Vert && recordI.m_r1.id != 0;

    if (!isPointOnSharedVert)
    {
      // the only way the point could be classified on an edge in both
      // triangles is if it is on the third edge.
      assert(!(recordI.m_r1.type == PointClassification::Edge && recordI.m_r2.type == PointClassification::Edge));

      // m_r1_i is the one classified on the edge
      bool isR1OnVertOrEdge =
          recordI.m_r1.type == PointClassification::Vert || recordI.m_r1.type == PointClassification::Edge;
      const PointRecordForTriangle& r1I = isR1OnVertOrEdge ? recordI.m_r1 : recordI.m_r2;
      const PointRecordForTriangle& r2I = isR1OnVertOrEdge ? recordI.m_r2 : recordI.m_r1;

      const PointRecordForTriangle& r1E = r1I.el == recordE.m_r1.el ? recordE.m_r1 : recordE.m_r2;
      const PointRecordForTriangle& r2E = r1I.el == recordE.m_r1.el ? recordE.m_r2 : recordE.m_r1;

      auto p = m_classifier.get_edge_intersection_xi(r1I, r1E, false);
      // if intersection found with non-interior edge
      if (p.record.id != -1 &&
          !(p.record.type == PointClassification::Edge && p.record.id == m_quadToTriangles.INTERIOR_EDGE))
      {
        return make_triangle_records_for_intersection_that_must_exist(recordI.el, recordI.m_r1.el, p);
      } else
      {
        // the line is between a point on one edge of the triangle and either
        //   1. passes through the interior edge and an exterior edge of the other
        //      triangle
        // or
        //   2. does not intersect the quad at all
        CornerRecordForTriangle c = m_classifier.get_edge_xi_corner(r2I, r2E);
        return get_exterior_quad_point_record(record1.el, c);
      }
    } else
    {
      // std::cout << "interior point is on shared vertex" << std::endl;
      //  interior point is on a shared vertex, need to call
      //  get_edge_intersection_xi() on both elements

      // check for intersections on record_i.r1.el
      const PointRecordForTriangle& rSame = recordI.m_r1.el == recordE.m_r1.el ? recordE.m_r1 : recordE.m_r2;
      auto p                              = m_classifier.get_edge_intersection_xi(recordI.m_r1, rSame, false);
      if (p.record.id != -1 &&
          !(p.record.type == PointClassification::Edge && p.record.id == m_quadToTriangles.INTERIOR_EDGE))
        return make_triangle_records_for_intersection_that_must_exist(recordI.el, recordI.m_r1.el, p);

      // else check for intersections on record_i.m_r2.el
      const PointRecordForTriangle& rOther = recordI.m_r2.el == recordE.m_r1.el ? recordE.m_r1 : recordE.m_r2;
      p                                    = m_classifier.get_edge_intersection_xi(recordI.m_r2, rOther, false);
      if (p.record.id != -1 &&
          !(p.record.type == PointClassification::Edge && p.record.id == m_quadToTriangles.INTERIOR_EDGE))
        return make_triangle_records_for_intersection_that_must_exist(recordI.el, recordI.m_r2.el, p);
      else
        return {PointRecord(PointClassification::Exterior, -1, record1.el), -1};
    }
  }
}

utils::Point PointClassifier::compute_xyz_coords(const PointRecord& record)
{
  assert(record.type != PointClassification::Exterior);
  assert(record.el->get_type() == mesh::MeshEntityType::Triangle ||
         record.el->get_type() == mesh::MeshEntityType::Quad);

  if (record.el->get_type() == mesh::MeshEntityType::Triangle)
  {
    return m_classifier.compute_xyz_coords(record.m_r1);
  } else // quad
  {
    m_quadToTriangles.set_triangles(record.el);
    if (record.type == PointClassification::Vert)
    {
      std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
      get_downward(record.el, 0, verts.data());
      return verts[record.id]->get_point_orig(0);
    } else
    {
      m_quadToTriangles.set_triangles(record.el);
      const PointRecordForTriangle recordTri =
          record.m_r1.type != PointClassification::Exterior ? record.m_r1 : record.m_r2;
      assert(recordTri.type != PointClassification::Exterior);
      return m_classifier.compute_xyz_coords(recordTri);
    }
  }
}

double PointClassifier::compute_orthogonal_dist(const PointRecord& record1, const int id)
{
  if (record1.el->get_type() == mesh::MeshEntityType::Triangle)
    return m_classifier.compute_orthogonal_dist(record1.m_r1, id);
  else // quad
  {
    int triId               = m_quadToTriangles.EDGEMAP_QUAD_TO_TRI[id].first;
    mesh::MeshEntityPtr tri = triId == 0 ? m_quadToTriangles.el1 : m_quadToTriangles.el2;
    int triEdgeId           = m_quadToTriangles.EDGEMAP_QUAD_TO_TRI[id].second;
    // const PointRecordForTriangle& r = record1.tri_id1 == tri_id ? record1.m_r1 : record1.m_r2;
    const PointRecordForTriangle& r = record1.m_r1.el == tri ? record1.m_r1 : record1.m_r2;
    return m_classifier.compute_orthogonal_dist(r, triEdgeId);
  }
}

double PointClassifier::get_edge_xi_orthogonal(const PointRecord& record1, const int id)
{
  if (record1.el->get_type() == mesh::MeshEntityType::Triangle)
    return m_classifier.get_edge_xi_orthogonal(record1.m_r1, id);
  else // quad
  {
    int triId               = m_quadToTriangles.EDGEMAP_QUAD_TO_TRI[id].first;
    mesh::MeshEntityPtr tri = triId == 0 ? m_quadToTriangles.el1 : m_quadToTriangles.el2;
    int triEdgeId           = m_quadToTriangles.EDGEMAP_QUAD_TO_TRI[id].second;
    // const PointRecordForTriangle& r = record1.tri_id1 == tri_id ? record1.m_r1 : record1.m_r2;
    const PointRecordForTriangle& r = record1.m_r1.el == tri ? record1.m_r1 : record1.m_r2;
    double xi                       = m_classifier.get_edge_xi_orthogonal(r, triEdgeId);
    return convert_edge_xi(xi, record1.el, id);
  }
}

double PointClassifier::convert_edge_xi(double triXi, mesh::MeshEntityPtr el, const int id)
{
  // the way m_mesh is set up the triangle edges are oriented the same
  // way as the quad reference element, so all we need to do here is
  // take into account the edge orientation

  double xi = triXi;
  if (el->get_down_orientation(id) == mesh::EntityOrientation::Reversed)
    xi = 1 - xi;

  return xi;
}

PointRecord PointClassifier::convert_point_record(const PointRecordForTriangle& record)
{
  return PointRecord(record.type, record.id, record.el, record);
}

CornerRecord PointClassifier::convert_tri_corners_to_quad_corner(mesh::MeshEntityPtr quadEl,
                                                                 const CornerRecordForTriangle& c1,
                                                                 const CornerRecordForTriangle& c2)
{
  assert(c1.hasIntersection || c2.hasIntersection);

  PossibleEdgeIntersection p1, p2;

  bool c1IntersectionOnInteriorEdge =
      (c1.record1.type == PointClassification::Edge && c1.record1.id == m_quadToTriangles.INTERIOR_EDGE) ||
      (c1.record2.type == PointClassification::Edge && c1.record2.id == m_quadToTriangles.INTERIOR_EDGE);
  // bool c2_intersection_on_interior_edge = c2.record1.type == PointClassification::Edge && c2.record1.id ==
  // m_quadToTriangles.INTERIOR_EDGE ||
  //                                         c2.record2.type == PointClassification::Edge && c2.record2.id ==
  //                                         m_quadToTriangles.INTERIOR_EDGE;
  bool c1HasTwoNonInteriorEdgeIntersections =
      (c1.record1.type != PointClassification::Exterior && c1.record2.type != PointClassification::Exterior) &&
      !(c1.record1.type == PointClassification::Edge && c1.record1.id == m_quadToTriangles.INTERIOR_EDGE) &&
      !(c1.record2.type == PointClassification::Edge && c1.record2.id == m_quadToTriangles.INTERIOR_EDGE);
  bool c2HasTwoNonInteriorEdgeIntersections =
      (c2.record1.type != PointClassification::Exterior && c2.record2.type != PointClassification::Exterior) &&
      !(c2.record1.type == PointClassification::Edge && c2.record1.id == m_quadToTriangles.INTERIOR_EDGE) &&
      !(c2.record2.type == PointClassification::Edge && c2.record2.id == m_quadToTriangles.INTERIOR_EDGE);
  if (c1HasTwoNonInteriorEdgeIntersections)
  {
    p1 = get_exterior_quad_point_record(quadEl, c1.record1, c1.xi1);
    p2 = get_exterior_quad_point_record(quadEl, c1.record2, c1.xi2);
  } else if (c2HasTwoNonInteriorEdgeIntersections)
  {
    p1 = get_exterior_quad_point_record(quadEl, c2.record1, c2.xi1);
    p2 = get_exterior_quad_point_record(quadEl, c2.record2, c2.xi2);

  } else if (c1IntersectionOnInteriorEdge) // TODO: or intersection is on a shared vertex
  {
    p1 = get_exterior_quad_point_record(quadEl, c1);
    p2 = get_exterior_quad_point_record(quadEl, c2);
  } else if (c1.hasIntersection) // the intersection is on the non-shared vert
  {
    p1 = get_exterior_quad_point_record(quadEl, c1.record1, c1.xi1);
  } else if (c2.hasIntersection)
  {
    p1 = get_exterior_quad_point_record(quadEl, c2.record1, c2.xi1);
  } else
    throw std::runtime_error("unreachable case");

  return CornerRecord(true, p1.record, p2.record, p1.edgeXi, p2.edgeXi);
}

PossibleEdgeIntersection PointClassifier::get_exterior_quad_point_record(mesh::MeshEntityPtr quadEl,
                                                                         const CornerRecordForTriangle& corner)
{
  if (!corner.hasIntersection)
    return {PointRecord(PointClassification::Exterior, -1, quadEl), -1};

  const PointRecordForTriangle* exteriorRecord;
  double xi = -1;
  if (corner.record1.type == PointClassification::Vert)
    exteriorRecord = &(corner.record1);
  else if (corner.record2.type == PointClassification::Vert)
    exteriorRecord = &(corner.record2);
  else
  {
    assert(corner.record1.id == m_quadToTriangles.INTERIOR_EDGE ||
           corner.record2.id == m_quadToTriangles.INTERIOR_EDGE);
    exteriorRecord = corner.record1.id == m_quadToTriangles.INTERIOR_EDGE ? &(corner.record2) : &(corner.record1);
    // Note: it would be ok to move this outside the if condition.  This would allow using branchless
    //       techniques to chose the exterior record
    xi = corner.record1.id == m_quadToTriangles.INTERIOR_EDGE ? corner.xi2 : corner.xi1;
  }

  return get_exterior_quad_point_record(quadEl, *exteriorRecord, xi);
}

PossibleEdgeIntersection PointClassifier::get_exterior_quad_point_record(mesh::MeshEntityPtr quadEl,
                                                                         const PointRecordForTriangle& record,
                                                                         double xi)
{
  assert(record.type == PointClassification::Vert || record.type == PointClassification::Edge);
  assert(!(record.type == PointClassification::Edge && record.id == m_quadToTriangles.INTERIOR_EDGE));

  mesh::MeshEntityPtr triEl = record.el;
  const int* vertmap =
      triEl == m_quadToTriangles.el1 ? m_quadToTriangles.VERTMAP_TRI1_TO_QUAD : m_quadToTriangles.VERTMAP_TRI2_TO_QUAD;
  const int* edgemap =
      triEl == m_quadToTriangles.el1 ? m_quadToTriangles.EDGEMAP_TRI1_TO_QUAD : m_quadToTriangles.EDGEMAP_TRI2_TO_QUAD;
  int entityId = record.type == PointClassification::Vert ? vertmap[record.id] : edgemap[record.id];

  PointRecord recordOut(record.type, entityId, quadEl, record);
  m_quadToTriangles.enforce_record_consistency(recordOut.m_r1, recordOut.m_r2);

  if (record.type == PointClassification::Edge)
    xi = convert_edge_xi(xi, quadEl, entityId);
  else
    xi = -1;

  return {recordOut, xi};
}

std::pair<PointRecordForTriangle, PointRecordForTriangle>
PointClassifier::make_triangle_records_for_quad_vert_intersection(mesh::MeshEntityPtr triEl, const int vertId)
{
  const int* vertmap1 =
      triEl == m_quadToTriangles.el1 ? m_quadToTriangles.VERTMAP_TRI1_TO_QUAD : m_quadToTriangles.VERTMAP_TRI2_TO_QUAD;
  mesh::MeshEntityPtr otherEl = triEl == m_quadToTriangles.el1 ? m_quadToTriangles.el2 : m_quadToTriangles.el1;

  auto pt =
      m_elemOps.get_projection()->project_plane_coords(m_quadToTriangles.verts[vertmap1[vertId]]->get_point_orig(0));
  auto ptXi1                = m_elemOps.compute_tri_xi_coords(triEl, pt);
  auto ptXi2                = m_elemOps.compute_tri_xi_coords(otherEl, pt);
  PointClassification type1 = PointClassification::Vert;
  PointClassification type2 = vertId == 0 ? PointClassification::Exterior : PointClassification::Vert;

  PointRecordForTriangle r1(type1, vertId, triEl, ptXi1);
  PointRecordForTriangle r2(type2, 1 - vertId, otherEl, ptXi2);

  return std::make_pair(r1, r2);
}

std::pair<PointRecordForTriangle, PointRecordForTriangle>
PointClassifier::make_triangle_records_for_quad_edge_intersection(mesh::MeshEntityPtr triEl, const int edgeId,
                                                                  const double xi)
{
  mesh::MeshEntityPtr otherEl = triEl == m_quadToTriangles.el1 ? m_quadToTriangles.el2 : m_quadToTriangles.el1;

  auto pt                   = m_elemOps.compute_edge_coords(triEl->get_down(edgeId), xi);
  auto ptXi1                = m_elemOps.compute_tri_xi_coords(triEl, pt);
  auto ptXi2                = m_elemOps.compute_tri_xi_coords(otherEl, pt);
  PointClassification type1 = PointClassification::Edge;
  PointClassification type2 = PointClassification::Exterior;

  PointRecordForTriangle r1(type1, edgeId, triEl, ptXi1);
  PointRecordForTriangle r2(type2, 0, otherEl, ptXi2);

  return std::make_pair(r1, r2);
}

PossibleEdgeIntersection PointClassifier::make_triangle_records_for_intersection_that_must_exist(
    mesh::MeshEntityPtr quadEl, mesh::MeshEntityPtr triEl, PossibleEdgeIntersectionForTriangle p)
{
  assert(p.record.id != -1);
  const int* vertmap1 =
      triEl == m_quadToTriangles.el1 ? m_quadToTriangles.VERTMAP_TRI1_TO_QUAD : m_quadToTriangles.VERTMAP_TRI2_TO_QUAD;
  const int* edgemap1 =
      triEl == m_quadToTriangles.el1 ? m_quadToTriangles.EDGEMAP_TRI1_TO_QUAD : m_quadToTriangles.EDGEMAP_TRI2_TO_QUAD;

  if (p.record.type == PointClassification::Vert)
  {
    auto pr = make_triangle_records_for_quad_vert_intersection(triEl, p.record.id);
    PointRecord r(PointClassification::Vert, vertmap1[p.record.id], quadEl);
    r.m_r1 = pr.first;
    r.m_r2 = pr.second;
    return {r, -1.0};
  } else // type must be edge
  {
    assert(p.record.type == PointClassification::Edge);
    assert(p.record.id != m_quadToTriangles.INTERIOR_EDGE);
    double xi = p.edgeXi;
    auto pr   = make_triangle_records_for_quad_edge_intersection(triEl, p.record.id, xi);
    PointRecord r(PointClassification::Edge, edgemap1[p.record.id], quadEl);
    r.m_r1 = pr.first;
    r.m_r2 = pr.second;

    xi = convert_edge_xi(xi, quadEl, edgemap1[p.record.id]);
    return {r, xi};
  }
}

} // namespace impl

} // namespace predicates
} // namespace middle_mesh
} // namespace stk
