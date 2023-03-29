#include "stk_middle_mesh/predicates/point_classifier.hpp"

#include <cassert>
#include <cmath>
#include <limits>

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

void PointClassifierForTriangle::set_projection(utils::impl::Projection* proj)
{
  m_elemOps.set_projection(proj);
}

PointRecordForTriangle PointClassifierForTriangle::classify(mesh::MeshEntityPtr face, const utils::Point& pt)
{
  utils::Point ptXi = m_elemOps.compute_xi_coords(face, pt);

  PointRecordForTriangle r;
  if (face->get_type() == stk::middle_mesh::mesh::MeshEntityType::Triangle)
    r = m_xiClassifier.classify_triangle(ptXi);
  else if (face->get_type() == stk::middle_mesh::mesh::MeshEntityType::Quad)
    r = m_xiClassifier.classify_quad(ptXi);
  else
    throw std::invalid_argument("face must be either a triangle and a quad");

  r.el = face;
  return r;
}

double PointClassifierForTriangle::get_edge_xi(const PointRecordForTriangle& record)
{
  return m_triangleCoordUtils.get_edge_xi(record);
}

double PointClassifierForTriangle::get_edge_xi(stk::middle_mesh::mesh::MeshEntityType type, const int id,
                                               const utils::Point& ptXi, mesh::EntityOrientation orient)
{
  return m_triangleCoordUtils.get_edge_xi(type, id, ptXi, orient);
}

PossibleEdgeIntersectionForTriangle PointClassifierForTriangle::get_edge_xi(const PointRecordForTriangle& record1,
                                                                            const PointRecordForTriangle& record2)
{
  assert((record1.type == PointClassification::Interior && record2.type == PointClassification::Exterior) ||
         (record1.type == PointClassification::Exterior && record2.type == PointClassification::Interior));
  assert(record1.el == record2.el);
  auto el = record1.el;
  assert(el->get_type() == stk::middle_mesh::mesh::MeshEntityType::Triangle ||
         el->get_type() == stk::middle_mesh::mesh::MeshEntityType::Quad);

  std::pair<int, double> edgeIntersection;
  if (el->get_type() == stk::middle_mesh::mesh::MeshEntityType::Triangle)
    edgeIntersection = get_edge_xi_triangle(el, record1, record2);
  else if (el->get_type() == stk::middle_mesh::mesh::MeshEntityType::Quad)
    edgeIntersection = get_edge_xi_quad(el, record1, record2);

  return classify_intersection_on_vert_if_possible(el, edgeIntersection);
}

CornerRecordForTriangle PointClassifierForTriangle::get_edge_xi_corner(const PointRecordForTriangle& record1,
                                                                       const PointRecordForTriangle& record2)
{
  assert(record1.el == record2.el);
  auto el = record1.el;

  assert(el->get_type() == stk::middle_mesh::mesh::MeshEntityType::Quad ||
         el->get_type() == stk::middle_mesh::mesh::MeshEntityType::Triangle);

  assert(record1.type == PointClassification::Exterior || record1.type == PointClassification::Edge);
  assert(record2.type == PointClassification::Exterior || record2.type == PointClassification::Edge);

  double dxi1 = record2.m_ptXi.x - record1.m_ptXi.x;
  double dxi2 = record2.m_ptXi.y - record1.m_ptXi.y;

  // for vertical and horizontal lines, which edges they might intersect
  const static int quadIds[4] = {0, 2, 1, 3};
  const static int triIds[4]  = {0, 1, 1, 2};
  const int* ids              = quadIds;
  if (el->get_type() == stk::middle_mesh::mesh::MeshEntityType::Triangle)
    ids = triIds;

  CornerRecordForTriangle recordOut;
  if (std::abs(dxi1) > m_smallSlopeTol && std::abs(dxi2) > m_smallSlopeTol)
  {
    compute_intersection_general_corner(el, record1, record2, recordOut);
  } else if (std::abs(dxi1) < m_smallSlopeTol)
  {
    // vertical line
    int id1 = ids[0], id2 = ids[1];
    compute_intersection_degenerate_line(el, record1, record2, id1, id2, recordOut);
  } else if (std::abs(dxi2) < m_smallSlopeTol)
  {
    // horizontal line
    int id1 = ids[2], id2 = ids[3];
    compute_intersection_degenerate_line(el, record1, record2, id1, id2, recordOut);
  } else
    throw std::invalid_argument("line is too small");

  return recordOut;
}

double PointClassifierForTriangle::get_edge_xi_orthogonal(const PointRecordForTriangle& record1, const int id)
{
  auto ptXi = m_triangleCoordUtils.orthogonal_proj(record1.el->get_type(), record1.m_ptXi, id);
  return get_edge_xi(record1.el->get_type(), id, ptXi, record1.el->get_down_orientation(id));
}

PossibleEdgeIntersectionForTriangle
PointClassifierForTriangle::get_edge_intersection_xi(const PointRecordForTriangle& record1,
                                                     const PointRecordForTriangle& record2, const bool excludeEdgeVerts)
{
  assert(record1.el == record2.el);
  auto el1 = record1.el;

  assert(el1->get_type() == stk::middle_mesh::mesh::MeshEntityType::Quad ||
         el1->get_type() == stk::middle_mesh::mesh::MeshEntityType::Triangle);
  assert(((record1.type == PointClassification::Vert || record1.type == PointClassification::Edge) &&
          record2.type == PointClassification::Exterior) ||
         ((record2.type == PointClassification::Vert || record2.type == PointClassification::Edge) &&
          record1.type == PointClassification::Exterior));

  const PointRecordForTriangle& recordI = record1.type == PointClassification::Exterior ? record2 : record1;
  const PointRecordForTriangle& recordE = record1.type == PointClassification::Exterior ? record1 : record2;

  double dxi1 = recordE.m_ptXi.x - recordI.m_ptXi.x;
  double dxi2 = recordE.m_ptXi.y - recordI.m_ptXi.y;

  bool isLineHorizontal = std::abs(dxi2) < m_eps;
  bool isLineVertical   = std::abs(dxi1) < m_eps;

  // setup masks and return early in some simple cases
  bool excludedEdges[mesh::MAX_DOWN], excludedVerts[mesh::MAX_DOWN];
  for (int i = 0; i < mesh::MAX_DOWN; ++i)
  {
    excludedEdges[i] = false;
    excludedVerts[i] = false;
  }

  int nedges = el1->count_down();
  int nverts = nedges; // true for quads and triangles
  PointRecordForTriangle recordOut(PointClassification::Exterior, -1, el1, utils::Point(0, 0));
  if (recordI.type == PointClassification::Vert)
  {
    // can't possibly intersect a non-excluded edge
    if (excludeEdgeVerts && (isLineHorizontal || isLineVertical))
      return {recordOut, recordOut.m_ptXi.x};

    excludedEdges[recordI.id]                         = true;
    excludedEdges[(recordI.id - 1 + nedges) % nedges] = true;

    excludedVerts[recordI.id] = true;
    if (excludeEdgeVerts)
    {
      excludedVerts[(recordI.id + 1) % nverts]          = true;
      excludedVerts[(recordI.id - 1 + nverts) % nverts] = true;
    }
  } else if (recordI.type == PointClassification::Edge)
  {
    bool isEdgeHorizontal, isEdgeVertical;
    if (el1->get_type() == stk::middle_mesh::mesh::MeshEntityType::Quad)
    {
      isEdgeHorizontal = (recordI.id % 2) == 0;
      isEdgeVertical   = !isEdgeHorizontal;
    } else // triangle
    {
      isEdgeHorizontal = recordI.id == 0;
      isEdgeVertical   = recordI.id == 2;
    }

    if (excludeEdgeVerts)
    {
      if (isLineHorizontal && isEdgeHorizontal)
        return {recordOut, recordOut.m_ptXi.x};

      if (isLineVertical && isEdgeVertical)
        return {recordOut, recordOut.m_ptXi.x};
    }

    excludedEdges[recordI.id] = true;
    if (excludeEdgeVerts)
    {
      excludedVerts[recordI.id]                = true;
      excludedVerts[(recordI.id + 1) % nverts] = true;
    }
  } else
    throw std::runtime_error("unreachable case");

  recordOut = get_edge_intersection_xi(el1, record1, record2, excludedVerts, excludedEdges);

  return {recordOut, recordOut.m_ptXi.x};
}

utils::Point PointClassifierForTriangle::compute_xyz_coords(const PointRecordForTriangle& record)
{
  return m_triangleCoordUtils.compute_xyz_coords(record);
}

double PointClassifierForTriangle::compute_orthogonal_dist(const PointRecordForTriangle& record1, const int id)
{
  return m_triangleCoordUtils.compute_orthogonal_dist(record1, id);
}

void PointClassifierForTriangle::compute_intersection_general_corner(mesh::MeshEntityPtr el,
                                                                     const PointRecordForTriangle& record1,
                                                                     const PointRecordForTriangle& record2,
                                                                     CornerRecordForTriangle& recordOut)
{
  double alphas[mesh::MAX_DOWN];
  double xis[mesh::MAX_DOWN];
  int nedges   = el->count_down();
  int nInRange = 0;
  for (int i = 0; i < nedges; ++i)
  {
    alphas[i] = compute_edge_intersection(el->get_type(), record1, record2, i);

    auto ptXi = compute_point_on_line(record1, record2, alphas[i]);
    xis[i]    = get_edge_xi(el->get_type(), i, ptXi, el->get_down_orientation(i));

    if (in_range(alphas[i], 0, 1, m_eps) && in_range(xis[i], 0, 1, m_eps))
      nInRange++;
  }

  int id1 = -1, id2 = -1;
  if (nInRange == 2)
  {
    // figure out the two edges intersected
    for (int i = 0; i < nedges; ++i)
      if (in_range(alphas[i], 0, 1, m_eps) && in_range(xis[i], 0, 1, m_eps))
        id1 == -1 ? id1 = i : id2 = i;

    // TODO: should be r_in_range > 0, to allow for several matches due to
    //       tolerances
  } else if (nInRange == 1)
  {
    // this shouldn't be possible, assume there are two intersections and
    // find the edge that is closest to being intersected

    // get the intersect that is in range
    double minDeviation = std::numeric_limits<double>::max();
    for (int i = 0; i < nedges; ++i)
      if (in_range(alphas[i], 0, 1, m_eps) && in_range(xis[i], 0, 1, m_eps))
        id1 = i;
      else
      {
        double deviationA  = std::min(std::abs(alphas[i]), std::abs(alphas[i] - 1));
        double deviationXi = std::min(std::abs(xis[i]), std::abs(xis[i] - 1));
        double deviation   = deviationA + deviationXi;
        if (deviation < minDeviation)
        {
          minDeviation = deviation;
          id2          = i;
        }
      }
  } else if (nInRange > 2)
  {
    // pick the two intersections that are most "central"
    double deviation1 = std::numeric_limits<double>::max();
    double deviation2 = std::numeric_limits<double>::max();
    for (int i = 0; i < nedges; ++i)
      if (in_range(alphas[i], 0, 1, m_eps) && in_range(xis[i], 0, 1, m_eps))
      {
        auto deviation = std::abs(alphas[i] - 0.5) + std::abs(xis[i] - 0.5);
        if (deviation < deviation1)
        {
          id2        = id1;
          deviation2 = deviation1;
          id1        = i;
          deviation1 = deviation;
        } else if (deviation < deviation2)
        {
          id2        = i;
          deviation2 = deviation;
        }
      }
  }

  if (nInRange >= 1)
  {
    assert(id1 != -1 && id2 != -1);
    assert(id1 != id2);
    set_corner_record(el, id1, id2, xis[id1], xis[id2], recordOut);
  } // else no intersection, no need to change record_out
}

void PointClassifierForTriangle::compute_intersection_degenerate_line(mesh::MeshEntityPtr el,
                                                                      const PointRecordForTriangle& record1,
                                                                      const PointRecordForTriangle& record2,
                                                                      const int id1, const int id2,
                                                                      CornerRecordForTriangle& recordOut)
{
  double alpha1 = compute_edge_intersection(el->get_type(), record1, record2, id1);
  double alpha2 = compute_edge_intersection(el->get_type(), record1, record2, id2);

  auto ptXi1 = compute_point_on_line(record1, record2, alpha1);
  auto ptXi2 = compute_point_on_line(record1, record2, alpha2);
  double xi1 = get_edge_xi(el->get_type(), id1, ptXi1, el->get_down_orientation(id1));
  double xi2 = get_edge_xi(el->get_type(), id2, ptXi2, el->get_down_orientation(id2));

  // if either one is in range, force them both to be
  if ((in_range(alpha1, 0, 1, m_eps) && in_range(xi1, 0, 1, m_eps)) ||
      (in_range(alpha2, 0, 1, m_eps) && in_range(xi2, 0, 1, m_eps)))
  {
    set_corner_record(el, id1, id2, xi1, xi2, recordOut);
  }
}

void PointClassifierForTriangle::set_corner_record(mesh::MeshEntityPtr el, int id1, int id2, double xi1, double xi2,
                                                   CornerRecordForTriangle& recordOut)
{
  xi1 = clamp(xi1, 0, 1);
  xi2 = clamp(xi2, 0, 1);

  auto p1                   = classify_intersection_on_vert_if_possible(el, std::make_pair(id1, xi1));
  auto p2                   = classify_intersection_on_vert_if_possible(el, std::make_pair(id2, xi2));
  recordOut.hasIntersection = true;
  recordOut.record1         = p1.record;
  recordOut.xi1             = p1.edgeXi;

  // TODO: make this branchless
  if (!(p1.record.type == p2.record.type && p1.record.id == p2.record.id))
  {
    recordOut.record2 = p2.record;
    recordOut.xi2     = p2.edgeXi;
  } else
  {
    recordOut.xi2 = -1;
  }
}

std::pair<int, double> PointClassifierForTriangle::get_edge_xi_triangle(mesh::MeshEntityPtr el,
                                                                        const PointRecordForTriangle& record1,
                                                                        const PointRecordForTriangle& record2)
{
  double dxi1             = record2.m_ptXi.x - record1.m_ptXi.x;
  double dxi2             = record2.m_ptXi.y - record1.m_ptXi.y;
  const static int ids[4] = {0, 1, 1, 2};

  if (std::abs(dxi1) > m_smallSlopeTol && std::abs(dxi2) > m_smallSlopeTol)
  {
    for (int i = 0; i < 3; ++i)
    {
      auto alpha = compute_edge_intersection_triangle(record1, record2, i);
      auto pt    = compute_point_on_line(record1, record2, alpha);
      double xi  = get_edge_xi(el->get_type(), i, pt, el->get_down_orientation(i));

      if (in_range(alpha, 0, 1, m_eps) && in_range(xi, 0, 1, m_eps))
      {
        return std::make_pair(i, xi);
      }
    }

  } else if ((std::abs(dxi1) < m_smallSlopeTol) || (std::abs(dxi2) < m_smallSlopeTol))
  {
    int offset = std::abs(dxi1) < m_smallSlopeTol ? 0 : 2;
    int id1 = ids[offset], id2 = ids[offset + 1];

    double alpha1 = compute_edge_intersection_triangle(record1, record2, id1);
    auto pt1      = compute_point_on_line(record1, record2, alpha1);
    double xi1    = get_edge_xi(el->get_type(), id1, pt1, el->get_down_orientation(id1));

    double alpha2 = compute_edge_intersection_triangle(record1, record2, id2);
    auto pt2      = compute_point_on_line(record1, record2, alpha2);
    double xi2    = get_edge_xi(el->get_type(), id2, pt2, el->get_down_orientation(id2));

    if (in_range(alpha1, 0, 1, m_eps) && in_range(xi1, 0, 1, m_eps))
    {
      return std::make_pair(id1, xi1);
    } else if (in_range(alpha2, 0, 1, m_eps) && in_range(xi2, 0, 1, m_eps))
    {
      return std::make_pair(id2, xi2);
    } else
    {
      return compute_xi_closest(el, record1, record2);
    }
  } // else
  //  {
  throw std::invalid_argument("line is too small");
  //  }
}

std::pair<int, double> PointClassifierForTriangle::get_edge_xi_quad(mesh::MeshEntityPtr el,
                                                                    const PointRecordForTriangle& record1,
                                                                    const PointRecordForTriangle& record2)
{
  double dxi1 = record2.m_ptXi.x - record1.m_ptXi.x;
  double dxi2 = record2.m_ptXi.y - record1.m_ptXi.y;

  // first two entries are the faces that a vertical line might intersect,
  // second two are for horizontal lines
  const static int ids[4] = {2, 0, 3, 1};

  if (std::abs(dxi1) > m_smallSlopeTol && std::abs(dxi2) > m_smallSlopeTol)
  {
    for (int i = 0; i < 4; ++i)
    {
      auto alpha = compute_edge_intersection_quad(record1, record2, i);
      auto pt    = compute_point_on_line(record1, record2, alpha);
      double xi  = get_edge_xi(el->get_type(), i, pt, el->get_down_orientation(i));

      if (in_range(alpha, 0, 1, m_eps) && in_range(xi, 0, 1, m_eps))
      {
        return std::make_pair(i, xi);
      }
    }

    throw std::runtime_error("did not find intersection");
  } else if (std::abs(dxi1) < m_smallSlopeTol || std::abs(dxi2) < m_smallSlopeTol)
  {
    // horizontal or vertical lines can possibly intersect two faces
    int offset = std::abs(dxi1) < m_smallSlopeTol ? 0 : 2;
    int id1 = ids[offset], id2 = ids[offset + 1];

    double alpha1 = compute_edge_intersection_quad(record1, record2, id1);
    auto pt1      = compute_point_on_line(record1, record2, alpha1);
    double xi1    = get_edge_xi(el->get_type(), id1, pt1, el->get_down_orientation(id1));

    double alpha2 = compute_edge_intersection_quad(record1, record2, id2);
    auto pt2      = compute_point_on_line(record1, record2, alpha2);
    double xi2    = get_edge_xi(el->get_type(), id2, pt2, el->get_down_orientation(id2));

    if (in_range(alpha1, 0, 1, m_eps) && in_range(xi1, 0, 1, m_eps))
    {
      return std::make_pair(id1, xi1);
    } else if (in_range(alpha2, 0, 1, m_eps) && in_range(xi2, 0, 1, m_eps))
    {
      return std::make_pair(id2, xi2);
    } else
      return compute_xi_closest(el, record1, record2);
  } else // both dxi1 and dxi2 are too small
    throw std::invalid_argument("line is too small");
}

PossibleEdgeIntersectionForTriangle
PointClassifierForTriangle::classify_intersection_on_vert_if_possible(mesh::MeshEntityPtr el,
                                                                      const std::pair<int, double>& edgeIntersection)
{
  int edgeId    = edgeIntersection.first;
  double edgeXi = edgeIntersection.second;

  int numVerts = el->count_down();
  if (std::abs(edgeXi) < get_eps())
  {
    int vertId = edgeId;
    if (el->get_down_orientation(edgeId) == mesh::EntityOrientation::Reversed)
      vertId = (edgeId + 1) % numVerts;

    return {PointRecordForTriangle(PointClassification::Vert, vertId, el), -1};
  } else if (std::abs(edgeXi - 1) < get_eps())
  {
    int vertId = (edgeId + 1) % numVerts;
    if (el->get_down_orientation(edgeId) == mesh::EntityOrientation::Reversed)
      vertId = edgeId;

    return {PointRecordForTriangle(PointClassification::Vert, vertId, el), -1};
  } else
  {
    return {PointRecordForTriangle(PointClassification::Edge, edgeId, el), edgeXi};
  }
}

utils::Point PointClassifierForTriangle::compute_point_on_line(const PointRecordForTriangle& record1,
                                                               const PointRecordForTriangle& record2,
                                                               const double alpha)
{
  double dxi1 = record2.m_ptXi.x - record1.m_ptXi.x;
  double dxi2 = record2.m_ptXi.y - record1.m_ptXi.y;
  double xi1a = record1.m_ptXi.x;
  double xi2a = record1.m_ptXi.y;

  return utils::Point(dxi1 * alpha + xi1a, dxi2 * alpha + xi2a);
}

double PointClassifierForTriangle::compute_edge_intersection(stk::middle_mesh::mesh::MeshEntityType type,
                                                             const PointRecordForTriangle& record1,
                                                             const PointRecordForTriangle& record2, const int id)
{
  if (type == stk::middle_mesh::mesh::MeshEntityType::Quad)
    return compute_edge_intersection_quad(record1, record2, id);
  else if (type == stk::middle_mesh::mesh::MeshEntityType::Triangle)
    return compute_edge_intersection_triangle(record1, record2, id);
  else
    throw std::invalid_argument("type must be Quad or Triangle");
}

double PointClassifierForTriangle::compute_edge_intersection_quad(const PointRecordForTriangle& record1,
                                                                  const PointRecordForTriangle& record2, const int id)
{
  double dxi1 = record2.m_ptXi.x - record1.m_ptXi.x;
  double dxi2 = record2.m_ptXi.y - record1.m_ptXi.y;
  double xi1a = record1.m_ptXi.x;
  double xi2a = record1.m_ptXi.y;

  switch (id)
  {
    case 0: {
      return -xi2a / dxi2;
    }
    case 1: {
      return (1 - xi1a) / dxi1;
    }
    case 2: {
      return (1 - xi2a) / dxi2;
    }
    case 3: {
      return -xi1a / dxi1;
    }
    default:
      throw std::invalid_argument("invalid id");
  }
}

double PointClassifierForTriangle::compute_edge_intersection_triangle(const PointRecordForTriangle& record1,
                                                                      const PointRecordForTriangle& record2,
                                                                      const int id)
{
  double xi1a = record1.m_ptXi.x;
  double xi2a = record1.m_ptXi.y;

  double dxi1 = record2.m_ptXi.x - record1.m_ptXi.x;
  double dxi2 = record2.m_ptXi.y - record1.m_ptXi.y;

  switch (id)
  {
    case 0: {
      return -xi2a / dxi2;
    }
    case 1: {
      return (1 - xi1a - xi2a) / (dxi1 + dxi2);
    }
    case 2: {
      return -xi1a / dxi1;
    }
    default:
      throw std::invalid_argument("invalid id");
  }
}

std::pair<int, double> PointClassifierForTriangle::compute_xi_closest(mesh::MeshEntityPtr el,
                                                                      const PointRecordForTriangle& record1,
                                                                      const PointRecordForTriangle& record2)
{
  double minDist = std::numeric_limits<double>::max();
  int id         = 0;
  for (int i = 0; i < el->count_down(); ++i)
  {
    double dist = 0;
    dist += m_triangleCoordUtils.compute_orthogonal_dist(el->get_type(), record1.m_ptXi, i);
    dist += m_triangleCoordUtils.compute_orthogonal_dist(el->get_type(), record2.m_ptXi, i);
    if (dist < minDist)
    {
      minDist = dist;
      id      = i;
    }
  }

  // compute midpoint and project onto edge
  auto ptMid = (record1.m_ptXi + record2.m_ptXi) / 2;
  auto pt    = m_triangleCoordUtils.orthogonal_proj(el->get_type(), ptMid, id);
  pt.x       = clamp(pt.x, 0, 1);
  pt.y       = clamp(pt.y, 0, 1);
  double xi  = get_edge_xi(el->get_type(), id, pt, el->get_down_orientation(id));

  return std::make_pair(id, xi);
}

PointRecordForTriangle PointClassifierForTriangle::get_edge_intersection_xi(mesh::MeshEntityPtr el1,
                                                                            const PointRecordForTriangle& inputRecordI,
                                                                            const PointRecordForTriangle& recordE,
                                                                            const bool excludedVerts[mesh::MAX_DOWN],
                                                                            const bool excludedEdges[mesh::MAX_DOWN])
{
  double alphas[mesh::MAX_DOWN], xis[mesh::MAX_DOWN];
  PointRecordForTriangle recordsOut[mesh::MAX_DOWN];
  for (int i = 0; i < mesh::MAX_DOWN; ++i)
    recordsOut[i].id = -1;

  int nedges = el1->count_down();
  for (int i = 0; i < nedges; ++i)
  {
    if (excludedEdges[i])
      continue;

    // TODO: why does computePointOnLine take el->get_type()?  it can
    //       retrieve the information from record1.el
    double alpha = compute_edge_intersection(el1->get_type(), inputRecordI, recordE, i);
    auto pt      = compute_point_on_line(inputRecordI, recordE, alpha);
    double xi    = get_edge_xi(el1->get_type(), i, pt, el1->get_down_orientation(i));

    if (!(in_range(alpha, 0, 1, m_eps) && in_range(xi, 0, 1, m_eps)))
      continue;

    PointRecordForTriangle recordI;
    if (el1->get_type() == stk::middle_mesh::mesh::MeshEntityType::Triangle)
      recordI = m_xiClassifier.classify_triangle(pt);
    else // if (el1->get_type() == stk::middle_mesh::mesh::MeshEntityType::Quad)
      recordI = m_xiClassifier.classify_quad(pt);
    recordI.el = el1;

    if (recordI.type == PointClassification::Vert && excludedVerts[recordI.id])
      continue;

    // if we get here, this must be a valid intersection
    if (recordI.type == PointClassification::Vert)
      recordI.m_ptXi = utils::Point(0, 0);
    else
    {
      recordI.type   = PointClassification::Edge;
      recordI.m_ptXi = utils::Point(xi, 0);
    }
    recordsOut[i] = recordI;
    alphas[i]     = alpha;
    xis[i]        = xi;
  }

  // If there are multiple intersections, prefer the one that is:
  //   1. classified on a vertex, or
  //   2. closest to alpha = 0,5, xi = 0.5
  double minDeviation = std::numeric_limits<double>::max();
  int minIdx          = -1;
  int vertIdx         = -1;
  for (int i = 0; i < nedges; ++i)
  {
    if (excludedEdges[i] || recordsOut[i].id == -1)
      continue;

    if (inputRecordI.type == PointClassification::Vert)
      vertIdx = i;
    else // interior of edge
    {
      double deviation = std::abs(alphas[i] - 0.5) + std::abs(xis[i] - 0.5);
      if (deviation < minDeviation)
      {
        minDeviation = deviation;
        minIdx       = i;
      }
    }
  }

  if (minIdx == -1 && vertIdx == -1)
    return PointRecordForTriangle{PointClassification::Exterior, -1, el1, utils::Point(0, 0)};
  else if (vertIdx != -1)
    return recordsOut[vertIdx];
  else // min_idx != -1
    return recordsOut[minIdx];
}

} // namespace impl
} // namespace predicates
} // namespace middle_mesh
} // namespace stk
