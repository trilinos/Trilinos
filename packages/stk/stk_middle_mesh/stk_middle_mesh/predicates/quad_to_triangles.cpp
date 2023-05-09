#include "stk_middle_mesh/predicates/quad_to_triangles.hpp"

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

// static variable defs
const int QuadToTriangles::VERTMAP_TRI1_TO_QUAD[3]                            = {0, 1, 3};
const int QuadToTriangles::VERTMAP_TRI2_TO_QUAD[3]                            = {2, 3, 1};
const int QuadToTriangles::EDGEMAP_TRI1_TO_QUAD[3]                            = {0, -1, 3};
const int QuadToTriangles::EDGEMAP_TRI2_TO_QUAD[3]                            = {2, -1, 1};

const std::array<std::pair<int, int>, 4> QuadToTriangles::VERTMAP_QUAD_TO_TRI = 
    {std::make_pair(0, 0), std::make_pair(0, 1), std::make_pair(1, 0), std::make_pair(0, 2)};

const std::array<std::pair<int, int>, 4> QuadToTriangles::EDGEMAP_QUAD_TO_TRI = {
    std::make_pair(0, 0), std::make_pair(1, 2), std::make_pair(1, 0), std::make_pair(0, 2)};
const int QuadToTriangles::INTERIOR_EDGE = 1;

PointRecord QuadToTriangles::get_quad_record(mesh::MeshEntityPtr quad, PointRecordForTriangle& r1Input,
                                             PointRecordForTriangle& r2Input)
{
  PointRecord r3;
  r3.el = quad;

  PointRecordForTriangle& r1 = r1Input.el == el1 ? r1Input : r2Input;
  PointRecordForTriangle& r2 = r1Input.el == el1 ? r2Input : r1Input;

  assert(r2.el == el2 || r2.el == nullptr);

  enforce_record_consistency(r1, r2);
  if (r1.type == PointClassification::Vert)
  {
    r3.type = r1.type;
    r3.id   = VERTMAP_TRI1_TO_QUAD[r1.id];
    r3.m_r1 = r1;
    r3.m_r2 = r2;

  } else if (r2.type == PointClassification::Vert)
  {
    r3.type = r2.type;
    r3.id   = VERTMAP_TRI2_TO_QUAD[r2.id];
    r3.m_r1 = r2;
    r3.m_r2 = r1;
  } else if (r1.type == PointClassification::Edge && r1.id != INTERIOR_EDGE)
  {
    r3.type = r1.type;
    r3.id   = EDGEMAP_TRI1_TO_QUAD[r1.id];
    r3.m_r1 = r1;
    r3.m_r2 = r2;
  } else if (r2.type == PointClassification::Edge && r2.id != INTERIOR_EDGE)
  {
    r3.type = r2.type;
    r3.id   = EDGEMAP_TRI2_TO_QUAD[r2.id];
    r3.m_r1 = r2;
    r3.m_r2 = r1;
  } else if (r1.type == PointClassification::Interior ||
             (r1.type == PointClassification::Edge && r1.id == INTERIOR_EDGE))
  {
    r3.type = PointClassification::Interior;
    r3.id   = 0;
    r3.m_r1 = r1;
    r3.m_r2 = r2;
  } else if (r2.type == PointClassification::Interior ||
             (r2.type == PointClassification::Edge && r2.id == INTERIOR_EDGE))
  {
    r3.type = PointClassification::Interior;
    r3.id   = 0;
    r3.m_r1 = r2;
    r3.m_r2 = r1;
  } else // exterior
  {
    r3.type = PointClassification::Exterior;
    r3.id   = 0;
    r3.m_r1 = r1;
    r3.m_r2 = r2;
  }

  return r3;
}

void QuadToTriangles::enforce_record_consistency(PointRecordForTriangle& r1, PointRecordForTriangle& r2)
{
  if (r1.type == PointClassification::Vert)
  {
    int quadVertId = VERTMAP_TRI1_TO_QUAD[r1.id];
    int r2Id       = find_common_vertex(VERTMAP_TRI2_TO_QUAD, quadVertId);

    if (r2Id != -1)
    {
      r2.type = PointClassification::Vert;
      r2.id   = r2Id;
    }

  } else if (r2.type == PointClassification::Vert)
  {
    int quadVertId = VERTMAP_TRI2_TO_QUAD[r2.id];
    int r1Id       = find_common_vertex(VERTMAP_TRI1_TO_QUAD, quadVertId);

    if (r1Id != -1)
    {
      r1.type = PointClassification::Vert;
      r1.id   = r1Id;
    }
  } else if (r1.type == PointClassification::Edge && r1.id == INTERIOR_EDGE)
  {
    r2.type = PointClassification::Edge;
    r2.id   = INTERIOR_EDGE;
    r2.m_ptXi = {1 - r1.m_ptXi.x, 1 - r1.m_ptXi.y};
  } else if (r2.type == PointClassification::Edge && r2.id == INTERIOR_EDGE)
  {
    r1.type = PointClassification::Edge;
    r1.id   = INTERIOR_EDGE;
    r1.m_ptXi = {1 - r2.m_ptXi.x, 1 - r2.m_ptXi.y};
  }
}

utils::Point QuadToTriangles::compute_xyz_coords(const PointRecord& quadRecord, bool allowExterior)
{
  const PointRecordForTriangle* r;
  if (quadRecord.m_r1.type != PointClassification::Exterior || quadRecord.m_r2.type != PointClassification::Exterior)
    r = quadRecord.m_r1.type != PointClassification::Exterior ? &(quadRecord.m_r1) : &(quadRecord.m_r2);
  else
  {
    double deviation1 = m_triangleUtils.compute_exterior_deviation(quadRecord.m_r1);
    double deviation2 = m_triangleUtils.compute_exterior_deviation(quadRecord.m_r2);
    r = deviation1 < deviation2 ? &(quadRecord.m_r1) : &(quadRecord.m_r2);
  }

  return m_triangleUtils.compute_xyz_coords(*r, allowExterior);
}

utils::Point QuadToTriangles::get_quad_xi_coords(const PointRecord& quadRecord, bool allowExterior)
{
  if (quadRecord.type == PointClassification::Vert)
  {
    switch (quadRecord.id)
    {
      case 0: return {0, 0};
      case 1: return {1, 0};
      case 2: return {1, 1};
      case 3: return {0, 1};
      default:
        throw std::runtime_error("unrecognized quadRecord vertex id");
    }
  } else if (quadRecord.type == PointClassification::Edge)
  {
    const PointRecordForTriangle& r1 = quadRecord.m_r1;


    if (quadRecord.id == 0)
    {
      assert(r1.el == el1);
      return {r1.m_ptXi.x, 0};
    } else if (quadRecord.id == 1)
    {
      assert(r1.el == el2);
      return {1, 1 - r1.m_ptXi.y};
    } else if (quadRecord.id == 2)
    {
      assert(r1.el == el2);
      return {1 - r1.m_ptXi.x, 1};
    } else if (quadRecord.id == 3)
    {
      assert(r1.el == el1);
      return {0, r1.m_ptXi.y};
    } else
      throw std::runtime_error("unrecognized quadRecord edge id");
  } else if (quadRecord.type == PointClassification::Interior || (quadRecord.type == PointClassification::Exterior && allowExterior))
  {
    const PointRecordForTriangle& r1 = quadRecord.m_r1;
    //const PointRecordForTriangle& r2 = quadRecord.m_r2;
    assert(r1.type == PointClassification::Interior || r1.type == PointClassification::Exterior || 
           (r1.type == PointClassification::Edge && r1.id == INTERIOR_EDGE));

    // use linear interpolation to compute quad xi coordinates
    // The algebra simplifies to the result below
    if (r1.el == el1)
    {
      return {r1.m_ptXi.x, r1.m_ptXi.y};
    } else
    {
      return {1 - r1.m_ptXi.x, 1 - r1.m_ptXi.y};
    }
  } else
    throw std::runtime_error("cannot compute xi coordinates on the exterior of a quad");
}


PointRecord QuadToTriangles::create_record(mesh::MeshEntityPtr quad, int edgeId, double edgeXi)
{
  assert(edgeId >= 0 && edgeId < 4);

  int triId     = EDGEMAP_QUAD_TO_TRI[edgeId].first;
  int triEdgeId = EDGEMAP_QUAD_TO_TRI[edgeId].second;

  mesh::MeshEntityPtr el = triId == 0 ? el1 : el2;
  // because of the way the reference elements are defined, the triangleEdgeXi
  // is always the same as the quad edge xi
  double triangleEdgeXi = edgeXi;

  PointRecordForTriangle r1 = m_triangleUtils.create_record(el, triEdgeId, triangleEdgeXi);
  PointRecord record(PointClassification::Edge, edgeId, quad, r1);
  enforce_record_consistency(r1, record.m_r2);

  return record;
}

PointRecord QuadToTriangles::create_record(mesh::MeshEntityPtr quad, int vertId)
{
  int elId      = VERTMAP_QUAD_TO_TRI[vertId].first;
  int triVertId = VERTMAP_QUAD_TO_TRI[vertId].second;

  mesh::MeshEntityPtr el = elId == 0 ? el1 : el2;
  PointRecordForTriangle r1 = m_triangleUtils.create_record(el, triVertId);
  PointRecordForTriangle r2;

  enforce_record_consistency(r1, r2);

  PointRecord record(PointClassification::Vert, vertId, quad, r1);
  record.m_r2 = r2;

  return record;
}


PointRecord QuadToTriangles::classify_onto(const PointRecord& record, mesh::MeshEntityPtr el)
{
  assert(record.type == PointClassification::Vert || record.type == PointClassification::Edge);
  assert(record.el->get_type() == mesh::MeshEntityType::Quad);
  assert(record.el->get_type() == el->get_type());

  if (record.el == el)
    return record;

  int dim = record.type == PointClassification::Vert ? 0 : 1;
  mesh::MeshEntityPtr entityOld = get_entity(record);
  std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> entitiesNew;
  int nentities = mesh::get_downward(el, dim, entitiesNew.data());

  int entityNewLocalId = -1;
  for (int i=0; i < nentities; ++i)
    if (entitiesNew[i] == entityOld)
    {
      entityNewLocalId = i;
      break;
    }


  if (record.type == PointClassification::Vert)
  {
    // get vertex on record.el
    // find corresponding vertex on el
    // use (create) VERTMAP_QUAD_TO_TRI to figure out triangle, vertex id
    // call TriangleCoordUtils::create_record(el, vertexId)
    // call enforceRecordConsistency() to create r2
    if (entityNewLocalId == -1)
      throw std::runtime_error("specified vertex is not shared between record.el and el");

    return create_record(el, entityNewLocalId);

  } else // record.type == Edge
  {
    if (entityNewLocalId == -1)
      throw std::runtime_error("specified edge is not shared between record.el and el");
    // get edge on record.el
    // get corresponding edge on el
    // get quad edgeXi on record.el
    // convert quad edgeXi to edge on el
    // call create_record(el, edgeId, quadEdgeXiNewElement)
    utils::Point ptXi = get_quad_xi_coords(record);
    double oldEdgeXi;
    switch (record.id)
    {
      case 0: {oldEdgeXi = ptXi.x;     break;}
      case 1: {oldEdgeXi = ptXi.y;     break;}
      case 2: {oldEdgeXi = 1 - ptXi.x; break; }
      case 3: {oldEdgeXi = 1 - ptXi.y; break;}
      default:
        throw std::runtime_error("a quad only has 4 edges");
    };

    // because of the way the reference element is defined, the
    // edge xi on adjacent elements is always reversed
    double newEntityEdgeXi = 1 - oldEdgeXi;

    return create_record(el, entityNewLocalId, newEntityEdgeXi);
  }
}


// returns the vertex index on the triangle that is on the given quad
// vertex id.  Returns -1 if not found
int QuadToTriangles::find_common_vertex(const int* vertmap, int quadId)
{
  for (int i = 0; i < 3; ++i)
    if (vertmap[i] == quadId)
      return i;

  return -1;
}

} // namespace impl

} // namespace predicates
} // namespace middle_mesh
} // namespace stk
