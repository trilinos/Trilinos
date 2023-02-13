#include "predicates/quad_to_triangles.h"

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

// static variable defs
const int QuadToTriangles::VERTMAP_TRI1_TO_QUAD[3]                            = {0, 1, 3};
const int QuadToTriangles::VERTMAP_TRI2_TO_QUAD[3]                            = {2, 3, 1};
const int QuadToTriangles::EDGEMAP_TRI1_TO_QUAD[3]                            = {0, -1, 3};
const int QuadToTriangles::EDGEMAP_TRI2_TO_QUAD[3]                            = {2, -1, 1};
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
  } else if (r2.type == PointClassification::Edge && r2.id == INTERIOR_EDGE)
  {
    r1.type = PointClassification::Edge;
    r1.id   = INTERIOR_EDGE;
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
