#include "active_vert_data.hpp"
#include "distortion_metric.hpp"
#include "mesh_entity.hpp"
#include "utils.hpp"

namespace stk {
namespace middle_mesh {
namespace opt {
namespace impl {

const std::array<int, 3> ActiveVertData::M_VERTMAP_QUAD_TO_TRI1 = {0, 1, 3};
const std::array<int, 3> ActiveVertData::M_VERTMAP_QUAD_TO_TRI2 = {1, 2, 3};
const std::array<int, 3> ActiveVertData::M_VERTMAP_QUAD_TO_TRI3 = {0, 1, 2};
const std::array<int, 3> ActiveVertData::M_VERTMAP_QUAD_TO_TRI4 = {0, 2, 3};
const std::array<int, 4> ActiveVertData::M_EXCLUDED_TRI         = {1, 3, 0, 2};

void ActiveVertData::set_sizes(mesh::MeshEntityPtr vert)
{
  std::vector<mesh::MeshEntityPtr> els;
  int nelems     = get_upward(vert, 2, els);
  int ntriangles = 0;
  for (int i = 0; i < nelems; ++i)
    ntriangles += count_num_triangles(vert, els[i]);

  // size data structures
  m_triVertIndices.resize(ntriangles, 3);
  m_currentEntityIdx.resize(ntriangles);
}

void ActiveVertData::set_triangles(mesh::MeshEntityPtr vert)
{
  std::unordered_map<mesh::MeshEntityPtr, int> vertIndices;
  std::vector<mesh::MeshEntityPtr> els;
  int nelems = get_upward(vert, 2, els);

  int triIdx = 0;
  std::array<mesh::MeshEntityPtr, 12> triVerts;
  // ensure the vertex the patch is centered on is the first vert
  m_localVertsUnique.push_back(vert);
  m_allPointsOrig.push_back(vert->get_point_orig(0));
  vertIndices[vert] = 0;
  for (int i = 0; i < nelems; ++i)
  {
    int ntris  = get_triangle_verts(els[i], vert, triVerts);
    int nverts = 3 * ntris;

    // update verts_unique
    for (int j = 0; j < nverts; ++j)
    {
      auto vertJ = triVerts[j];
      if (vertIndices.count(vertJ) == 0)
      {
        m_localVertsUnique.push_back(vertJ);
        m_allPointsOrig.push_back(vertJ->get_point_orig(0));
        vertIndices[vertJ] = m_localVertsUnique.size() - 1;
      }
    }

    // update m_tri_verts
    for (int t = 0; t < ntris; ++t)
    {
      for (int j = 0; j < 3; ++j)
        m_triVertIndices(triIdx, j) = vertIndices[triVerts[3 * t + j]];

      triIdx++;
    }
  }

  m_allPoints = m_allPointsOrig;
}

// gets the local index of vert in each triangle and writes it to
// m_current_entity
void ActiveVertData::set_current_indices(mesh::MeshEntityPtr vert)
{
  for (unsigned int i = 0; i < m_currentEntityIdx.size(); ++i)
  {
    int idx = -1;
    for (int j = 0; j < 3; ++j)
      if (m_localVertsUnique[m_triVertIndices(i, j)] == vert)
      {
        idx = j;
        break;
      }

    assert(idx != -1);
    m_currentEntityIdx[i] = idx;
  }
}

int ActiveVertData::count_num_triangles(mesh::MeshEntityPtr /*vert*/, mesh::MeshEntityPtr el)
{
  if (el->get_type() == mesh::MeshEntityType::Triangle)
    return 1;
  else
  {
    return 3;
    /*
    // the convention is that the quad is split into triangles
    // along line 1-3
    int idx_local = getLocalIdx(el, vert);
    if (idx_local == 0 || idx_local == 2)
      return 1;
    else
      return 2;
    */
  }
}

// gets the local index of vert within el's verts
int ActiveVertData::get_local_idx(mesh::MeshEntityPtr el, mesh::MeshEntityPtr vert)
{
  mesh::MeshEntityPtr verts[mesh::MAX_DOWN];
  int dim   = get_type_dimension(vert->get_type());
  int ndown = mesh::get_downward(el, dim, verts);
  for (int i = 0; i < ndown; ++i)
  {  
    if (verts[i] == vert)
      return i;
  }

  throw std::invalid_argument("vert is not downward adjacent to el");
}

// returns the vertices of the triangles in tri_verts, and returns
// the number of triangles.  If there is only 1 triangle, only
// the first 3 entries of tri_verts are populated
int ActiveVertData::get_triangle_verts(mesh::MeshEntityPtr el, mesh::MeshEntityPtr vert,
                                       std::array<mesh::MeshEntityPtr, 12>& triVerts)
{
  mesh::MeshEntityPtr verts[mesh::MAX_DOWN];
  mesh::get_downward(el, 0, verts);
  if (el->get_type() == mesh::MeshEntityType::Triangle)
  {
    for (int i = 0; i < 3; ++i)
      triVerts[i] = verts[i];
    return 1;
  } else // quad
  {
    const std::array<const std::array<int, 3>, 4> vertmaps = {M_VERTMAP_QUAD_TO_TRI1, M_VERTMAP_QUAD_TO_TRI2,
                                                              M_VERTMAP_QUAD_TO_TRI3, M_VERTMAP_QUAD_TO_TRI4};

    int localIdx = get_local_idx(el, vert);
    int idx      = 0;
    for (int i = 0; i < 4; ++i)
    {
      if (i == M_EXCLUDED_TRI[localIdx])
        continue;

      for (int j = 0; j < 3; ++j)
      {
        triVerts[idx] = verts[vertmaps[i][j]];
        idx++;
      }
    }

    return 3;

    /*
    int local_idx = getLocalIdx(el, vert);
    if (local_idx == 0)
    {
      for (int i=0; i < 3; ++i)
        tri_verts[i] = verts[M_VERTMAP_QUAD_TO_TRI1[i]];
      return 1;
    } else if (local_idx == 2)
    {
      for (int i=0; i < 3; ++i)
        tri_verts[i] = verts[M_VERTMAP_QUAD_TO_TRI2[i]];
      return 1;
    } else
    {
      for (int i=0; i < 3; ++i)
      {
        tri_verts[i]     = verts[M_VERTMAP_QUAD_TO_TRI1[i]];
        tri_verts[i + 3] = verts[M_VERTMAP_QUAD_TO_TRI2[i]];
      }

      return 2;
    }
    */
  }
}

void ActiveVertData::add_local_vert(mesh::MeshEntityPtr vert)
{
  m_localVertsUnique.push_back(vert);
  m_allPointsOrig.push_back(vert->get_point_orig(0));
  m_allPoints.push_back(vert->get_point_orig(0));  
}

void ActiveVertData::add_remote_vert(const mesh::RemoteSharedEntity& owner)
{
  m_remoteVertOwners.push_back(owner);
  m_allPoints.emplace_back(0, 0, 0);
  m_allPointsOrig.emplace_back(0, 0, 0);
}

void ActiveVertData::add_remote_element(const std::array<int, 3>& triVertIndices)
{
  int nrows = m_triVertIndices.extent0();
  m_triVertIndices.resize(nrows+1, 3);
  int centerVertIdx = -1;
  for (int i=0; i < 3; ++i)
  {
    m_triVertIndices(nrows, i) = triVertIndices[i];
    if (triVertIndices[i] == 0)
      centerVertIdx = i;
  }

  assert(centerVertIdx >= 0 && centerVertIdx < 3);
  m_currentEntityIdx.push_back(centerVertIdx);
}

} // namespace impl

} // namespace opt
} // namespace middle_mesh
} // namespace stk
