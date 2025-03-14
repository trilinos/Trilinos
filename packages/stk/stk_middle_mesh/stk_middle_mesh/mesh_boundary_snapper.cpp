#include "mesh_boundary_snapper.hpp"
#include "field.hpp"
#include <cassert>
#include <stdexcept>

#include "mesh_io.hpp" //TODO: DEBUGGING

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

namespace {
using VertEdgePair = MeshBoundarySnapper::VertEdgePair;
using MEPair       = std::pair<MeshEntityPtr, MeshEntityPtr>;
} // namespace

void MeshBoundarySnapper::snap(std::shared_ptr<Mesh> mesh1, std::shared_ptr<Mesh> mesh2, MPI_Comm unionComm)
{
  auto comm = unionComm;
  auto nprocs = utils::impl::comm_size(comm);
  auto rank = utils::impl::comm_rank(comm);
  auto root = 0;

  if (nprocs == 1) {
    start_local_snap(mesh1, mesh2);
  } else {
    MeshExchangeBoundaryEdges exchange1(mesh1, unionComm, root);
    auto meshToSnap1 = exchange1.get_boundary_edge_mesh();

    MeshExchangeBoundaryEdges exchange2(mesh2, unionComm, root);
    auto meshToSnap2 = exchange2.get_boundary_edge_mesh();

    if(rank == root) {
      start_local_snap(meshToSnap1, meshToSnap2);
    }
    
    exchange1.update_remote_vertices();
    exchange2.update_remote_vertices();
  }
}

void MeshBoundarySnapper::start_local_snap(std::shared_ptr<Mesh> mesh1, std::shared_ptr<Mesh> mesh2)
{
  m_seenEdges1  = create_field<FakeBool>(mesh1, FieldShape(0, 1, 0), 1, false);
  m_seenEdges2  = create_field<FakeBool>(mesh2, FieldShape(0, 1, 0), 1, false);
  int nedges1   = count_boundary_edges(mesh1);
  m_nedgesSeen1 = 0;
  int nedges2   = count_boundary_edges(mesh2);
  m_nedgesSeen2 = 0;

  VertEdgePair p1, p2;
  while (m_nedgesSeen1 != nedges1 && m_nedgesSeen2 != nedges2)
  {
    double maxSnapDist = 0;
    get_start_entity(mesh1, mesh2, p1, p2);
    snap_verts(p1, p2, maxSnapDist);
    check_zero_length_edges(p1);
    check_zero_length_edges(p2);
  }

  // free memory
  m_seenEdges1 = nullptr;
  m_seenEdges2 = nullptr;
}

void MeshBoundarySnapper::get_start_entity(std::shared_ptr<Mesh> mesh1, std::shared_ptr<Mesh> mesh2, VertEdgePair& p1,
                                           VertEdgePair& p2)
{
  // find coincident verts
  auto pv = find_coincident_verts(mesh1, mesh2);
  auto pe = get_corresponding_edges(pv.first, pv.second);

  p1.first  = pv.first;
  p1.second = pe.first;
  p2.first  = pv.second;
  p2.second = pe.second;
}

MEPair MeshBoundarySnapper::get_corresponding_edges(MeshEntityPtr v1, MeshEntityPtr v2)
{
  std::vector<MeshEntityPtr> candidates1, candidates2;
  get_boundary_edges(v1, candidates1);
  get_boundary_edges(v2, candidates2);

  assert(candidates1.size() > 0);
  assert(candidates2.size() > 0);

  MeshEntityPtr edge1 = candidates1[0];
  // std::cout << "edge1 = " << edge1 << std::endl;

  // get the edge that is closest to parallel with edge1
  MeshEntityPtr otherVert = get_other_vert(v1, edge1);
  utils::Point b1         = otherVert->get_point_orig(0) - v1->get_point_orig(0);
  b1                      = b1 / std::sqrt(dot(b1, b1)); // normalize

  double minDev         = std::numeric_limits<double>::max();
  MeshEntityPtr edgeMin = nullptr;
  for (auto& edge2 : candidates2)
  {
    // std::cout << "edge2 = " << edge2 << std::endl;
    otherVert       = get_other_vert(v2, edge2);
    utils::Point b2 = otherVert->get_point_orig(0) - v2->get_point_orig(0);
    b2              = b2 / std::sqrt(dot(b2, b2)); // normalize

    double val = dot(b1, b2); // cos(theta)
    // std::cout << "val = " << val << std::endl;
    if (std::abs(val - 1) < minDev)
    {
      // std::cout << "found minimum" << std::endl;
      minDev  = std::abs(val - 1);
      edgeMin = edge2;
    }
  }

  assert(edgeMin);
  return std::make_pair(edge1, edgeMin);
}

void MeshBoundarySnapper::snap_verts(VertEdgePair p1, VertEdgePair p2, double& maxSnapDist)
{
  auto p1In = p1;
  // MeshEntityPtr v1_next = get_other_vert(p1.first, p1.second);
  utils::Point nextVert = get_other_vert(p1.first, p1.second)->get_point_orig(0);
  if (M_OUTPUT)
  {
    std::cout << "p1.first = " << p1.first << std::endl;
    std::cout << "p1.second = " << p1.second << std::endl;
    std::cout << "v1_next = " << nextVert << std::endl;
    std::cout << "first mesh2 vert " << p2.first->get_id() << std::endl;
  }

  // TODO: in parallel, we the local part may not have the entire boundary
  //       curve
  do {
    // because we choose the further away point for the next target point,
    // it is always guaranteed that we can advance the starting point by 1.
    // This also helps deal with non-monotonicity causing the algorithm
    // to get stuck
    VertEdgePair p1SearchStart = get_next_pair(p1);
    VertEdgePair p2SearchStart = get_next_pair(p2);
    VertEdgePair p2Next        = traverse_next(p2SearchStart, nextVert);
    VertEdgePair p1Next        = traverse_next(p1SearchStart, p2Next.first->get_point_orig(0));

    if (would_create_zero_length_edge(p2Next.first, p1Next.first->get_point_orig(0)))
    {
      if (M_OUTPUT)
        std::cout << "skipping" << std::endl;
      nextVert = get_next_point(p1, p1Next, p2, p2Next);
    } else
    {
      // snap mesh2 vert to mesh1
      if (M_OUTPUT)
      {
        std::cout << "snapping boundary vert " << p2Next.first << " to " << p1Next.first << std::endl;
        std::cout << "snapping mesh2 vertex " << p2Next.first->get_id() << " to mesh1 vertex " << p1Next.first->get_id()
                  << std::endl;
      }

      maxSnapDist = std::max(maxSnapDist, std::sqrt(compute_dist(p1Next.first, p2Next.first->get_point_orig(0))));
      p2Next.first->set_point_orig(0, p1Next.first->get_point_orig(0));

      if (p1.first != p1Next.first)
        snap_intermediate_verts(p1, p1Next, true, maxSnapDist);

      if (p2.first != p2Next.first)
        snap_intermediate_verts(p2, p2Next, false, maxSnapDist);

      p1 = p1Next;
      p2 = p2Next;
      // v1_next = get_other_vert(p1.first, p1.second);
      nextVert = get_next_point(p1Next, p2Next);
      // p1 = traverse_next(p1, v1_next->get_point_orig(0));
    }

    if (M_OUTPUT)
      std::cout << "next_vert = " << nextVert << std::endl;
  } while (p1.first != p1In.first);
}

void MeshBoundarySnapper::check_zero_length_edges(VertEdgePair p)
{
  // std::cout << "\nChecking for zero length edges" << std::endl;
  double len = get_edge_length(p.second);
  // std::cout << "vert = " << p.first->get_id() << ", has coords " << p.first->get_point_orig(0) << std::endl;
  // std::cout << "edge has vert ids " << p.second->get_down(0)->get_id() << ", " << p.second->get_down(1)->get_id() <<
  // std::endl; std::cout << "len = " << len << std::endl;
  if (len < 1e-13)
    throw std::runtime_error("found zero length edge");

  MeshEntityPtr startingVert = p.first;
  // MeshEntityPtr v2 = get_other_vert(p.first, p.second);
  // std::cout << "other vert id = " << v2->get_id() << std::endl;
  p = get_next_pair(p);
  while (p.first != startingVert)
  {
    len = get_edge_length(p.second);

    // std::cout << "vert = " << p.first->get_id() << ", has coords " << p.first->get_point_orig(0) << std::endl;
    // std::cout << "edge has vert ids " << p.second->get_down(0)->get_id() << ", " << p.second->get_down(1)->get_id()
    // << std::endl; std::cout << "len = " << len << std::endl;
    if (len < 1e-13)
      throw std::runtime_error("found zero length edge");

    // v2 = get_other_vert(p.first, p.second);
    // std::cout << "other vert id = " << v2->get_id() << std::endl;
    // p = traverse_next(p, v2->get_point_orig(0));
    p = get_next_pair(p);
  }
}

utils::Point MeshBoundarySnapper::get_next_point(VertEdgePair p1Start, VertEdgePair p1Current, VertEdgePair /*p2Start*/,
                                                 VertEdgePair /*p2Current*/)
{
  utils::Point p1Beginning = p1Start.first->get_point_orig(0);
  utils::Point p1End       = get_other_vert(p1Current.first, p1Current.second)->get_point_orig(0);
  utils::Point disp1       = p1End - p1Beginning;
  double dist1             = dot(disp1, disp1);

  utils::Point p2Beginning = p1Start.first->get_point_orig(0);
  utils::Point p2End       = get_other_vert(p1Current.first, p1Current.second)->get_point_orig(0);
  utils::Point disp2       = p2End - p2Beginning;
  double dist2             = dot(disp2, disp2);

  return dist1 > dist2 ? p1End : p2End;
}

utils::Point MeshBoundarySnapper::get_next_point(VertEdgePair p1, VertEdgePair p2)
{
  double edge1Len = get_edge_length(p1.second);
  double edge2Len = get_edge_length(p2.second);

  // MeshEntityPtr v1 = get_other_vert(p1.first, p1.second);
  // MeshEntityPtr v2 = get_other_vert(p2.first, p2.second);
  // std::cout << "mesh1 candidate = " << v1->get_point_orig(0) << ", id = " << v1->get_id() << std::endl;
  // std::cout << "mesh2 candidate = " << v2->get_point_orig(0) << ", id = " << v2->get_id() << std::endl;

  // std::cout << "edge1_len = " << edge1_len << ", edge2_len = " << edge2_len << std::endl;

  MeshEntityPtr v = edge1Len > edge2Len ? get_other_vert(p1.first, p1.second) : get_other_vert(p2.first, p2.second);
  return v->get_point_orig(0);
}

double MeshBoundarySnapper::get_edge_length(MeshEntityPtr edge)
{
  assert(edge->get_type() == stk::middle_mesh::mesh::MeshEntityType::Edge);
  utils::Point v1   = edge->get_down(0)->get_point_orig(0);
  utils::Point v2   = edge->get_down(1)->get_point_orig(0);
  utils::Point disp = v2 - v1;
  return std::sqrt(dot(disp, disp));
}

bool MeshBoundarySnapper::would_create_zero_length_edge(MeshEntityPtr vert, const utils::Point& newPt)
{
  assert(vert->get_type() == stk::middle_mesh::mesh::MeshEntityType::Vertex);

  for (int i = 0; i < vert->count_up(); ++i)
  {
    MeshEntityPtr edge = vert->get_up(i);
    if (is_boundary_edge(edge))
    {
      utils::Point otherPt = get_other_vert(vert, edge)->get_point_orig(0);
      utils::Point disp    = newPt - otherPt;
      double dist          = std::sqrt(dot(disp, disp));
      if (dist < 1e-13)
        return true;
    }
  }

  return false;
}

// given a range of vertices, snaps all the vertices bewteen them
// to be on the line
void MeshBoundarySnapper::snap_intermediate_verts(VertEdgePair pStart, VertEdgePair pEnd, bool isMesh1,
                                                  double& maxSnapDist)
{
  if (M_OUTPUT)
    std::cout << "\nEntered snapIntermediateVerts" << std::endl;
  utils::Point pt1 = pStart.first->get_point_orig(0);
  utils::Point pt2 = pEnd.first->get_point_orig(0);

  if (M_OUTPUT)
  {
    std::cout << "pt1 = " << pt1 << std::endl;
    std::cout << "pt2 = " << pt2 << std::endl;
  }

  isMesh1 ? mark_edge1_seen(pStart.second) : mark_edge2_seen(pStart.second);
  VertEdgePair pI = get_next_pair(pStart);
  while (pI.first != pEnd.first)
  {
    auto ptNew = compute_closest_point(pt1, pt2, pI.first->get_point_orig(0));
    if (M_OUTPUT)
      std::cout << "snapping intermediate vert " << pI.first << " to " << ptNew << std::endl;

    if (!isMesh1 && M_OUTPUT)
      std::cout << "snapping intermediate mesh2 vertex " << pI.first->get_id() << std::endl;
    maxSnapDist = std::max(maxSnapDist, std::sqrt(compute_dist(pI.first, ptNew)));
    pI.first->set_point_orig(0, ptNew);
    isMesh1 ? mark_edge1_seen(pI.second) : mark_edge2_seen(pI.second);

    pI = get_next_pair(pI);
  }
}

VertEdgePair MeshBoundarySnapper::traverse_next(VertEdgePair p, const utils::Point& pt)
{
  // std::cout << "advancing from pt " << p.first->get_point_orig(0) << " to closest point to " << pt << std::endl;
  double distPrev = compute_dist(p.first, pt);
  double distI    = 0;
  while (true) // ick
  {
    VertEdgePair pNext = get_next_pair(p);
    distI              = compute_dist(pNext.first, pt);
    // std::cout << "candidate point " << p_next.first->get_point_orig(0) << " has distance " << dist_i << std::endl;

    if (distI > distPrev)
    {
      // std::cout << "returning" << std::endl;
      return p;
    }

    p        = pNext;
    distPrev = distI;
  }
}

// gets the boundary edge adjacent to the "other" vert of the
// VertEdgePair
VertEdgePair MeshBoundarySnapper::get_next_pair(VertEdgePair p)
{
  MeshEntityPtr vOther = get_other_vert(p.first, p.second);
  for (int i = 0; i < vOther->count_up(); ++i)
  {
    auto edge = vOther->get_up(i);

    if (is_boundary_edge(edge) && edge != p.second)
      return std::make_pair(vOther, edge);
  }

  throw std::runtime_error("could not find next edge");
}

// compute distance squared between two verts
double MeshBoundarySnapper::compute_dist(MeshEntityPtr v1, const utils::Point& pt2)
{
  auto pt1 = v1->get_point_orig(0);
  auto d   = pt2 - pt1;
  return dot(d, d);
}

// given the line from pt1 to pt2, computes the closest point to pt3 that
// is on the line
utils::Point MeshBoundarySnapper::compute_closest_point(const utils::Point& pt1, const utils::Point& pt2,
                                                        const utils::Point& pt3)
{
  /*
  double num = dot(pt1 - pt3, pt2 - pt1);
  double den = dot(pt2 - pt1, pt2 - pt1);
  double xi = -num/den;

  double eps = 1e-13;
  assert(xi > -m_eps && xi - 1 < m_eps);
  */
  auto b1    = pt3 - pt1;
  auto b2    = pt2 - pt1;
  auto b2Mag = std::sqrt(dot(b2, b2));
  b2         = b2 / b2Mag;

  auto b1Parallel = dot(b1, b2) * b2;

  return b1Parallel + pt1;

  // return (pt2 - pt1)*xi + pt1;
}

MEPair MeshBoundarySnapper::find_coincident_verts(std::shared_ptr<Mesh> mesh1, std::shared_ptr<Mesh> mesh2)
{
  double eps     = 1e-12;
  double minDist = std::numeric_limits<double>::max();
  for (auto& edge1 : mesh1->get_edges())
    if (edge1 && is_boundary_edge(edge1) && !((*m_seenEdges1)(edge1, 0, 0)))
      for (auto& edge2 : mesh2->get_edges())
        if (edge2 && is_boundary_edge(edge2) && !(*m_seenEdges2)(edge2, 0, 0))
        {
          // std::cout << "checking edges " << edge1 << " and " << edge2 << std::endl;
          auto v1   = edge1->get_down(0);
          auto v2   = edge1->get_down(1);
          auto v3   = edge2->get_down(0);
          auto v4   = edge2->get_down(1);
          double d1 = std::sqrt(compute_dist(v1, v3->get_point_orig(0)));
          double d2 = std::sqrt(compute_dist(v1, v4->get_point_orig(0)));
          double d3 = std::sqrt(compute_dist(v2, v3->get_point_orig(0)));
          double d4 = std::sqrt(compute_dist(v2, v4->get_point_orig(0)));

          // std::cout << "d1 = " << d1 << std::endl;
          // std::cout << "d2 = " << d2 << std::endl;
          // std::cout << "d3 = " << d3 << std::endl;
          // std::cout << "d4 = " << d4 << std::endl;
          //
          // std::cout << "v1 = " << v1 << std::endl;
          // std::cout << "v2 = " << v2 << std::endl;
          // std::cout << "v3 = " << v3 << std::endl;
          // std::cout << "v4 = " << v4 << std::endl;
          // std::cout << "v1_orig = " << v1->get_point_orig(0) << std::endl;
          // std::cout << "v3_orig = " << v3->get_point_orig(0) << std::endl;

          if (d1 < eps)
          {
            // std::cout << "case 1" << std::endl;
            return std::make_pair(v1, v3);
          } else if (d2 < eps)
          {
            // std::cout << "case 2" << std::endl;
            return std::make_pair(v1, v4);
          } else if (d3 < eps)
          {
            // std::cout << "case 3" << std::endl;
            return std::make_pair(v2, v3);
          } else if (d4 < eps)
          {
            // std::cout << "case 4" << std::endl;
            return std::make_pair(v2, v4);
          }

          minDist = std::min(minDist, d1);
          minDist = std::min(minDist, d2);
          minDist = std::min(minDist, d3);
          minDist = std::min(minDist, d4);
        }

  throw std::runtime_error("could not find coincident verts, closest found = " + std::to_string(minDist));
}

void MeshBoundarySnapper::get_boundary_edges(MeshEntityPtr v1, std::vector<MeshEntityPtr>& candidates)
{
  candidates.clear();
  for (int i = 0; i < v1->count_up(); ++i)
  {
    auto edge = v1->get_up(i);
    if (is_boundary_edge(edge))
      candidates.push_back(edge);
  }
}

bool MeshBoundarySnapper::is_boundary_edge(MeshEntityPtr edge)
{
  return (edge && ((edge->count_up() == 1 && edge->count_remote_shared_entities() == 0) || edge->count_up() == 0));
}

int MeshBoundarySnapper::count_boundary_edges(std::shared_ptr<Mesh> mesh)
{
  int count = 0;
  for (auto& edge : mesh->get_edges())
    if (is_boundary_edge(edge))
      count++;
  return count;
}

void MeshBoundarySnapper::mark_edge1_seen(MeshEntityPtr edge)
{
  assert(!(*m_seenEdges1)(edge, 0, 0));
  (*m_seenEdges1)(edge, 0, 0) = true;
  m_nedgesSeen1++;
}

void MeshBoundarySnapper::mark_edge2_seen(MeshEntityPtr edge)
{
  assert(!(*m_seenEdges2)(edge, 0, 0));
  (*m_seenEdges2)(edge, 0, 0) = true;
  m_nedgesSeen2++;
}

} // namespace impl
} // namespace mesh
} // namespace middle_mesh
} // namespace stk
