#ifndef MESH_BOUNDARY_SNAPPER_H
#define MESH_BOUNDARY_SNAPPER_H

#include "mesh.hpp"
#include "mesh_exchange_boundary_edges.hpp"
#include <unordered_map>

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

class MeshBoundarySnapper
{
  public:
    using VertEdgePair = std::pair<MeshEntityPtr, MeshEntityPtr>;

    void snap(std::shared_ptr<Mesh> mesh1, std::shared_ptr<Mesh> mesh2, MPI_Comm unionComm);

    virtual ~MeshBoundarySnapper(){}; // TODO: why is this here?  There are no derived classes

  protected:
    virtual void get_start_entity(std::shared_ptr<Mesh> mesh1, std::shared_ptr<Mesh> mesh2, VertEdgePair& p1,
                                  VertEdgePair& p2);

    virtual std::pair<MeshEntityPtr, MeshEntityPtr> get_corresponding_edges(MeshEntityPtr v1, MeshEntityPtr v2);

    static constexpr bool M_OUTPUT = false;

  private:
    // Field<bool> doesn't work because std::vector<bool> is non-standard
    using FakeBool = int_least8_t;

    void start_local_snap(std::shared_ptr<Mesh> mesh1, std::shared_ptr<Mesh> mesh2);

    void snap_verts(VertEdgePair p1, VertEdgePair p2, double& maxSnapDist);

    void check_zero_length_edges(VertEdgePair p);

    utils::Point get_next_point(VertEdgePair p1Start, VertEdgePair p1Current, VertEdgePair p2Start,
                                VertEdgePair p2Current);

    utils::Point get_next_point(VertEdgePair p1, VertEdgePair p2);

    double get_edge_length(MeshEntityPtr edge);

    // returns true if moving the vert to the given coordinates would make any
    // boundary edge have zero length
    bool would_create_zero_length_edge(MeshEntityPtr vert, const utils::Point& newPt);

    // given a range of vertices, snaps all the vertices bewteen them
    // to be on the line
    void snap_intermediate_verts(VertEdgePair pStart, VertEdgePair p1End, bool isMesh1, double& maxSnapDist);

    VertEdgePair traverse_next(VertEdgePair p, const utils::Point& pt);

    // gets the boundary edge adjacent to the "other" vert of the
    // VertEdgePair
    VertEdgePair get_next_pair(VertEdgePair p);

    // compute distance squared between two verts
    double compute_dist(MeshEntityPtr v1, const utils::Point& pt2);

    // given the line from pt1 to pt2, computes the closest point to pt3 that
    // is on the line
    utils::Point compute_closest_point(const utils::Point& pt1, const utils::Point& pt2, const utils::Point& pt3);

    std::pair<MeshEntityPtr, MeshEntityPtr> find_coincident_verts(std::shared_ptr<Mesh> mesh1,
                                                                  std::shared_ptr<Mesh> mesh2);

    void get_boundary_edges(MeshEntityPtr v1, std::vector<MeshEntityPtr>& candidates);

    bool is_boundary_edge(MeshEntityPtr edge);

    int count_boundary_edges(std::shared_ptr<Mesh> mesh);

    void mark_edge1_seen(MeshEntityPtr edge);

    void mark_edge2_seen(MeshEntityPtr edge);

    std::shared_ptr<Field<FakeBool>> m_seenEdges1;
    std::shared_ptr<Field<FakeBool>> m_seenEdges2;
    int m_nedgesSeen1 = 0;
    int m_nedgesSeen2 = 0;
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
