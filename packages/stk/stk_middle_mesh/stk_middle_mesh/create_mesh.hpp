#ifndef CREATE_MESH_H
#define CREATE_MESH_H

#include "mesh.hpp"
#include "utils.hpp"
#include <iostream>

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

struct MeshSpec
{
    int numelX     = -1;
    int numelY     = -1;
    double xmin    = 0;
    double xmax    = 0;
    double ymin    = 0;
    double ymax    = 0;
    bool xPeriodic = false; // true if the verts at y = ymax should be the
                            // same as y = ymin
    bool yPeriodic = false; // true if the verts at x = xmax should be the
                            // same as x = xmin
};

inline std::ostream& operator<<(std::ostream& os, const MeshSpec& spec)
{
  os << "MeshSpec: grid of " << spec.numelX << " by " << spec.numelY << " elements" << std::endl;
  os << "          domain [" << spec.xmin << ", " << spec.xmax << "] by [" << spec.ymin << ", " << spec.ymax << "]"
     << std::endl;
  os << "          x_periodic: " << spec.xPeriodic << ", yperiodic = " << spec.yPeriodic << std::endl;

  return os;
}

class EdgeNeighborData;

class MeshGenerator
{
  public:
    MeshGenerator(MPI_Comm comm)
      : m_comm(comm)
      , m_mesh(make_empty_mesh(comm))
    {
      MPI_Comm_rank(comm, &m_myrank);
    }

    template <typename Func>
    std::shared_ptr<Mesh> create_mesh(const MeshSpec& spec, Func func, const bool triangles = false,
                                      bool reverseXEdges = false)
    {
      assert(spec.numelX > 0);
      assert(spec.numelY > 0);
      m_meshspecGlobal     = spec;
      m_processorGridShape = get_process_grid_shape(utils::impl::comm_size(m_comm));
      m_processIndices     = compute_local_indices(m_processorGridShape, m_myrank);
      m_meshspecLocal      = compute_local_mesh_spec();

      auto vertCoords = compute_vertex_coordinates(m_meshspecLocal, func);

      auto verts  = create_vertices(m_meshspecLocal, vertCoords);
      Edges edges = create_edges(m_meshspecLocal, verts, triangles, reverseXEdges);
      create_mesh_elements(m_meshspecLocal, edges, triangles);
      set_remote_shared_entities(m_mesh, verts);

      return m_mesh;
    }

  private:
    class Edges
    {
      public:
        Edges(const MeshSpec& spec)
          : m_spec(spec)
          , m_edgesX(spec.numelX * (spec.numelY + 1), nullptr)
          , m_edgesOrientationX(spec.numelX * (spec.numelY + 1), EntityOrientation::Standard)
          , m_edgesY((spec.numelX + 1) * spec.numelY, nullptr)
          , m_edgesDiagonal(spec.numelX * spec.numelY, nullptr)
        {}

        MeshEntityPtr& get_x_edge(int i, int j) { return m_edgesX[get_idx_x(i, j)]; }

        EntityOrientation& get_x_edge_orientation(int i, int j) { return m_edgesOrientationX[get_idx_x(i, j)]; }

        MeshEntityPtr& get_y_edge(int i, int j) { return m_edgesY[get_idx_y(i, j)]; }

        MeshEntityPtr& get_diag_edge(int i, int j) { return m_edgesDiagonal[get_idx_diag(i, j)]; }

      private:
        int get_idx_x(int i, int j)
        {
          assert(i >= 0 && i < m_spec.numelX);
          assert(j >= 0 && j < m_spec.numelY + 1);
          return i + j * m_spec.numelX;
        }

        int get_idx_y(int i, int j)
        {
          assert(i >= 0 && i < m_spec.numelX + 1);
          assert(j >= 0 && j < m_spec.numelY);
          return i + j * (m_spec.numelX + 1);
        }

        int get_idx_diag(int i, int j)
        {
          assert(i >= 0 && i < m_spec.numelX);
          assert(j >= 0 && j < m_spec.numelY);
          return i + j * m_spec.numelX;
        }

        MeshSpec m_spec;
        std::vector<MeshEntityPtr> m_edgesX;
        std::vector<EntityOrientation> m_edgesOrientationX;
        std::vector<MeshEntityPtr> m_edgesY;
        std::vector<MeshEntityPtr> m_edgesDiagonal;
    };

    int get_idx(const MeshSpec& spec, const int i, const int j);

    std::array<int, 2> get_process_grid_shape(int commsize);

    std::array<int, 2> compute_local_indices(const std::array<int, 2>& gridShape, int myrank);

    int wrap_processor_index(int idx, int size);

    std::array<int, 2> wrap_processor_indices(const std::array<int, 2>& gridIndices);

    bool does_block_exist(const std::array<int, 2>& gridIndices);

    int get_process_rank(const std::array<int, 2>& gridIndices);

    bool is_x_periodic_local();

    bool is_y_periodic_local();

    MeshSpec compute_local_mesh_spec();

    std::pair<int, double> get_num_el_in_direction(int numelGlobal, double coordMin, double coordMax, int numProcs,
                                                   int procIdx);

    template <typename Func>
    std::vector<utils::Point> compute_vertex_coordinates(const MeshSpec& spec, Func func);

    std::vector<MeshEntityPtr> create_vertices(const MeshSpec& spec, const std::vector<utils::Point>& vertCoords);

    Edges create_edges(const MeshSpec& spec, const std::vector<MeshEntityPtr>& verts, bool triangles,
                       bool reverseXEdges);

    // create mesh of quads from a vector of utils::Points from createVertices
    void create_mesh_elements(const MeshSpec& spec, Edges& edges, const bool triangles = false);

    void set_remote_shared_entities(std::shared_ptr<Mesh> mesh, std::vector<MeshEntityPtr>& verts);

    void set_edge_adjacent_block_shared_entities(std::shared_ptr<Mesh> mesh, std::vector<MeshEntityPtr>& verts);

    void get_y_shared_entities(EdgeNeighborData& leftEdge, EdgeNeighborData& rightEdge,
                               std::vector<MeshEntityPtr>& verts);

    void get_x_shared_entities(EdgeNeighborData& bottomEdge, EdgeNeighborData& topEdge,
                               std::vector<MeshEntityPtr>& verts);

    void set_corner_adjacent_block_shared_entities(std::shared_ptr<Mesh> mesh, std::vector<MeshEntityPtr>& verts);

    MPI_Comm m_comm;
    MeshSpec m_meshspecGlobal;
    int m_myrank;
    std::array<int, 2> m_processorGridShape;
    std::array<int, 2> m_processIndices;
    MeshSpec m_meshspecLocal;
    std::shared_ptr<Mesh> m_mesh;
};

// create vector of utils::Points on a grid, using the maping function
// utils::Point = func(utils::Point);
template <typename Func>
std::vector<utils::Point> MeshGenerator::compute_vertex_coordinates(const MeshSpec& spec, Func func)
{
  assert(spec.numelX > 0);
  assert(spec.numelY > 0);
  std::vector<utils::Point> verts((spec.numelX + 1) * (spec.numelY + 1));

  double dx = (spec.xmax - spec.xmin) / spec.numelX;
  double dy = (spec.ymax - spec.ymin) / spec.numelY;
  for (int i = 0; i < spec.numelX + 1; ++i)
    for (int j = 0; j < spec.numelY + 1; ++j)
    {
      utils::Point pt(spec.xmin + dx * i, spec.ymin + dy * j);
      verts[get_idx(spec, i, j)] = func(pt);
    }

  return verts;
}

template <typename Func>
std::shared_ptr<Mesh> create_mesh(const MeshSpec& spec, Func func, MPI_Comm comm = parallel_machine_world(),
                                  const bool triangles = false, bool reverseXEdges = false)
{
  assert(spec.numelX > 0);
  assert(spec.numelY > 0);

  MeshGenerator gen(comm);
  return gen.create_mesh(spec, func, triangles, reverseXEdges);
}

// create a sphere with radius r_outer.  This function works by creating a quarter annulus
// and then projecting the points up onto the surface of a sphere with radius r_outer.
// r_inner gives the inner radius of the annulus
std::shared_ptr<Mesh> create_eigth_sphere(int numelR, int numelTheta, double rInner, double rOuter, MPI_Comm comm=MPI_COMM_WORLD, bool createTriangles=false); // CHECK: ALLOW MPI_COMM_WORLD

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
