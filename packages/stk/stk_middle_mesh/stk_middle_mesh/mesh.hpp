#ifndef MESH_H
#define MESH_H

#include "field_manager.hpp"
#include "geo_classification.hpp"
#include "mesh_entity.hpp"
#include "stk_util/parallel/Parallel.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {

template <typename T>
class Field;

class Mesh
{
    using TStorage = std::vector<MeshEntityPtr>;

  public:
    ~Mesh();

    MPI_Comm get_comm() { return m_comm; }

    MeshEntityPtr create_vertex(const double x, const double y, const double z = 0);

    MeshEntityPtr create_vertex(const stk::middle_mesh::utils::Point& pt)
    {
      return create_vertex(pt.get_x(), pt.get_y(), pt.get_z());
    }

    MeshEntityPtr create_edge(MeshEntityPtr vert1, MeshEntityPtr vert2);

    // edge1_orient determines the orientation of edge1.  The orientation of
    // the other edges can be computed from this using the ordering convention
    MeshEntityPtr create_triangle(MeshEntityPtr edge1, MeshEntityPtr edge2, MeshEntityPtr edge3,
                                  EntityOrientation edge1Orient = EntityOrientation::Standard);

    // create element from vertices, creating edges as neeeded
    MeshEntityPtr create_triangle_from_verts(MeshEntityPtr vert1, MeshEntityPtr vert2, MeshEntityPtr vert3);

    MeshEntityPtr create_quad(MeshEntityPtr edge1, MeshEntityPtr edge2, MeshEntityPtr edge3, MeshEntityPtr edge4,
                              EntityOrientation edge1Orient = EntityOrientation::Standard);

    MeshEntityPtr create_quad_from_verts(MeshEntityPtr vert1, MeshEntityPtr vert2, MeshEntityPtr vert3,
                                         MeshEntityPtr vert4);

    const impl::GeoClassification& get_geo_class(MeshEntityPtr e);

    void set_geo_class(MeshEntityPtr e, const impl::GeoClassification& classI);

    // create a new vertex at specified location in the edge
    // This leaves the mesh in a topologically invalid state because
    // it does not split faces that use the edge
    // The new vertex is returned.  The new edges are returned in the array
    // new_edges
    MeshEntityPtr split_edge_broken(MeshEntityPtr e, double xi, MeshEntityPtr* edgesNew = nullptr);

    // splits an edge at the specified xi coordinate.  The edge must be bounded
    // by two triangles.  The result is that the edge and the two triangles
    // are deleted, and four new edges and four new triangles are created,
    // If an array is supplied to the faces argument, it will be populated with
    // the newly created faces
    void split_edge(MeshEntityPtr e, double xi, MeshEntityPtr* faces = nullptr);

    // reassigns entity IDs so that they are contiguous, and shrinks the underlying arrays
    // This does not currently invalidate MeshEntity pointers, but it might in the future
    void condense_arrays();

    // get all entities of a given type.  Might contain nullptrs
    const TStorage& get_vertices() const { return m_verts; }

    const TStorage& get_edges() const { return m_edges; }

    const TStorage& get_elements() const { return m_elements; }

    const TStorage& get_mesh_entities(const int dim);

    void attach_field(std::shared_ptr<impl::FieldBase> field) { m_fieldManager.add_field(field); }

    // deletes a mesh face, and all downward adjacent entities that would
    // be orphaned by deleting it.  The mesh remains topologically valid
    void delete_face(MeshEntityPtr face);

  private:
    explicit Mesh(MPI_Comm comm = parallel_machine_world())
      : m_comm(comm)
    {
    }

    void set_geo_classification(std::shared_ptr<Field<impl::GeoClassification>> field);

    void condense_array(std::vector<MeshEntityPtr>& entities);

    std::vector<EntityOrientation> get_edges_orientation(std::vector<MeshEntityPtr> edges,
                                                         EntityOrientation edge1Orient);

    // removes entity from storage.  Does not remove it from all topological
    // data structures
    void p_delete_face(MeshEntityPtr face)
    {
      m_elements[face->get_id()] = nullptr;
      delete face;
    }

    void delete_edge(MeshEntityPtr edge)
    {
      m_edges[edge->get_id()] = nullptr;
      delete edge;
    }

    void delete_vert(MeshEntityPtr vert)
    {
      m_verts[vert->get_id()] = nullptr;
      delete vert;
    }

    bool does_edge_exist(MeshEntityPtr vert1, MeshEntityPtr vert2);

    bool does_face_exist(MeshEntityPtr edge1, MeshEntityPtr edge2, MeshEntityPtr edge3);

    bool does_face_exist(MeshEntityPtr edge1, MeshEntityPtr edge2, MeshEntityPtr edge3, MeshEntityPtr edge4);

    MPI_Comm m_comm;
    TStorage m_verts;
    TStorage m_edges;
    TStorage m_elements;
    impl::FieldManager m_fieldManager;
    std::shared_ptr<Field<impl::GeoClassification>> m_geoClass;

    friend std::shared_ptr<Mesh> make_empty_mesh(MPI_Comm comm);
};

int count_valid(const std::vector<MeshEntityPtr>& entities);

// count number of elements that have angles outside the range
// [theta_min, theta_max], angles in degrees
int check_angles(std::shared_ptr<Mesh> mesh, const double thetaMin, const double thetaMax);

void check_topology(std::shared_ptr<Mesh> mesh, int maxDim=2);

void check_topology_down(const std::vector<MeshEntityPtr>& entities, const std::vector<MeshEntityPtr>& entitiesDown);

void check_topology_up(const std::vector<MeshEntityPtr>& entities, const std::vector<MeshEntityPtr>& entitiesUp,
                       const int minUp);

void check_edge_orientation_parallel(std::shared_ptr<Mesh> mesh);

void check_remotes_symmetric(std::shared_ptr<Mesh> mesh);

void check_remotes_unique(std::shared_ptr<Mesh> mesh);

void check_coordinate_field(std::shared_ptr<Mesh> mesh);

void reverse_edge(MeshEntityPtr edge);

void apply_orientation(EntityOrientation flag, MeshEntityPtr* down, const int n);

int get_vertices(MeshEntityPtr tri, MeshEntityPtr*);

int get_downward(MeshEntityPtr e, int dim, MeshEntityPtr*);

int get_upward(MeshEntityPtr e, int dim, std::vector<MeshEntityPtr>& up);

int get_bridge_adjacent(MeshEntityPtr e, const int viaDim, const int targetDim, std::vector<MeshEntityPtr>& entities);

int get_owner(std::shared_ptr<Mesh> mesh, MeshEntityPtr entity);

int get_local_id(MeshEntityPtr higherDimensionEntity, MeshEntityPtr lowerDimensionEntity);

RemoteSharedEntity get_owner_remote(std::shared_ptr<Mesh> mesh, MeshEntityPtr entity);

bool check_is_entity_owner(std::shared_ptr<Mesh> mesh, MeshEntityPtr entity);

RemoteSharedEntity get_remote_shared_entity(MeshEntityPtr entity, int rank);

MeshEntityPtr get_other_vert(MeshEntityPtr v, MeshEntityPtr edge);

// returns the edge connecting vert1 and vert2, or nullptr if no
// such edge exists
MeshEntityPtr get_common_edge(MeshEntityPtr vert1, MeshEntityPtr vert2);

// uses original coordinates
utils::Point compute_edge_coords_orig(MeshEntityPtr edge, const double xi);

utils::Point compute_edge_coords(const utils::Point& n1, const utils::Point& n2, const double xi);

utils::Point compute_quad_centroid(MeshEntityPtr quad);

utils::Point compute_quad_centroid_3d(MeshEntityPtr quad);

void compute_lagrange_vals(const double xi, double vals[2]);

void compute_lagrange_derivs(const double xi, double derivs[2]);

// input is xi coords, output is x-y coords
utils::Point compute_quad_coords_from_xi(MeshEntityPtr quad, const utils::Point& pt);

utils::Point compute_quad_coords_from_xi(const double xi[2], const utils::Point& pt1, const utils::Point& pt2,
                                         const utils::Point& pt3, const utils::Point& pt4);

utils::Point compute_quad_coords_from_xi_3d(MeshEntityPtr quad, const utils::Point& ptXi);

utils::Point compute_quad_coords_from_xi_3d(const double xi[2], const utils::Point& pt1, const utils::Point& pt2,
                                            const utils::Point& pt3, const utils::Point& pt4);

utils::Point compute_tri_coords_from_xi(MeshEntityPtr tri, const utils::Point& ptXi);

utils::Point compute_tri_coords_from_xi_3d(MeshEntityPtr tri, const utils::Point& ptXi);

utils::Point compute_tri_centroid(MeshEntityPtr tri);

utils::Point compute_tri_centroid_3d(MeshEntityPtr tri);

utils::Point compute_centroid(MeshEntityPtr el);

utils::Point compute_centroid_3d(MeshEntityPtr el);

utils::Point compute_coords_from_xi(MeshEntityPtr el, const utils::Point& pt);

utils::Point compute_coords_from_xi_3d(MeshEntityPtr el, const utils::Point& pt);

double convert_xi_coords_from_range(double xiStart, double xiEnd, double xi);

utils::Point convert_xi_coords_from_range(double xiStart, double xiEnd, const utils::Point& ptXi);

double convert_xi_coords_to_range(double xiStart, double xiEnd, double xi);

utils::Point convert_xi_coords_to_range(double xiStart, double xiEnd, const utils::Point& ptXi);

template <typename T>
T interpolate_triangle(const T* vertVals, const utils::Point& ptXi)
{
  double xi0 = 1 - ptXi.x - ptXi.y;
  return vertVals[0]*xi0 + vertVals[1] * ptXi.x + vertVals[2] * ptXi.y;
}

template <typename T>
T interpolate_quad(const T* vertVals, const utils::Point& ptXi)
{
  double lagXi[2], lagEta[2];
  compute_lagrange_vals(ptXi.x, lagXi);
  compute_lagrange_vals(ptXi.y, lagEta);

  double l00 = lagXi[0] * lagEta[0],
         l10 = lagXi[1] * lagEta[0],
         l11 = lagXi[1] * lagEta[1],
         l01 = lagXi[0] * lagEta[1];  

  return l00 * vertVals[0] + l10 * vertVals[1] + l11 * vertVals[2] + l01 * vertVals[3];
}

template <typename T>
T interpolate(MeshEntityType type, const T* vertVals, const utils::Point& ptXi)
{
  switch (type)
  {
    case MeshEntityType::Triangle: { return interpolate_triangle(vertVals, ptXi); }
    case MeshEntityType::Quad:     { return interpolate_quad(vertVals, ptXi); }
    default:
      throw std::runtime_error("unsupported MeshEntityType");
  }
}

int compute_angles(MeshEntityPtr el, double angles[MAX_DOWN]);

bool is_unique(std::vector<MeshEntityPtr> entities);

bool is_null(std::vector<MeshEntityPtr> entities);

void check_vertices_null(std::vector<MeshEntityPtr> entities);

int count_entities_of_type(std::shared_ptr<mesh::Mesh> mesh, MeshEntityType type);


std::shared_ptr<Mesh> make_empty_mesh(MPI_Comm comm = parallel_machine_world());

} // namespace mesh

} // namespace middle_mesh
} // namespace stk
#endif
