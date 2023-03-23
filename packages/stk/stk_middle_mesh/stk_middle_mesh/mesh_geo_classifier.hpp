#ifndef MESH_GEO_CLASSIFIER_H
#define MESH_GEO_CLASSIFIER_H

#include "mesh.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

// figure out a geometric classification for the mesh
// this only works for all quad meshes

class MeshGeoClassifier
{
  public:
    void classify(std::shared_ptr<Mesh> mesh);

  private:
    // classify all mesh entities on geometric face 0
    void set_all_interior(std::shared_ptr<Mesh> mesh);

    // sets the dimension (but not id) of boundary edges
    void set_boundary_edge_dim(std::shared_ptr<Mesh> mesh);

    // classify vertices on geometric vertices
    // returns vector of mesh verticies classified on geometric vertices
    std::vector<MeshEntityPtr> set_vertex_verts(std::shared_ptr<Mesh> mesh);

    // assign geometric ids to boundary edges
    void set_boundary_edge_ids(std::shared_ptr<Mesh> mesh, std::vector<MeshEntityPtr>& verts);

    // gets a vertex at which to start the edge iteration
    MeshEntityPtr get_start_vertex(std::shared_ptr<Mesh> mesh, std::vector<MeshEntityPtr>& verts);

    // returns a mesh edge that is:
    //   1. adjacent to v
    //   2. is classified on a geometric edge
    //   3. has geometric edge id greater than iter
    //   4. is closest parallel with the vector <1, 0, 0>
    MeshEntityPtr get_next_edge(std::shared_ptr<Mesh> mesh, MeshEntityPtr v, const int iter);

    // classifies verts on geometric edges
    void set_boundary_vert_ids(std::shared_ptr<Mesh> mesh);
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
