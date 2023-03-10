#ifndef MESH_IO_H
#define MESH_IO_H

#include "mesh.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

void print_tin(const std::string& fname, std::shared_ptr<Mesh> mesh);

void print_tin(const std::string& fname, std::vector<MeshEntityPtr>& els);

void print_vert_edges(const std::string& fname, std::shared_ptr<Mesh> meshIn);

void print_vert_edges(const std::string& fname, const std::vector<MeshEntityPtr>& els);

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
