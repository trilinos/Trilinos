#ifndef MESH_SNAPPER_OPTS
#define MESH_SNAPPER_OPTS

namespace stk {
namespace middle_mesh {

struct MeshSnapperOpts
{
    double tol = 1e-8; // relative tolerance for snapping (relative to each element edge length)
};

} // namespace middle_mesh
} // namespace stk
#endif
