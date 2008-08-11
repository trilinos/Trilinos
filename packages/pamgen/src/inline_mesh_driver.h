#ifndef inline_mesh_driverH
#define inline_mesh_driverH
/* class Inline_Mesh_Desc; */
namespace ms_rw {
  class Mesh_Specification;
}
ms_rw::Mesh_Specification * buildMeshSpecification(Inline_Mesh_Desc *,int rank, int num_procs);

#endif //inline_mesh_driverH
