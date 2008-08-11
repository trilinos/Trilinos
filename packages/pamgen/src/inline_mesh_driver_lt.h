#ifndef inline_mesh_driver_LTH
#define inline_mesh_driver_LTH
/* class Inline_Mesh_Desc; */
namespace ms_lt {
  class Mesh_Specification;
}
ms_lt::Mesh_Specification * buildMeshSpecification_LT(PAMGEN_NEVADA::Inline_Mesh_Desc *,int rank, int num_procs);

#endif //inline_mesh_driver_LTH
