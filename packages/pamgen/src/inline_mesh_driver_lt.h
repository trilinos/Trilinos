// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef inline_mesh_driver_LTH
#define inline_mesh_driver_LTH
/* class Inline_Mesh_Desc; */
namespace ms_lt {
  class Mesh_Specification;
}
ms_lt::Mesh_Specification * buildMeshSpecification_LT(PAMGEN_NEVADA::Inline_Mesh_Desc *,long long rank, long long num_procs);
ms_lt::Mesh_Specification * consolidateMeshSpecification_LT(ms_lt::Mesh_Specification *);

#endif //inline_mesh_driver_LTH
