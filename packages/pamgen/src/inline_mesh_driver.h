// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef inline_mesh_driverH
#define inline_mesh_driverH
/* class Inline_Mesh_Desc; */
namespace ms_rw {
  class Mesh_Specification;
}
ms_rw::Mesh_Specification * buildMeshSpecification(Inline_Mesh_Desc *,long long rank, long long num_procs);

#endif //inline_mesh_driverH
