// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef GETNGPMESH_HPP
#define GETNGPMESH_HPP

#include "stk_mesh/base/Ngp.hpp"
#include "stk_mesh/base/NgpMesh.hpp"
#include "stk_mesh/base/BulkData.hpp"

namespace stk {
namespace mesh {

template<typename NgpMemSpace = NgpMeshDefaultMemSpace>
inline NgpMeshT<NgpMemSpace> & get_updated_ngp_mesh(const BulkData & bulk)
{
  static_assert(std::is_same_v<NgpMeshT<NgpMemSpace>,HostMeshT<NgpMemSpace>> ||
                Kokkos::SpaceAccessibility<NgpMemSpace,NgpMeshDefaultMemSpace>::accessible,
                "In a GPU-enabled build, get_updated_ngp_mesh requires a device-accessible memory-space.");

  STK_ThrowRequireMsg(!bulk.in_modifiable_state(), "NgpMesh cannot be updated during a mesh modification.");

  NgpMeshBase * ngpMeshBase = impl::get_ngp_mesh(bulk);

  if (ngpMeshBase == nullptr) {
    ngpMeshBase = new NgpMeshT<NgpMemSpace>(bulk);
    impl::set_ngp_mesh(bulk, ngpMeshBase);
  }
  else {
    ngpMeshBase->update_mesh();
  }
  return dynamic_cast<NgpMeshT<NgpMemSpace>&>(*ngpMeshBase);
}

inline NgpMesh & get_updated_ngp_mesh(const BulkData & bulk)
{
  STK_ThrowRequireMsg(!bulk.in_modifiable_state(), "NgpMesh cannot be updated during a mesh modification.");

  return get_updated_ngp_mesh<NgpMeshDefaultMemSpace>(bulk);
}

}}

#endif // GETNGPMESH_HPP
