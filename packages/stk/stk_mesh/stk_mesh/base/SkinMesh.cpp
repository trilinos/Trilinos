// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
//

#include <cassert>
#include <stk_mesh/base/SkinMesh.hpp>
#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/MetaData.hpp>   // for BulkData, etc
#include <stk_mesh/base/Selector.hpp>   // for Selector, operator&
#include <vector>                       // for vector
#include "stk_mesh/base/Types.hpp"      // for EntityVector, etc
#include "stk_mesh/base/ElemElemGraph.hpp"
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowErrorMsgIf

namespace stk { namespace mesh {

void skin_mesh( BulkData & mesh, PartVector const& skin_parts)
{
  skin_mesh(mesh, mesh.mesh_meta_data().universal_part(), skin_parts);
}

void skin_mesh( BulkData & mesh, Selector const& element_selector, PartVector const& skin_parts,
                const Selector * secondary_selector)
{
    ThrowErrorMsgIf( mesh.in_modifiable_state(), "mesh is not SYNCHRONIZED" );

    Selector *air = nullptr;
    Selector tmp;

    if(secondary_selector != nullptr)
    {
        tmp = !(*secondary_selector);
        air = &tmp;
    }

    stk::mesh::ElemElemGraph elem_elem_graph(mesh, element_selector, air);
    elem_elem_graph.skin_mesh(skin_parts);
}

}} // namespace stk::mesh
