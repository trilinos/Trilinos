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

#ifndef stk_mesh_SkinBoundary_hpp
#define stk_mesh_SkinBoundary_hpp

#include <ostream>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { class Part; } }

namespace stk { namespace mesh {
/*
 *
 * API interface to skin the boundary.
 *
 */
void create_exposed_block_boundary_sides(BulkData &bulkData, const Selector& blocksToSkin, const stk::mesh::PartVector& partToPutSidesInto);
#ifndef STK_HIDE_DEPRECATED_CODE // Delete after March 1st 2018
STK_DEPRECATED void create_exposed_block_boundary_sides(BulkData &bulkData, const Selector& blocksToSkin, const stk::mesh::PartVector& partToPutSidesInto, const stk::mesh::Selector* air);
#endif // STK_HIDE_DEPRECATED_CODE
void create_exposed_block_boundary_sides(BulkData &bulkData, const Selector& blocksToSkin, const stk::mesh::PartVector& partToPutSidesInto, const stk::mesh::Selector& air);
void create_interior_block_boundary_sides(BulkData&, const Selector& blocksToConsider, const stk::mesh::PartVector& partToPutSidesInto);
void create_all_sides(stk::mesh::BulkData &bulkData, const stk::mesh::Selector &blocksToConsider, const stk::mesh::PartVector& partToPutSidesInto, bool connect_faces_to_edges);

bool check_exposed_block_boundary_sides(BulkData &bulkData, const Selector& skinnedBlock, stk::mesh::Part& skinnedPart);
bool check_exposed_block_boundary_sides(BulkData &bulkData, const Selector& skinnedBlock, stk::mesh::Part& skinnedPart, std::ostream &stream);
bool check_interior_block_boundary_sides(stk::mesh::BulkData &bulkData, const stk::mesh::Selector &skinnedBlock, stk::mesh::Part &skinnedPart);
bool check_interior_block_boundary_sides(stk::mesh::BulkData &bulkData, const stk::mesh::Selector &skinnedBlock, stk::mesh::Part &skinnedPart, std::ostream &stream);
bool check_all_sides(stk::mesh::BulkData &bulkData, const stk::mesh::Selector &skinnedBlock, stk::mesh::Part& skinnedPart);
bool check_all_sides(stk::mesh::BulkData &bulkData, const stk::mesh::Selector &skinnedBlock, stk::mesh::Part& skinnedPart, std::ostream &stream);

}} // namespace stk::mesh
#endif
