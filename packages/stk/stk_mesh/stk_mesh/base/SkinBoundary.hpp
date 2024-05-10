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
void create_all_block_boundary_sides(BulkData &bulkData, const Selector& blocksToSkin, const PartVector& partToPutSidesInto, const PartVector* separatePartForInteriorSides = nullptr);
void create_exposed_block_boundary_sides(BulkData &bulkData, const Selector& blocksToSkin, const PartVector& partToPutSidesInto);
void create_exposed_block_boundary_sides(BulkData &bulkData, const Selector& blocksToSkin, const PartVector& partToPutSidesInto, const Selector& air);
void create_interior_block_boundary_sides(BulkData&, const Selector& blocksToConsider, const PartVector& partToPutSidesInto);
void create_all_sides(BulkData &bulkData, const Selector &blocksToConsider, const PartVector& partToPutSidesInto=PartVector(), bool connect_faces_to_edges=true);

bool check_exposed_block_boundary_sides(BulkData &bulkData, const Selector& skinnedBlock, Part& skinnedPart);
bool check_exposed_block_boundary_sides(BulkData &bulkData, const Selector& skinnedBlock, Part& skinnedPart, std::ostream &stream);
bool check_interior_block_boundary_sides(BulkData &bulkData, const Selector &skinnedBlock, Part &skinnedPart);
bool check_interior_block_boundary_sides(BulkData &bulkData, const Selector &skinnedBlock, Part &skinnedPart, std::ostream &stream);
bool check_all_sides(BulkData &bulkData, const Selector &skinnedBlock, Part& skinnedPart);
bool check_all_sides(BulkData &bulkData, const Selector &skinnedBlock, Part& skinnedPart, std::ostream &stream);

}} // namespace stk::mesh
#endif
