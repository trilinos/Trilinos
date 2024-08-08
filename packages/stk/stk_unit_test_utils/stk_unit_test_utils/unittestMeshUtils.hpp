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

#ifndef unittestMeshUtils_hpp
#define unittestMeshUtils_hpp

#include <stk_mesh/base/Types.hpp>      // for EntityVector, PartVector
#include <string>                       // for string
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Part; } }

namespace stk
{

namespace unit_test_util
{

void put_mesh_into_part(stk::mesh::BulkData& bulkData, stk::mesh::Part& part);

std::string get_name_of_generated_mesh(int xdim, int ydim, int zdim, const std::string &options);

void move_killed_elements_out_of_parts(stk::mesh::BulkData& bulkData,
                                  const stk::mesh::EntityVector& killedElements,
                                  const stk::mesh::PartVector& removeParts);

struct ElementAndPart {
  stk::mesh::EntityId id;
  std::string partName;
};

void put_elements_into_part(stk::mesh::BulkData& bulkData, const std::vector<ElementAndPart> & entries);


namespace simple_fields {

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void put_mesh_into_part(stk::mesh::BulkData& bulkData, stk::mesh::Part& part);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
std::string get_name_of_generated_mesh(int xdim, int ydim, int zdim, const std::string &options);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void move_killed_elements_out_of_parts(stk::mesh::BulkData& bulkData,
                                  const stk::mesh::EntityVector& killedElements,
                                  const stk::mesh::PartVector& removeParts);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void put_elements_into_part(stk::mesh::BulkData& bulkData, const std::vector<ElementAndPart> & entries);

} // namespace simple_fields

} // namespace unit_test_util
} // namespace stk


#endif // unittestMeshUtils_hpp
