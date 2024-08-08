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

#ifndef GENERATE_ALEFRA_MESH_HPP
#define GENERATE_ALEFRA_MESH_HPP

#include <stk_util/stk_config.h>
#include <string>

namespace stk {

namespace mesh {
  class BulkData;
  class Part;
}

namespace unit_test_util {

enum SidesetDirection {LEFT = 0, RIGHT, DOUBLE};
enum ElementOrdering {INCREASING, DECREASING};

void create_AA_mesh(stk::mesh::BulkData &bulk, ElementOrdering elemOrdering = INCREASING);
void create_AB_mesh(stk::mesh::BulkData &bulk, ElementOrdering elemOrdering = INCREASING);

stk::mesh::Part* create_AA_mesh_with_sideset(stk::mesh::BulkData &bulk, SidesetDirection direction,
                                             ElementOrdering elemOrdering = INCREASING);
stk::mesh::Part* create_AB_mesh_with_sideset(stk::mesh::BulkData &bulk, SidesetDirection direction,
                                             ElementOrdering elemOrdering = INCREASING);

stk::mesh::Part* create_AA_mesh_with_sideset_and_field(stk::mesh::BulkData &bulk, SidesetDirection direction,
                                                       ElementOrdering elemOrdering, const std::string & fieldName);
stk::mesh::Part* create_AB_mesh_with_sideset_and_field(stk::mesh::BulkData &bulk, SidesetDirection direction,
                                                       ElementOrdering elemOrdering, const std::string & fieldName);

stk::mesh::Part* create_AB_mesh_with_sideset_and_distribution_factors(stk::mesh::BulkData &bulk,
                                                                      SidesetDirection direction,
                                                                      ElementOrdering elemOrdering,
                                                                      const std::string & fieldName,
                                                                      const double initValue);

namespace simple_fields {

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void create_AA_mesh(stk::mesh::BulkData &bulk, ElementOrdering elemOrdering = INCREASING);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void create_AB_mesh(stk::mesh::BulkData &bulk, ElementOrdering elemOrdering = INCREASING);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
stk::mesh::Part* create_AA_mesh_with_sideset(stk::mesh::BulkData &bulk, SidesetDirection direction,
                                             ElementOrdering elemOrdering = INCREASING);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
stk::mesh::Part* create_AB_mesh_with_sideset(stk::mesh::BulkData &bulk, SidesetDirection direction,
                                             ElementOrdering elemOrdering = INCREASING);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
stk::mesh::Part* create_AA_mesh_with_sideset_and_field(stk::mesh::BulkData &bulk, SidesetDirection direction,
                                                       ElementOrdering elemOrdering, const std::string & fieldName);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
stk::mesh::Part* create_AB_mesh_with_sideset_and_field(stk::mesh::BulkData &bulk, SidesetDirection direction,
                                                       ElementOrdering elemOrdering, const std::string & fieldName);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
stk::mesh::Part* create_AB_mesh_with_sideset_and_distribution_factors(stk::mesh::BulkData &bulk,
                                                                      SidesetDirection direction,
                                                                      ElementOrdering elemOrdering,
                                                                      const std::string & fieldName,
                                                                      const double initValue);

} // namespace simple_fields

}
}

#endif
