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

#ifndef STK_SEARCH_UTIL_OBJECT_COORDINATES_HPP
#define STK_SEARCH_UTIL_OBJECT_COORDINATES_HPP

#include <map>
#include <ostream>                                   // for operator<<, etc
#include <stk_mesh/base/Types.hpp>                   // for EntityRank
#include <string>                                    // for allocator, etc
#include <utility>                                   // for pair, make_pair
#include <vector>                                    // for vector

#include "stk_mesh/base/MetaData.hpp"                // for MetaData
#include "stk_mesh/base/BulkData.hpp"                // for BulkData
#include "stk_mesh/base/Entity.hpp"                  // for Entity
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_topology/topology.hpp"                 // for topology, etc
#include "stk_util/util/ReportHandler.hpp"           // for ThrowRequireMsg
#include "stk_search_util/MasterElementProvider.hpp"

namespace stk {
namespace search {

void determine_centroid(const unsigned spatial_dim, stk::mesh::Entity element,
                        const stk::mesh::FieldBase& nodal_coord, double* real_loc);

void determine_centroid(const unsigned spatial_dim, stk::mesh::Entity element,
                        const stk::mesh::FieldBase& nodal_coord, std::vector<double>& real_loc);

void determine_gauss_points(const stk::mesh::BulkData& recvBulk, stk::mesh::Entity element,
                            const stk::search::MasterElementProviderInterface& masterElemProvider,
                            const stk::mesh::FieldBase& coordinateField, std::vector<double>& location);

double distance_from_nearest_entity_node(const stk::mesh::BulkData& bulk, stk::mesh::Entity entity,
                                         const stk::mesh::FieldBase* coord_field, const double* point);

double distance_from_nearest_entity_node(const stk::mesh::BulkData& bulk, stk::mesh::Entity entity,
                                         const stk::mesh::FieldBase* coord_field, const std::vector<double>& point);


double distance_squared_from_nearest_entity_node(const stk::mesh::BulkData& bulk, stk::mesh::Entity entity,
                                                 const stk::mesh::FieldBase* coord_field, const double* point);

double distance_squared_from_nearest_entity_node(const stk::mesh::BulkData& bulk, stk::mesh::Entity entity,
                                                 const stk::mesh::FieldBase* coord_field, const std::vector<double>& point);

} // end namespace search
} // end namespace stk

#endif

