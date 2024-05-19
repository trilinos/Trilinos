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

#ifndef STKIOUTILS_HPP_
#define STKIOUTILS_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <cstddef>                          // for size_t
#include <map>                              // for map, map<>::value_compare
#include <string>                           // for string
#include <utility>                          // for pair
#include <vector>                           // for vector
#include "Ioss_Field.h"                     // for Field, Field::BasicType
#include "stk_mesh/base/Types.hpp"          // for EntityRank, EntityVector
#include "stk_util/util/ParameterList.hpp"  // for Type
namespace stk { namespace io { class StkMeshIoBroker; } }
namespace stk { namespace io { struct OutputParams; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Selector; } }

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace io {

namespace impl {

stk::mesh::Selector internal_build_selector(const stk::mesh::Selector *subset_selector,
                                            const stk::mesh::Selector *output_selector,
                                            const stk::mesh::Selector *shared_selector,
                                            const stk::mesh::Part &part,
                                            bool include_shared);

} // namespace impl


size_t get_entities(OutputParams &params,
                    const stk::mesh::Part &part,
                    stk::mesh::EntityRank type,
                    stk::mesh::EntityVector &entities,
                    bool include_shared);

size_t get_entities_for_nodeblock(OutputParams &params,
                    const stk::mesh::Part &part,
                    stk::mesh::EntityRank type,
                    stk::mesh::EntityVector &entities,
                    bool include_shared);

stk::mesh::EntityRank part_primary_entity_rank(const stk::mesh::Part &part);

typedef std::map<std::string, std::vector<std::string>> IossBlockMembership;

IossBlockMembership get_block_memberships(stk::io::StkMeshIoBroker& stkIo);

void fill_block_parts_given_names(const std::vector<std::string>& side_block_names,
                                              stk::mesh::MetaData& meta,
                                              std::vector<const stk::mesh::Part*>& blocks);

void throw_if_any_elem_block_has_invalid_topology(const stk::mesh::MetaData& meta,
                                                  const std::string& msgRegionName);

stk::mesh::FieldVector get_fields_with_role(const stk::mesh::MetaData &meta, const Ioss::Field::RoleType role);
stk::mesh::FieldVector get_fields_with_role(const stk::mesh::MetaData &meta, const stk::mesh::EntityRank rank,
                                            const Ioss::Field::RoleType role);

stk::mesh::FieldVector get_transient_fields(const stk::mesh::MetaData &meta);
stk::mesh::FieldVector get_transient_fields(const stk::mesh::MetaData &meta, const stk::mesh::EntityRank rank);

template<typename DATA_TYPE>
void write_global_to_stk_io(stk::io::StkMeshIoBroker& stkIo, size_t dbIndex,
                            const std::string& externalName,
                            size_t component_count, const void* ptr);

std::pair<size_t, stk::util::ParameterType::Type>
get_parameter_type_from_storage(const std::string &storage,
                                stk::util::ParameterType::Type scalar,
                                stk::util::ParameterType::Type vector);

std::pair<size_t, stk::util::ParameterType::Type>
get_parameter_type_from_field_representation(const std::string &storage,
                                             Ioss::Field::BasicType dataType,
                                             int copies = 1);

std::pair<size_t, Ioss::Field::BasicType>
get_io_parameter_size_and_type(const stk::util::ParameterType::Type type,
                               const std::any &value);

void superset_mesh_parts(const stk::mesh::Part& part, stk::mesh::PartVector& supersetParts);

stk::mesh::Selector construct_sideset_selector(stk::io::OutputParams &params);

std::string construct_parallel_filename(const std::string &baseFilename, int numSubdomains, int subdomainIndex);

std::string construct_filename_for_serial_or_parallel(const std::string &baseFilename, int numSubdomains, int subdomainIndex);

}}


#endif /* STKIOUTILS_HPP_ */
