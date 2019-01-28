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

#ifndef STK_IO_IOHELPERS_HPP
#define STK_IO_IOHELPERS_HPP
// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <Ioss_Field.h>                            // for Field, etc
#include <Ioss_PropertyManager.h>                  // for PropertyManager
#include <stddef.h>                                // for size_t
#include <Teuchos_RCP.hpp>                         // for RCP::RCP<T>, etc
#include <algorithm>                               // for swap
#include <stk_io/DatabasePurpose.hpp>              // for DatabasePurpose
#include <stk_io/IossBridge.hpp>
#include <stk_io/MeshField.hpp>                    // for MeshField, etc
#include <stk_mesh/base/BulkData.hpp>              // for BulkData
#include <stk_mesh/base/Selector.hpp>              // for Selector
#include <stk_util/parallel/Parallel.hpp>          // for ParallelMachine
#include <stk_util/util/ParameterList.hpp>         // for Type
#include <string>                                  // for string
#include <vector>                                  // for vector
#include "Teuchos_RCPDecl.hpp"                     // for RCP
#include "mpi.h"                                   // for MPI_Comm, etc
#include "stk_mesh/base/Types.hpp"                 // for FieldVector
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssert, etc
namespace Ioss { class Property; }
namespace Ioss { class Region; }
namespace boost { class any; }
namespace stk { namespace io { class InputFile; } }
namespace stk { namespace io { class SidesetUpdater; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################
namespace stk { namespace mesh { class BulkData; } }


namespace stk {
namespace io {

inline bool fieldOrdinalSort(const stk::io::FieldAndName& f1, const stk::io::FieldAndName &f2) {
  return f1.field()->mesh_meta_data_ordinal() < f2.field()->mesh_meta_data_ordinal();
}

template <typename DataType>
void internal_write_global(Teuchos::RCP<Ioss::Region> output_region, const std::string &globalVarName,
                           DataType globalVarData);

template <typename DataType>
void internal_write_global(Teuchos::RCP<Ioss::Region> output_region, const std::string &globalVarName,
                           std::vector<DataType> &globalVarData);

bool internal_has_global(Teuchos::RCP<Ioss::Region> input_region, const std::string &globalVarName);

template <typename DataType>
bool internal_read_global(Teuchos::RCP<Ioss::Region> input_region, const std::string &globalVarName,
                          DataType &globalVarData, Ioss::Field::BasicType iossType,
                          bool abort_if_not_found);

template <typename DataType>
bool internal_read_global(Teuchos::RCP<Ioss::Region> input_region, const std::string &globalVarName,
                          std::vector<DataType> &globalVarData, Ioss::Field::BasicType iossType,
                          bool abort_if_not_found);

void internal_write_parameter(Teuchos::RCP<Ioss::Region> output_region,
                              const std::string &name, const boost::any &any_value,
                              stk::util::ParameterType::Type type);

void write_defined_global_any_fields(Teuchos::RCP<Ioss::Region> region,
                                     std::vector<stk::io::GlobalAnyVariable> &global_any_fields);

bool internal_read_parameter(Teuchos::RCP<Ioss::Region> input_region,
                             const std::string &globalVarName,
                             boost::any &any_value, stk::util::ParameterType::Type type,
                             bool abort_if_not_found);

void internal_add_global(Teuchos::RCP<Ioss::Region> region,
                         const std::string &globalVarName,
                         const std::string &storage,
                         Ioss::Field::BasicType dataType,
                         int copies = 1,
                         Ioss::Field::RoleType role = Ioss::Field::TRANSIENT);

void internal_add_global(Teuchos::RCP<Ioss::Region> region,
                         const std::string &globalVarName,
                         int component_count,
                         Ioss::Field::BasicType dataType,
                         int copies = 1,
                         Ioss::Field::RoleType role = Ioss::Field::TRANSIENT);

size_t get_entities(const stk::mesh::Part &part,
                    const stk::mesh::BulkData &bulk,
                    std::vector<stk::mesh::Entity> &entities,
                    const stk::mesh::Selector *subset_selector);

bool is_skipped_attribute_field(const std::string &name, size_t numAttrFields);

void process_surface_entity_df(const Ioss::SideSet* sset, stk::mesh::BulkData & bulk);

template <typename INT>
void process_node_coords_and_attributes(Ioss::Region &region, stk::mesh::BulkData &bulk);

template <typename INT>
void process_elem_attributes_and_implicit_ids(Ioss::Region &region, stk::mesh::BulkData &bulk, const bool shouldAutoLoadAttributes);

template <typename INT>
void process_nodesets_df(Ioss::Region &region, stk::mesh::BulkData &bulk);

void process_sidesets_df(Ioss::Region &region, stk::mesh::BulkData &bulk);

void internal_fill_output_entities(Ioss::GroupingEntity *io_entity,
                                 stk::mesh::Part *part,
                                 stk::mesh::EntityRank part_type,
                                 OutputParams &params,
                                 std::vector<stk::mesh::Entity> &entities);

void put_field_data(OutputParams &params,
                  stk::mesh::Part &part,
                  stk::mesh::EntityRank part_type,
                  Ioss::GroupingEntity *io_entity,
                  const std::vector<stk::io::FieldAndName> &namedFields,
                  Ioss::Field::RoleType filter_role,
                  const stk::mesh::FieldState *state);

void put_field_data(OutputParams &params,
                  stk::mesh::Part &part,
                  stk::mesh::EntityRank part_type,
                  Ioss::GroupingEntity *io_entity,
                  const std::vector<stk::io::FieldAndName> &namedFields,
                  const stk::mesh::FieldState *state=nullptr);


}
}
#endif
