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

#ifndef STK_STK_IO_STK_IO_INPUTQUERY_HPP_
#define STK_STK_IO_STK_IO_INPUTQUERY_HPP_

#include <stk_mesh/base/Types.hpp>
#include <stk_io/DatabasePurpose.hpp>   // for DatabasePurpose
#include <stk_io/MeshField.hpp>
#include <stk_io/IossBridge.hpp>
#include "Ioss_EntityType.h"

namespace Ioss {
class PropertyManager;
class GroupingEntity;
class Region;
class DatabaseIO;
}

namespace stk {
namespace mesh {
class MetaData;
class BulkData;
class Part;
}

namespace io {
class StkMeshIoBroker;

using MissingFieldMap = std::map<stk::mesh::FieldBase *, const stk::io::MeshField *>;

class InputQuery
{
public:
  InputQuery(const Ioss::Region& region,
             const stk::mesh::MetaData& meta,
             const DatabasePurpose dbPurpose,
             const std::vector<std::string>* multiStateSuffixes = nullptr);

  ~InputQuery() { }

  int build_field_part_associations(stk::io::MeshField& mf,
                                    std::vector<stk::io::MeshField> *missingFields = nullptr,
                                    const bool throwOnErrorMessage = true);

  int build_field_part_associations_from_grouping_entity(stk::io::MeshField& mf,
                                                         std::vector<stk::io::MeshField> *missingFields = nullptr,
                                                         const bool throwOnErrorMessage = true);

  void build_field_part_associations_for_part(stk::io::MeshField &mf, const stk::mesh::Part * part);

  bool process_fields_for_grouping_entity(stk::io::MeshField &mf,
                                          const stk::mesh::Part &part,
                                          Ioss::GroupingEntity *ioEntity,
                                          MissingFieldMap *missingFieldsCollectorPtr = nullptr);

  bool build_field_part_associations(stk::io::MeshField &mesh_field,
                                     const stk::mesh::Part &part,
                                     const stk::mesh::EntityRank rank,
                                     Ioss::GroupingEntity *ioEntity,
                                     MissingFieldMap *missingFields = nullptr);

private:
  const Ioss::Region& m_region;
  const stk::mesh::MetaData& m_meta;
  DatabasePurpose m_dbPurpose;
  const std::vector<std::string>* m_multiStateSuffixes = nullptr;
};

bool verify_field_request(const StkMeshIoBroker &broker,
                          const stk::io::MeshField &meshField,
                          bool printWarning = true);

bool verify_field_request(const Ioss::Region& region,
                          const stk::mesh::MetaData& meta,
                          const DatabasePurpose dbPurpose,
                          const std::vector<std::string>& multiStateSuffixes,
                          const stk::io::MeshField &meshField,
                          bool printWarning = true);

}
}

#endif /* STK_STK_IO_STK_IO_INPUTQUERY_HPP_ */
