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

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stk_io/InputQuery.hpp>
#include <exception>                    // for exception
#include <algorithm>                           // for copy, sort, max, find
#include <cmath>                               // for fmod
#include <cstddef>                             // for size_t
#include <iostream>                            // for operator<<, basic_ostream
#include <limits>                              // for numeric_limits
#include <stdexcept>                           // for runtime_error
#include <stk_io/DatabasePurpose.hpp>          // for READ_RESTART, Database...
#include <stk_io/DbStepTimeInterval.hpp>       // for DBStepTimeInterval
#include <stk_io/InputFile.hpp>
#include <stk_io/IossBridge.hpp>               // for is_part_io_part, all_f...
#include <stk_io/MeshField.hpp>                // for MeshField, MeshField::...
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/FieldBase.hpp>         // for FieldBase, FieldBase::...
#include <stk_mesh/base/FindRestriction.hpp>   // for find_restriction
#include <stk_mesh/base/MetaData.hpp>          // for MetaData
#include <stk_util/environment/FileUtils.hpp>  // for filename_substitution
#include "stk_util/environment/RuntimeWarning.hpp"  // for RuntimeWarningAdHoc
#include <stk_util/util/ReportHandler.hpp>     // for ThrowErrorMsgIf, Throw...
#include <utility>                             // for move, pair
#include "Ioss_DBUsage.h"                      // for DatabaseUsage, READ_MODEL
#include "Ioss_DatabaseIO.h"                   // for DatabaseIO
#include "Ioss_EntityType.h"                   // for SIDESET, EntityType
#include "Ioss_Field.h"                        // for Field, Field::TRANSIENT
#include "Ioss_GroupingEntity.h"               // for GroupingEntity
#include "Ioss_IOFactory.h"                    // for IOFactory
#include "Ioss_MeshType.h"                     // for MeshType, MeshType::UN...
#include "Ioss_NodeBlock.h"                    // for NodeBlock
#include "Ioss_NodeSet.h"                      // for NodeSet
#include "Ioss_Property.h"                     // for Property
#include "Ioss_Region.h"                       // for Region, NodeBlockConta...
#include "Ioss_SideBlock.h"                    // for SideBlock
#include "Ioss_SideSet.h"                      // for SideSet
#include "StkIoUtils.hpp"                      // for part_primary_entity_rank
#include "stk_mesh/base/BulkData.hpp"          // for BulkData
#include "stk_mesh/base/FieldState.hpp"        // for FieldState
#include "stk_mesh/base/Part.hpp"              // for Part
#include "stk_mesh/base/Types.hpp"             // for PartVector, EntityRank
#include "stk_topology/topology.hpp"           // for topology, topology::NO...
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace {
void add_missing_fields(std::vector<stk::io::MeshField> *missingFields,
                        stk::io::MissingFieldMap& missingFieldsCollector)
{
  if (nullptr != missingFields) {
    std::vector<stk::io::MeshField> discoveredMissingFields;
    for (auto missingStatedFieldIter : missingFieldsCollector)
    {
      discoveredMissingFields.push_back(stk::io::MeshField(missingStatedFieldIter.first,
                                                           missingStatedFieldIter.second->db_name()));
    }
    std::sort(discoveredMissingFields.begin(), discoveredMissingFields.end(),
              [](const stk::io::MeshField &a, const stk::io::MeshField &b) {
                   return (a.db_name() < b.db_name())
                           || ((a.db_name() == b.db_name()) && (a.field()->name() < b.field()->name())); });

    for(stk::io::MeshField &missingField : *missingFields) {
      std::vector<stk::io::MeshField>::iterator iter = std::find(discoveredMissingFields.begin(), discoveredMissingFields.end(), missingField);
      if(iter != discoveredMissingFields.end()) {
        discoveredMissingFields.erase(iter);
      }
    }

    missingFields->insert(missingFields->end(), discoveredMissingFields.begin(), discoveredMissingFields.end());
  }
}
}

namespace stk {
namespace io {

  InputQuery::InputQuery(const Ioss::Region& region,
                         const stk::mesh::MetaData& meta,
                         const DatabasePurpose dbPurpose,
                         const std::vector<std::string>* multiStateSuffixes)
    : m_region(region),
      m_meta(meta),
      m_dbPurpose(dbPurpose),
      m_multiStateSuffixes(multiStateSuffixes)
  {
  }

  bool InputQuery::build_field_part_associations(stk::io::MeshField &meshField,
                                                 const stk::mesh::Part &part,
                                                 const stk::mesh::EntityRank rank,
                                                 Ioss::GroupingEntity *ioEntity,
                                                 MissingFieldMap *missingFieldsCollector)
  {
    bool fieldIsMissing = false;
    stk::mesh::FieldBase *f = meshField.field();
    // Only add TRANSIENT Fields -- check role; if not present assume transient...
    const Ioss::Field::RoleType *role = stk::io::get_field_role(*f);
    if (role == nullptr || *role == Ioss::Field::TRANSIENT) {
      if (stk::io::is_field_on_part(f, rank, part)) {
        const stk::mesh::FieldBase::Restriction &res = stk::mesh::find_restriction(*f, rank, part);
        FieldType fieldType;
        stk::io::get_io_field_type(f, res, &fieldType);

        if (fieldType.type != Ioss::Field::INVALID) {
          const std::string &dbName = meshField.db_name();
          unsigned numStates = f->number_of_states();
          std::vector<stk::mesh::FieldState> missingStates;
          if (numStates > 1) {
            bool hasAllStates = all_field_states_exist_on_io_entity(dbName, f, ioEntity, missingStates, m_multiStateSuffixes);
            if(hasAllStates == false) {
              fieldIsMissing = true;
              if (missingFieldsCollector) {
                for (stk::mesh::FieldState missingState : missingStates)
                  (*missingFieldsCollector)[f->field_state(missingState)] = &meshField;
              }
            }
          }

          bool fieldExists = ioEntity->field_exists(dbName);
          if (!fieldExists) {
            fieldIsMissing = true;
            if (missingFieldsCollector) {
              (*missingFieldsCollector)[f] = &meshField;
            }
          }

          // See if field with that name exists on ioEntity...
          if (fieldExists) {
            meshField.add_part(rank, part, ioEntity);
            meshField.set_single_state((m_dbPurpose == stk::io::READ_RESTART) ? false : true);
            meshField.set_active();
          }
        }
      }
    }
    return fieldIsMissing;
  }

  int InputQuery::build_field_part_associations(stk::io::MeshField& mf,
                                                std::vector<stk::io::MeshField> *missingFields,
                                                const bool throwOnErrorMessage)
  {
    MissingFieldMap missingFieldsCollector;
    MissingFieldMap *missingFieldsCollectorPtr = (missingFields ? &missingFieldsCollector : nullptr);

    // Each input field will have a list of the Parts that the field exists on...
    // Create this list.
    int numMissingFields = 0;
    // First handle any fields that are sub-setted (restricted to a specified list of parts)

    const stk::mesh::FieldBase *f = mf.field();

    for (const stk::mesh::Part *part : mf.m_subsetParts) {
      stk::mesh::EntityRank rank = part_primary_entity_rank(*part);
      bool fieldIsMissing = false;

      if (f->entity_rank() == rank) {
        Ioss::GroupingEntity *ioEntity = m_region.get_entity(part->name());
        STK_ThrowErrorMsgIf( ioEntity == nullptr,
                             "ERROR: For field '" <<
                             mf.field()->name() <<
                             "' Could not find database entity corresponding to the part named '" <<
                             part->name() << "'.");
        fieldIsMissing |= build_field_part_associations(mf, *part, rank, ioEntity, missingFieldsCollectorPtr);
      }

      // If rank is != NODE_RANK, then see if field is defined on the nodes of this part
      if (rank != stk::topology::NODE_RANK && f->entity_rank() == stk::topology::NODE_RANK) {
        Ioss::GroupingEntity *nodeEntity = nullptr;
        std::string nodesName = part->name() + "_nodes";
        nodeEntity = m_region.get_entity(nodesName);
        if (nodeEntity == nullptr) {
          nodeEntity = m_region.get_entity("nodeblock_1");
        }
        if (nodeEntity != nullptr) {
          fieldIsMissing |= build_field_part_associations(mf, *part, stk::topology::NODE_RANK, nodeEntity,
                                                          missingFieldsCollectorPtr);
        }
      }

      if (fieldIsMissing) {
        ++numMissingFields;
      }
    }


    // Now handle the non-subsetted fields...

    // Check universal_part() NODE_RANK first...
    if (mf.m_subsetParts.empty()) {
      if (f->entity_rank() == stk::topology::NODE_RANK) {
        Ioss::GroupingEntity *nodeEntity = m_region.get_node_blocks()[0];
        bool fieldIsMissing = build_field_part_associations(mf, m_meta.universal_part(), stk::topology::NODE_RANK,
                                                            nodeEntity, missingFieldsCollectorPtr);
        if (fieldIsMissing) {
          ++numMissingFields;
        }
      }
    }

    // Now handle all non-nodeblock parts...
    for ( stk::mesh::Part * const part : m_meta.get_parts()) {
      // Check whether this part is an input part...
      if (stk::io::is_part_io_part(*part)) {
        stk::mesh::EntityRank rank = part_primary_entity_rank(*part);
        // Get Ioss::GroupingEntity corresponding to this part...
        Ioss::GroupingEntity *entity = m_region.get_entity(part->name());
        if (entity != nullptr && entity->type() != Ioss::SIDESET) {

          if (mf.m_subsetParts.empty()) {
            f = mf.field();
            bool fieldIsMissing = false;
            if (f->entity_rank() == rank) {
              fieldIsMissing |= build_field_part_associations(mf, *part, rank, entity, missingFieldsCollectorPtr);
            }

            // If rank is != NODE_RANK, then see if field is defined on the nodes of this part
            if (rank != stk::topology::NODE_RANK && f->entity_rank() == stk::topology::NODE_RANK) {
              Ioss::GroupingEntity *nodeEntity = nullptr;
              std::string nodesName = part->name() + "_nodes";
              nodeEntity = m_region.get_entity(nodesName);
              if (nodeEntity == nullptr) {
                nodeEntity = m_region.get_entity("nodeblock_1");
              }
              if (nodeEntity != nullptr) {
                fieldIsMissing |= build_field_part_associations(mf, *part, stk::topology::NODE_RANK, nodeEntity,
                                                                missingFieldsCollectorPtr);
              }
            }

            if (fieldIsMissing) {
              ++numMissingFields;
            }
          }
        }
      }
    }

    if (numMissingFields > 0 && missingFields==nullptr && throwOnErrorMessage) {
        std::ostringstream msg;
        msg << "ERROR: Input field processing could not find " << numMissingFields << " fields.\n";
        throw std::runtime_error( msg.str() );
    }

    add_missing_fields(missingFields, missingFieldsCollector);

    return numMissingFields;
  }

  bool InputQuery::process_fields_for_grouping_entity(stk::io::MeshField &mf,
                                                      const stk::mesh::Part &part,
                                                      Ioss::GroupingEntity *ioEntity,
                                                      MissingFieldMap *missingFieldsCollectorPtr)
  {
    STK_ThrowRequireMsg(ioEntity != nullptr, "Null IO entity");

    bool doesFieldExist = false;

    stk::mesh::FieldBase *f = mf.field();

    stk::mesh::EntityRank rank = part_primary_entity_rank(part);
    if(f->entity_rank() == rank) {
      const std::string &dbName = mf.db_name();
      unsigned numStates = f->number_of_states();
      std::vector<stk::mesh::FieldState> missingStates;
      if (numStates > 1) {
        bool hasAllStates = all_field_states_exist_on_io_entity(dbName, f, ioEntity, missingStates, m_multiStateSuffixes);
        if(hasAllStates == false) {
          if (missingFieldsCollectorPtr) {
            for (stk::mesh::FieldState missingState : missingStates) {
              (*missingFieldsCollectorPtr)[f->field_state(missingState)] = &mf;
            }
          }
        } else {
          doesFieldExist = true;
        }
      }

      if(doesFieldExist == false) {
        doesFieldExist = ioEntity->field_exists(dbName);
        if (!doesFieldExist) {
          if (missingFieldsCollectorPtr) {
            (*missingFieldsCollectorPtr)[f] = &mf;
          }
        }
      }

      // See if field with that name exists on ioEntity...
      if (doesFieldExist) {
        mf.add_part(f->entity_rank(), part, ioEntity);
        mf.set_single_state((m_dbPurpose == stk::io::READ_RESTART) ? false : true);
        mf.set_active();
      }
    }

    return doesFieldExist;
  }

  int InputQuery::build_field_part_associations_from_grouping_entity(stk::io::MeshField& mf,
                                                                     std::vector<stk::io::MeshField> *missingFields,
                                                                     const bool throwOnErrorMessage)
  {
    int numMissingFields = 0;

    if(mf.is_active()) {
      return numMissingFields;
    }

    MissingFieldMap missingFieldCollector;
    bool doesFieldExist = false;
    stk::mesh::Part &universalPart = m_meta.universal_part();
    Ioss::GroupingEntity * universalNodeEntity = m_region.get_entity("nodeblock_1");

    doesFieldExist |= process_fields_for_grouping_entity(mf, universalPart, universalNodeEntity, &missingFieldCollector);

    for ( stk::mesh::Part * const part : m_meta.get_parts() ) {
      // Check whether this part is an input part...
      if (stk::io::is_part_io_part(*part)) {
        // Get Ioss::GroupingEntity corresponding to this part...
        Ioss::GroupingEntity *ioEntity = m_region.get_entity(part->name());

        if(ioEntity == nullptr) {
          continue;
        }

        doesFieldExist |= process_fields_for_grouping_entity(mf, *part, ioEntity, &missingFieldCollector);

        if(ioEntity->type() == Ioss::SIDEBLOCK || ioEntity->type() == Ioss::SIDESET) {
          static const std::string s_nodeset_suffix("_n");

          std::string nsName = part->name();
          nsName += s_nodeset_suffix;
          Ioss::NodeSet *ioNodeSet = m_region.get_nodeset(nsName);
          if(ioNodeSet != nullptr) {
            // Process hidden nodesets
            doesFieldExist |= process_fields_for_grouping_entity(mf, *part, ioNodeSet, &missingFieldCollector);
          }
        }

        if(ioEntity->type() == Ioss::SIDESET) {
          Ioss::SideSet* sideSet = dynamic_cast<Ioss::SideSet*>(ioEntity);
          auto faceBlocks = sideSet->get_side_blocks();
          for (auto faceBlock : faceBlocks) {
            doesFieldExist |= process_fields_for_grouping_entity(mf, *part, faceBlock, &missingFieldCollector);
          }
        }
      }
    }

    if (!doesFieldExist) {
      numMissingFields += missingFieldCollector.size();
      if (nullptr != missingFields) {
        add_missing_fields(missingFields, missingFieldCollector);
      }
      else {
        for (auto missingField : missingFieldCollector) {
          std::cout << "Missing field: " << missingField.second->db_name() << std::endl;
        }
      }
    }

    if (numMissingFields > 0 && missingFields==nullptr && throwOnErrorMessage) {
      std::ostringstream msg;
      msg << "ERROR: Input field processing could not find " << numMissingFields << " fields.\n";
      throw std::runtime_error( msg.str() );
    }

    return numMissingFields;
  }

  void InputQuery::build_field_part_associations_for_part(stk::io::MeshField &mf, const stk::mesh::Part * part)
  {
    stk::mesh::FieldBase *f = mf.field();
    stk::mesh::EntityRank rank = part_primary_entity_rank(*part);
    // Get Ioss::GroupingEntity corresponding to this part...
    Ioss::GroupingEntity *entity = m_region.get_entity(part->name());

    if (entity != nullptr) {
      if (f->entity_rank() == rank) {
        build_field_part_associations(mf, *part, rank, entity);
        process_fields_for_grouping_entity(mf, *part, entity);

        if(entity->type() == Ioss::SIDESET) {
          auto io_side_set = dynamic_cast<Ioss::SideSet*>(entity);
          STK_ThrowRequire(io_side_set != nullptr);
          auto fbs = io_side_set->get_side_blocks();

          for(auto& io_fblock : fbs) {
            build_field_part_associations(mf, *part, rank, io_fblock);
            process_fields_for_grouping_entity(mf, *part, io_fblock);
          }
        }
      }

      // If rank is != NODE_RANK, then see if field is defined on the nodes of this part
      if (rank != stk::topology::NODE_RANK && f->entity_rank() == stk::topology::NODE_RANK) {
        Ioss::GroupingEntity *nodeEntity = nullptr;
        std::string nodes_name = part->name() + "_nodes";

        nodeEntity = m_region.get_entity(nodes_name);

        if (nodeEntity == nullptr) {
          nodes_name = part->name() + "_n";
          nodeEntity = m_region.get_entity(nodes_name);
        }

        if (nodeEntity == nullptr) {
          nodeEntity = m_region.get_entity("nodeblock_1");
        }
        if (nodeEntity != nullptr) {
          build_field_part_associations(mf, *part, stk::topology::NODE_RANK, nodeEntity);
          process_fields_for_grouping_entity(mf, *part, nodeEntity);
        }
      }
    }
  }

  bool verify_field_request(const Ioss::Region& region,  const stk::mesh::MetaData& meta,
                            const DatabasePurpose dbPurpose, const std::vector<std::string>& multiStateSuffixes,
                            const stk::io::MeshField &meshField, bool /*printWarning*/)
  {
    stk::io::InputQuery iq(region, meta, dbPurpose, (multiStateSuffixes.empty() ? nullptr : &multiStateSuffixes));

    stk::io::MeshField mf(meshField.field(), meshField.db_name());
    std::vector<stk::io::MeshField> missingFields;

    iq.build_field_part_associations(mf, &missingFields, false);
    iq.build_field_part_associations_from_grouping_entity(mf, &missingFields, false);

    if(missingFields.size() > 0) {
      std::ostringstream oss;
      oss << "For input IO field: "
          << meshField.db_name()
          << " the following associated fields for the requested STK field: "
          << meshField.field()->name()
          << " of rank: "
          << meshField.field()->entity_rank()
          << ", are missing in database: "
          << region.get_database()->get_filename()
          << std::endl;

      for(auto & missingField : missingFields) {
        oss << "\t" << missingField.field()->name() << std::endl;
      }

      stk::RuntimeWarning() << oss.str();
    }

    return mf.is_active();
  }

  bool verify_field_request(const StkMeshIoBroker &broker, const MeshField &meshField, bool printWarning)
  {
    auto region = broker.get_input_ioss_region();
    if(!region) {
      if(printWarning) {
        stk::RuntimeWarning() << "Broker has no input Ioss::Region" << std::endl;
      }

      return false;
    }

    if(broker.is_meta_data_null()) {
      if(printWarning) {
        stk::RuntimeWarning() << "Broker has no stk::mesh::MetaData defined" << std::endl;
      }

      return false;
    }

    const stk::mesh::MetaData &meta = broker.meta_data();
    InputFile& inputFile = broker.get_mesh_database(broker.get_active_mesh());

    return verify_field_request(*region, meta, inputFile.get_database_purpose(),
                                inputFile.get_multistate_suffixes(),
                                meshField, printWarning);
  }
}
}

