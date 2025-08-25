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
#include <stk_io/MeshField.hpp>
#include <Ioss_GroupingEntity.h>            // for GroupingEntity
#include <Ioss_VariableType.h>              // for VariableType
#include <cassert>                          // for assert
#include <ostream>                          // for operator<<, basic_ostream
#include <stk_io/DbStepTimeInterval.hpp>    // for DBStepTimeInterval
#include <stk_io/IossBridge.hpp>            // for STKIORequire, get_input_e...
#include <stk_mesh/base/FieldBase.hpp>      // for FieldBase, field_data
#include <stk_mesh/base/MetaData.hpp>       // for MetaData
#include <stk_mesh/base/Part.hpp>           // for Part
#include <stk_util/util/ReportHandler.hpp>  // for ThrowErrorMsgIf
#include <string>                           // for string, operator<<, opera...
#include <utility>                          // for swap
#include "Ioss_EntityType.h"                // for NODEBLOCK
#include "Ioss_Field.h"                     // for Field
#include "Ioss_Region.h"                    // for Region
#include "stk_mesh/base/BulkData.hpp"       // for BulkData
#include "stk_mesh/base/Entity.hpp"         // for Entity
#include "stk_mesh/base/EntityKey.hpp"      // for operator<<
#include "stk_mesh/base/FieldState.hpp"     // for FieldState, StateNew
#include "stk_topology/topology.hpp"        // for operator<<, topology, top...
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
  namespace io {

MeshField::MeshField(stk::mesh::FieldBase *field,
		     const std::string &db_name,
		     TimeMatchOption tmo)
  : m_field(field),
    m_dbName(db_name),
    m_timeToRead(0.0),
    m_timeMatch(tmo),
    m_oneTimeOnly(false),
    m_singleState(true),
    m_isActive(false),
    m_fieldRestored(false),
    m_timeRestored(std::numeric_limits<double>::max())
{
  if (db_name == "") {
    m_dbName = field->name();
  }

  STK_ThrowErrorMsgIf(m_timeMatch == LINEAR_INTERPOLATION && m_field->type_is<int>(),
		  "ERROR: Input interpolation field '" << m_field->name()
		  << "' is an integer field.  Only double fields can be interpolated.");
}
    
MeshField::MeshField(stk::mesh::FieldBase &field,
		     const std::string &db_name,
		     TimeMatchOption tmo)
  : m_field(&field),
    m_dbName(db_name),
    m_timeToRead(0.0),
    m_timeMatch(tmo),
    m_oneTimeOnly(false),
    m_singleState(true),
    m_isActive(false),
    m_fieldRestored(false),
    m_timeRestored(std::numeric_limits<double>::max())
{
  if (db_name == "") {
    m_dbName = field.name();
  }

  STK_ThrowErrorMsgIf(m_timeMatch == LINEAR_INTERPOLATION && m_field->type_is<int>(),
		  "ERROR: Input interpolation field '" << m_field->name()
		  << "' is an integer field.  Only double fields can be interpolated.");
}

MeshField::~MeshField()
{}

MeshField& MeshField::set_read_time(double time_to_read)
{
  m_timeToRead = time_to_read;
  m_timeMatch = SPECIFIED;
  m_oneTimeOnly = true;
  return *this;
}

MeshField& MeshField::set_active()
{
  m_isActive = true;
  return *this;
}
MeshField& MeshField::set_inactive()
{
  m_isActive = false;
  return *this;
}

MeshField& MeshField::set_single_state(bool yesno)
{
  m_singleState = yesno;
  return *this;
}

MeshField& MeshField::set_read_once(bool yesno)
{
  m_oneTimeOnly = yesno;
  return *this;
}

MeshField& MeshField::add_subset(const stk::mesh::Part &part)
{
  m_subsetParts.push_back(&part);
  return *this;
}

void MeshField::add_part(const stk::mesh::EntityRank rank,
			 const stk::mesh::Part &part,
			 Ioss::GroupingEntity *io_entity)
{
  m_fieldParts.emplace_back(rank, &part, io_entity, m_dbName);
}

bool MeshField::operator==(const MeshField &other) const
{
  // NOTE: Do not check 'm_dbName'.  The behavior is that
  // if the user attempts to add to MeshFields that only differ by
  // the database name, the name is updated to the most recent
  // MeshField database name.
  return m_field       == other.m_field &&
         m_subsetParts == other.m_subsetParts;
}

void MeshField::fill_entity_list_cache(const stk::mesh::BulkData &bulk)
{
  for (auto & meshFieldPart : m_fieldParts) {
    const stk::mesh::EntityRank rank = meshFieldPart.get_entity_rank();
    Ioss::GroupingEntity *io_entity = meshFieldPart.get_io_entity();
    meshFieldPart.set_cached_entity_list(stk::io::get_input_entity_list(io_entity, rank, bulk));
  }
}

double MeshField::restore_field_data_at_step(Ioss::Region *region,
                                             stk::mesh::BulkData &bulk,
                                             int step,
                                             bool ignore_missing_fields,
                                             std::vector<std::string>* multiStateSuffixes,
                                             bool useEntityCacheList)
{
  STKIORequire(step > 0);

  double time_read = region->begin_state(step);

  for (auto & meshFieldPart : m_fieldParts) {
    const stk::mesh::EntityRank rank = meshFieldPart.get_entity_rank();
    Ioss::GroupingEntity *io_entity = meshFieldPart.get_io_entity();
    const std::vector<stk::mesh::Entity> & entity_list = (useEntityCacheList) ?
                                                         meshFieldPart.get_cached_entity_list() :
                                                         stk::io::get_input_entity_list(io_entity, rank, bulk);
    const stk::mesh::Part *stk_part = meshFieldPart.get_stk_part();

    // If the field being restored is a nodal field stored on the
    // Ioss::Nodeblock on the database, but is not being applied to
    // the stk universal part, then we need to transfer a subset of
    // the data to the stk field. The subset will be defined as a
    // selector of the stk part.
    bool subsetted = rank == stk::topology::NODE_RANK &&
        io_entity->type() == Ioss::NODEBLOCK &&
        *stk_part != bulk.mesh_meta_data().universal_part();

    size_t state_count = m_field->number_of_states();
    stk::mesh::FieldState state = m_field->state();
    // If the multi-state field is not "set" at the newest state, then the user has
    // registered the field at a specific state and only that state should be input.
    if(m_singleState || state_count == 1 || state != stk::mesh::StateNew) {
      if (subsetted) {
        stk::io::subsetted_field_data_from_ioss(bulk, m_field, entity_list,
                                                io_entity, stk_part, m_dbName);
      } else {
        stk::io::field_data_from_ioss(bulk, m_field, entity_list,
                                      io_entity, m_dbName);
      }
    } else {
      if (subsetted) {
        stk::io::subsetted_multistate_field_data_from_ioss(bulk, m_field, entity_list,
                                                           io_entity, stk_part, m_dbName, state_count,
                                                           ignore_missing_fields, multiStateSuffixes);
      } else {
        stk::io::multistate_field_data_from_ioss(bulk, m_field, entity_list,
                                                 io_entity, m_dbName, state_count,
                                                 ignore_missing_fields, multiStateSuffixes);
      }
    }

    if (m_oneTimeOnly) {
      meshFieldPart.release_field_data();
    }
  }
  region->end_state(step);

  return time_read;
}

double MeshField::restore_field_data(stk::mesh::BulkData &bulk,
                                     const stk::io::DBStepTimeInterval &sti,
                                     bool ignore_missing_fields,
                                     std::vector<std::string>* multiStateSuffixes)
{
  m_fieldRestored = false;
  m_timeRestored = std::numeric_limits<double>::max();

  double time_read = -1.0;
  if (!is_active())
    return time_read;

  if (m_timeMatch == CLOSEST || m_timeMatch == SPECIFIED) {
    int step = 0;
    if (m_timeMatch == CLOSEST) {
      step = sti.get_closest_step();
    }
    else if (m_timeMatch == SPECIFIED) {
      DBStepTimeInterval sti2(sti.region, m_timeToRead);
      step = sti2.get_closest_step();
    }

    time_read = restore_field_data_at_step(sti.region, bulk, step, ignore_missing_fields, multiStateSuffixes);
  }
  else if (m_timeMatch == LINEAR_INTERPOLATION) {
    // Interpolation only handles single-state fields with state StateNew
    size_t state_count = m_field->number_of_states();
    stk::mesh::FieldState state = m_field->state();
    STKIORequire(m_singleState || state_count == 1 || state != stk::mesh::StateNew);

    stk::mesh::field_data_execute<double, stk::mesh::OverwriteAll>(*m_field,
      [&](auto& fieldData) {
        for (auto &field_part : m_fieldParts) {
          // Get data at beginning of interval...
          std::vector<double> values;
          field_part.get_interpolated_field_data(sti, values);

          Ioss::GroupingEntity *io_entity = field_part.get_io_entity();
          const Ioss::Field &io_field = io_entity->get_fieldref(m_dbName);
          int field_component_count = io_field.transformed_storage()->component_count();

          const stk::mesh::EntityRank rank = field_part.get_entity_rank();
          std::vector<stk::mesh::Entity> entity_list = stk::io::get_input_entity_list(io_entity, rank, bulk);

          for (size_t i=0; i < entity_list.size(); ++i) {
            if (bulk.is_valid(entity_list[i]) && m_field->defined_on(entity_list[i])) {
              auto fldValues = fieldData.entity_values(entity_list[i]);
              STK_ThrowRequireMsg(fldValues.num_scalars() == field_component_count, "Size mis-match for field "<<m_field->name()<<", num-scalars = "<<fldValues.num_scalars()<<" but IO-field reports num-components = "<<field_component_count);

              for(stk::mesh::ScalarIdx idx : fldValues.scalars()) {
                fldValues(idx) = values[i*field_component_count+idx];
              }
            }
          }
          if (m_oneTimeOnly) {
            field_part.release_field_data();
          }
        }
      }
    );

    time_read = sti.t_analysis;
  }
  if (m_oneTimeOnly) {
    set_inactive();
  }

  m_fieldRestored = true;
  m_timeRestored = time_read;

  return time_read;
}

void MeshField::clear_field_parts()
{
  m_fieldParts.clear();
}

void MeshFieldPart::release_field_data()
{
  m_preStep = 0;
  m_postStep = 0;
  std::vector<double>().swap(m_preData);
  std::vector<double>().swap(m_postData);
}

void MeshFieldPart::load_field_data(const DBStepTimeInterval &sti)
{
  // Use cached data if possible; avoid reading from disk...
  
  if (sti.exists_before && m_preStep != sti.s_before) {
    assert(sti.s_before > 0);

    if (sti.s_before == m_postStep) {
      m_preData.swap(m_postData);
      std::swap(m_preStep, m_postStep);
    }
    else {
      // See if postStep can use my current data...
      if (sti.exists_after && sti.s_after == m_preStep) {
	m_postData.swap(m_preData);
	std::swap(m_postStep, m_preStep);
      }
      m_preStep = sti.s_before;
      sti.region->begin_state(m_preStep);
      m_ioEntity->get_field_data(m_dbName, m_preData);
      sti.region->end_state(m_preStep);
    }
  }

  if (sti.exists_after && m_postStep != sti.s_after) {
    m_postStep = sti.s_after;
    assert(m_postStep > 0);

    if (m_preStep == m_postStep) {
      m_postData.resize(0);
      m_postData.reserve(m_preData.size());
      m_postData.insert(m_postData.end(), m_preData.begin(), m_preData.end());
    }
    else {
      sti.region->begin_state(m_postStep);
      m_ioEntity->get_field_data(m_dbName, m_postData);
      sti.region->end_state(m_postStep);
    }
  }
}

void MeshFieldPart::get_interpolated_field_data(const DBStepTimeInterval &sti, std::vector<double> &values)
{
  load_field_data(sti);
  size_t values_size = sti.exists_before ? m_preData.size() : m_postData.size();
  values.resize(0);
  values.reserve(values_size);

  if (sti.exists_before && !sti.exists_after) {
    values.insert(values.end(), m_preData.begin(), m_preData.end()); 
  }
  else if (!sti.exists_before && sti.exists_after) {
    values.insert(values.end(), m_postData.begin(), m_postData.end()); 
  }
  else {
    assert(sti.exists_before && sti.exists_after);
    if (sti.s_after == sti.s_before) {
      // No interpolation. preData and postData contain the same step.
      values.insert(values.end(), m_preData.begin(), m_preData.end()); 
    }
    else {
      // Interpolate
      double tb = sti.region->get_state_time(sti.s_before);
      double ta = sti.region->get_state_time(sti.s_after);
      double delta =  ta - tb;
      double frac  = (sti.t_analysis - tb) / delta;

      assert(m_preData.size() == m_postData.size());
      for (size_t i=0; i < m_preData.size(); i++) {
	values.push_back((1.0 - frac) * m_preData[i] + frac * m_postData[i]);
      }
    }
  }
}

}}  // Close namespaces
