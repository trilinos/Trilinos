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

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stk_util/environment/Env.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <Ionit_Initializer.h>                       // for Initializer
#include <assert.h>                                  // for assert
#include <stdlib.h>                                  // for exit, etc
#include <string.h>                                  // for memcpy
#include <cstdint>                                   // for int64_t
#include <iostream>                                  // for operator<<, etc
#include <iterator>
#include <limits>                                    // for numeric_limits
#include <map>
#include <stdexcept>                                 // for runtime_error
#include <stk_io/IossBridge.hpp>                     // for FieldAndName, etc
#include <stk_io/Heartbeat.hpp>                      // for Heartbeat
#include <stk_io/IOHelpers.hpp>
#include <stk_mesh/base/BulkData.hpp>                // for BulkData, etc
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Field.hpp>                   // for Field
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>                // for MetaData, etc
#include <stk_util/environment/FileUtils.hpp>
#include <stk_util/util/ReportHandler.hpp>    // for ThrowErrorMsgIf, etc
#include <utility>                                   // for pair, make_pair
#include "Ioss_CodeTypes.h"                          // for NameList
#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h"                         // for DatabaseIO
#include "Ioss_ElementBlock.h"                       // for ElementBlock
#include "Ioss_ElementTopology.h"                    // for ElementTopology
#include "Ioss_EntityType.h"
#include "Ioss_Field.h"
#include "Ioss_GroupingEntity.h"                     // for GroupingEntity
#include "Ioss_IOFactory.h"                          // for IOFactory
#include "Ioss_NodeBlock.h"                          // for NodeBlock
#include "Ioss_NodeSet.h"                            // for NodeSet
#include "Ioss_ParallelUtils.h"                      // for ParallelUtils
#include "Ioss_Property.h"                           // for Property
#include "Ioss_PropertyManager.h"                    // for PropertyManager
#include "Ioss_Region.h"                             // for Region, etc
#include "Ioss_SideBlock.h"                          // for SideBlock
#include "Ioss_SideSet.h"                            // for SideSet
#include "Ioss_State.h"
#include "Ioss_VariableType.h"                       // for VariableType
#include "ProcessSetsOrBlocks.hpp"
#include "SidesetTranslator.hpp"
#include "StkIoUtils.hpp"
#include "Teuchos_RCP.hpp"                           // for RCP::operator->, etc
#include "boost/any.hpp"                             // for any_cast, any
#include "stk_io/DatabasePurpose.hpp"                // for DatabasePurpose, etc
#include "stk_io/MeshField.hpp"                      // for MeshField, etc
#include "stk_io/SidesetUpdater.hpp"
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/Entity.hpp"                  // for Entity
#include "stk_mesh/base/FieldBase.hpp"               // for FieldBase
#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_mesh/base/FieldState.hpp"              // for FieldState
#include "stk_mesh/base/Part.hpp"                    // for Part
#include "stk_mesh/base/Selector.hpp"                // for Selector, etc
#include "stk_mesh/base/Types.hpp"                   // for FieldVector, etc
#include "stk_topology/topology.hpp"                 // for topology, etc
#include "stk_util/parallel/Parallel.hpp"            // for ParallelMachine, etc
#include "stk_util/util/ParameterList.hpp"           // for Type, etc
#include "stk_util/diag/StringUtil.hpp"           // for Type, etc
#include "stk_util/util/string_case_compare.hpp"

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################


namespace stk {
namespace io {


impl::Heartbeat::Heartbeat(const std::string &filename, HeartbeatType hb_type,
                           Ioss::PropertyManager properties, stk::ParallelMachine comm,
                           bool openFileImmediately)
: m_currentStep(0), m_processor(0)
{
    if (comm != MPI_COMM_NULL) {
        m_processor = stk::parallel_machine_rank(comm);
    }

    std::string db_io_type = "exodusII";
    Ioss::DatabaseUsage db_usage = Ioss::WRITE_HISTORY;

    if (hb_type != BINARY) {
        db_io_type = "heartbeat";
        db_usage = Ioss::WRITE_HEARTBEAT;

        // Always add the "time" field to all heartbeat outputs...
        if (!properties.exists("SHOW_TIME_FIELD")) {
            properties.add(Ioss::Property("SHOW_TIME_FIELD", true));
        }

        if (hb_type == SPYHIS) {
            if (!properties.exists("FILE_FORMAT")) {
                properties.add(Ioss::Property("FILE_FORMAT", "spyhis"));
            }
        }
        else if (hb_type == CSV) {
            if (!properties.exists("SHOW_TIME_STAMP")) {
                properties.add(Ioss::Property("SHOW_TIME_STAMP", false));
            }
            if (!properties.exists("FIELD_SEPARATOR")) {
                properties.add(Ioss::Property("FIELD_SEPARATOR", ", "));
            }
        }
        else if (hb_type == TS_CSV) {
            if (!properties.exists("SHOW_TIME_STAMP")) {
                properties.add(Ioss::Property("SHOW_TIME_STAMP", true));
            }
            if (!properties.exists("FIELD_SEPARATOR")) {
                properties.add(Ioss::Property("FIELD_SEPARATOR", ", "));
            }
        }
        else if (hb_type == TEXT) {
            if (!properties.exists("SHOW_TIME_STAMP")) {
                properties.add(Ioss::Property("SHOW_TIME_STAMP", false));
            }
            if (!properties.exists("FIELD_SEPARATOR")) {
                properties.add(Ioss::Property("FIELD_SEPARATOR", "\t"));
            }
        }
        else if (hb_type == TS_TEXT) {
            if (!properties.exists("SHOW_TIME_STAMP")) {
                properties.add(Ioss::Property("SHOW_TIME_STAMP", true));
            }
            if (!properties.exists("FIELD_SEPARATOR")) {
                properties.add(Ioss::Property("FIELD_SEPARATOR", "\t"));
            }
        }
    }

    Ioss::DatabaseIO *db = Ioss::IOFactory::create(db_io_type, filename,
                                                   db_usage, comm, properties);
    if (db == nullptr || (openFileImmediately && !db->ok())) {
        std::cerr << "ERROR: Could not open history/heartbeat database '" << filename << "'\n";
        return;
    }

    // NOTE: 'region' owns 'db' pointer at this time...
    m_region = Teuchos::rcp(new Ioss::Region(db, filename));

}

void impl::Heartbeat::begin_define_transient()
{
    if (m_processor == 0) {
        ThrowErrorMsgIf (m_currentStep != 0,
                         "At least one output step has been written to the history/heartbeat file. "
                         "Variables cannot be added anymore.");

        Ioss::State currentState = m_region->get_state();
        if(currentState != Ioss::STATE_DEFINE_TRANSIENT) {
            m_region->begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
        }
    }
}

void impl::Heartbeat::end_define_transient()
{
    if (m_processor == 0) {
        Ioss::State currentState = m_region->get_state();
        if(currentState == Ioss::STATE_DEFINE_TRANSIENT) {
            m_region->end_mode(Ioss::STATE_DEFINE_TRANSIENT);
        }
    }
}

bool impl::Heartbeat::has_global(const std::string &name)
{
    return m_region->field_exists(name);
}

void impl::Heartbeat::define_global_ref(const std::string &name,
                                        const boost::any *value,
                                        stk::util::ParameterType::Type type,
                                        int copies,
                                        Ioss::Field::RoleType role)
{
    if (m_processor == 0) {
        ThrowErrorMsgIf (m_currentStep != 0,
                         "At least one output step has been written to the history/heartbeat file. "
                         "Variables cannot be added anymore.");

        // Determine name and type of parameter...
        std::pair<size_t, Ioss::Field::BasicType> parameter_type = get_io_parameter_size_and_type(type, *value);
        internal_add_global(m_region, name, parameter_type.first, parameter_type.second, copies, role);
        m_fields.emplace_back(name, value, type);
    }
}

void impl::Heartbeat::add_global_ref(const std::string &name,
                                     const boost::any *value,
                                     stk::util::ParameterType::Type type,
                                     int copies,
                                     Ioss::Field::RoleType role)
{
    if (m_processor == 0) {
        ThrowErrorMsgIf (m_currentStep != 0,
                         "At least one output step has been written to the history/heartbeat file. "
                         "Variables cannot be added anymore.");

        Ioss::State currentState = m_region->get_state();
        if(currentState != Ioss::STATE_DEFINE_TRANSIENT) {
            m_region->begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
        }

        define_global_ref(name, value, type, copies, role);
    }
}

void impl::Heartbeat::define_global_ref(const std::string &name,
                                        const boost::any *value,
                                        const std::string &storage,
                                        Ioss::Field::BasicType dataType,
                                        int copies,
                                        Ioss::Field::RoleType role)
{
    if (m_processor == 0) {
        ThrowErrorMsgIf (m_currentStep != 0,
                         "At least one output step has been written to the history/heartbeat file. "
                         "Variables cannot be added anymore.");

        std::pair<size_t, stk::util::ParameterType::Type> type = get_parameter_type_from_field_representation(storage, dataType, copies);

        // Determine name and type of parameter...
        std::pair<size_t, Ioss::Field::BasicType> parameter_type = get_io_parameter_size_and_type(type.second, *value);
        ThrowRequireMsg(dataType == parameter_type.second, "data type must be consistent");
        internal_add_global(m_region, name, storage, dataType, copies, role);
        m_fields.emplace_back(name, value, type.second);
    }
}

void impl::Heartbeat::add_global_ref(const std::string &name,
                                     const boost::any *value,
                                     const std::string &storage,
                                     Ioss::Field::BasicType dataType,
                                     int copies,
                                     Ioss::Field::RoleType role)
{
    if (m_processor == 0) {
        ThrowErrorMsgIf (m_currentStep != 0,
                         "At least one output step has been written to the history/heartbeat file. "
                         "Variables cannot be added anymore.");

        Ioss::State currentState = m_region->get_state();
        if(currentState != Ioss::STATE_DEFINE_TRANSIENT) {
            m_region->begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
        }

        define_global_ref(name, value, storage, dataType, copies, role);
    }
}

void impl::Heartbeat::process_output_pre_write(int step, double time)
{
    if (m_processor == 0) {
        Ioss::State currentState = m_region->get_state();
        if(currentState == Ioss::STATE_DEFINE_TRANSIENT) {
            m_region->end_mode(Ioss::STATE_DEFINE_TRANSIENT);
        }

        m_region->begin_mode(Ioss::STATE_TRANSIENT);
        m_currentStep = m_region->add_state(time);
        m_region->begin_state(m_currentStep);
    }
}

void impl::Heartbeat::process_output_write(int step, double time)
{
    if (m_processor == 0) {
        write_defined_global_any_fields(m_region, m_fields);
    }
}

void impl::Heartbeat::process_output_post_write(int step, double time)
{
    if (m_processor == 0) {
        m_region->end_state(m_currentStep);
        m_region->end_mode(Ioss::STATE_TRANSIENT);
    }
}

void impl::Heartbeat::process_output(int step, double time)
{
    process_output_pre_write(step, time);
    process_output_write(step, time);
    process_output_post_write(step, time);
}

void impl::Heartbeat::flush_output() const
{
    if (m_processor == 0) {
        m_region->get_database()->flush_database();
    }
}

} // namespace io
} // namespace stk
