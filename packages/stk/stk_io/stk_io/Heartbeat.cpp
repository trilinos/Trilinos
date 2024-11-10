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
#include <stk_io/Heartbeat.hpp>
#include <exception>                 // for exception
#include <cstddef>                          // for size_t
#include <iostream>                         // for operator<<, basic_ostream
#include <stk_io/IOHelpers.hpp>             // for impl::add_global, writ...
#include <stk_io/IossBridge.hpp>            // for GlobalAnyVariable
#include <stk_util/util/ReportHandler.hpp>  // for ThrowErrorMsgIf, ThrowReq...
#include <utility>                          // for pair, move
#include "Ioss_DBUsage.h"                   // for DatabaseUsage, WRITE_HEAR...
#include "Ioss_DatabaseIO.h"                // for DatabaseIO
#include "Ioss_Field.h"                     // for Field, Field::BasicType
#include "Ioss_IOFactory.h"                 // for IOFactory
#include "Ioss_Property.h"                  // for Property
#include "Ioss_PropertyManager.h"           // for PropertyManager
#include "Ioss_Region.h"                    // for Region
#include "Ioss_State.h"                     // for STATE_DEFINE_TRANSIENT
#include "StkIoUtils.hpp"                   // for get_io_parameter_size_and...
#include "stk_util/parallel/Parallel.hpp"   // for parallel_machine_rank
#include "stk_util/util/ParameterList.hpp"  // for Parameter, Type, ...

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
            if (!properties.exists("FILE_FORMAT")) {
                properties.add(Ioss::Property("FILE_FORMAT", "csv"));
            }
        }
        else if (hb_type == TS_CSV) {
            if (!properties.exists("FILE_FORMAT")) {
                properties.add(Ioss::Property("FILE_FORMAT", "ts_csv"));
            }
        }
        else if (hb_type == TEXT) {
            if (!properties.exists("FILE_FORMAT")) {
                properties.add(Ioss::Property("FILE_FORMAT", "text"));
            }
        }
        else if (hb_type == TS_TEXT) {
            if (!properties.exists("FILE_FORMAT")) {
                properties.add(Ioss::Property("FILE_FORMAT", "ts_text"));
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
    m_region = std::make_shared<Ioss::Region>(db, filename);

}

void impl::Heartbeat::begin_define_transient()
{
    if (m_processor == 0) {
        STK_ThrowErrorMsgIf(m_currentStep != 0,
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

void impl::Heartbeat::internal_define_global_ref(const std::string &name,
                                        const std::any *value,
                                        stk::util::ParameterType::Type type,
                                        int copies,
                                        Ioss::Field::RoleType role)
{
    if (m_processor == 0) {
        STK_ThrowErrorMsgIf(m_currentStep != 0,
                         "At least one output step has been written to the history/heartbeat file. "
                         "Variables cannot be added anymore.");

        // Determine name and type of parameter...
        std::pair<size_t, Ioss::Field::BasicType> parameter_type = get_io_parameter_size_and_type(type, *value);
        impl::add_global(m_region, name, parameter_type.first, parameter_type.second, copies, role);
        m_fields.emplace_back(name, value, type);
    }
}

void impl::Heartbeat::define_global_ref(const std::string &name,
                                        const stk::util::Parameter &param,
                                        int copies,
                                        Ioss::Field::RoleType role)
{
  internal_define_global_ref(name, &param.value, param.type, copies, role);
}

void impl::Heartbeat::add_global_ref(const std::string &name,
                                     const stk::util::Parameter &param,
                                     int copies,
                                     Ioss::Field::RoleType role)
{
    if (m_processor == 0) {
        STK_ThrowErrorMsgIf(m_currentStep != 0,
                         "At least one output step has been written to the history/heartbeat file. "
                         "Variables cannot be added anymore.");

        Ioss::State currentState = m_region->get_state();
        if(currentState != Ioss::STATE_DEFINE_TRANSIENT) {
            m_region->begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
        }

        internal_define_global_ref(name, &param.value, param.type, copies, role);
    }
}

void impl::Heartbeat::internal_define_global_ref(const std::string &name,
                                        const std::any *value,
                                        const std::string &storage,
                                        Ioss::Field::BasicType dataType,
                                        int copies,
                                        Ioss::Field::RoleType role)
{
    if (m_processor == 0) {
        STK_ThrowErrorMsgIf(m_currentStep != 0,
                         "At least one output step has been written to the history/heartbeat file. "
                         "Variables cannot be added anymore.");

        std::pair<size_t, stk::util::ParameterType::Type> type = get_parameter_type_from_field_representation(storage, dataType, copies);

        // Determine name and type of parameter...
        std::pair<size_t, Ioss::Field::BasicType> parameter_type = get_io_parameter_size_and_type(type.second, *value);
        STK_ThrowRequireMsg(dataType == parameter_type.second, "data type must be consistent");
        impl::add_global(m_region, name, storage, dataType, copies, role);
        m_fields.emplace_back(name, value, type.second);
    }
}

void impl::Heartbeat::define_global_ref(const std::string &name,
                                        const stk::util::Parameter &param,
                                        const std::string &storage,
                                        Ioss::Field::BasicType dataType,
                                        int copies,
                                        Ioss::Field::RoleType role)
{
  internal_define_global_ref(name, &param.value, storage, dataType, copies, role);
}

void impl::Heartbeat::add_global_ref(const std::string &name,
                                     const stk::util::Parameter &param,
                                     const std::string &storage,
                                     Ioss::Field::BasicType dataType,
                                     int copies,
                                     Ioss::Field::RoleType role)
{
    if (m_processor == 0) {
        STK_ThrowErrorMsgIf(m_currentStep != 0,
                         "At least one output step has been written to the history/heartbeat file. "
                         "Variables cannot be added anymore.");

        Ioss::State currentState = m_region->get_state();
        if(currentState != Ioss::STATE_DEFINE_TRANSIENT) {
            m_region->begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
        }

        internal_define_global_ref(name, &param.value, storage, dataType, copies, role);
    }
}

void impl::Heartbeat::process_output_pre_write(int step, double time)
{
    if (m_processor == 0) {
        Ioss::State currentState = m_region->get_state();
        if(currentState == Ioss::STATE_DEFINE_TRANSIENT) {
            m_region->end_mode(Ioss::STATE_DEFINE_TRANSIENT);
        }
	if (currentState != Ioss::STATE_TRANSIENT) {
	    m_region->begin_mode(Ioss::STATE_TRANSIENT);
	}

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
