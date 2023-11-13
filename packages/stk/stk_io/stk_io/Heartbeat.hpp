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

#ifndef STK_IO_HEARTBEAT_HPP
#define STK_IO_HEARTBEAT_HPP
// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <Ioss_Field.h>                     // for Field, Field::TRANSIENT
#include <stk_io/IossBridge.hpp>            // for GlobalAnyVariable
#include <stk_util/util/ParameterList.hpp>  // for Type
#include <string>                           // for string
#include <vector>                           // for vector
#include "mpi.h"                            // for MPI_Comm
namespace Ioss { class PropertyManager; }
namespace Ioss { class Region; }
namespace Ioss { class Property; }
namespace stk { namespace io { class InputFile; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################
namespace stk { namespace mesh { class BulkData; } }

namespace Ioss { class DatabaseIO; }

namespace stk {
namespace io {

enum HeartbeatType {
  BINARY = 1, /* Exodus (history file) */
  CSV,        /* Comma-seperated values */
  TS_CSV,     /* same as CSV except lines preceded by timestamp*/
  TEXT,       /* Same as CSV except fields separated by tab (by default) */
  TS_TEXT,    /* same as TEXT except lines preceded by timestamp*/
  SPYHIS,     /* Format for use in spyhis plotter */
  NONE        /* Ignored in this class, can be used by apps */
};

namespace impl
{
class Heartbeat {
public:
    Heartbeat(const std::string &filename, HeartbeatType db_type,
              Ioss::PropertyManager properties, MPI_Comm comm,
              bool openFileImmediately = true);
    ~Heartbeat() {};

    void define_global_ref(const std::string &variableName,
                           const stk::util::Parameter &param,
                           int copies = 1,
                           Ioss::Field::RoleType role = Ioss::Field::TRANSIENT);

    void define_global_ref(const std::string &name,
                           const stk::util::Parameter &param,
                           const std::string &storage,
                           Ioss::Field::BasicType dataType,
                           int copies = 1,
                           Ioss::Field::RoleType role = Ioss::Field::TRANSIENT);

    void add_global_ref(const std::string &variableName,
                        const stk::util::Parameter &param,
                        int copies = 1,
                        Ioss::Field::RoleType role = Ioss::Field::TRANSIENT);

    void add_global_ref(const std::string &name,
                        const stk::util::Parameter &param,
                        const std::string &storage,
                        Ioss::Field::BasicType dataType,
                        int copies = 1,
                        Ioss::Field::RoleType role = Ioss::Field::TRANSIENT);

    void process_output(int step, double time);
    void process_output_pre_write(int step, double time);
    void process_output_write(int step, double time);
    void process_output_post_write(int step, double time);

    void flush_output() const;

    std::shared_ptr<Ioss::Region> get_heartbeat_ioss_region() {
        return m_region;
    }

    void begin_define_transient();
    void end_define_transient();

    bool has_global(const std::string &name);

private:
    void internal_define_global_ref(const std::string &variableName,
                           const std::any *value,
                           stk::util::ParameterType::Type type,
                           int copies = 1,
                           Ioss::Field::RoleType role = Ioss::Field::TRANSIENT);

    void internal_define_global_ref(const std::string &name,
                           const std::any *value,
                           const std::string &storage,
                           Ioss::Field::BasicType dataType,
                           int copies = 1,
                           Ioss::Field::RoleType role = Ioss::Field::TRANSIENT);

    std::vector<GlobalAnyVariable> m_fields;
    std::shared_ptr<Ioss::Region> m_region;

    int m_currentStep;
    int m_processor;
};

}
}
}
#endif
