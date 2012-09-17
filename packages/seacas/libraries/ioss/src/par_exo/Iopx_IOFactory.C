// Copyright(C) 2012
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
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

#include <exodusII/Ioex_DatabaseIO.h>   // for Ioex DatabaseIO
#include <par_exo/Iopx_DatabaseIO.h>    // for Iopx DatabaseIO

#include <par_exo/Iopx_IOFactory.h>    // for IOFactory

#include <stddef.h>                     // for NULL
#include <string>                       // for string

#include "Ioss_CodeTypes.h"             // for MPI_Comm
#include "Ioss_DBUsage.h"               // for DatabaseUsage
#include "Ioss_IOFactory.h"             // for IOFactory

namespace Ioss { class DatabaseIO; }

namespace Iopx {

  const IOFactory* IOFactory::factory()
  {
    static IOFactory registerThis;
    return &registerThis;
  }

  IOFactory::IOFactory()
    : Ioss::IOFactory("parallel_exodus")
  {
    Ioss::IOFactory::alias("parallel_exodus", "dof_exodus");
    Ioss::IOFactory::alias("parallel_exodus", "dof");
  }

  Ioss::DatabaseIO* IOFactory::make_IO(const std::string& filename,
				       Ioss::DatabaseUsage db_usage,
				       MPI_Comm communicator,
				       const Ioss::PropertyManager &properties) const
  {
    int proc_count = 1;
    MPI_Comm_size(communicator, &proc_count);
    // READ_RESTART not officially supported.  Here just so gdsjaar can experiment
    if (proc_count > 1 && (db_usage == Ioss::READ_MODEL || db_usage == Ioss::READ_RESTART)) {
      return new Iopx::DatabaseIO(NULL, filename, db_usage, communicator, properties);
    } else {
      return new Ioex::DatabaseIO(NULL, filename, db_usage, communicator, properties);
    }
  }

}
