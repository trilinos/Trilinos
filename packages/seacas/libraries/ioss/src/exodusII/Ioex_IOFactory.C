// Copyright(C) 1999-2010
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

#include <exodusII/Ioex_IOFactory.h>    // for Ioex IOFactory

#include <exodusII/Ioex_DatabaseIO.h>   // for Ioex DatabaseIO
#include <par_exo/Iopx_DatabaseIO.h>    // for Iopx DatabaseIO
#include <tokenize.h>

#include <stddef.h>                     // for NULL
#include <string>                       // for string

#include "Ioss_CodeTypes.h"             // for MPI_Comm
#include "Ioss_DBUsage.h"               // for DatabaseUsage
#include "Ioss_IOFactory.h"             // for IOFactory

namespace Ioss { class DatabaseIO; }

namespace {
  std::string check_external_decomposition_property(MPI_Comm comm);
}

namespace Ioex {

  const IOFactory* IOFactory::factory()
  {
    static IOFactory registerThis;
    return &registerThis;
  }

  IOFactory::IOFactory()
    : Ioss::IOFactory("exodusII")
  {
    Ioss::IOFactory::alias("exodusII", "exodusii");
    Ioss::IOFactory::alias("exodusII", "exodus");
    Ioss::IOFactory::alias("exodusII", "genesis");
  }

  Ioss::DatabaseIO* IOFactory::make_IO(const std::string& filename,
				       Ioss::DatabaseUsage db_usage,
				       MPI_Comm communicator,
				       const Ioss::PropertyManager &properties) const
  {
    // The "exodus" and "parallel_exodus" databases can both be accessed
    // from this factory.  The "parallel_exodus" is returned only if the following
    // are true:
    // 0. The db_usage is 'READ_MODEL' (not suppported for READ_RESTART yet)
    // 1. Parallel run with >1 processor
    // 2. There is a DECOMPOSITION_METHOD specified in 'properties'
    // 3. The decomposition method is not "EXTERNAL"

    bool decompose = false;
    if (db_usage == Ioss::READ_MODEL) {
      int proc_count = 1;
      if (communicator != MPI_COMM_NULL) {
	MPI_Comm_size(communicator, &proc_count);

	if (proc_count > 1) {
	  // Check for property...
	  if (properties.exists("DECOMPOSITION_METHOD")) {
	    std::string method = properties.get("DECOMPOSITION_METHOD").get_string();
	    if (method != "EXTERNAL") {
	      decompose = true;
	    }
	  } else {
	    std::string method = check_external_decomposition_property(communicator);
	    if (!method.empty() && method != "EXTERNAL") {
	      decompose = true;
	    }
	  }
	}
      }

      if (decompose)
	return new Iopx::DatabaseIO(NULL, filename, db_usage, communicator, properties);
      else
	return new Ioex::DatabaseIO(NULL, filename, db_usage, communicator, properties);
    }
  }
}

namespace {
  std::string check_external_decomposition_property(MPI_Comm comm)
  {
    // Check environment variable IOSS_PROPERTIES. If it exists, parse
    // the contents and see if it specifies a decomposition method.
  
    std::string decomp_method;
    
    Ioss::ParallelUtils util(comm);
    std::string env_props;
    if (util.get_environment("IOSS_PROPERTIES", env_props, true)) {
      // env_props string should be of the form
      // "PROP1=VALUE1:PROP2=VALUE2:..."
      std::vector<std::string> prop_val;
      Ioss::tokenize(env_props, ":", prop_val);
    
      for (size_t i=0; i < prop_val.size(); i++) {
	std::vector<std::string> property;
	Ioss::tokenize(prop_val[i], "=", property);
	if (property.size() != 2) {
	  std::ostringstream errmsg;
	  errmsg << "ERROR: Invalid property specification found in IOSS_PROPERTIES environment variable\n"
		 << "       Found '" << prop_val[i] << "' which is not of the correct PROPERTY=VALUE form";
	  IOSS_ERROR(errmsg);
	}
	std::string prop = Ioss::Utils::uppercase(property[0]);
	if (prop == "DECOMPOSITION_METHOD") {
	  std::string value = property[1];
	  decomp_method = Ioss::Utils::uppercase(value);
	}
      }
    }
    return decomp_method;
  }
}
