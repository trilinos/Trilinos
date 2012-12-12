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
#if defined(HAVE_MPI) && !defined(NO_DOF_EXODUS_SUPPORT)
#include <par_exo/Iopx_IOFactory.h>    // for Iopx DatabaseIO
#endif
#include <tokenize.h>

#include <stddef.h>                     // for NULL
#include <string>                       // for string

#include "Ioss_CodeTypes.h"             // for MPI_Comm
#include "Ioss_DBUsage.h"               // for DatabaseUsage
#include "Ioss_IOFactory.h"             // for IOFactory

namespace Ioss { class DatabaseIO; }

#if defined(HAVE_MPI) && !defined(NO_DOF_EXODUS_SUPPORT)
namespace {
  std::string check_external_decomposition_property(MPI_Comm comm);
  bool check_external_composition_property(MPI_Comm comm, Ioss::DatabaseUsage db_usage);
}
#endif

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
#if defined(HAVE_MPI) && !defined(NO_DOF_EXODUS_SUPPORT)
    // The "exodus" and "parallel_exodus" databases can both be accessed
    // from this factory.  The "parallel_exodus" is returned only if the following
    // are true:
    // 0. The db_usage is 'READ_MODEL' (not suppported for READ_RESTART yet)
    // 1. Parallel run with >1 processor
    // 2. There is a DECOMPOSITION_METHOD specified in 'properties'
    // 3. The decomposition method is not "EXTERNAL"

    int proc_count = 1;
    if (communicator != MPI_COMM_NULL) {
      MPI_Comm_size(communicator, &proc_count);
    }
    
    bool decompose = false;
    if (proc_count > 1) {
      if (db_usage == Ioss::READ_MODEL) {
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
      else if (db_usage == Ioss::WRITE_RESULTS) {
	if (properties.exists("COMPOSE_RESULTS_FILE")) {
	  decompose = true;
	} else if (check_external_composition_property(communicator, db_usage)) {
	  decompose = true;
	}
	decompose = true;
      }
      else if (db_usage == Ioss::WRITE_RESTART) {
	if (properties.exists("COMPOSE_RESTART_FILE")) {
	  decompose = true;
	} else if (check_external_composition_property(communicator, db_usage)) {
	  decompose = true;
	}
      }
    }

    // Could call Iopx::DatabaseIO constructor directly, but that leads to some circular
    // dependencies and other yuks.
    if (decompose)
      return Ioss::IOFactory::create("dof_exodus", filename, db_usage, communicator, properties);
    else
#endif
      return new Ioex::DatabaseIO(NULL, filename, db_usage, communicator, properties);
  }
}

#if defined(HAVE_MPI) && !defined(NO_DOF_EXODUS_SUPPORT)
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

  bool check_external_composition_property(MPI_Comm comm, Ioss::DatabaseUsage db_usage)
  {
    // Check environment variable IOSS_PROPERTIES. If it exists, parse
    // the contents and see if it specifies the use of a single file for output...

    bool compose = false;
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
	std::string prop = Ioss::Utils::uppercase(property[0]);
	if (db_usage == Ioss::WRITE_RESULTS && prop == "COMPOSE_RESULTS_FILE") {
	  compose = true;
	  break;
	}
	else if (db_usage == Ioss::WRITE_RESTART && prop == "COMPOSE_RESTART_FILE") {
	  compose = true;
	  break;
	}
      }
    }
    return compose;
  }
}
#endif
