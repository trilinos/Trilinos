// Copyright(C) 1999-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
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
//     * Neither the name of NTESS nor the names of its
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

#include <cgns/Iocgns_DatabaseIO.h> // for DatabaseIO -- serial
#include <cgns/Iocgns_IOFactory.h>
#include <cstddef> // for nullptr
#ifdef SEACAS_HAVE_MPI
#include <cgns/Iocgns_ParallelDatabaseIO.h> // for DatabaseIO -- parallel
#endif
#include "Ioss_DBUsage.h"   // for DatabaseUsage
#include "Ioss_IOFactory.h" // for IOFactory
#include <string>           // for string
#include <tokenize.h>

namespace Ioss {
  class PropertyManager;
}

#if defined(SEACAS_HAVE_MPI)
namespace {
  std::string check_decomposition_property(MPI_Comm comm, const Ioss::PropertyManager &properties,
                                           Ioss::DatabaseUsage db_usage);
  bool        check_composition_property(MPI_Comm comm, const Ioss::PropertyManager &properties,
                                         Ioss::DatabaseUsage db_usage);
} // namespace
#endif

namespace Iocgns {

  const IOFactory *IOFactory::factory()
  {
    static IOFactory registerThis;
    return &registerThis;
  }

  IOFactory::IOFactory() : Ioss::IOFactory("cgns")
  {
#if defined(SEACAS_HAVE_MPI)
    Ioss::IOFactory::alias("cgns", "dof_cgns");
    Ioss::IOFactory::alias("cgns", "par_cgns");
#endif
  }

  Ioss::DatabaseIO *IOFactory::make_IO(const std::string &filename, Ioss::DatabaseUsage db_usage,
                                       MPI_Comm                     communicator,
                                       const Ioss::PropertyManager &properties) const
  {
// The "cgns" and "parallel_cgns" databases can both be accessed from
// this factory.  The "parallel_cgns" is returned if being run on more
// than 1 processor unless the decomposition property is set and the
// value is "external" or the composition property is set with value "external"
#if defined(SEACAS_HAVE_MPI)
    int proc_count = 1;
    if (communicator != MPI_COMM_NULL) {
      MPI_Comm_size(communicator, &proc_count);
    }

    bool decompose = false;

    if (proc_count > 1) {
      decompose = true; // Default to decompose instead of file-per-processor if parallel.
      if (db_usage == Ioss::READ_MODEL || db_usage == Ioss::READ_RESTART) {
        std::string method = check_decomposition_property(communicator, properties, db_usage);
        if (!method.empty() && method == "EXTERNAL") {
          decompose = false;
        }
      }
      else if (db_usage == Ioss::WRITE_RESULTS || db_usage == Ioss::WRITE_RESTART) {
        decompose = check_composition_property(communicator, properties, db_usage);
      }
    }

    if (decompose)
      return new Iocgns::ParallelDatabaseIO(nullptr, filename, db_usage, communicator, properties);
    else
#endif
      return new Iocgns::DatabaseIO(nullptr, filename, db_usage, communicator, properties);
  }
} // namespace Iocgns

#if defined(SEACAS_HAVE_MPI)
namespace {
  std::string check_decomposition_property(MPI_Comm comm, const Ioss::PropertyManager &properties,
                                           Ioss::DatabaseUsage db_usage)
  {
    std::string decomp_method;
    std::string decomp_property;
    if (db_usage == Ioss::READ_MODEL) {
      decomp_property = "MODEL_DECOMPOSITION_METHOD";
    }
    else if (db_usage == Ioss::READ_RESTART) {
      decomp_property = "RESTART_DECOMPOSITION_METHOD";
    }

    // Applies to either read_model or read_restart
    if (properties.exists("DECOMPOSITION_METHOD")) {
      std::string method = properties.get("DECOMPOSITION_METHOD").get_string();
      return Ioss::Utils::uppercase(method);
    }

    // Check for property...
    if (properties.exists(decomp_property)) {
      std::string method = properties.get(decomp_property).get_string();
      return Ioss::Utils::uppercase(method);
    }
    return decomp_method;
  }

  bool check_composition_property(MPI_Comm comm, const Ioss::PropertyManager &properties,
                                  Ioss::DatabaseUsage db_usage)
  {
    bool        compose          = true;
    std::string compose_property = "COMPOSE_INVALID";
    if (db_usage == Ioss::WRITE_RESULTS) {
      compose_property = "COMPOSE_RESULTS";
    }
    else if (db_usage == Ioss::WRITE_RESTART) {
      compose_property = "COMPOSE_RESTART";
    }

    Ioss::Utils::check_set_bool_property(properties, compose_property, compose);
    return compose;
  }
} // namespace
#endif
