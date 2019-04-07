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

#include <Ionit_Initializer.h>
#include <Ioss_CodeTypes.h>

#if defined(SEACAS_HAVE_EXODUS)
#include <exo_fac/Ioex_IOFactory.h>
#endif

#include <gen_struc/Iogs_DatabaseIO.h>
#include <generated/Iogn_DatabaseIO.h>
#include <heartbeat/Iohb_DatabaseIO.h>

#ifdef HAVE_SEACASIOSS_ADIOS2
#include <adios/Ioad_Initializer.h>
#endif

#if defined(SEACAS_HAVE_PAMGEN)
#include <pamgen/Iopg_DatabaseIO.h>
#endif

#if defined(SEACAS_HAVE_DATAWAREHOUSE)
#include <data_warehouse/Iodw_DatabaseIO.h>
#endif

#if defined(SEACAS_HAVE_CGNS)
#include <cgns/Iocgns_IOFactory.h>
#endif

#include <Ioss_ConcreteVariableType.h>
#include <Ioss_Initializer.h>
#include <transform/Iotr_Initializer.h>
#include <visualization/Iovs_IOFactory.h>

namespace {
#if defined(IOSS_THREADSAFE)
  std::mutex m_;
#endif
} // namespace

namespace Ioss {
  namespace Init {
    Initializer &Initializer::initialize_ioss()
    {
      static Initializer ionit;
      return ionit;
    }

    /** \brief Initialize the Ioss library.
     *
     *  Calls appropriate internal functions and methods to
     *  initialize the Ioss library. Initializes all database
     *  types.
     */
    Initializer::Initializer()
    {
      IOSS_FUNC_ENTER(m_);

#if defined(SEACAS_HAVE_EXODUS)
      Ioex::IOFactory::factory(); // Exodus
#endif
#if defined(SEACAS_HAVE_PAMGEN)
      Iopg::IOFactory::factory(); // Pamgen
#endif
#if defined(SEACAS_HAVE_DATAWAREHOUSE)
      Iodw::IOFactory::factory(); // DataWarehouse
#endif
#if defined(SEACAS_HAVE_CGNS)
      Iocgns::IOFactory::factory();
#endif

      Iovs::IOFactory::factory(); // Visualization
      Iohb::IOFactory::factory(); // HeartBeat
      Iogn::IOFactory::factory(); // Generated
      Iogs::IOFactory::factory(); // Structured Mesh Generator
      Ioss::StorageInitializer();
      Ioss::Initializer();
      Iotr::Initializer();
      #ifdef HAVE_SEACASIOSS_ADIOS2
      Ioad::Initializer(); // ADIOS2
      #endif
    }

    Initializer::~Initializer()
    {
      try {
        IOSS_FUNC_ENTER(m_);
        Ioss::IOFactory::clean();
        // Put code here that should run after sierra is finished
        // executing...
      }
      catch (...) {
      }
    }
  } // namespace Init
} // namespace Ioss
