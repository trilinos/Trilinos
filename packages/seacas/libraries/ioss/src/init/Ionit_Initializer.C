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

#include <init/Ionit_Initializer.h>

#if !defined(NO_EXODUS_SUPPORT)
#include <exodusII/Ioex_IOFactory.h>
#if defined(HAVE_MPI) && !defined(NO_DOF_EXODUS_SUPPORT)
#include <par_exo/Iopx_IOFactory.h>
#endif
#endif
#include <heartbeat/Iohb_DatabaseIO.h>
#include <generated/Iogn_DatabaseIO.h>
#if !defined(NO_PAMGEN_SUPPORT)
#include <pamgen/Iopg_DatabaseIO.h>
#endif

// with introduction of paraview sierra catalyst plugin, the Iovs stuff is
// always included and NO_PARAVIEWMESH_SUPPORT is never defined.  With the
// plugin architecture, there is no overhead for sierra when the plugin is
// not loaded.  The #define test is left here for now in case developers
// need to use it.

// NOTE: (gdsjaar) -- Do *not* remove the NO_PARAVIEWIMESH_SUPPORT define.
//                    The Ioss is used in more products than just Sierra,
//                    so we cannot always rely on the paraview catalyst
//                    plugin being available.

#if !defined(NO_PARAVIEWIMESH_SUPPORT)
#include <visualization/Iovs_IOFactory.h>
#endif

#include <Ioss_ConcreteVariableType.h>
#include <Ioss_Initializer.h>
#include <transform/Iotr_Initializer.h>

namespace Ioss {
  namespace Init {
    Initializer::Initializer()
    {
#if !defined(NO_EXODUS_SUPPORT)
      Ioex::IOFactory::factory();    // ExodusII
#if defined(HAVE_MPI) && !defined(NO_DOF_EXODUS_SUPPORT)
      Iopx::IOFactory::factory();    // ExodusII
#endif
#endif
      Iohb::IOFactory::factory();   // HeartBeat
      Iogn::IOFactory::factory();  // Generated
#if !defined(NO_PAMGEN_SUPPORT)
      Iopg::IOFactory::factory(); // Pamgen
#endif
#if !defined(NO_PARAVIEWIMESH_SUPPORT)
      Iovs::IOFactory::factory(); // Visualization
#endif
      
      Ioss::StorageInitializer();
      Ioss::Initializer();
      Iotr::Initializer();
    }

    Initializer::~Initializer()
    {
      try {
	Ioss::IOFactory::clean();
	// Put code here that should run after sierra is finished
	// executing...
      } catch (...) {}
    }
  }
}
