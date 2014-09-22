// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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

#include <Ioss_ConcreteVariableType.h>  // for StorageInitializer
#include <Ioss_ElementTopology.h>       // for ElementTopology
#include <Ioss_Initializer.h>           // for Initializer
#include <stdlib.h>                     // for NULL, EXIT_FAILURE, etc
#include <iomanip>                      // for operator<<, setw
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_io/IossBridge.hpp>        // for map_ioss_topology_to_stk
#include <stk_topology/topology.hpp>    // for topology, operator++, etc
#include <string>                       // for operator<<, operator!=, etc
#include "Ioss_CodeTypes.h"             // for HAVE_MPI
#include "Ioss_VariableType.h"          // for NameList
#include <stk_util/parallel/Parallel.hpp>

#define OUTPUT std::cerr

// ========================================================================
static int convert_ioss_to_stk_topology();
static int convert_stk_to_ioss_topology();
// ========================================================================


// TODO: Check that stk::topology and Ioss::ElementTopology give similar results
//       for all queries (num_nodes, ...)

int main(int argc, char *argv[])
{
  stk::parallel_machine_init(&argc, &argv);

  Ioss::StorageInitializer initialize_storage;
  Ioss::Initializer  initialize_topologies;

  int err_count = convert_ioss_to_stk_topology();
  err_count += convert_stk_to_ioss_topology();

  stk::parallel_machine_finalize();
  OUTPUT << "\n" << argv[0];;
  if (err_count == 0) {
    OUTPUT << "\nSIERRA execution successful." << std::endl;
    return EXIT_SUCCESS;
  } else {
    OUTPUT << "\nSIERRA execution failed." << std::endl;
    return EXIT_FAILURE;
  }
}

int convert_ioss_to_stk_topology()
{
  int err_count = 0;

  Ioss::NameList topologies;
  int topology_count = Ioss::ElementTopology::describe(&topologies);

  OUTPUT.setf(std::ios::left);
  for (int i=0; i < topology_count; i++) {
    Ioss::ElementTopology *topo = Ioss::ElementTopology::factory(topologies[i], false);
    if (topologies[i] != topo->name())
      continue; // Alias

    OUTPUT << "Testing ioss topology: " << std::setw(20) << topologies[i] << "\n";
    Ioss::ElementTopology *ioss_topo = Ioss::ElementTopology::factory(topologies[i], false);
    stk::topology stk_topo = stk::io::map_ioss_topology_to_stk(ioss_topo);

    if (stk_topo == stk::topology::INVALID_TOPOLOGY && topologies[i] != "unknown") {
      OUTPUT << "ERROR: IOSS topology '" << topologies[i] << "' could not be converted to STK topology.\n";
      err_count++;
      continue;
    }

    // Convert back to Ioss::Topology and see if we get the same type...
    Ioss::ElementTopology *new_topo = Ioss::ElementTopology::factory(stk_topo.name(), true);
    if (new_topo == NULL) {
      OUTPUT << "ERROR: STK Topology '" << stk_topo.name() << "', created from IOSS topology '" << topologies[i]
             << "' could not be converted back to IOSS topology.\n";
      err_count++;
      continue;
    }
    if (new_topo->name() != ioss_topo->name()) {
      if (new_topo->name() == "edge2" || new_topo->name() == "edge3") {
	OUTPUT << "ERROR: Mismatch in topology names. Expected '" << ioss_topo->name()
	       << "' Got '" << new_topo->name() << "' (OK FOR NOW)\n";
      } else {
	OUTPUT << "ERROR: Mismatch in topology names. Expected '" << ioss_topo->name()
	       << "' Got '" << new_topo->name() << "'\n";
	err_count++;
      }
    }
  }
  return err_count;
}

int convert_stk_to_ioss_topology()
{
  int err_count = 0;

  for (stk::topology topo = stk::topology::BEGIN_TOPOLOGY; topo < stk::topology::END_TOPOLOGY; ++topo) {
    OUTPUT << "Testing stk topology: " << std::setw(20) << topo.name() << "\n";

    Ioss::ElementTopology *ioss_topo = Ioss::ElementTopology::factory(topo.name(), true);
    if (ioss_topo == NULL) {
      OUTPUT << "ERROR: STK Topology '" << topo.name() << "' could not be converted to IOSS topology.\n";
      err_count++;
      continue;
    }

    // See if get the same type back...
    stk::topology stk_topo = stk::io::map_ioss_topology_to_stk(ioss_topo);
    if (stk_topo == stk::topology::INVALID_TOPOLOGY && ioss_topo->name() != "unknown") {
      OUTPUT << "ERROR: IOSS topology '" << ioss_topo->name() << "' created from stk topology '"
          << topo.name() << "' could not be converted to back STK topology.\n";
      err_count++;
    }
  }
  return err_count;
}
