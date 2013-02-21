#include <stk_topology/topology.hpp>
#include <stk_io/IossBridge.hpp>
#include <Ioss_Initializer.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_ConcreteVariableType.h>

using namespace Ioss;
#define OUTPUT std::cerr

// ========================================================================
static int convert_ioss_to_stk_topology();
static int convert_stk_to_ioss_topology();
// ========================================================================


// TODO: Check that stk::topology and Ioss::ElementTopology give similar results
//       for all queries (num_nodes, ...)

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  Ioss::StorageInitializer initialize_storage;
  Ioss::Initializer  initialize_topologies;

  int err_count = convert_ioss_to_stk_topology();
  err_count += convert_stk_to_ioss_topology();

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
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

  NameList topologies;
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
