#include "KokkosCompat_ClassicNodeAPI_Wrapper.hpp"
#include "Kokkos_Core.hpp"
#include <iostream>
#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE

// Test that Node does not initializes Kokkos if it is already
// initialized, and does not attempt to double-finalize Kokkos
// in this case.

namespace { // (anonymous)

// Is Kokkos initialized?
bool isKokkosInitialized ()
{
  // When https://github.com/kokkos/kokkos/issues/1060 is fixed, use
  // Kokkos::is_initialized here, instead of asking the default
  // execution space whether it is initialized.
  typedef Kokkos::DefaultExecutionSpace::execution_space
    execution_space;
  // The merge of Trilinos' Pull Request #1911 into Trilinos:develop
  // fixes Kokkos issue https://github.com/kokkos/kokkos/issues/1184.
  // Thus, this function should now return meaningful results, even if
  // Kokkos::Serial is the default execution space.
  return execution_space::is_initialized ();
}

} // namespace (anonymous)

int
main (int argc, char* argv[])
{
  using std::cerr;
  using std::endl;
  typedef Kokkos::DefaultExecutionSpace execution_space;
  typedef Kokkos::Compat::KokkosDeviceWrapperNode<execution_space>
    node_type;

  // The merge of Trilinos' Pull Request #1911 into Trilinos:develop
  // fixes Kokkos issue https://github.com/kokkos/kokkos/issues/1184.
  // Thus, Kokkos::Serial::is_initialized() should return false if
  // Kokkos has not yet been initialized.
  if (isKokkosInitialized ()) {
    cerr << "FAILED: Before calling Kokkos::initialize, "
      "Kokkos reports that it has been initialized!" << endl;
    return EXIT_FAILURE;
  }
  Kokkos::initialize (argc, argv);
  if (! isKokkosInitialized ()) {
    cerr << "FAILED: After calling Kokkos::initialize, "
      "Kokkos is still not initialized!" << endl;
    return EXIT_FAILURE;
  }

  // Create a node instance.  Make sure that this does not make Kokkos
  // raise an exception, e.g., due to double initialization.
  node_type node1;

  if (! isKokkosInitialized ()) {
    cerr << "FAILED: After calling Kokkos::initialize, "
      "and after calling Node's constructor once, "
      "Kokkos is still not initialized!" << endl;
    return EXIT_FAILURE;
  }
  {
    node_type node2; // inner scope, so its destructor gets invoked
    if (! isKokkosInitialized ()) {
      cerr << "FAILED: After calling Kokkos::initialize, "
        "and after calling Node's constructor twice, "
        "Kokkos is still not initialized!" << endl;
      return EXIT_FAILURE;
    }
  }
  if (! isKokkosInitialized ()) {
    cerr << "FAILED: After calling Kokkos::initialize, "
      "after calling Node's constructor twice, "
      "and after calling Node's destructor once, "
      "Kokkos is still not initialized!" << endl;
    return EXIT_FAILURE;
  }

  Kokkos::finalize ();
  // Node's destructor will be called after Kokkos::finalize.  We thus
  // implicitly test whether this makes Kokkos raise an exception.
  // Don't print PASSED!  That will make the test think that it
  // passed.  We use FAIL_REGULAR_EXPRESSION (see CMakeLists.txt) so
  // that the test fails if and only if it prints "FAILED:".
  return EXIT_SUCCESS;
}

