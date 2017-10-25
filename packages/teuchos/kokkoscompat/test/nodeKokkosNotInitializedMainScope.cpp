#include "KokkosCompat_ClassicNodeAPI_Wrapper.hpp"
#include "Kokkos_Core.hpp"
#include <iostream>
#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE
#ifdef KOKKOS_ENABLE_SERIAL
#  include <type_traits>
#endif // KOKKOS_ENABLE_SERIAL

// Test that Node initializes Kokkos if it is not already initialized.
// It must also finalize Kokkos at exit, but the only portable way to
// test that would be to run Valgrind and ensure no memory leaks.

namespace { // (anonymous)

// Is Kokkos initialized?
//
// https://github.com/kokkos/kokkos/issues/1184 means that this
// function doesn't return something meaningful with Kokkos::Serial.
bool isKokkosInitialized ()
{
  // When https://github.com/kokkos/kokkos/issues/1060 is fixed, use
  // Kokkos::is_initialized here, instead of asking the default
  // execution space whether it is initialized.
  typedef Kokkos::DefaultExecutionSpace::execution_space
    execution_space;
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
#ifdef KOKKOS_ENABLE_SERIAL
  constexpr bool defaultExecSpaceIsSerial =
    std::is_same<execution_space, Kokkos::Serial>::value;
#else
  constexpr bool defaultExecSpaceIsSerial = false;
#endif // KOKKOS_ENABLE_SERIAL

  if (! defaultExecSpaceIsSerial) {
    // https://github.com/kokkos/kokkos/issues/1184 means that this
    // function doesn't return something meaningful with
    // Kokkos::Serial.
    if (isKokkosInitialized ()) {
      cerr << "FAILED: "
        "Before calling Node's constructor, "
        "and without having called Kokkos::initialize ourselves, "
        "Kokkos claims to be initialized!" << endl;
      return EXIT_FAILURE;
    }
  }

  // Create a node instance.  Make sure that this does not make Kokkos
  // raise an exception, e.g., due to double initialization.  Do this
  // at main() scope, so we can test that case.  We'll have another
  // test for where node instances are always in an inner scope.
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

  // Don't print PASSED!  That will make the test think that it
  // passed, even though things are still happening (due to atexit)
  // after main() exits.  We use FAIL_REGULAR_EXPRESSION (see
  // CMakeLists.txt) so that the test fails if and only if it prints
  // "FAILED:".
  return EXIT_SUCCESS;
}

