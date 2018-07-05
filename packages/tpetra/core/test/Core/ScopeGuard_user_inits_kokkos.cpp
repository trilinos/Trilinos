#include <cstdlib>
#include <iostream>
#include "Tpetra_Core.hpp"
#include "Kokkos_Core.hpp"

namespace { // (anonymous)

// NOTE TO TEST AUTHORS: The code that calls this function captures
// std::cerr, so don't write to std::cerr on purpose in this function.
void testMain (bool& success, int argc, char* argv[])
{
  using std::cout;
  using std::endl;

  if (Kokkos::is_initialized ()) {
    success = false;
    cout << "Kokkos::is_initialized() is true, "
      "even before Kokkos::initialize was called." << endl;
    return;
  }
  Kokkos::initialize (argc, argv);
  if (! Kokkos::is_initialized ()) {
    success = false;
    cout << "Kokkos::is_initialized() is false, "
      "even after Kokkos::initialize was called." << endl;
    return;
  }

  // In this example, the "user" has called Kokkos::initialize before
  // Tpetra::ScopeGuard's constructor is called.  ScopeGuard's
  // constructor must not try to call it again.
  {
    Tpetra::ScopeGuard tpetraScope (&argc, &argv);
    if (! Tpetra::isInitialized ()) {
      success = false;
      cout << "Tpetra::isInitialized() is false, "
        "even after Tpetra::ScopeGuard::ScopeGuard was called." << endl;
      return;
    }
    if (! Kokkos::is_initialized ()) {
      success = false;
      cout << "Kokkos::is_initialized() is false, even after "
        "Kokkos::initialize and Tpetra::ScopeGuard::ScopeGuard "
        "were called." << endl;
      // Let the program keep going, so MPI (if applicable) gets finalized.
    }

    cout << "About to leave Tpetra scope" << endl;
  }
  cout << "Left Tpetra scope" << endl;
  if (Tpetra::isInitialized ()) {
    success = false;
    cout << "Tpetra::isInitialized() is true, "
      "even after Tpetra::ScopeGuard::~ScopeGuard was called." << endl;
    // Let the program keep going, so Kokkos gets finalized.
  }

  // Since the "user" is responsible for calling Kokkos::finalize,
  // Tpetra::ScopeGuard's destructor should NOT have called
  // Kokkos::finalize.
  if (! Kokkos::is_initialized ()) {
    success = false;
    cout << "Kokkos::is_initialized() is false, "
      "after Tpetra::ScopeGuard::~ScopeGuard was called." << endl;
  }
  Kokkos::finalize ();
  if (Kokkos::is_initialized ()) {
    success = false;
    cout << "Kokkos::is_initialized() is true, "
      "even after Kokkos::initialize was called." << endl;
  }
}

class CaptureOstream {
public:
  CaptureOstream (std::ostream& stream) :
    originalStream_ (stream),
    originalBuffer_ (stream.rdbuf ())
  {
    originalStream_.rdbuf (tempStream_.rdbuf ());
  }

  std::string getCapturedOutput () const {
    return tempStream_.str ();
  }

  ~CaptureOstream () {
    originalStream_.rdbuf (originalBuffer_);
  }
private:
  std::ostream& originalStream_;
  std::ostringstream tempStream_;
  using buf_ptr_type = decltype (originalStream_.rdbuf ());
  buf_ptr_type originalBuffer_;
};

} // namespace (anonymous)

int main (int argc, char* argv[])
{
  using std::cout;
  using std::endl;

  bool success = true;
  {
    // Capture std::cerr output, so we can tell if Tpetra::initialize
    // printed a warning message.
    CaptureOstream captureCerr (std::cerr);
    testMain (success, argc, argv);
    const std::string capturedOutput = captureCerr.getCapturedOutput ();
    cout << "Captured output: " << capturedOutput << endl;
    if (capturedOutput.size () != 0) {
      success = false; // should NOT have printed in this case
      cout << "Captured output is empty!" << endl;
    }
  }

  cout << "End Result: TEST " << (success ? "PASSED" : "FAILED") << endl;
  return EXIT_SUCCESS;
}
