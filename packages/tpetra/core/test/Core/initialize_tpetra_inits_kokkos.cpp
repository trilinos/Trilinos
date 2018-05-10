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

  // In this example, Tpetra::initialize is responsible for calling
  // Kokkos::initialize and Kokkos::finalize.
  if (Kokkos::is_initialized ()) {
    success = false;
    cout << "Kokkos::is_initialized() is true, "
      "before Tpetra::initialize was called." << endl;
    return;
  }
  Tpetra::initialize (&argc, &argv);

  if (! Kokkos::is_initialized ()) {
    success = false;
    cout << "Kokkos::is_initialized() is false, "
      "after Tpetra::initialize was called." << endl;
  }
  if (! Tpetra::isInitialized ()) {
    success = false;
    cout << "Tpetra::isInitialized() is false, "
      "even after Tpetra::initialize was called." << endl;
  }

  auto comm = Tpetra::getDefaultComm ();
  if (comm.is_null ()) {
    success = false;
    cout << "Tpetra::getDefaultComm() is null." << endl;
  }

  cout << "About to call Tpetra::finalize." << endl;
  Tpetra::finalize ();
  cout << "Called Tpetra::finalize." << endl;

  // Kokkos is like Tpetra; Kokkos::is_initialized() means "was
  // initialized and was not finalized."  That differs from MPI, where
  // MPI_Initialized only refers to MPI_Init and MPI_Finalized only
  // refers to MPI_Finalize.
  if (Kokkos::is_initialized ()) {
    success = false;
    cout << "Tpetra::finalize did not call Kokkos::finalize." << endl;
    return;
  }

  // MPI is no longer initialized, so we can't all-reduce on this.
  if (Tpetra::isInitialized ()) {
    success = false;
    cout << "Tpetra::isInitialized() is true, "
      "even after Tpetra::finalize has been called" << endl;
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
