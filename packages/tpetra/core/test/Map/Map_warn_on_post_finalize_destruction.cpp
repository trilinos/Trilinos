// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <cstdlib>
#include <iostream>
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"

namespace { // (anonymous)

// NOTE TO TEST AUTHORS: The code that calls this function captures
// std::cerr, so don't write to std::cerr on purpose in this function.
void testMain (bool& success, int argc, char* argv[])
{
  using std::cout;
  using std::endl;

  if (Tpetra::isInitialized ()) {
    success = false;
    cout << "Hm, Tpetra::isInitialized() returns true "
      "even before Tpetra::initialize() has been called" << endl;
    return;
  }
  
  Tpetra::initialize (&argc, &argv);

  auto comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  if (myRank == 0) {
    cout << "Tpetra::initialize has been called; about to create Map" << endl;
  }

  // NOTE: This is INCORRECT BEHAVIOR, because the Map outlives
  // Tpetra::finalize.  DO NOT WRITE CODE LIKE THIS!  The point is to
  // test whether Map's destructor prints an error message per GitHub
  // Issue #2372.  Map's destructor must print to std::cerr in this
  // case; the code that calls testMain will capture std::cerr output
  // in order to check correct behavior.
  Tpetra::Map<> map (numProcs, 0, comm);

  if (myRank == 0) {
    cout << "Created Map; about to call Tpetra::finalize" << endl;
  }

  Tpetra::finalize ();
  if (myRank == 0) {
    cout << "Called Tpetra::finalize" << endl;
  }

  if (Tpetra::isInitialized ()) {
    success = false;
    if (myRank == 0) {
      cout << "Hm, Tpetra::isInitialized() returns true "
	"even after Tpetra::finalize() has been called" << endl;
    }
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
    // Capture std::cerr output, so we can tell if Map's destructor
    // printed a warning message.
    CaptureOstream captureCerr (std::cerr);
    testMain (success, argc, argv);
    const std::string capturedOutput = captureCerr.getCapturedOutput ();
    cout << "Captured output: " << capturedOutput << endl;
    if (capturedOutput.size () == 0) {
      success = false; // SHOULD have printed in this case
      cout << "Captured output is empty!" << endl;
    }
  }
  
  cout << "End Result: TEST " << (success ? "PASSED" : "FAILED") << endl;
  return EXIT_SUCCESS;
}
