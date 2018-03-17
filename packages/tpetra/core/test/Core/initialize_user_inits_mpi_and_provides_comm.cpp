#include <cstdlib>
#include <iostream>
#include "Tpetra_Core.hpp"

#if ! defined(HAVE_TPETRACORE_MPI)
#  error "Building and testing this example requires MPI."
#endif // ! defined(HAVE_TPETRACORE_MPI)
#include "mpi.h"
#include "Tpetra_Details_extractMpiCommFromTeuchos.hpp"

namespace { // (anonymous)

bool isMpiInitialized ()
{
  int mpiInitializedInt = 0;
  (void) MPI_Initialized (&mpiInitializedInt);
  return mpiInitializedInt != 0;
}

int getRankInComm (MPI_Comm comm)
{
  int myRank = 0;
  (void) MPI_Comm_rank (comm, &myRank);
  return myRank;
}
  
bool allTrueInComm (const bool lclTruth, MPI_Comm comm)
{
  int lclTruthInt = lclTruth ? 1 : 0;
  int gblTruthInt = 0;
  (void) MPI_Allreduce (&lclTruthInt, &gblTruthInt, 1,
			MPI_INT, MPI_MIN, comm);
  return gblTruthInt != 0;
}

bool
tpetraCommIsLocallyLegit (const Teuchos::Comm<int>* wrappedTpetraComm,
			  MPI_Comm expectedMpiComm)
{
  if (wrappedTpetraComm == nullptr) {
    return false;
  }
  MPI_Comm tpetraComm;
  try {
    using Tpetra::Details::extractMpiCommFromTeuchos;
    tpetraComm = extractMpiCommFromTeuchos (*wrappedTpetraComm);
  }
  catch (...) {
    return false;
  }
  if (tpetraComm == MPI_COMM_NULL) {
    return false;
  }
  int result = MPI_UNEQUAL;
  (void) MPI_Comm_compare (expectedMpiComm, tpetraComm, &result);
  // Tpetra reserves the right to MPI_Comm_dup the input comm.
  return result == MPI_IDENT || result == MPI_CONGRUENT;
}
  

// NOTE TO TEST AUTHORS: The code that calls this function captures
// std::cerr, so don't write to std::cerr on purpose in this function.
void testMain (bool& success, int argc, char* argv[])
{
  using std::cout;
  using std::endl;

  if (isMpiInitialized ()) {
    success = false;
    cout << "MPI_Initialized claims MPI is initialized, "
      "before MPI_Init was called" << endl;
    return;
  }
  (void) MPI_Init (&argc, &argv);
  if (! isMpiInitialized ()) {
    success = false;
    cout << "MPI_Initialized claims MPI is not initialized, "
      "even after MPI_Init was called" << endl;
    return;
  }

  // Split off Process 0 into a comm by itself.  Only invoke
  // Tpetra::initialize and Tpetra::finalize on the remaining
  // processes.
  const int color = (myRank == 0) ? 0 : 1;  
  MPI_Comm splitComm;
  {
    const int key = 0;
    const int errCode =
      MPI_Comm_split (MPI_COMM_WORLD, color, key, &splitComm);
    if (errCode != MPI_SUCCESS) {
      success = false;
      cout << "MPI_Comm_split failed!" << endl;
      (void) MPI_Abort (MPI_COMM_WORLD, EXIT_FAILURE);
    }
  }

  if (color == 0) { // this subcomm doesn't participate
    MPI_Comm_free (&splitComm);
    MPI_Finalize ();
    return;
  }

  Tpetra::initialize (&argc, &argv, splitComm);
  if (! isMpiInitialized ()) {
    success = false;
    cout << "MPI_Initialized claims MPI was not initialized, "
      "even after MPI_Init and Tpetra::initialize were called" << endl;
    Tpetra::finalize (); // just for completeness
    return;
  }

  // MPI is initialized, so we can check whether all processes report
  // Tpetra as initialized.
  const bool tpetraIsNowInitialized =
    allTrueInComm (Tpetra::isInitialized (), splitComm);
  if (! tpetraIsNowInitialized) {
    success = false;
    if (getRankInComm (splitComm) == 0) {
      cout << "Tpetra::isInitialized() is false on at least one process"
	", even after Tpetra::initialize has been called." << endl;
    }
    (void) MPI_Abort (MPI_COMM_WORLD, EXIT_FAILURE);
  }

  auto comm = Tpetra::getDefaultComm ();
  const bool tpetraCommGloballyValid =
    allTrueInComm (tpetraCommIsLocallyLegit (comm.get (), splitComm));
  if (! tpetraCommGloballyValid) {
    success = false;
    if (getRankInComm (splitComm) == 0) {
      cout << "Tpetra::getDefaultComm() returns an invalid comm "
	"on at least one process." << endl;
    }
  }

  const int myTpetraRank = comm.is_null () ? 0 : comm->getRank ();
  const bool ranksSame =
    allTrueInCommWorld (getRankInComm (splitComm) == myTpetraRank);
  if (! ranksSame) {
    success = false;
    if (myRank == 0) {
      cout << "MPI rank does not match Tpetra rank "
	"on at least one process" << endl;
    }
  }

  if (myRank == 0) {
    cout << "About to call Tpetra::finalize" << endl;
  }
  Tpetra::finalize ();
  if (myRank == 0) {
    cout << "Called Tpetra::finalize" << endl;
  }
  // Since the "user" is responsible for calling MPI_Finalize,
  // Tpetra::finalize should NOT have called MPI_Finalize.
  if (! isMpiInitialized ()) {
    success = false;
    cout << "Tpetra::finalize() seems to have called MPI_Finalize, "
      "even though the user was responsible for initializing and "
      "finalizing MPI." << endl;
    return;
  }

  // MPI is still initialized, so we can check whether processes are
  // consistent.
  const bool tpetraGloballyFinalized =
    allTrueInCommWorld (! Tpetra::isInitialized ());
  if (! tpetraGloballyFinalized) {
    success = false;
    if (myRank == 0) {
      cout << "Tpetra::isInitialized() returns true on some process, "
	"even after Tpetra::finalize() has been called" << endl;
    }
  }

  (void) MPI_Finalize ();
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
