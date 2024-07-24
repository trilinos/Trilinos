// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*!
\example lesson06_custom_operator.cpp
\brief Custom subclass of Tpetra::Operator that uses Tpetra::Import to
  handle MPI communication.

\ref Tpetra_Lesson06 explains this example in detail.
*/

#include <Tpetra_Core.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Operator.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE
#include <iostream>

// Define a class for our user-defined operator.
// In this case, it is the tridiagonal matrix [-1,2,-1].
// You may define it to be whatever you like.
//
// In general, Trilinos does NOT require the user to deal with MPI
// communication explicitly.  If you want to define your own operator
// though, there's no getting around it.  Fortunately, Trilinos makes
// this relatively straightforward with the use of Map and Import
// objects.  All you have to do is define your initial data
// distribution (which is a block row distribution here), and the data
// distribution you need to perform the operations of your
// matrix-vector multiply.  For instance, when performing a
// matrix-vector multiply with a tridiagonal matrix (with a block row
// distribution), each process needs to know the last element owned by
// the previous process and the first element owned by the next
// process.
//
// If you are only interested in running the code sequentially, you
// may safely ignore everything here regarding Map and Import objects.
class MyOp : public Tpetra::Operator<> {
public:
  // Tpetra::Operator subclasses should always define these four typedefs.
  typedef Tpetra::Operator<>::scalar_type scalar_type;
  typedef Tpetra::Operator<>::local_ordinal_type local_ordinal_type;
  typedef Tpetra::Operator<>::global_ordinal_type global_ordinal_type;
  typedef Tpetra::Operator<>::node_type node_type;

  // The type of the input and output arguments of apply().
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type,
                              global_ordinal_type, node_type> MV;
  // The Map specialization used by this class.
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

private:
  // This is an implementation detail; users don't need to see it.
  typedef Tpetra::Import<local_ordinal_type, global_ordinal_type,
                         node_type> import_type;

public:
  // Constructor
  //
  // n: Global number of rows and columns in the operator.
  // comm: The communicator over which to distribute those rows and columns.
  MyOp (const global_ordinal_type n,
        const Teuchos::RCP<const Teuchos::Comm<int> > comm)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::cout;
    using std::endl;

    TEUCHOS_TEST_FOR_EXCEPTION(
      comm.is_null (), std::invalid_argument,
      "MyOp constructor: The input Comm object must be nonnull.");

    //
    // Get the rank of this process and the number of processes
    // We're going to have to do something special with the first and last processes
    //
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();
    if (myRank == 0) {
      cout << "MyOp constructor" << endl;
    }

    //
    // Construct a map for our block row distribution
    //
    const global_ordinal_type indexBase = 0;
    opMap_ = rcp (new map_type (n, indexBase, comm));

    //
    // Get the local number of rows
    //
    local_ordinal_type nlocal = opMap_->getLocalNumElements ();

    //
    // Define the distribution that you need for the matvec.  When you
    // define this for your own operator, it is helpful to draw
    // pictures on a sheet of paper to keep track of who needs to
    // receive which entries of the source vector.  Here, each process
    // needs to receive one entry from each of its neighbors.
    //

    // All processes but the first will receive one element from the
    // previous process.
    if (myRank > 0) {
      ++nlocal;
    }
    // All processes but the last will receive one element from the
    // next process.
    if (myRank < numProcs - 1) {
      ++nlocal;
    }
    // Construct a list of columns where this process has nonzero
    // elements.  For our tridiagonal matrix, this is
    // firstRowItOwns-1:lastRowItOwns+1.
    std::vector<global_ordinal_type> indices;
    indices.reserve (nlocal);
    // The first process is a special case...
    if (myRank > 0) {
      indices.push_back (opMap_->getMinGlobalIndex () - 1);
    }
    for (global_ordinal_type i = opMap_->getMinGlobalIndex ();
         i <= opMap_->getMaxGlobalIndex (); ++i) {
      indices.push_back (i);
    }
    // So is the last process...
    if (myRank < numProcs - 1) {
      indices.push_back (opMap_->getMaxGlobalIndex () + 1);
    }

    // Wrap our vector in an array view, which is like a pointer
    Teuchos::ArrayView<const global_ordinal_type> elementList (indices);

    // Make a column Map for handling the redistribution.
    //
    // There will be some redundancies (i.e., some of the entries will
    // be owned by multiple processes).  Those redundancies will help
    // express the communication pattern for the sparse mat-vec.
    const global_ordinal_type numGlobalElements = n + 2*(numProcs - 1);
    redistMap_ = rcp (new map_type (numGlobalElements, elementList, indexBase, comm));

    // Make an Import object that describes how data will be
    // redistributed.  It takes a Map describing who owns what
    // originally, and a Map that describes who you WANT to own what.
    importer_= rcp (new import_type (opMap_, redistMap_));
  };

  //
  // These functions are required since we inherit from Tpetra::Operator
  //

  // Destructor
  virtual ~MyOp () {}

  // Get the domain Map of this Operator subclass.
  Teuchos::RCP<const map_type> getDomainMap() const { return opMap_; }

  // Get the range Map of this Operator subclass.
  Teuchos::RCP<const map_type> getRangeMap() const { return opMap_; }

  // Compute Y := alpha Op X + beta Y.
  //
  // We ignore the cases alpha != 1 and beta != 0 for simplicity.
  void
  apply (const MV& X,
         MV& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one (),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero ()) const
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::cout;
    using std::endl;
    typedef Teuchos::ScalarTraits<scalar_type> STS;

    RCP<const Teuchos::Comm<int> > comm = opMap_->getComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    if (myRank == 0) {
      cout << "MyOp::apply" << endl;
    }

    // We're writing the Operator subclass, so we are responsible for
    // error handling.  You can decide how much error checking you
    // want to do.  Just remember that checking things like Map
    // sameness or compatibility are expensive.
    TEUCHOS_TEST_FOR_EXCEPTION(
      X.getNumVectors () != Y.getNumVectors (), std::invalid_argument,
      "X and Y do not have the same numbers of vectors (columns).");

    // Let's make sure alpha is 1 and beta is 0...
    // This will throw an exception if that is not the case.
    TEUCHOS_TEST_FOR_EXCEPTION(
      alpha != STS::one() || beta != STS::zero(), std::logic_error,
      "MyOp::apply was given alpha != 1 or beta != 0. "
      "These cases are not implemented.");

    // Get the number of vectors (columns) in X (and Y).
    const size_t numVecs = X.getNumVectors ();

    // Make a temporary multivector for holding the redistributed
    // data.  You could also create this in the constructor and reuse
    // it across different apply() calls, but you would need to be
    // careful to reallocate if it has a different number of vectors
    // than X.  The number of vectors in X can vary across different
    // apply() calls.
    RCP<MV> redistData = rcp (new MV (redistMap_, numVecs));

    // Redistribute the data.
    // This will do all the necessary communication for you.
    // All processes now own enough data to do the matvec.
    redistData->doImport (X, *importer_, Tpetra::INSERT);

    // Get the number of local rows in X, on the calling process.
    const local_ordinal_type nlocRows =
      static_cast<local_ordinal_type> (X.getLocalLength ());

    // Perform the matvec with the data we now locally own.
    //
    // For each column...
    for (size_t c = 0; c < numVecs; ++c) {
      // Get a view of the desired column
      Teuchos::ArrayRCP<scalar_type> colView = redistData->getDataNonConst (c);

      local_ordinal_type offset;
      // Y[0,c] = -colView[0] + 2*colView[1] - colView[2] (using local indices)
      if (myRank > 0) {
        Y.replaceLocalValue (0, c, -colView[0] + 2*colView[1] - colView[2]);
        offset = 0;
      }
      // Y[0,c] = 2*colView[1] - colView[2] (using local indices)
      else {
        Y.replaceLocalValue (0, c, 2*colView[0] - colView[1]);
        offset = 1;
      }

      // Y[r,c] = -colView[r-offset] + 2*colView[r+1-offset] - colView[r+2-offset]
      for (local_ordinal_type r = 1; r < nlocRows - 1; ++r) {
        const scalar_type newVal =
          -colView[r-offset] + 2*colView[r+1-offset] - colView[r+2-offset];
        Y.replaceLocalValue (r, c, newVal);
      }
      // Y[nlocRows-1,c] = -colView[nlocRows-1-offset] + 2*colView[nlocRows-offset]
      //                   - colView[nlocRows+1-offset]
      if (myRank < numProcs - 1) {
        const scalar_type newVal =
          -colView[nlocRows-1-offset] + 2*colView[nlocRows-offset]
          - colView[nlocRows+1-offset];
        Y.replaceLocalValue (nlocRows-1, c, newVal);
      }
      // Y[nlocRows-1,c] = -colView[nlocRows-1-offset] + 2*colView[nlocRows-offset]
      else {
        const scalar_type newVal =
          -colView[nlocRows-1-offset] + 2*colView[nlocRows-offset];
        Y.replaceLocalValue (nlocRows-1, c, newVal);
      }
    }
  }

private:
  Teuchos::RCP<const map_type> opMap_, redistMap_;
  Teuchos::RCP<const import_type> importer_;
};


int
main (int argc, char *argv[])
{
  using Teuchos::RCP;
  using std::cout;
  using std::endl;

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  bool success = true;
  {
    auto comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();

    //
    // Get parameters from command-line processor
    //
    MyOp::global_ordinal_type n = 100;
    Teuchos::CommandLineProcessor cmdp (false, true);
    cmdp.setOption ("n", &n, "Number of rows of our operator.");
    if (cmdp.parse (argc, argv) !=
	Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return EXIT_FAILURE;
    }

    // Construct the operator.  Note that the operator does not have
    // to be an explicitly stored matrix.  Here, we are using our
    // user-defined operator.
    MyOp K (n, comm);

    // Construct a Vector of all ones, using the above Operator's
    // domain Map.
    typedef Tpetra::Vector<MyOp::scalar_type,
			   MyOp::local_ordinal_type,
			   MyOp::global_ordinal_type,
			   MyOp::node_type> vec_type;
    vec_type x (K.getDomainMap ());
    x.putScalar (1.0);
    // Construct an output Vector for K*x.
    vec_type y (K.getRangeMap ());
    K.apply (x, y); // Compute y := K*x.

    // The operator has a stencil (-1, 2, -1), except for the
    // boundaries.  At the left boundary (global row 0), the stencil is
    // (2, -1), and at the right boundary (global row n-1), the stencil
    // is (-1, 2).  Thus, we know that if all entries of the input
    // Vector are 1, then all entries of the output Vector are 0, except
    // for the boundary entries, which are both 1.
    //
    // To test this, construct the expected output vector y_expected,
    // and compare y to y_expected using the max norm.  Even in single
    // precision, the max norm of y - y_expected should be exactly zero.

    using map_type = MyOp::map_type;
    RCP<const map_type> rangeMap = K.getRangeMap ();

    vec_type y_expected (rangeMap);
    y_expected.putScalar (0.0);

    if (rangeMap->isNodeGlobalElement (0)) {
      y_expected.replaceGlobalValue (0, 1.0);
    }
    if (rangeMap->isNodeGlobalElement (n - 1)) {
      y_expected.replaceGlobalValue (n - 1, 1.0);
    }

    y_expected.update (1.0, y, -1.0); // y_expected := y - y_expected
    typedef vec_type::mag_type mag_type; // type of a norm of vec_type
    const mag_type diffMaxNorm = y_expected.normInf ();

    if (myRank == 0) {
      if (diffMaxNorm == 0.0) {
	// This tells the Trilinos test framework that the test passed.
	cout << "Yay!  ||y - y_expected||_inf = 0." << endl
	     << "End Result: TEST PASSED" << endl;
      }
      else {
	success = false;
	// This tells the Trilinos test framework that the test passed.
	cout << "Oops!  ||y - y_expected||_inf = " << diffMaxNorm
	     << " != 0." << endl
	     << "End Result: TEST FAILED" << endl;
      }
    }
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
