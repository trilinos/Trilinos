#include <string>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Tpetra_DefaultPlatform.hpp>
#include <Epetra_CrsMatrix.h>

#include "Amesos2_EpetraCrsMatrixAdapter.hpp"
#include "Amesos2_Util_is_same.hpp"

namespace {

using std::cout;
using std::endl;
using std::string;

using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::rcp;
using Teuchos::Comm;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::tuple;
using Teuchos::ScalarTraits;
using Teuchos::OrdinalTraits;
using Teuchos::FancyOStream;
using Teuchos::VerboseObjectBase;

using Amesos::MatrixAdapter;

using Amesos::Util::is_same;

typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;

bool testMpi = false;

// Where to look for input files
string filedir;

TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.setOption("filedir",&filedir,"Directory of matrix files.");
  clp.addOutputSetupOptions(true);
  clp.setOption("test-mpi", "test-serial", &testMpi,
    "Test Serial by default (for now) or force MPI test.  In a serial build,"
    " this option is ignored and a serial comm is always used." );
}

const Epetra_Comm& getDefaultComm()
{
#ifdef EPETRA_MPI
  return Epetra_MpiComm( MPI_COMM_WORLD );
#else
  return Epetra_SerialComm();
#endif
}

RCP<FancyOStream> getDefaultOStream()
{
  return( VerboseObjectBase::getDefaultOStream() );
}


/*
 * UNIT TESTS
 */
TEUCHOS_UNIT_TEST( CrsMatrixAdapter, Initialization )
{
  /* Test correct initialization of the MatrixAdapter
   *
   * - All Constructors
   * - Correct initialization of class members
   * - Correct typedefs ( using Amesos::is_same<> )
   */
  typedef Epetra_CrsMatrix MAT;
  typedef MatrixAdapter<MAT> ADAPT;

  Epetra_SerialComm comm;
  //const size_t numprocs = comm.NumProc();
  //  const int rank     = comm.MyPID();
  // create a Map
  const int num_eqn = 100;

  Epetra_Map Map(num_eqn, 0, comm);
  
  // Get update list and number of local equations from newly created Map.

  int NumMyElements = Map.NumMyElements();

  std::vector<int> MyGlobalElements(NumMyElements);
  Map.MyGlobalElements(&MyGlobalElements[0]);

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
  // on this processor

  std::vector<int> NumNz(NumMyElements);

  for (int i = 0; i < NumMyElements; ++i){
    NumNz[i] = 1;
  }

  // Create a Epetra_Matrix
  Teuchos::RCP<Epetra_CrsMatrix> eye = rcp(new Epetra_CrsMatrix(Copy, Map, &NumNz[0]));
  
  std::vector<double> Values(2);
  std::vector<int> Indices(2);
  double one = 1.0;
  
  for (int i = 0; i < NumMyElements; ++i)
  {
    int ierr = eye->InsertGlobalValues(MyGlobalElements[i], 1, &one, &MyGlobalElements[i]);
    TEST_FOR_EXCEPTION(ierr != 0,std::runtime_error,"Error inserting value into matrix");
  }

  eye->FillComplete();
  // Constructor from RCP
  RCP<ADAPT> adapter  = rcp(new MatrixAdapter<MAT>( eye ));
  // Copy constructor
  RCP<ADAPT> adapter2 = rcp(new MatrixAdapter<MAT>(*adapter));

  // The following should all pass at compile time
  TEST_ASSERT( (is_same<double,ADAPT::scalar_type>::value) );
  TEST_ASSERT( (is_same<int,ADAPT::local_ordinal_type>::value) );
  TEST_ASSERT( (is_same<int,ADAPT::global_ordinal_type>::value) );
  TEST_ASSERT( (is_same<size_t,ADAPT::global_size_type>::value) );
  TEST_ASSERT( (is_same<MAT,ADAPT::matrix_type>::value) );

}

TEUCHOS_UNIT_TEST( CrsMatrixAdapter, Dimensions )
{
  // Test that the dimensions reported by the adapter match those as reported
  // by the Tpetra::CrsMatrix
  // Check equality of mapped method calls
  typedef Epetra_CrsMatrix MAT;
  typedef MatrixAdapter<MAT> ADAPT;

  Epetra_SerialComm comm;
  //const size_t numprocs = comm.NumProc();
  //  const int rank     = comm.MyPID();
  // create a Map
  const int num_eqn = 100;

  Epetra_Map Map(num_eqn, 0, comm);
  
  // Get update list and number of local equations from newly created Map.

  int NumMyElements = Map.NumMyElements();

  std::vector<int> MyGlobalElements(NumMyElements);
  Map.MyGlobalElements(&MyGlobalElements[0]);

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
  // on this processor

  std::vector<int> NumNz(NumMyElements);

  for (int i = 0; i < NumMyElements; ++i){
    NumNz[i] = 1;
  }

  // Create a Epetra_Matrix
  Teuchos::RCP<Epetra_CrsMatrix> eye = rcp(new Epetra_CrsMatrix(Copy, Map, &NumNz[0]));
  
  std::vector<double> Values(2);
  std::vector<int> Indices(2);
  double one = 1.0;
  
  for (int i = 0; i < NumMyElements; ++i)
  {
    int ierr = eye->InsertGlobalValues(MyGlobalElements[i], 1, &one, &MyGlobalElements[i]);
    TEST_FOR_EXCEPTION(ierr != 0,std::runtime_error,"Error inserting value into matrix");
  }

  eye->FillComplete();
  // Constructor from RCP
  RCP<ADAPT> adapter  = rcp(new MatrixAdapter<MAT>( eye ));

  TEST_EQUALITY(Teuchos::as<ADAPT::global_size_type>(eye->NumGlobalNonzeros()), adapter->getGlobalNNZ());
  TEST_EQUALITY(Teuchos::as<size_t>(eye->NumMyNonzeros()), adapter->getLocalNNZ());
  TEST_EQUALITY(Teuchos::as<ADAPT::global_size_type>(eye->NumGlobalRows()), adapter->getGlobalNumRows());
  TEST_EQUALITY(Teuchos::as<size_t>(eye->NumMyRows()), adapter->getLocalNumRows());
  TEST_EQUALITY(Teuchos::as<size_t>(eye->NumMyCols()), adapter->getLocalNumCols());
  TEST_EQUALITY(Teuchos::as<ADAPT::global_size_type>(eye->NumGlobalCols()), adapter->getGlobalNumCols());
  TEST_EQUALITY(Teuchos::as<size_t>(eye->MaxNumEntries()), adapter->getMaxNNZ());
}


TEUCHOS_UNIT_TEST( CrsMatrixAdapter, CRS )
{
  /* Test the getCrs() method of MatrixAdapter.  We check against a simple
   * test matrix that we construct on the fly.
   */
  typedef Epetra_CrsMatrix MAT;
  typedef MatrixAdapter<MAT> ADAPT;

  /* We will be using the following matrix for this test:
   * 
   * [ [ 7,  0,  -3, 0,  -1, 0 ]
   *   [ 2,  8,  0,  0,  0,  0 ]
   *   [ 0,  0,  1,  0,  0,  0 ]
   *   [ -3, 0,  0,  5,  0,  0 ]
   *   [ 0,  -1, 0,  0,  4,  0 ]
   *   [ 0,  0,  0,  -2, 0,  6 ] ]
   */
  Epetra_SerialComm comm;

  // create a Map
  const int num_eqn = 6;

  Epetra_Map Map(num_eqn, 0, comm);
  
  // Get update list and number of local equations from newly created Map.

  int NumMyElements = Map.NumMyElements();

  std::vector<int> MyGlobalElements(NumMyElements);
  Map.MyGlobalElements(&MyGlobalElements[0]);

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
  // on this processor

  int NumNz[] = {3, 2, 1, 2, 2, 2};

   // Create a Epetra_Matrix
  Teuchos::RCP<Epetra_CrsMatrix> mat = rcp(new Epetra_CrsMatrix(Copy, Map, NumNz));
  
  // Construct matrix
  mat->InsertGlobalValues(0,NumNz[0],
    tuple<ADAPT::scalar_type>(7,-3,-1).getRawPtr(),
    tuple<ADAPT::global_ordinal_type>(0,2,4).getRawPtr());
  mat->InsertGlobalValues(1,NumNz[1],
    tuple<ADAPT::scalar_type>(2,8).getRawPtr(),
    tuple<ADAPT::global_ordinal_type>(0,1).getRawPtr());
  mat->InsertGlobalValues(2,NumNz[2],
    tuple<ADAPT::scalar_type>(1).getRawPtr(),
    tuple<ADAPT::global_ordinal_type>(2).getRawPtr());
  mat->InsertGlobalValues(3,NumNz[3],
    tuple<ADAPT::scalar_type>(-3,5).getRawPtr(),
    tuple<ADAPT::global_ordinal_type>(0,3).getRawPtr());
  mat->InsertGlobalValues(4,NumNz[4],
    tuple<ADAPT::scalar_type>(-1,4).getRawPtr(),
    tuple<ADAPT::global_ordinal_type>(1,4).getRawPtr());
  mat->InsertGlobalValues(5,NumNz[5],
    tuple<ADAPT::scalar_type>(-2,6).getRawPtr(),
    tuple<ADAPT::global_ordinal_type>(3,5).getRawPtr());
  mat->FillComplete();

  // Print for sanity sake
  // RCP<FancyOStream> os = getDefaultOStream();
  // mat->describe(*os,Teuchos::VERB_EXTREME);

  RCP<ADAPT> adapter = rcp(new ADAPT(mat));

  Array<ADAPT::scalar_type> nzvals_test(tuple<ADAPT::scalar_type>(7,-3,-1,2,8,1,-3,5,-1,4,-2,6));
  Array<ADAPT::global_ordinal_type> colind_test(tuple<ADAPT::global_ordinal_type>(0,2,4,0,1,2,0,3,1,4,3,5));
  Array<ADAPT::global_size_type> rowptr_test(tuple<ADAPT::global_size_type>(0,3,5,6,8,10,12));

  Array<ADAPT::scalar_type> nzvals(adapter->getGlobalNNZ());
  Array<ADAPT::global_ordinal_type> colind(adapter->getGlobalNNZ());
  Array<ADAPT::global_size_type> rowptr(adapter->getGlobalNumRows() + 1);
  size_t nnz;

  adapter->getCrs(nzvals,colind,rowptr,nnz);

  TEST_COMPARE_ARRAYS(nzvals,nzvals_test);
  TEST_COMPARE_ARRAYS(colind,colind_test);
  TEST_COMPARE_ARRAYS(rowptr,rowptr_test);
  TEST_EQUALITY(nnz,12);
}

TEUCHOS_UNIT_TEST( CrsMatrixAdapter, CCS )
{
  /* Test the getCcs() method of MatrixAdapter.  Again, check against a known
   * matrix
   */
  typedef Epetra_CrsMatrix MAT;
  typedef MatrixAdapter<MAT> ADAPT;

  /* We will be using the following matrix for this test:
   * 
   * [ [ 7,  0,  -3, 0,  -1, 0 ]
   *   [ 2,  8,  0,  0,  0,  0 ]
   *   [ 0,  0,  1,  0,  0,  0 ]
   *   [ -3, 0,  0,  5,  0,  0 ]
   *   [ 0,  -1, 0,  0,  4,  0 ]
   *   [ 0,  0,  0,  -2, 0,  6 ] ]
   */
  Epetra_SerialComm comm;

  // create a Map
  const int num_eqn = 6;

  Epetra_Map Map(num_eqn, 0, comm);
  
  // Get update list and number of local equations from newly created Map.

  int NumMyElements = Map.NumMyElements();

  std::vector<int> MyGlobalElements(NumMyElements);
  Map.MyGlobalElements(&MyGlobalElements[0]);

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
  // on this processor

  int NumNz[] = {3, 2, 1, 2, 2, 2};

   // Create a Epetra_Matrix
  Teuchos::RCP<Epetra_CrsMatrix> mat = rcp(new Epetra_CrsMatrix(Copy, Map, NumNz));
  
  // Construct matrix
  mat->InsertGlobalValues(0,NumNz[0],
    tuple<ADAPT::scalar_type>(7,-3,-1).getRawPtr(),
    tuple<ADAPT::global_ordinal_type>(0,2,4).getRawPtr());
  mat->InsertGlobalValues(1,NumNz[1],
    tuple<ADAPT::scalar_type>(2,8).getRawPtr(),
    tuple<ADAPT::global_ordinal_type>(0,1).getRawPtr());
  mat->InsertGlobalValues(2,NumNz[2],
    tuple<ADAPT::scalar_type>(1).getRawPtr(),
    tuple<ADAPT::global_ordinal_type>(2).getRawPtr());
  mat->InsertGlobalValues(3,NumNz[3],
    tuple<ADAPT::scalar_type>(-3,5).getRawPtr(),
    tuple<ADAPT::global_ordinal_type>(0,3).getRawPtr());
  mat->InsertGlobalValues(4,NumNz[4],
    tuple<ADAPT::scalar_type>(-1,4).getRawPtr(),
    tuple<ADAPT::global_ordinal_type>(1,4).getRawPtr());
  mat->InsertGlobalValues(5,NumNz[5],
    tuple<ADAPT::scalar_type>(-2,6).getRawPtr(),
    tuple<ADAPT::global_ordinal_type>(3,5).getRawPtr());
  mat->FillComplete();

  // Print for sanity sake
  // RCP<FancyOStream> os = getDefaultOStream();
  // mat->describe(*os,Teuchos::VERB_EXTREME);

  RCP<ADAPT> adapter = rcp(new ADAPT(mat));

  Array<ADAPT::scalar_type> nzvals_test(tuple<ADAPT::scalar_type>(7,2,-3,8,-1,-3,1,5,-2,-1,4,6));
  Array<ADAPT::global_ordinal_type> rowind_test(tuple<ADAPT::global_ordinal_type>(0,1,3,1,4,0,2,3,5,0,4,5));
  Array<ADAPT::global_size_type> colptr_test(tuple<ADAPT::global_size_type>(0,3,5,7,9,11,12));

  Array<ADAPT::scalar_type> nzvals(adapter->getGlobalNNZ());
  Array<ADAPT::global_ordinal_type> rowind(adapter->getGlobalNNZ());
  Array<ADAPT::global_size_type> colptr(adapter->getGlobalNumRows() + 1);
  size_t nnz;

  adapter->getCcs(nzvals,rowind,colptr,nnz);

  TEST_COMPARE_ARRAYS(nzvals,nzvals_test);
  TEST_COMPARE_ARRAYS(rowind,rowind_test);
  TEST_COMPARE_ARRAYS(colptr,colptr_test);
  TEST_EQUALITY(nnz,12);
}

/* Also need to test the updateValues[Crs,Ccs]() methods, to make sure that
 * those changes are persistent in the underlying matrices (once implemented).
 */

} // end anonymous namespace
