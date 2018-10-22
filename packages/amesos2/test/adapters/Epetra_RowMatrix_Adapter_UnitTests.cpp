// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#include <string>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>

#include <Epetra_CrsMatrix.h>
#ifdef HAVE_MPI
#  include <Epetra_MpiComm.h>
#else
#  include <Epetra_SerialComm.h>
#endif

#include "Amesos2_MatrixAdapter_def.hpp"
#include "Amesos2_Util.hpp"

namespace {

  using std::cout;
  using std::endl;
  using std::string;

  using Teuchos::as;
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

  using Amesos2::MatrixAdapter;

  using Amesos2::Meta::is_same;

  using Amesos2::ROOTED;
  using Amesos2::GLOBALLY_REPLICATED;

  using Amesos2::Util::to_teuchos_comm;

  typedef Tpetra::Map<>::node_type Node;

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

  const RCP<Epetra_Comm> getDefaultComm()
  {
#ifdef EPETRA_MPI
    return rcp(new Epetra_MpiComm( MPI_COMM_WORLD ));
#else
    return rcp(new Epetra_SerialComm());
#endif
  }

  /*
  RCP<FancyOStream> getDefaultOStream()
  {
    return( VerboseObjectBase::getDefaultOStream() );
  }
  */

  template<typename T1, typename T2>
  const RCP<Array<std::pair<T1,T2> > >
  zip(const ArrayView<T1>& a, const ArrayView<T2>& b)
  {
    typedef std::pair<T1,T2> result_type;
    size_t size = std::min(a.size(), b.size());
    RCP<Array<result_type > > r = rcp( new Array<result_type>(size) );
    for( size_t i = 0; i < size; ++i){
      (*r)[i] = std::make_pair(a[i], b[i]);
    }
    return(r);
  }
  
  template<typename T>
  bool contains(const ArrayView<T> a, T t)
  {
    typedef typename ArrayView<T>::iterator it;
    it first = a.begin();
    it last  = a.end();
    return( std::find(first, last, t) != last );
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
     * - Correct typedefs ( using Amesos2::is_same<> )
     */
    typedef Epetra_CrsMatrix MAT;
    typedef MatrixAdapter<MAT> ADAPT;

    RCP<Epetra_Comm> comm = getDefaultComm();
    //const size_t numprocs = comm->NumProc();
    //  const int rank     = comm->MyPID();
    // create a Map
    const int num_eqn = 100;

    Epetra_Map Map(num_eqn, 0, *comm);
  
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
	TEUCHOS_TEST_FOR_EXCEPTION(ierr != 0,std::runtime_error,"Error inserting value into matrix");
      }

    eye->FillComplete();
    RCP<ADAPT> adapter = Amesos2::createMatrixAdapter<MAT>(eye);

    // The following should all pass at compile time
    TEST_ASSERT( (is_same<double,ADAPT::scalar_t>::value) );
    TEST_ASSERT( (is_same<int,ADAPT::local_ordinal_t>::value) );
    TEST_ASSERT( (is_same<int,ADAPT::global_ordinal_t>::value) );
    TEST_ASSERT( (is_same<size_t,ADAPT::global_size_t>::value) );
    TEST_ASSERT( (is_same<MAT,ADAPT::matrix_t>::value) );

  }

  TEUCHOS_UNIT_TEST( CrsMatrixAdapter, Dimensions )
  {
    // Test that the dimensions reported by the adapter match those as reported
    // by the Tpetra::CrsMatrix
    // Check equality of mapped method calls
    typedef Epetra_CrsMatrix MAT;
    typedef MatrixAdapter<MAT> ADAPT;

    RCP<Epetra_Comm> comm = getDefaultComm();
    //const size_t numprocs = comm->NumProc();
    //  const int rank     = comm->MyPID();
    // create a Map
    const int num_eqn = 100;

    Epetra_Map Map(num_eqn, 0, *comm);
  
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
	TEUCHOS_TEST_FOR_EXCEPTION(ierr != 0,std::runtime_error,"Error inserting value into matrix");
      }

    eye->FillComplete();
    // Constructor from RCP
    RCP<ADAPT> adapter  = Amesos2::createMatrixAdapter<MAT>( eye );

    TEST_EQUALITY(Teuchos::as<ADAPT::global_size_t>(eye->NumGlobalNonzeros()), adapter->getGlobalNNZ());
    TEST_EQUALITY(Teuchos::as<ADAPT::global_size_t>(eye->NumGlobalRows()), adapter->getGlobalNumRows());
    TEST_EQUALITY(Teuchos::as<ADAPT::global_size_t>(eye->NumGlobalCols()), adapter->getGlobalNumCols());
  }


  TEUCHOS_UNIT_TEST( CrsMatrixAdapter, CRS_Serial )
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
    RCP<Epetra_Comm> comm = getDefaultComm();
    int rank = comm->MyPID();

    // create a Map
    const int num_eqn = 6;

    Epetra_Map Map(num_eqn, 0, *comm);
  
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
			    tuple<ADAPT::scalar_t>(7,-3,-1).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(0,2,4).getRawPtr());
    mat->InsertGlobalValues(1,NumNz[1],
			    tuple<ADAPT::scalar_t>(2,8).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(0,1).getRawPtr());
    mat->InsertGlobalValues(2,NumNz[2],
			    tuple<ADAPT::scalar_t>(1).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(2).getRawPtr());
    mat->InsertGlobalValues(3,NumNz[3],
			    tuple<ADAPT::scalar_t>(-3,5).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(0,3).getRawPtr());
    mat->InsertGlobalValues(4,NumNz[4],
			    tuple<ADAPT::scalar_t>(-1,4).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(1,4).getRawPtr());
    mat->InsertGlobalValues(5,NumNz[5],
			    tuple<ADAPT::scalar_t>(-2,6).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(3,5).getRawPtr());
    mat->FillComplete();

    // Print for sanity sake
    // RCP<FancyOStream> os = getDefaultOStream();
    // mat->describe(*os,Teuchos::VERB_EXTREME);

    RCP<ADAPT> adapter = Amesos2::createMatrixAdapter<MAT>(mat);

    Array<ADAPT::scalar_t> nzvals_test(tuple<ADAPT::scalar_t>(7,-3,-1,2,8,1,-3,5,-1,4,-2,6));
    Array<ADAPT::global_ordinal_t> colind_test(tuple<ADAPT::global_ordinal_t>(0,2,4,0,1,2,0,3,1,4,3,5));
    Array<ADAPT::global_size_t> rowptr_test(tuple<ADAPT::global_size_t>(0,3,5,6,8,10,12));

    Array<ADAPT::scalar_t> nzvals(adapter->getGlobalNNZ());
    Array<ADAPT::global_ordinal_t> colind(adapter->getGlobalNNZ());
    Array<ADAPT::global_size_t> rowptr(adapter->getGlobalNumRows() + 1);
    size_t nnz;

    adapter->getCrs(nzvals,colind,rowptr,nnz,ROOTED);

    // getCrs does not guarantee the sorted-ness of the column
    // indices, so this test might fail

    if(rank == 0)
    {
      TEST_COMPARE_ARRAYS(nzvals, nzvals_test);
      TEST_COMPARE_ARRAYS(colind, colind_test);
      TEST_COMPARE_ARRAYS(rowptr, rowptr_test);
      TEST_EQUALITY_CONST(nnz, 12);
    }

    /////////////////////////////////////////////
    // Check now a rooted, sorted-indices repr //
    /////////////////////////////////////////////

    adapter->getCrs(nzvals,colind,rowptr,nnz,ROOTED,Amesos2::SORTED_INDICES);

    if ( rank == 0 ){
      // Now the arrays should compare directly
      TEST_COMPARE_ARRAYS(nzvals, nzvals_test);
      TEST_COMPARE_ARRAYS(colind, colind_test);
      TEST_COMPARE_ARRAYS(rowptr, rowptr_test);
      TEST_EQUALITY_CONST(nnz, 12);
    }
  }

  TEUCHOS_UNIT_TEST( CrsMatrixAdapter, CRS_Replicated )
  {
    /* Test the getCrs() method of MatrixAdapter.  We check against a simple
     * test matrix that we construct on the fly.
     */
    typedef Epetra_CrsMatrix MAT;
    typedef MatrixAdapter<MAT> ADAPT;
    typedef int GO;
    typedef int global_size_t;
    typedef std::pair<double,GO> my_pair_t;

    /* We will be using the following matrix for this test:
     * 
     * [ [ 7,  0,  -3, 0,  -1, 0 ]
     *   [ 2,  8,  0,  0,  0,  0 ]
     *   [ 0,  0,  1,  0,  0,  0 ]
     *   [ -3, 0,  0,  5,  0,  0 ]
     *   [ 0,  -1, 0,  0,  4,  0 ]
     *   [ 0,  0,  0,  -2, 0,  6 ] ]
     */
    RCP<Epetra_Comm> comm = getDefaultComm();

    // create a Map
    const int num_eqn = 6;

    Epetra_Map Map(num_eqn, 0, *comm);
  
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
			    tuple<ADAPT::scalar_t>(7,-3,-1).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(0,2,4).getRawPtr());
    mat->InsertGlobalValues(1,NumNz[1],
			    tuple<ADAPT::scalar_t>(2,8).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(0,1).getRawPtr());
    mat->InsertGlobalValues(2,NumNz[2],
			    tuple<ADAPT::scalar_t>(1).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(2).getRawPtr());
    mat->InsertGlobalValues(3,NumNz[3],
			    tuple<ADAPT::scalar_t>(-3,5).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(0,3).getRawPtr());
    mat->InsertGlobalValues(4,NumNz[4],
			    tuple<ADAPT::scalar_t>(-1,4).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(1,4).getRawPtr());
    mat->InsertGlobalValues(5,NumNz[5],
			    tuple<ADAPT::scalar_t>(-2,6).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(3,5).getRawPtr());
    mat->FillComplete();

    // Print for sanity sake
    // RCP<FancyOStream> os = getDefaultOStream();
    // mat->describe(*os,Teuchos::VERB_EXTREME);

    RCP<ADAPT> adapter = Amesos2::createMatrixAdapter<MAT>(mat);

    Array<ADAPT::scalar_t> nzvals_test(tuple<ADAPT::scalar_t>(7,-3,-1,2,8,1,-3,5,-1,4,-2,6));
    Array<ADAPT::global_ordinal_t> colind_test(tuple<ADAPT::global_ordinal_t>(0,2,4,0,1,2,0,3,1,4,3,5));
    Array<ADAPT::global_size_t> rowptr_test(tuple<ADAPT::global_size_t>(0,3,5,6,8,10,12));

    Array<ADAPT::scalar_t> nzvals(adapter->getGlobalNNZ());
    Array<ADAPT::global_ordinal_t> colind(adapter->getGlobalNNZ());
    Array<ADAPT::global_size_t> rowptr(adapter->getGlobalNumRows() + 1);
    size_t nnz;

    ////////////////////////////////////////////////////
    // Now check a globally replicated representation //
    ////////////////////////////////////////////////////

    adapter->getCrs(nzvals,colind,rowptr,nnz,GLOBALLY_REPLICATED);

    // All processes check
    
    // getCrs() does not guarantee a permutation of the non-zero
    // values and the column indices in the Arbitrary case, we just
    // know that they need to match up with what is expected.
    GO maxRow = Map.MaxAllGID();
    for ( GO row = Map.MinAllGID(); row <= maxRow; ++row ){
      global_size_t rp  = rowptr[row];
      global_size_t nrp = rowptr[row+1];
      global_size_t row_nnz = nrp - rp;
      TEST_ASSERT( rp < as<global_size_t>(nzvals.size()) );
      TEST_ASSERT( rp < as<global_size_t>(colind.size()) );
      const RCP<Array<my_pair_t> > expected_pairs
	= zip(nzvals_test.view(rp,row_nnz), colind_test.view(rp,row_nnz));
      const RCP<Array<my_pair_t> > got_pairs
	= zip(nzvals.view(rp,row_nnz), colind.view(rp,row_nnz));
      for ( global_size_t i = 0; i < row_nnz; ++i ){
	TEST_ASSERT( contains((*got_pairs)(), (*expected_pairs)[i]) );
      }
    }
    TEST_COMPARE_ARRAYS(rowptr, rowptr_test);
    TEST_EQUALITY_CONST(nnz, 12);

    ///////////////////////////////////////////
    // Check globally-repl, with sorted indx //
    ///////////////////////////////////////////

    adapter->getCrs(nzvals,colind,rowptr,nnz,
		    GLOBALLY_REPLICATED,Amesos2::SORTED_INDICES);

    TEST_COMPARE_ARRAYS(nzvals, nzvals_test);
    TEST_COMPARE_ARRAYS(colind, colind_test);
    TEST_COMPARE_ARRAYS(rowptr, rowptr_test);
    TEST_EQUALITY_CONST(nnz, 12);
  }

  TEUCHOS_UNIT_TEST( CrsMatrixAdapter, CRS_Map )
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
    RCP<Epetra_Comm> comm = getDefaultComm();
    int numprocs = comm->NumProc();
    int rank = comm->MyPID();

    // create a Map
    const int num_eqn = 6;

    Epetra_Map Map(num_eqn, 0, *comm);
  
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
			    tuple<ADAPT::scalar_t>(7,-3,-1).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(0,2,4).getRawPtr());
    mat->InsertGlobalValues(1,NumNz[1],
			    tuple<ADAPT::scalar_t>(2,8).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(0,1).getRawPtr());
    mat->InsertGlobalValues(2,NumNz[2],
			    tuple<ADAPT::scalar_t>(1).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(2).getRawPtr());
    mat->InsertGlobalValues(3,NumNz[3],
			    tuple<ADAPT::scalar_t>(-3,5).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(0,3).getRawPtr());
    mat->InsertGlobalValues(4,NumNz[4],
			    tuple<ADAPT::scalar_t>(-1,4).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(1,4).getRawPtr());
    mat->InsertGlobalValues(5,NumNz[5],
			    tuple<ADAPT::scalar_t>(-2,6).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(3,5).getRawPtr());
    mat->FillComplete();

    // Print for sanity sake
    // RCP<FancyOStream> os = getDefaultOStream();
    // mat->describe(*os,Teuchos::VERB_EXTREME);

    RCP<ADAPT> adapter = Amesos2::createMatrixAdapter<MAT>(mat);

    Array<ADAPT::scalar_t> nzvals_test(tuple<ADAPT::scalar_t>(7,-3,-1,2,8,1,-3,5,-1,4,-2,6));
    Array<ADAPT::global_ordinal_t> colind_test(tuple<ADAPT::global_ordinal_t>(0,2,4,0,1,2,0,3,1,4,3,5));
    Array<ADAPT::global_size_t> rowptr_test(tuple<ADAPT::global_size_t>(0,3,5,6,8,10,12));

    Array<ADAPT::scalar_t> nzvals(adapter->getGlobalNNZ());
    Array<ADAPT::global_ordinal_t> colind(adapter->getGlobalNNZ());
    Array<ADAPT::global_size_t> rowptr(adapter->getGlobalNumRows() + 1);
    size_t nnz;

    /**
     * Check the getCrs overload that accepts a row-distribution map
     * as input.  Divide the 6-row matrix in two, give the top half to
     * rank 0, and the bottom half to rank 1.  Then check the results.
     */
    size_t my_num_rows = OrdinalTraits<size_t>::zero();
    if ( numprocs > 1 ){
      if ( rank < 2 ){
	my_num_rows = 3;		// total num_rows is 6
      }
    } else {			// We only have 1 proc, then she just takes it all
      my_num_rows = 6;
    }
    const Tpetra::Map<int,int> half_map(6, my_num_rows, 0,
					to_teuchos_comm(comm));

    adapter->getCrs(nzvals,colind,rowptr,nnz, Teuchos::ptrInArg(half_map), Amesos2::SORTED_INDICES, Amesos2::DISTRIBUTED); // ROOTED = default distribution

    /*
     * Check that you got the entries you'd expect
     *
     * It's convenient that exactly half of the non-zero entries are
     * found in the top half of the rows, and the other half are found
     * in the bottom half of the rows.
     */
    if(rank == 0) {
      TEST_COMPARE_ARRAYS(nzvals.view(0,6), nzvals_test.view(0,6));
      TEST_COMPARE_ARRAYS(colind.view(0,6), colind_test.view(0,6));
      for(int i = 0; i < 4; i++) {
        TEST_EQUALITY_CONST(rowptr[i], rowptr_test[i]);
      }
    } else if(rank == 1) {
      TEST_COMPARE_ARRAYS(nzvals.view(0,6), nzvals_test.view(6,6));
      TEST_COMPARE_ARRAYS(colind.view(0,6), colind_test.view(6,6));
      for(int i = 0; i < 4; i++) {
        TEST_EQUALITY_CONST(rowptr[i], rowptr_test[3 + i] - rowptr_test[3]);
      }
    }
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
    RCP<Epetra_Comm> comm = getDefaultComm();

    // create a Map
    const int num_eqn = 6;

    Epetra_Map Map(num_eqn, 0, *comm);
  
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
			    tuple<ADAPT::scalar_t>(7,-3,-1).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(0,2,4).getRawPtr());
    mat->InsertGlobalValues(1,NumNz[1],
			    tuple<ADAPT::scalar_t>(2,8).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(0,1).getRawPtr());
    mat->InsertGlobalValues(2,NumNz[2],
			    tuple<ADAPT::scalar_t>(1).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(2).getRawPtr());
    mat->InsertGlobalValues(3,NumNz[3],
			    tuple<ADAPT::scalar_t>(-3,5).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(0,3).getRawPtr());
    mat->InsertGlobalValues(4,NumNz[4],
			    tuple<ADAPT::scalar_t>(-1,4).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(1,4).getRawPtr());
    mat->InsertGlobalValues(5,NumNz[5],
			    tuple<ADAPT::scalar_t>(-2,6).getRawPtr(),
			    tuple<ADAPT::global_ordinal_t>(3,5).getRawPtr());
    mat->FillComplete();

    // Print for sanity sake
    // RCP<FancyOStream> os = getDefaultOStream();
    // mat->describe(*os,Teuchos::VERB_EXTREME);

    RCP<ADAPT> adapter = Amesos2::createMatrixAdapter<MAT>(mat);

    Array<ADAPT::scalar_t> nzvals_test(tuple<ADAPT::scalar_t>(7,2,-3,8,-1,-3,1,5,-2,-1,4,6));
    Array<ADAPT::global_ordinal_t> rowind_test(tuple<ADAPT::global_ordinal_t>(0,1,3,1,4,0,2,3,5,0,4,5));
    Array<ADAPT::global_size_t> colptr_test(tuple<ADAPT::global_size_t>(0,3,5,7,9,11,12));

    Array<ADAPT::scalar_t> nzvals(adapter->getGlobalNNZ());
    Array<ADAPT::global_ordinal_t> rowind(adapter->getGlobalNNZ());
    Array<ADAPT::global_size_t> colptr(adapter->getGlobalNumRows() + 1);
    size_t nnz;

    adapter->getCcs(nzvals,rowind,colptr,nnz,GLOBALLY_REPLICATED);

    TEST_COMPARE_ARRAYS(nzvals,nzvals_test);
    TEST_COMPARE_ARRAYS(rowind,rowind_test);
    TEST_COMPARE_ARRAYS(colptr,colptr_test);
    TEST_EQUALITY(nnz,12);
  }

  /* Also need to test the updateValues[Crs,Ccs]() methods, to make sure that
   * those changes are persistent in the underlying matrices (once implemented).
   */

} // end anonymous namespace
