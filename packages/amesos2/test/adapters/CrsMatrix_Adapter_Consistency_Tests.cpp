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

/**
 * \file   CrsMatrix_Adapter_Consistency_Tests.cpp
 * \author Eric Bavier <etbavie@muenster.srn.sandia.gov>
 * \date   Fri Jul 15 07:36:37 2011
 * 
 * \brief Checks consistency between the Epetra_CrsMatrix and
 * Tpetra::CrsMatrix adapters.
 * 
 * For the same input matrices, the adapters should return the same
 * output values.  These tests only check double scalar and int
 * ordinal.
 */


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
#include <Tpetra_CrsMatrix.hpp>

#include <Epetra_CrsMatrix.h>
#ifdef HAVE_MPI
#  include <Epetra_MpiComm.h>
#else
#  include <Epetra_SerialComm.h>
#endif

#include <MatrixMarket_Tpetra.hpp> // reading matrixmarket files for Tpetra
#include <EpetraExt_CrsMatrixIn.h> // reading matrixmarket files for Epetra

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
  
  using Amesos2::Util::to_teuchos_comm;

  using Amesos2::ROOTED;
  using Amesos2::GLOBALLY_REPLICATED;
  using Amesos2::SORTED_INDICES;

  // Where to look for input files
  string filedir = "../matrices/";
  string test_mat_file = "arc130.mtx";

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.setOption("filedir",&filedir,"Directory of matrix files.");
    clp.setOption("testfile",&test_mat_file,"Matrix-Market file to use for tests");
    clp.addOutputSetupOptions(true);
  }

  const RCP<const Epetra_Comm> getDefaultEComm()
  {
#ifdef HAVE_MPI
    return rcp(new Epetra_MpiComm( MPI_COMM_WORLD ));
#else
    return rcp(new Epetra_SerialComm());
#endif
  }

  RCP<const Comm<int> > getDefaultTComm()
  {
#ifdef HAVE_MPI
    return rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));
#else
    return rcp(new Teuchos::SerialComm<int>());
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


  TEUCHOS_UNIT_TEST( CrsMatrixAdapter, Dimensions )
  {
    // Test that the dimensions reported by the adapter match those as reported
    // by the Tpetra::CrsMatrix
    // Check equality of mapped method calls
    typedef Epetra_CrsMatrix EMAT;
    typedef Tpetra::CrsMatrix<double,int,int> TMAT;
    typedef MatrixAdapter<EMAT> EADAPT;
    typedef MatrixAdapter<TMAT> TADAPT;

    RCP<const Epetra_Comm> ecomm = getDefaultEComm();
    RCP<const Comm<int> > tcomm = getDefaultTComm();

    // Read in test matrix file
    std::string mat_pathname = filedir + test_mat_file;
    RCP<TMAT> TA = Tpetra::MatrixMarket::Reader<TMAT>::readSparseFile(mat_pathname,tcomm);
    EMAT* EA;
    EpetraExt::MatrixMarketFileToCrsMatrix(mat_pathname.c_str(), *ecomm, EA, false, false);

    // Create adapters
    RCP<TADAPT> tadapter = Amesos2::createMatrixAdapter<TMAT>( TA );
    RCP<EADAPT> eadapter = Amesos2::createMatrixAdapter<EMAT>( rcp(EA) ); // this rcp becomes owning

    // Check that adapters return the same dimensions and stats

    TEST_EQUALITY(tadapter->getGlobalNNZ(), eadapter->getGlobalNNZ());
    TEST_EQUALITY(tadapter->getGlobalNumRows(), eadapter->getGlobalNumRows());
    TEST_EQUALITY(tadapter->getGlobalNumCols(), eadapter->getGlobalNumCols());

    // The column and row maps could have been initialized differently
    // in the input readers, so it would be meaningless to check their
    // equality or the equality of the local statistics.
  }


  TEUCHOS_UNIT_TEST( CrsMatrixAdapter, CRS_Serial )
  {
    /* Test the getCrs() method of MatrixAdapter.
     */
    typedef Epetra_CrsMatrix EMAT;
    typedef Tpetra::CrsMatrix<double,int,int> TMAT;
    typedef MatrixAdapter<EMAT> EADAPT;
    typedef MatrixAdapter<TMAT> TADAPT;

    RCP<const Epetra_Comm> ecomm = getDefaultEComm();
    RCP<const Comm<int> > tcomm = getDefaultTComm();

    // Read in test matrix file
    std::string mat_pathname = filedir + test_mat_file;
    RCP<TMAT> TA = Tpetra::MatrixMarket::Reader<TMAT>::readSparseFile(mat_pathname,tcomm);
    EMAT* EA;
    EpetraExt::MatrixMarketFileToCrsMatrix(mat_pathname.c_str(), *ecomm, EA, false, false);

    // Create adapters
    RCP<TADAPT> tadapter = Amesos2::createMatrixAdapter<TMAT>( TA );
    RCP<EADAPT> eadapter = Amesos2::createMatrixAdapter<EMAT>( rcp(EA) );

    Array<TADAPT::scalar_t> tnzvals(tadapter->getGlobalNNZ());
    Array<TADAPT::global_ordinal_t> tcolind(tadapter->getGlobalNNZ());
    Array<TADAPT::global_size_t> trowptr(tadapter->getGlobalNumRows() + 1);
    size_t tnnz = 0;

    Array<EADAPT::scalar_t> enzvals(eadapter->getGlobalNNZ());
    Array<EADAPT::global_ordinal_t> ecolind(eadapter->getGlobalNNZ());
    Array<EADAPT::global_size_t> erowptr(eadapter->getGlobalNumRows() + 1);
    size_t ennz = 0;

    tadapter->getCrs(tnzvals, tcolind, trowptr, tnnz, ROOTED, SORTED_INDICES);
    eadapter->getCrs(enzvals, ecolind, erowptr, ennz, ROOTED, SORTED_INDICES);

    // The nzvals, colind, rowptr, and nnz values for each adapter
    // should be the same at the root.  It is ok to to these tests at
    // non-root rank because the arrays will have been initialized to
    // zero and the nnz as well.

    TEST_COMPARE_ARRAYS(tnzvals, enzvals);
    TEST_COMPARE_ARRAYS(tcolind, ecolind);
    TEST_COMPARE_ARRAYS(trowptr, erowptr);
    TEST_EQUALITY(tnnz, ennz);
  }

  TEUCHOS_UNIT_TEST( CrsMatrixAdapter, CRS_Replicated )
  {
    /* Test the getCrs() method of MatrixAdapter.  This unit test is
     * essentially the same as the above, but the global matrix
     * representation should be valid on all processors
     */
    typedef Epetra_CrsMatrix EMAT;
    typedef Tpetra::CrsMatrix<double,int,int> TMAT;
    typedef MatrixAdapter<EMAT> EADAPT;
    typedef MatrixAdapter<TMAT> TADAPT;

    RCP<const Epetra_Comm> ecomm = getDefaultEComm();
    RCP<const Comm<int> > tcomm = getDefaultTComm();

    // Read in test matrix file
    std::string mat_pathname = filedir + test_mat_file;
    RCP<TMAT> TA = Tpetra::MatrixMarket::Reader<TMAT>::readSparseFile(mat_pathname,tcomm);
    EMAT* EA;
    EpetraExt::MatrixMarketFileToCrsMatrix(mat_pathname.c_str(), *ecomm, EA, false, false);

    // Create adapters
    RCP<TADAPT> tadapter = Amesos2::createMatrixAdapter<TMAT>( TA );
    RCP<EADAPT> eadapter = Amesos2::createMatrixAdapter<EMAT>( rcp(EA) );

    Array<TADAPT::scalar_t> tnzvals(tadapter->getGlobalNNZ());
    Array<TADAPT::global_ordinal_t> tcolind(tadapter->getGlobalNNZ());
    Array<TADAPT::global_size_t> trowptr(tadapter->getGlobalNumRows() + 1);
    size_t tnnz = 0;

    Array<EADAPT::scalar_t> enzvals(eadapter->getGlobalNNZ());
    Array<EADAPT::global_ordinal_t> ecolind(eadapter->getGlobalNNZ());
    Array<EADAPT::global_size_t> erowptr(eadapter->getGlobalNumRows() + 1);
    size_t ennz = 0;

    tadapter->getCrs(tnzvals, tcolind, trowptr, tnnz, GLOBALLY_REPLICATED, SORTED_INDICES);
    eadapter->getCrs(enzvals, ecolind, erowptr, ennz, GLOBALLY_REPLICATED, SORTED_INDICES);

    // The nzvals, colind, rowptr, and nnz values for each adapter
    // should be the same at the root.  It is ok to to these tests at
    // non-root rank because the arrays will have been initialized to
    // zero and the nnz as well.

    TEST_COMPARE_ARRAYS(tnzvals, enzvals);
    TEST_COMPARE_ARRAYS(tcolind, ecolind);
    TEST_COMPARE_ARRAYS(trowptr, erowptr);
    TEST_EQUALITY(tnnz, ennz);
  }

  TEUCHOS_UNIT_TEST( CrsMatrixAdapter, CRS_Map )
  {
    /* Test the getCrs() method of MatrixAdapter.  We check against a simple
     * test matrix that we construct on the fly.
     */
    typedef Epetra_CrsMatrix EMAT;
    typedef Tpetra::CrsMatrix<double,int,int> TMAT;
    typedef MatrixAdapter<EMAT> EADAPT;
    typedef MatrixAdapter<TMAT> TADAPT;

    RCP<const Epetra_Comm> ecomm = getDefaultEComm();
    RCP<const Comm<int> > tcomm = getDefaultTComm();
    const int numprocs = tcomm->getSize();
    const int rank = tcomm->getRank();

    // Read in test matrix file
    std::string mat_pathname = filedir + test_mat_file;
    RCP<TMAT> TA = Tpetra::MatrixMarket::Reader<TMAT>::readSparseFile(mat_pathname,tcomm);
    EMAT* EA;
    EpetraExt::MatrixMarketFileToCrsMatrix(mat_pathname.c_str(), *ecomm, EA, false, false);

    // Create adapters
    RCP<TADAPT> tadapter = Amesos2::createMatrixAdapter<TMAT>( TA );
    RCP<EADAPT> eadapter = Amesos2::createMatrixAdapter<EMAT>( rcp(EA) );

    TADAPT::global_size_t g_num_rows = tadapter->getGlobalNumRows();

    Array<TADAPT::scalar_t> tnzvals(tadapter->getGlobalNNZ());
    Array<TADAPT::global_ordinal_t> tcolind(tadapter->getGlobalNNZ());
    Array<TADAPT::global_size_t> trowptr(g_num_rows + 1);
    size_t tnnz = 0;

    Array<EADAPT::scalar_t> enzvals(eadapter->getGlobalNNZ());
    Array<EADAPT::global_ordinal_t> ecolind(eadapter->getGlobalNNZ());
    Array<EADAPT::global_size_t> erowptr(g_num_rows + 1);
    size_t ennz = 0;

    /**
     * Check the getCrs overload that accepts a row-distribution map
     * as input.  Divide the 6-row matrix in two, give the top half to
     * rank 0, and the bottom half to rank 1.  Then check the results.
     */
    TADAPT::global_size_t my_num_rows = OrdinalTraits<TADAPT::global_size_t>::zero();
    if ( numprocs > 1 ){
      if ( rank < 2 ){
	my_num_rows = g_num_rows / 2;
      }
      // If we have an odd number of rows, rank=0 gets the remainder
      if ( rank == 0 ) my_num_rows += g_num_rows % 2;
    } else {			// We only have 1 proc, then she just takes it all
      my_num_rows = g_num_rows;
    }
    const Tpetra::Map<int,int> half_map(g_num_rows, my_num_rows, 0, tcomm);

    tadapter->getCrs(tnzvals, tcolind, trowptr, tnnz, Teuchos::ptrInArg(half_map), SORTED_INDICES, Amesos2::DISTRIBUTED); // ROOTED = default distribution
    eadapter->getCrs(enzvals, ecolind, erowptr, ennz, Teuchos::ptrInArg(half_map), SORTED_INDICES, Amesos2::DISTRIBUTED);

    /*
     * Check that you got the entries you'd expect
     *
     * It's convenient that exactly half of the non-zero entries are
     * found in the top half of the rows, and the other half are found
     * in the bottom half of the rows.
     */
    if ( rank == 0 ){
      TEST_EQUALITY(tnnz, ennz);
      TEST_EQUALITY_CONST(trowptr[my_num_rows], tnnz);
      TEST_EQUALITY_CONST(erowptr[my_num_rows], ennz);
	
      TEST_COMPARE_ARRAYS(tnzvals.view(0,tnnz), enzvals.view(0,ennz));
      TEST_COMPARE_ARRAYS(tcolind.view(0,tnnz), ecolind.view(0,ennz));
      TEST_COMPARE_ARRAYS(trowptr.view(0,my_num_rows), erowptr.view(0,my_num_rows));
    } else if ( rank == 1 ){
      TEST_EQUALITY(tnnz, ennz);
      TEST_EQUALITY_CONST(trowptr[my_num_rows], tnnz);
      TEST_EQUALITY_CONST(erowptr[my_num_rows], ennz);
	
      TEST_COMPARE_ARRAYS(tnzvals.view(0,tnnz), enzvals.view(0,ennz));
      TEST_COMPARE_ARRAYS(tcolind.view(0,tnnz), ecolind.view(0,ennz));
      TEST_COMPARE_ARRAYS(trowptr.view(0,my_num_rows), erowptr.view(0,my_num_rows));
    }
  }
} // end anonymous namespace

  
