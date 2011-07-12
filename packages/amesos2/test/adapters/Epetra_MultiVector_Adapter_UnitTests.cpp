// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2010 Sandia Corporation
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

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_as.hpp>

#include <Tpetra_DefaultPlatform.hpp>

#include "Amesos2_EpetraMultiVecAdapter.hpp"
#include "Amesos2_Util.hpp"
#include "Amesos2_Util_is_same.hpp"

namespace {

  using std::cout;
  using std::endl;
  using std::string;

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::ptrInArg;
  using Teuchos::outArg;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::tuple;
  using Teuchos::ScalarTraits;
  using Teuchos::OrdinalTraits;
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObjectBase;

  using Amesos::MultiVecAdapter;
  using Amesos::createMultiVecAdapter;

  using Amesos::Util::getDistributionMap;
  using Amesos::Util::to_teuchos_comm;
  using Amesos::Util::is_same;
  using Amesos::Util::get_1d_copy_helper;
  using Amesos::Util::put_1d_data_helper;

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

  const Epetra_Comm* getDefaultComm()
  {
#ifdef EPETRA_MPI
    return new Epetra_MpiComm( MPI_COMM_WORLD );
#else
    return new Epetra_SerialComm();
#endif
  }

  RCP<FancyOStream> getDefaultOStream()
  {
    return( VerboseObjectBase::getDefaultOStream() );
  }

  /*
   * UNIT TESTS
   */

  TEUCHOS_UNIT_TEST( MultiVecAdapter, Initialization )
  {
    /* Test correct initialization of the MultiVecAdapter
     *
     * - All Constructors
     * - Correct initialization of class members
     * - Correct typedefs ( using Amesos::is_same<> )
     */
    typedef ScalarTraits<double> ST;
    typedef Epetra_MultiVector MV;
    typedef MultiVecAdapter<MV> ADAPT;

    Epetra_SerialComm comm;
    // create a Map
    const size_t numLocal = 10;
    Epetra_Map map(-1,numLocal,0,comm);

    RCP<MV> mv = rcp(new MV(map,11));
    mv->Random();

    RCP<ADAPT> adapter = createMultiVecAdapter(mv);

    // The following should all pass at compile time
    TEST_ASSERT( (is_same<double, ADAPT::scalar_t>::value) );
    TEST_ASSERT( (is_same<int,    ADAPT::local_ordinal_t>::value) );
    TEST_ASSERT( (is_same<int,    ADAPT::global_ordinal_t>::value) );
    TEST_ASSERT( (is_same<Node,   ADAPT::node_t>::value) );
    TEST_ASSERT( (is_same<size_t, ADAPT::global_size_t>::value) );
    TEST_ASSERT( (is_same<MV,     ADAPT::multivec_t>::value) );

  }

  TEUCHOS_UNIT_TEST( MultiVecAdapter, Dimensions )
  {
    // Test that the dimensions reported by the adapter match those as reported
    // by the Tpetra::MultiVector
    typedef ScalarTraits<double> ST;
    typedef Epetra_MultiVector MV;
    typedef MultiVecAdapter<MV> ADAPT;

    const Epetra_Comm* comm = getDefaultComm();
    // create a Map
    const size_t numLocal = 10;
    Epetra_Map map(-1,numLocal,0,*comm);

    RCP<MV> mv = rcp(new MV(map,11));
    mv->Random();

    RCP<ADAPT> adapter = createMultiVecAdapter(mv);

    TEST_EQUALITY( mv->MyLength(),     as<int>(adapter->getLocalLength())     );
    TEST_EQUALITY( mv->NumVectors(),   as<int>(adapter->getLocalNumVectors()) );
    TEST_EQUALITY( mv->NumVectors(),   as<int>(adapter->getGlobalNumVectors()));
    TEST_EQUALITY( mv->GlobalLength(), as<int>(adapter->getGlobalLength())    );
    TEST_EQUALITY( mv->Stride(),       as<int>(adapter->getStride())          );

    delete comm;
  }

  TEUCHOS_UNIT_TEST( MultiVecAdapter, Copy )
  {
    /* Test the get1dCopy() method of MultiVecAdapter.  We can check against a
     * known multivector and also check against what is returned by the
     * Tpetra::MultiVector.
     */
    typedef ScalarTraits<double> ST;
    typedef Epetra_MultiVector MV;
    typedef MultiVecAdapter<MV> ADAPT;

    const Epetra_Comm* comm = getDefaultComm();
    int numprocs = comm->NumProc();
    int rank = comm->MyPID();

    // create a Map
    const size_t numVectors = 3;
    const size_t numLocal = 13;
    Epetra_Map map(-1,numLocal,0,*comm);

    RCP<MV> mv = rcp(new MV(map,numVectors));
    mv->Random();

    // mv->Print(std::cout);

    RCP<ADAPT> adapter = createMultiVecAdapter(mv);
    Array<double> original(numVectors*numLocal*numprocs);
    Array<double> copy(numVectors*numLocal*numprocs);

    get_1d_copy_helper<ADAPT,double>::do_get(ptrInArg(*adapter), copy(),
					     numLocal*numprocs, Amesos::Util::Rooted);

    // Just rank==0 process has global copy of mv data, check against an import
    int my_elements = 0;
    if( rank == 0 ){
      my_elements = numLocal*numprocs;
    }
    Epetra_Map root_map(-1,my_elements,0,*comm);
    MV root_mv(root_map, numVectors);
    Epetra_Import importer(root_map, map);
    root_mv.Import(*mv, importer, Insert);

    // root_mv.Print(std::cout);

    root_mv.ExtractCopy(original.getRawPtr(), numLocal*numprocs);

    TEST_EQUALITY( original, copy );


    // Check getting copy of just local data
    original.clear();
    original.resize(numVectors*numLocal);
    copy.clear();
    copy.resize(numVectors*numLocal);
    mv->Random();
    
    mv->ExtractCopy(original.getRawPtr(),mv->MyLength());
    get_1d_copy_helper<ADAPT,double>::do_get(ptrInArg(*adapter), copy(),
					     numLocal, Amesos::Util::Distributed);
    
    // Check that the values remain the same
    TEST_EQUALITY( original, copy );

    delete comm;
  }

  // Do not check Views, since their use is deprecated already

  TEUCHOS_UNIT_TEST( MultiVecAdapter, Globalize )
  {
    typedef ScalarTraits<double> ST;
    typedef Epetra_MultiVector MV;
    typedef MultiVecAdapter<MV> ADAPT;

    const Epetra_Comm* comm = getDefaultComm();
    int numprocs = comm->NumProc();
    int rank = comm->MyPID();

    // create a Map
    const size_t numVectors = 7;
    const size_t numLocal = 13;
    Epetra_Map map(-1,numLocal,0,*comm);

    RCP<MV> mv = rcp(new MV(map,numVectors));
    mv->Random();

    RCP<ADAPT> adapter = createMultiVecAdapter(mv);
    Array<double> original(numVectors*numLocal*numprocs);
    Array<double> copy(numVectors*numLocal*numprocs);

    if( rank == 0 ){
      std::fill(original.begin(), original.end(), 1.9);
    }

    // distribute rank 0's data
    put_1d_data_helper<ADAPT,double>::do_put(outArg(*adapter), original(),
					     numLocal*numprocs,
					     Amesos::Util::Rooted);

    // Send rank 0's array to everyone else
    comm->Broadcast(original.getRawPtr(), original.size(), 0);
    
    // Now have everyone get a copy from the multivector adapter
    get_1d_copy_helper<ADAPT,double>::do_get(ptrInArg(*adapter), copy(),
					     numLocal*numprocs,
					     Amesos::Util::Globally_Replicated);

    // Check that the values are the same
    TEST_EQUALITY( original, copy );

    delete comm;
  }

} // end anonymous namespace
