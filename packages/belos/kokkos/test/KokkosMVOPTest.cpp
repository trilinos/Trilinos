// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  This test uses the MVOPTester.hpp functions to test the Belos adapters
//  to Kokkos.
//


#include "BelosConfigDefs.hpp"
#include "BelosKokkosMVOPTester.hpp"
#include "BelosOutputManager.hpp"
#include "BelosKokkosAdapter.hpp"
#include "KokkosSparse_IOUtils.hpp"

#include "Teuchos_StandardCatchMacros.hpp"
#ifdef HAVE_MPI
  #include <mpi.h>
#endif

using std::cout;
using std::endl;
using Teuchos::RCP;
using Belos::Warnings;
using Belos::OutputManager;

template <class ScalarType>
bool TestKokkosMultiVecOneScalar(const Teuchos::RCP<OutputManager<ScalarType> >& );

template <class ScalarType1, class ScalarType2>
bool TestKokkosMultiVecTwoScalar(const Teuchos::RCP<OutputManager<ScalarType1> > & outputMgr);

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif
  bool isPassed;
  bool success = true;
  Kokkos::initialize();
  {
  bool verbose = false;
  if (argc>1) {
    if (argv[1][0]=='-' && argv[1][1]=='v') {
      verbose = true;
    }
  }
    typedef double ScalarType;
    typedef float ScalarType2;
    typedef Belos::MultiVec<ScalarType> KMV;
    typedef Belos::Operator<ScalarType> KOP; 
    typedef Belos::MultiVec<ScalarType2> KMV2;
    typedef Belos::Operator<ScalarType2> KOP2; 
    typedef Kokkos::DefaultExecutionSpace     EXSP;

    // Create an output manager to handle the I/O from the solver (defaults to std::cout).
    Teuchos::RCP<Belos::OutputManager<ScalarType> > myOutputMgr = Teuchos::rcp( new Belos::OutputManager<ScalarType>() );
    Teuchos::RCP<Belos::OutputManager<ScalarType2> > myOutputMgr2 = Teuchos::rcp( new Belos::OutputManager<ScalarType2>() );
    if (verbose) {
      myOutputMgr->setVerbosity( Warnings );
      myOutputMgr2->setVerbosity( Warnings );
    }
try {
    // number of global elements
    int dim = 10;
    int blockSize = 5;
    std::vector<ScalarType> norms(blockSize);


    //*****************************************************************************************************************************
    //*****************************************************************************************************************************
    // Run native Kokkos::Multivec test (double).
    isPassed = TestKokkosMultiVecOneScalar<ScalarType>(myOutputMgr);
    success = isPassed && success;
    if (isPassed) {
      myOutputMgr->print(Belos::Warnings,"*** KokkosAdapter PASSED TestKokkosMultiVecOneScalar() double. \n");
    }
    else {
      myOutputMgr->print(Belos::Warnings,"*** KokkosAdapter FAILED TestKokkosMultiVecOneScalar() double. ***\n\n");
    }

    // Run native Kokkos::Multivec test (single).
    isPassed = TestKokkosMultiVecOneScalar<ScalarType2>(myOutputMgr2);
    success = isPassed && success;
    if (isPassed) {
      myOutputMgr->print(Belos::Warnings,"*** KokkosAdapter PASSED TestKokkosMultiVecOneScalar() single. \n");
    }
    else {
      myOutputMgr->print(Belos::Warnings,"*** KokkosAdapter FAILED TestKokkosMultiVecOneScalar() single. ***\n\n");
    }

    //*****************************************************************************************************************************
    //*****************************************************************************************************************************
    // Run native Kokkos::Multivec tests with mixed single-double precisions:
    isPassed = TestKokkosMultiVecTwoScalar<ScalarType, ScalarType2>(myOutputMgr);
    success = isPassed && success;
    if (isPassed) {
      myOutputMgr->print(Belos::Warnings,"*** KokkosAdapter PASSED TestKokkosMultiVecTwoScalar() \n");
    }
    else {
      myOutputMgr->print(Belos::Warnings,"*** KokkosAdapter FAILED TestKokkosMultiVecTwoScalar() ***\n\n");
    }

    //*****************************************************************************************************************************
    //*****************************************************************************************************************************
    // Run Belos tests on Kokkos MultiVec (double precision):
    Teuchos::RCP<Belos::KokkosMultiVec<ScalarType>> ivec = Teuchos::rcp( new Belos::KokkosMultiVec<ScalarType>(dim, blockSize) );
    ivec->MvRandom();
    isPassed = Belos::TestKokkosMultiVecTraits<ScalarType,KMV>(myOutputMgr,ivec);
    success = isPassed && success;
    if (isPassed) {
      myOutputMgr->print(Belos::Warnings,"*** KokkosAdapter PASSED TestMultiVecTraits() double.\n");
    }
    else {
      myOutputMgr->print(Belos::Warnings,"*** KokkosAdapter FAILED TestMultiVecTraits() double. ***\n\n");
    }

    // Run Belos tests on Kokkos MultiVec (single precision):
    Teuchos::RCP<Belos::KokkosMultiVec<ScalarType2>> ivec_2 = Teuchos::rcp( new Belos::KokkosMultiVec<ScalarType2>(dim, blockSize) );
    ivec_2->MvRandom();
    isPassed = Belos::TestKokkosMultiVecTraits<ScalarType2,KMV2>(myOutputMgr2,ivec_2);
    success = isPassed && success;
    if (isPassed) {
      myOutputMgr->print(Belos::Warnings,"*** KokkosAdapter PASSED TestMultiVecTraits() single. \n");
    }
    else {
      myOutputMgr->print(Belos::Warnings,"*** KokkosAdapter FAILED TestMultiVecTraits() single. ***\n\n");
    }
    //*****************************************************************************************************************************
    //*****************************************************************************************************************************
    // Read in a matrix Market file and use it to test the Kokkos Operator (double precision).
    KokkosSparse::CrsMatrix<ScalarType, int, EXSP> crsMat = 
            KokkosSparse::Impl::read_kokkos_crst_matrix<KokkosSparse::CrsMatrix<ScalarType, int, EXSP>>("bcsstk13.mtx"); 
    Teuchos::RCP<Belos::KokkosCrsOperator<ScalarType, int, EXSP>> myOp = 
            Teuchos::rcp(new Belos::KokkosCrsOperator<ScalarType,int,EXSP>(crsMat));
    
    Teuchos::RCP<Belos::KokkosMultiVec<ScalarType>> ivec3 = Teuchos::rcp( new Belos::KokkosMultiVec<ScalarType>(2003, 2) );

    isPassed = Belos::TestKokkosOperatorTraits<ScalarType,KMV,KOP>(myOutputMgr,ivec3,myOp);
    success = isPassed && success;
    if (isPassed) {
      myOutputMgr->print(Belos::Warnings,"*** KokkosAdapter PASSED TestOperatorTraits() double. \n");
    }
    else {
      myOutputMgr->print(Belos::Warnings,"*** KokkosAdapter FAILED TestOperatorTraits() double ***\n\n");
    }
    // Read in a matrix Market file and use it to test the Kokkos Operator (single precision).
    KokkosSparse::CrsMatrix<ScalarType2, int, EXSP> crsMat2 = 
            KokkosSparse::Impl::read_kokkos_crst_matrix<KokkosSparse::CrsMatrix<ScalarType2, int, EXSP>>("bcsstk13.mtx"); 
    Teuchos::RCP<Belos::KokkosCrsOperator<ScalarType2, int, EXSP>> myOp2 = 
            Teuchos::rcp(new Belos::KokkosCrsOperator<ScalarType2,int,EXSP>(crsMat2));
    
    Teuchos::RCP<Belos::KokkosMultiVec<ScalarType2>> ivec3_2 = Teuchos::rcp( new Belos::KokkosMultiVec<ScalarType2>(2003, 2) );

    isPassed = Belos::TestKokkosOperatorTraits<ScalarType2,KMV2,KOP2>(myOutputMgr2,ivec3_2,myOp2);
    success = isPassed && success;
    if (isPassed) {
      myOutputMgr->print(Belos::Warnings,"*** KokkosAdapter PASSED TestOperatorTraits() single\n");
    }
    else {
      myOutputMgr->print(Belos::Warnings,"*** KokkosAdapter FAILED TestOperatorTraits() single ***\n\n");
    }
    //*****************************************************************************************************************************
    //*****************************************************************************************************************************
    
    if (!success) { 
      myOutputMgr->print(Belos::Warnings,"End Result: TEST FAILED\n");
    } else {
      myOutputMgr->print(Belos::Warnings,"End Result: TEST PASSED\n");
    }
    }//end try block


    TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,success);
  }
  Kokkos::finalize();
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}

template <class ScalarType>
bool TestKokkosMultiVecOneScalar(const Teuchos::RCP<OutputManager<ScalarType> > & outputMgr){
  int dim = 10;
  int blockSize = 5;
  std::vector<ScalarType> norms(blockSize);

  /// Test KokkosMultiVec constructors:
  // Test constructor #1:
  Belos::KokkosMultiVec<ScalarType> myVec1("myLabel", dim, blockSize);
  if ( myVec1.GetNumberVecs() != blockSize ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec constructor 1 returned wrong value "
      << "for GetNumberVecs()." << endl;
    return false;
  }
  if ( myVec1.GetGlobalLength() != dim ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec constructor 1 returned wrong value "
      << "for GetGlobalLength()." << endl;
    return false;
  }
  myVec1.MvNorm(norms); 
  for(int i = 0; i < blockSize; i++){
    if( norms[i] != 0 ){
      outputMgr->stream(Warnings)
        << "*** ERROR *** KokkosMultiVec constructor 1 returned wrong nrm2 value. "
        << "Vector was not initialized to zeros." << endl;
      return false;
    }
  }
  // Test constructor #2:
  Belos::KokkosMultiVec<ScalarType> myVec2(dim, blockSize);
  if ( myVec2.GetNumberVecs() != blockSize ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec constructor 2 returned wrong value "
      << "for GetNumberVecs()." << endl;
    return false;
  }
  if ( myVec2.GetGlobalLength() != dim ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec constructor 2 returned wrong value "
      << "for GetGlobalLength()." << endl;
    return false;
  }
  myVec2.MvNorm(norms); 
  for(int i = 0; i < blockSize; i++){
    if( norms[i] != 0 ){
      outputMgr->stream(Warnings)
        << "*** ERROR *** KokkosMultiVec constructor 2 returned wrong nrm2 value. "
        << "Vector was not initialized to zeros." << endl;
      return false;
    }
  }
  // Test constructor #3:
  Belos::KokkosMultiVec<ScalarType> myVec3(2*dim);
  if ( myVec3.GetNumberVecs() != 1 ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec constructor 3 returned wrong value "
      << "for GetNumberVecs()." << endl;
    return false;
  }
  if ( myVec3.GetGlobalLength() != 2*dim ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec constructor 3 returned wrong value "
      << "for GetGlobalLength()." << endl;
    return false;
  }
  myVec3.MvNorm(norms); 
  if( norms[0] != 0 ){
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec constructor 3 returned wrong nrm2 value. "
      << "Vector was not initialized to zeros." << endl;
    return false;
  }
  // Test copy constructor (should deep copy).
  Belos::KokkosMultiVec<ScalarType> myVecCopy(myVec3);
  if ( myVecCopy.GetNumberVecs() != 1 ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec copy constructor returned wrong value "
      << "for GetNumberVecs()." << endl;
    return false;
  }
  if ( myVecCopy.GetGlobalLength() != 2*dim ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec copy constructor returned wrong value "
      << "for GetGlobalLength()." << endl;
    return false;
  }
  myVecCopy.MvRandom();
  myVecCopy.MvNorm(norms);
  if( norms[0] == 0 ){
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec MvRandom did not fill with random values. " << endl;
    return false;
  }
  myVec3.MvNorm(norms);
  if( norms[0] != 0 ){
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec copy constructor did not deep copy. " << endl;
    return false;
  }
  // Test assignment operator (should also deep copy). 
  myVecCopy = myVec2;
  if ( myVecCopy.GetNumberVecs() != blockSize ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec assignment = returned wrong value "
      << "for GetNumberVecs()." << endl;
    return false;
  }
  if ( myVecCopy.GetGlobalLength() != dim ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec assignment = returned wrong value "
      << "for GetGlobalLength()." << endl;
    return false;
  }
  myVec2.MvInit(3.0);
  myVecCopy.MvNorm(norms); 
  for(int i = 0; i < blockSize; i++){
    if( norms[i] != 0 ){
      outputMgr->stream(Warnings)
        << "*** ERROR *** KokkosMultiVec assignment = returned wrong nrm2 value. "
        << "Vector was not deep copied." << endl;
      return false;
    }
  }
  // Test view to multivec:
  int numCols2 = 4;
  int numRows2 = 60;
  Kokkos::View<ScalarType**, Kokkos::LayoutLeft> myView("View2MV", numRows2, numCols2);
  typename Kokkos::View<ScalarType**>::HostMirror myView_h("View2MV_host", numRows2, numCols2);
  Kokkos::deep_copy(myView, 42);
  Belos::KokkosMultiVec<ScalarType> myVec4( myView );
  if ( myVec4.GetNumberVecs() != (int)myView.extent(1) ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec view to multivec returned wrong value "
      << "for GetNumberVecs()." << endl;
    return false;
  }
  if ( myVec4.GetGlobalLength() != (int)myView.extent(0) ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec view to multivec returned wrong value "
      << "for GetGlobalLength()." << endl;
    return false;
  }
  Kokkos::deep_copy(myView, 55);
  Kokkos::deep_copy(myView_h, myVec4.GetInternalViewConst());
  if ( myView_h(5,1) != 42 ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultivec view to multivec did not make a deep copy!" << endl;
    return false;
  }
  // Tesst view to multivec with shallow copy:
  Kokkos::deep_copy(myView, 100);
  Belos::KokkosMultiVec<ScalarType> myVec5( myView, false );
  if ( myVec5.GetNumberVecs() != (int)myView.extent(1) ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec view to multivec shallow returned wrong value "
      << "for GetNumberVecs()." << endl;
    return false;
  }
  if ( myVec5.GetGlobalLength() != (int)myView.extent(0) ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec view to multivec shallow returned wrong value "
      << "for GetGlobalLength()." << endl;
    return false;
  }
  Kokkos::deep_copy(myView, 500);
  Kokkos::deep_copy(myView_h, myVec5.GetInternalViewConst());
  if ( myView_h(5,1) != 500 ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultivec view to multivec shallow made a deep copy!" << endl;
    return false;
  }
  // Test GetInternalViewNonConst:
  auto myView2 = myVec5.GetInternalViewNonConst();
  Kokkos::deep_copy(myView2, 0);
  std::vector<ScalarType> norms2(4);
  myVec5.MvNorm(norms2);
  for(int i = 0; i < (int)myView2.extent(1); i++){
    if( norms[i] != 0 ){
      outputMgr->stream(Warnings)
        << "*** ERROR *** KokkosMultiVec GetInternalViewNonConst returned wrong nrm2 value. "
        << "Vector was not editable." << endl;
      return false;
    }
  }

  return true;
}


template <class ScalarType1, class ScalarType2>
bool TestKokkosMultiVecTwoScalar(const Teuchos::RCP<OutputManager<ScalarType1> > & outputMgr){
  // Test view to multivec ST1 to ST2:
  int numCols = 4;
  int numRows = 60;
  Kokkos::View<ScalarType1**, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace> myViewST1("View2MV1", numRows, numCols);
  typename Kokkos::View<ScalarType2**>::HostMirror myViewST2_h("View2MV1_host", numRows, numCols);
  Kokkos::deep_copy(myViewST1, 42);
  Belos::KokkosMultiVec<ScalarType2> myVecST2A( myViewST1 );
  if ( myVecST2A.GetNumberVecs() != (int)myViewST1.extent(1) ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec view to multivec 2 scalar returned wrong value "
      << "for GetNumberVecs()." << endl;
    return false;
  }
  if ( myVecST2A.GetGlobalLength() != (int)myViewST1.extent(0) ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec view to multivec 2 scalar returned wrong value "
      << "for GetGlobalLength()." << endl;
    return false;
  }
  Kokkos::deep_copy(myViewST1, 55);
  std::cout << "At line 424. " << std::endl;
  Kokkos::deep_copy(myViewST2_h, myVecST2A.GetInternalViewConst());
  if ( myViewST2_h(5,1) != 42 ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultivec view to multivec 2 scalar did not make a deep copy!" << endl;
    return false;
  }
  // Test view to multivec ST2 to ST1:
  Kokkos::View<ScalarType2**, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace> myViewST2("View2MV2", numRows, numCols);
  typename Kokkos::View<ScalarType1**>::HostMirror myViewST1_h("View2MV2_host", numRows, numCols);
  Kokkos::deep_copy(myViewST2, 56);
  Belos::KokkosMultiVec<ScalarType1> myVecST1A( myViewST2 );
  if ( myVecST1A.GetNumberVecs() != (int)myViewST2.extent(1) ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec view to multivec 2 scalar returned wrong value "
      << "for GetNumberVecs()." << endl;
    return false;
  }
  if ( myVecST1A.GetGlobalLength() != (int)myViewST2.extent(0) ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec view to multivec 2 scalar returned wrong value "
      << "for GetGlobalLength()." << endl;
    return false;
  }
  Kokkos::deep_copy(myViewST2, 22);
  std::cout << "At line 449. " << std::endl;
  Kokkos::deep_copy(myViewST1_h, myVecST1A.GetInternalViewConst());
  if ( myViewST1_h(5,1) != 56 ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultivec view to multivec 2 scalar did not make a deep copy!" << endl;
    return false;
  }
  //Test multivec copy constructor ST1 to ST2:
  Belos::KokkosMultiVec<ScalarType2> myVecST2B( myVecST1A );
  if ( myVecST1A.GetNumberVecs() != myVecST2B.GetNumberVecs() ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec templated copy constructor returned wrong value "
      << "for GetNumberVecs()." << endl;
    return false;
  }
  if ( myVecST1A.GetGlobalLength() != myVecST2B.GetGlobalLength() ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec templated copy constructor returned wrong value "
      << "for GetGlobalLength()." << endl;
    return false;
  }
  myVecST1A.MvInit(3.0);
  std::cout << "At line 471. " << std::endl;
  Kokkos::deep_copy(myViewST2_h, myVecST2B.GetInternalViewConst());
  if ( myViewST2_h(5,1) != 56 ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultivec templated copy constructor did not make a deep copy!" << endl;
    return false;
  }
  //Test multivec copy constructor ST2 to ST1:
  Belos::KokkosMultiVec<ScalarType1> myVecST1B( myVecST2A );
  if ( myVecST1B.GetNumberVecs() != myVecST2A.GetNumberVecs() ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec templated copy constructor returned wrong value "
      << "for GetNumberVecs()." << endl;
    return false;
  }
  if ( myVecST1B.GetGlobalLength() != myVecST2A.GetGlobalLength() ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec templated copy constructor returned wrong value "
      << "for GetGlobalLength()." << endl;
    return false;
  }
  myVecST2A.MvInit(3.0);
  Kokkos::deep_copy(myViewST1_h, myVecST1B.GetInternalViewConst());
  if ( myViewST1_h(5,1) != 42 ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultivec templated copy constructor did not make a deep copy!" << endl;
    return false;
  }
  //Test assignment operator ST1 to ST2:
  Belos::KokkosMultiVec<ScalarType2> myVecST2C  = myVecST1A; 
  if ( myVecST1A.GetNumberVecs() != myVecST2C.GetNumberVecs() ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec templated assignment operator returned wrong value "
      << "for GetNumberVecs()." << endl;
    return false;
  }
  if ( myVecST1A.GetGlobalLength() != myVecST2C.GetGlobalLength() ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec templated copy constructor returned wrong value "
      << "*** ERROR *** KokkosMultiVec templated assignment operator returned wrong value "
      << "for GetGlobalLength()." << endl;
    return false;
  }
  myVecST1A.MvInit(1.0);
  std::cout << "At line 515. " << std::endl;
  Kokkos::deep_copy(myViewST2_h, myVecST2C.GetInternalViewConst());
  if ( myViewST2_h(5,1) != 3 ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultivec templated assign op did not make a deep copy!" << endl;
    return false;
  }
  //Test assignement operator ST2 to ST1:
  Belos::KokkosMultiVec<ScalarType1> myVecST1C  = myVecST2A; 
  if ( myVecST1C.GetNumberVecs() != myVecST2A.GetNumberVecs() ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec templated assignment operator returned wrong value "
      << "for GetNumberVecs()." << endl;
    return false;
  }
  if ( myVecST1C.GetGlobalLength() != myVecST2A.GetGlobalLength() ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultiVec templated assignment operator returned wrong value "
      << "for GetGlobalLength()." << endl;
    return false;
  }
  myVecST2A.MvInit(2.0);
  Kokkos::deep_copy(myViewST1_h, myVecST1C.GetInternalViewConst());
  if ( myViewST1_h(5,1) != 3 ) {
    outputMgr->stream(Warnings)
      << "*** ERROR *** KokkosMultivec templated assign op did not make a deep copy!" << endl;
    return false;
  }

  return true;
}
