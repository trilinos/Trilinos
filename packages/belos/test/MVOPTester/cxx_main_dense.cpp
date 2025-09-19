// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
#include "BelosDenseMatTraits.hpp"
#include "BelosTeuchosDenseAdapter.hpp"
#include "BelosDenseMatTester.hpp"
#include "BelosOutputManager.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_RCP.hpp"

#include "Teuchos_StandardCatchMacros.hpp"

#define OTYPE int

#if defined(HAVE_TEUCHOS_COMPLEX)
#define STYPE std::complex<double>
#else
#define STYPE double
#endif

template<typename TYPE>
int PrintTestResults(std::string, Teuchos::RCP<TYPE>, Teuchos::RCP<TYPE>, bool);

template<typename TYPE>
int PrintTestResults2(std::string, TYPE, TYPE, bool);

int main(int argc, char* argv[])
{
  typedef Teuchos::SerialDenseMatrix<int, STYPE> DMatrix;
  typedef Belos::DenseMatTraits<STYPE, Teuchos::SerialDenseMatrix<int,STYPE> > DMT;

  int i;
  bool verbose = 0;
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  if (verbose)
    std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  bool success = true;
  try {
    bool ierr=false;
    int numberFailedTests = 0;
    std::string testName = "";
    STYPE zero = Teuchos::ScalarTraits<STYPE>::zero();
    STYPE one = Teuchos::ScalarTraits<STYPE>::one();

    // Create an output manager to handle the I/O from the solver
    Teuchos::RCP<Belos::OutputManager<double> > MyOM = Teuchos::rcp( new Belos::OutputManager<double>() );
    Teuchos::RCP<Belos::OutputManager<float> > MyOMFloat = Teuchos::rcp( new Belos::OutputManager<float>() );
#if defined(HAVE_TEUCHOS_COMPLEX)
    Teuchos::RCP<Belos::OutputManager<std::complex<double>> > MyOMCplxDouble = Teuchos::rcp( new Belos::OutputManager<std::complex<double>>() );
    Teuchos::RCP<Belos::OutputManager<std::complex<float>> > MyOMCplxFloat = Teuchos::rcp( new Belos::OutputManager<std::complex<float>>() );
#endif
    if (verbose) {
      MyOM->setVerbosity( Belos::Warnings );
      MyOMFloat->setVerbosity( Belos::Warnings );
#if defined(HAVE_TEUCHOS_COMPLEX)
      MyOMCplxDouble->setVerbosity( Belos::Warnings );
      MyOMCplxFloat->setVerbosity( Belos::Warnings );
#endif
    }

    //*********************************************************************
    // Teuchos SerialDense MatrixTraits impl testing.
    //*********************************************************************
    ierr = Belos::TestDenseMatTraits<double,Teuchos::SerialDenseMatrix<int,double>>(MyOM);
    if (ierr) {
      MyOM->print(Belos::Warnings,"*** TeuchosAdapter PASSED TestDenseMatTraits() scalar double \n");
    }
    else {
      MyOM->print(Belos::Warnings,"*** TeuchosAdapter FAILED TestDenseMatTraits() scalar double ***\n\n");
      numberFailedTests++;
    }

    //*********************************************************************
    // Teuchos SerialDense MatrixTraits impl testing.
    //*********************************************************************
    ierr = Belos::TestDenseMatTraits<float,Teuchos::SerialDenseMatrix<int,float>>(MyOMFloat);
    if (ierr) {
      MyOMFloat->print(Belos::Warnings,"*** TeuchosAdapter PASSED TestDenseMatTraits() scalar float \n");
    }
    else {
      MyOMFloat->print(Belos::Warnings,"*** TeuchosAdapter FAILED TestDenseMatTraits() scalar float ***\n\n");
      numberFailedTests++;
    }

#if defined(HAVE_TEUCHOS_COMPLEX)
    //*********************************************************************
    // Teuchos SerialDense MatrixTraits impl testing.
    //*********************************************************************
    ierr = Belos::TestDenseMatTraits<std::complex<double>,Teuchos::SerialDenseMatrix<int,std::complex<double>>>(MyOMCplxDouble);
    if (ierr) {
      MyOMCplxDouble->print(Belos::Warnings,"*** TeuchosAdapter PASSED TestDenseMatTraits() scalar complex double \n");
    }
    else {
      MyOMCplxDouble->print(Belos::Warnings,"*** TeuchosAdapter FAILED TestDenseMatTraits() scalar complex double ***\n\n");
      numberFailedTests++;
    }

    //*********************************************************************
    // Teuchos SerialDense MatrixTraits impl testing.
    //*********************************************************************
    ierr = Belos::TestDenseMatTraits<std::complex<float>,Teuchos::SerialDenseMatrix<int,std::complex<float>>>(MyOMCplxFloat);
    if (ierr) {
      MyOMCplxFloat->print(Belos::Warnings,"*** TeuchosAdapter PASSED TestDenseMatTraits() complex float\n");
    }
    else {
      MyOMCplxFloat->print(Belos::Warnings,"*** TeuchosAdapter FAILED TestDenseMatTraits() complex float ***\n\n");
      numberFailedTests++;
    }
#endif

    if (verbose) std::cout<<std::endl<<"********** CHECKING TEUCHOS SERIAL DENSE MATRIX **********"<<std::endl<<std::endl;

    // default constructor test
    Teuchos::RCP<DMatrix> DefConTest = DMT::Create();
    if (verbose) std::cout <<"default constructor -- construct empty matrix ";
    if ( DMT::GetRawHostPtr( *DefConTest )!=NULL || DMT::GetNumRows(*DefConTest)!=0 || DMT::GetNumCols(*DefConTest)!=0 || DMT::GetStride(*DefConTest)!=0 ) {
	if (verbose) std::cout << "[ unsuccessful ]"<<std::endl;
	numberFailedTests++;
    } else {
	if (verbose) std::cout << "[ successful ]"<<std::endl;
    }

    // constructor 1 (matrix w/ dimension but empty)

    Teuchos::RCP<DMatrix> Con1Test = DMT::Create( 3, 4 );
    if (verbose) std::cout <<"constructor 1 -- empty matrix with given dimensions ";
    if ( DMT::GetNumRows(*Con1Test)!=3 || DMT::GetNumCols(*Con1Test)!=4 || DMT::Value(*Con1Test, 1, 2 )!=zero || DMT::ValueConst(*Con1Test, 1, 2 )!=zero ) {
	if (verbose) std::cout << "[ unsuccessful ]"<<std::endl;
        numberFailedTests++;
    } else {
        if (verbose) std::cout << "[ successful ]"<<std::endl;
    }
	
    // constructor 2 (from array) tests

    Teuchos::RCP<DMatrix> A33 = DMT::Create(3, 3);
    STYPE* A33data = DMT::GetRawHostPtr(*A33);
    for(i = 0; i < 9; i++)
    {
      A33data[i] = i;
    }
    DMT::SyncHostToDevice(*A33);

    Teuchos::RCP<DMatrix> Con2Test1ExpRes = DMT::Create();
    DMT::Reshape(*Con2Test1ExpRes, 2, 3);
    DMT::Value(*Con2Test1ExpRes, 0, 0) = 0;  
    DMT::Value(*Con2Test1ExpRes, 0, 1) = 3; 
    DMT::Value(*Con2Test1ExpRes, 0, 2) = 6;
    DMT::Value(*Con2Test1ExpRes, 1, 0) = 1;  
    DMT::Value(*Con2Test1ExpRes, 1, 1) = 4; 
    DMT::Value(*Con2Test1ExpRes, 1, 2) = 7;
    DMT::SyncHostToDevice(*Con2Test1ExpRes);

    Teuchos::RCP<DMatrix> Con2Test1 = DMT::SubviewCopy(*A33, 2, 3);
    numberFailedTests += PrintTestResults("constructor 2 -- construct matrix from array subrange", Con2Test1, Con2Test1ExpRes, verbose);

    // constructor 3 (submatrix)

    Teuchos::RCP<DMatrix> Con4TestOrig = DMT::SubviewCopy(*A33, 3, 3);
    Teuchos::RCP<DMatrix> Con4TestSubmatrix = DMT::Create(2, 2);
    DMT::Value(*Con4TestSubmatrix, 0, 0) = 4;
    DMT::Value(*Con4TestSubmatrix, 0, 1) = 7;
    DMT::Value(*Con4TestSubmatrix, 1, 0) = 5;
    DMT::Value(*Con4TestSubmatrix, 1, 1) = 8;
    DMT::SyncHostToDevice(*Con4TestSubmatrix);
    Teuchos::RCP<DMatrix> Con4TestCopy1 = DMT::SubviewCopy(*Con4TestOrig, 2, 2, 1, 1);
    numberFailedTests += PrintTestResults("constructor 4 -- submatrix copy", Con4TestCopy1, Con4TestSubmatrix, verbose);
    Teuchos::RCP<DMatrix> Con4TestCopy2 = DMT::SubviewCopy(*Con4TestOrig, 3, 3, 0, 0);
    numberFailedTests += PrintTestResults("constructor 4 -- full matrix copy", Con4TestCopy2, Con4TestOrig, verbose);
    Teuchos::RCP<DMatrix> Con4TestView1 = DMT::Subview(*Con4TestOrig, 2, 2, 1, 1);
    numberFailedTests += PrintTestResults("constructor 4 -- full matrix view", Con4TestView1, Con4TestSubmatrix, verbose);
    Teuchos::RCP<DMatrix> Con4TestView2 = DMT::Subview(*Con4TestOrig, 3, 3, 0, 0);
    numberFailedTests += PrintTestResults("constructor 4 -- submatrix view", Con4TestView2, Con4TestOrig, verbose);

    // Norm Tests

    Teuchos::RCP<DMatrix> AAA = DMT::Create(3, 3);
    DMT::Value(*AAA, 0, 0) = 1; 
    DMT::Value(*AAA, 0, 1) = 2; 
    DMT::Value(*AAA, 0, 2) = 3;
    DMT::Value(*AAA, 1, 0) = 4; 
    DMT::Value(*AAA, 1, 1) = 5; 
    DMT::Value(*AAA, 1, 2) = 6;
    DMT::Value(*AAA, 2, 0) = 7; 
    DMT::Value(*AAA, 2, 1) = 8; 
    DMT::Value(*AAA, 2, 2) = 9;
    DMT::SyncHostToDevice(*AAA);

    Teuchos::RCP<DMatrix> BBB = DMT::Create();
    numberFailedTests += PrintTestResults2("normOne of a 3x3", DMT::NormOne(*AAA), 18.0, verbose);
    DMT::PutScalar(*AAA, one);
    numberFailedTests += PrintTestResults2("normOne of a 3x3 (ones)", DMT::NormOne(*AAA), 3.0, verbose);
    numberFailedTests += PrintTestResults2("normFrobenius of a 3x3", DMT::NormFrobenius(*AAA), 3.0, verbose);
    numberFailedTests += PrintTestResults2("normOne of a 0x0", DMT::NormOne(*BBB), 0.0, verbose);
    numberFailedTests += PrintTestResults2("normFrobenius of a 0x0", DMT::NormFrobenius(*BBB), 0.0, verbose);

    //  Check scale methods.
  
    Teuchos::RCP<DMatrix> ScalTest = DMT::Create( 8, 8 );
    DMT::PutScalar(*ScalTest, one);

    //  Scale the entries by 8, it should be 8.
    if (verbose) std::cout << "scale() -- scale matrix by some number ";
    DMT::Scale(*ScalTest, 8.0);
    if (DMT::Value(*ScalTest, 2, 3) == 8.0) {
	if (verbose) std::cout<< "[ successful ]" <<std::endl;
    } else {
	if (verbose) std::cout<< "[ unsuccessful ]" <<std::endl;
	numberFailedTests++;
    }
    //  Pointwise scale the entries by zero, they all should be zero.
    if (verbose) std::cout << "scale() -- scale matrix by zero ";
    DMT::Scale(*ScalTest, zero);
    if (DMT::NormOne(*ScalTest) == 0.0) {
	if (verbose) std::cout<< "[ successful ]" <<std::endl;
    } else {
	if (verbose) std::cout<< "[ unsuccessful ]" <<std::endl;
	numberFailedTests++;
    }
    //
    //  Check set methods.
    //
    Teuchos::RCP<DMatrix> CCC = DMT::Create( 5, 5 );
    //  Randomize the entries in CCC.
    if (verbose) std::cout << "random() -- enter random entries into matrix ";
    DMT::Randomize(*CCC);
    if (DMT::NormOne(*CCC) != 0.0) {
	if (verbose) std::cout<< "[ successful ]" <<std::endl;
    } else {
	if (verbose) std::cout<< "[ unsuccessful ]" <<std::endl;
	numberFailedTests++;
    }
    //  Set the entries of CCC to 1.0.
    if (verbose) std::cout << "putScalar() -- set every entry of this matrix to 1.0 ";
    DMT::PutScalar(*CCC, one);
    if ( DMT::Value(*CCC, 3, 4) == one ) {
      if (verbose) std::cout<< "[ successful ]" <<std::endl;
    } else {
      if (verbose) std::cout<< "[ unsuccessful ]" <<std::endl;
      numberFailedTests++;
    }
    //  Check assignment operator.
    Teuchos::RCP<DMatrix> CCC2 = DMT::Create( 5, 5 );
    DMT::Assign(*CCC2, *CCC);
    if (verbose) std::cout <<  "assign() -- copy the values of an input matrix ";
    if ( DMT::Value(*CCC2, 3, 4) == one ) {
      if (verbose) std::cout<< "[ successful ]" <<std::endl;
    } else {
      if (verbose) std::cout<< "[ unsuccessful ]" <<std::endl;
      numberFailedTests++;
    }
    //  Create a view into a submatrix of CCC, add them together
    Teuchos::RCP<DMatrix> CCCview = DMT::Subview( *CCC, 3, 3 );
    Teuchos::RCP<DMatrix> CCCview2 = DMT::Subview( *CCC, 3, 3, 1, 1 );

    if (verbose) std::cout << "view copy large(orig) -> small(view) ";
    if (DMT::GetNumRows(*CCCview)==3 && (DMT::GetRawHostPtr(*CCCview)==DMT::GetRawHostPtr(*CCC))) {
      if (verbose) std::cout<< "[ successful ]" <<std::endl;
    } else {
      if (verbose) std::cout<< "[ unsuccessful ]" <<std::endl;
      numberFailedTests++;
    }
    if (verbose) std::cout << "view copy addition ";
    STYPE sum = DMT::Value(*CCCview, 2, 2) + DMT::Value(*CCCview2, 2, 2);
    DMT::Add(*CCCview, *CCCview2);
    if ((DMT::Value(*CCCview, 2, 2) == DMT::Value(*CCC, 2, 2)) && (DMT::Value(*CCC, 2, 2) == sum)) {
      if (verbose) std::cout<< "[ successful ]"<<std::endl;
    } else {
      if (verbose) std::cout<< "[ unsuccessful ]"<<std::endl;
      numberFailedTests++;
    }

/*  Commenting out here!
  DMatrix CCCtest2( 2, 2 );
  CCCtest2 = 3.0;
  CCCtest1 = CCCtest2;
  if (verbose) std::cout << "operator= -- large(copy) = small(copy) ";
  if (CCCtest1.numRows()==2 ) {
    if (verbose) std::cout<< "[ successful ]"<<std::endl;
  } else {
    if (verbose) std::cout<< "[ unsuccessful ]"<<std::endl;
    numberFailedTests++;
  }
  CCCtest1 = CCCview;
  if (verbose) std::cout << "operator= -- large(copy) = small(view) ";
  if (CCCtest1.numRows()==3 && CCCtest1.stride()==5) {
    if(verbose) std::cout<<"[ successful ]" <<std::endl;
  } else {
    if (verbose) std::cout<<"[ unsuccessful ]"<<std::endl;
    numberFailedTests++;
  }

  DMatrix CCCtest3( CCCview );
  CCCtest1 += CCCtest3;
  if (verbose) std::cout << "operator+= -- add two matrices of the same size, but different leading dimension ";
  if (CCCtest1(1,1)==2.0) {
    if(verbose) std::cout<<"[ successful ]" <<std::endl;
  } else {
    if (verbose) std::cout<<"[ unsuccessful ]"<<std::endl;
    numberFailedTests++;
  }
  if (verbose) std::cout << "operator+= -- add two matrices of different size (nothing should change) ";
  CCCtest1 += CCC;
  if (CCCtest1(1,1)==2.0) {
    if(verbose) std::cout<<"[ successful ]" <<std::endl;
  } else {
    if (verbose) std::cout<<"[ unsuccessful ]"<<std::endl;
    numberFailedTests++;
  }

  if (verbose) std::cout<<std::endl<<"********** CHECKING TEUCHOS SERIAL DENSE VECTOR **********"<<std::endl<<std::endl;

  DVector DefConTestV;
  if (verbose) std::cout <<"default constructor -- construct empty std::vector ";
  if ( DefConTestV.values()!=NULL || DefConTestV.length()!=0 || DefConTestV.numRows()!=0 ||DefConTestV.stride()!=0 ) {
        if (verbose) std::cout << "[ unsuccessful ]"<<std::endl;
        numberFailedTests++;
  } else {
        if (verbose) std::cout << "[ successful ]"<<std::endl;
  }

  // constructor 1 (matrix w/ dimension but empty)

  DVector Con1TestV( 3 );
  if (verbose) std::cout <<"constructor 1 -- empty std::vector with given dimensions ";
  if ( Con1TestV.length()!=3 || Con1TestV.numCols()!=1 || Con1TestV( 1 )!=0.0 ) {
	if (verbose) std::cout << "[ unsuccessful ]"<<std::endl;
        numberFailedTests++;
  } else {
      	if (verbose) std::cout << "[ successful ]"<<std::endl;
  }
	
  // constructor 2 (from array) tests

  DVector Con2Test1V(Teuchos::Copy, a, 4);
  if (verbose) std::cout <<"constructor 2 -- construct std::vector from array subrange ";
  if ( Con2Test1V.numRows()!=4 || Con2Test1V.numCols()!=1 || Con2Test1V[ 2 ]!=2.0 ) {
        if (verbose) std::cout << "[ unsuccessful ]"<<std::endl;
        numberFailedTests++;
  } else {
        if (verbose) std::cout << "[ successful ]"<<std::endl;
  }

  // constructor 3 (copy constructor)

  DVector Con3TestCopyV( Con2Test1V );
  if(verbose) std::cout <<"constructor 3 -- copy constructor ";
  if ( Con3TestCopyV != Con2Test1V ) {
        if (verbose) std::cout << "[ unsuccessful ]"<<std::endl;
        numberFailedTests++;
  } else {
        if (verbose) std::cout << "[ successful ]"<<std::endl;
  }

  // non-member helper function (construct vector view of matrix column)

  OTYPE col = Teuchos::OrdinalTraits<OTYPE>::one();
  DVector ColViewTestV = Teuchos::getCol<OTYPE,STYPE>( Teuchos::View, AAA, col );
  if (verbose) std::cout <<"non-method helper function -- construct vector view of second column of matrix ";
  if ( ColViewTestV.normInf() != 1.0 || ColViewTestV.normOne() != 3.0 ) {
        if (verbose) std::cout << "[ unsuccessful ]"<<std::endl;
        numberFailedTests++;
  } else {
        if (verbose) std::cout << "[ successful ]"<<std::endl;
  }

  // checking norms

  numberFailedTests += PrintTestResults("normOne of a 3x1 std::vector", Con2Test1V.normOne(), 6.0, verbose);
  numberFailedTests += PrintTestResults("normInf of a 3x1 std::vector", Con2Test1V.normInf(), 3.0, verbose);
  Con2Test1V = Teuchos::ScalarTraits<STYPE>::one();
  numberFailedTests += PrintTestResults("normFrobenius of a 3x1 std::vector", Con2Test1V.normFrobenius(), 2.0, verbose);

  // check size/resize

  DVector SizeTestV1;
  SizeTestV1.size( 5 );
  if(verbose) std::cout <<"size() -- test ";
  if (SizeTestV1( 4 )!= 0.0) {
    if (verbose) std::cout << "[ unsuccessful ]"<<std::endl;
    numberFailedTests++;
  } else {
    if (verbose) std::cout << "[ successful ]"<<std::endl;
  }
  SizeTestV1 = 2.0*Teuchos::ScalarTraits<STYPE>::one();
  SizeTestV1.resize( 10 );
  if(verbose) std::cout <<"resize() -- test small --> large ";
  if (SizeTestV1[ 4 ]!= 2.0 || SizeTestV1[ 8 ]!=0.0 ) {
    if (verbose) std::cout << "[ unsuccessful ]"<<std::endl;
    numberFailedTests++;
  } else {
    if (verbose) std::cout << "[ successful ]"<<std::endl;
  }
  SizeTestV1.resize( 3 );
  if(verbose) std::cout <<"resize() -- test large --> small ";
  if (SizeTestV1( 2 )!= 2.0) {
    if (verbose) std::cout << "[ unsuccessful ]"<<std::endl;
    numberFailedTests++;
  } else {
    if (verbose) std::cout << "[ successful ]"<<std::endl;
  }

  DVector OpEqTestV1( 10 ); OpEqTestV1 = 3.0*Teuchos::ScalarTraits<STYPE>::one();
  DVector OpEqTestV2( Teuchos::View, OpEqTestV1.values(), 3 );
  DVector OpEqTestV3( 2 );
  OpEqTestV3 = OpEqTestV2;
  if (verbose) std::cout << "operator= -- small(empty) = large(view) ";
  if (OpEqTestV3.length()==3 && OpEqTestV3.values()==OpEqTestV2.values()) {
    if (verbose) std::cout<< "successful"<<std::endl;
  } else {
    if (verbose) std::cout<< "unsuccessful"<<std::endl;
    numberFailedTests++;
  }
  OpEqTestV3 = OpEqTestV1;
  if (verbose) std::cout << "operator= -- small(view) = large(copy) ";
  if (OpEqTestV3.length()==10 && OpEqTestV3.values()!=OpEqTestV1.values()) {
    if (verbose) std::cout<< "successful"<<std::endl;
  } else {
    if (verbose) std::cout<< "unsuccessful"<<std::endl;
    numberFailedTests++;
  }
  OpEqTestV3.size(5);
  OpEqTestV3 = OpEqTestV1;
  if (verbose) std::cout << "operator= -- small(copy) = large(copy) ";
  if (OpEqTestV3.length()==10 && OpEqTestV3.values()!=OpEqTestV1.values() && OpEqTestV3[ 9 ]==3.0) {
    if (verbose) std::cout<< "successful"<<std::endl;
  } else {
    if (verbose) std::cout<< "unsuccessful"<<std::endl;
    numberFailedTests++;
  }

  DVector OpSumTestV1( OpEqTestV2 );
  OpSumTestV1 += OpEqTestV2;
  if (verbose) std::cout << "operator+= -- add two vectors of the same size, but different leading dimension ";
  if (OpSumTestV1( 1 )==6.0) {
    if (verbose) std::cout<<"successful" <<std::endl;
  } else {
    if (verbose) std::cout<<"unsuccessful"<<std::endl;
    numberFailedTests++;
  }
  if (verbose) std::cout << "operator+= -- add two vectors of different size (nothing should change) ";
  OpSumTestV1 += OpEqTestV1;
  if (OpSumTestV1( 1 )==6.0) {
    if (verbose) std::cout<<"successful" <<std::endl;
  } else {
    if (verbose) std::cout<<"unsuccessful"<<std::endl;
    numberFailedTests++;
  }

  DVector OpCompTestV1( 5 );
  OpCompTestV1 = 2.0*Teuchos::ScalarTraits<STYPE>::one();
  if(verbose) std::cout <<"operator== -- test large == small ";
  if (OpCompTestV1 == SizeTestV1) {
    if (verbose) std::cout << "[ unsuccessful ]"<<std::endl;
    numberFailedTests++;
  } else {
    if (verbose) std::cout << "[ successful ]"<<std::endl;
  }
  if(verbose) std::cout <<"operator!= -- test large != small ";
  if (OpCompTestV1 != SizeTestV1) {
    if (verbose) std::cout << "[ successful ]"<<std::endl;
  } else {
    if (verbose) std::cout << "[ successful ]"<<std::endl;
    numberFailedTests++;
  }

  DVector ColSetTestV( AAA.numRows() );
  ColSetTestV.putScalar( 2.0 );
  bool ret = Teuchos::setCol<OTYPE,STYPE>( ColSetTestV, col, AAA );
  if (verbose) std::cout <<"non-method helper function -- set second column of matrix with vector ";
  if ( ColViewTestV.normInf() != 2.0 || ColViewTestV.normOne() != 6.0 || ret == false ) {
        if (verbose) std::cout << "[ unsuccessful ]"<<std::endl;
        numberFailedTests++;
  } else {
        if (verbose) std::cout << "[ successful ]"<<std::endl;
  }
*/
    //
    // If a test failed output the number of failed tests.
    //
    if (numberFailedTests>0) {
      success = false;
      MyOM->print(Belos::Warnings,"End Result: TEST FAILED\n");
    } else {
      success = true;
      MyOM->print(Belos::Warnings,"End Result: TEST PASSED\n");
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,success);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}

template<typename TYPE>
int PrintTestResults(std::string testName, Teuchos::RCP<TYPE> calculatedResult, Teuchos::RCP<TYPE> expectedResult, bool verbose)
{
  int result;
  if((*calculatedResult) == (*expectedResult))
    {
      if(verbose) std::cout << testName << " [ successful ]" << std::endl;
      result = 0;
    }
  else
    {
      if(verbose) std::cout << testName << " [ unsuccessful ]" << std::endl;
      result = 1;
    }
  return result;
}

template<typename TYPE>
int PrintTestResults2(std::string testName, TYPE calculatedResult, TYPE expectedResult, bool verbose)
{
  int result;
  if( calculatedResult == expectedResult )
    {
      if(verbose) std::cout << testName << " [ successful ]" << std::endl;
      result = 0;
    }
  else
    {
      if(verbose) std::cout << testName << " [ unsuccessful ]" << std::endl;
      result = 1;
    }
  return result;
}

