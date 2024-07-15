// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Objective.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_PartitionedVector.hpp"
#include "ROL_ProfiledVector.hpp"
#include "ROL_VectorClone.hpp"




using RealT = double;
using size_type = typename std::vector<RealT>::size_type;

template<>
ROL::VectorFunctionCalls<size_type> 
ROL::ProfiledVector<size_type,RealT>::functionCalls_ = ROL::VectorFunctionCalls<size_type>();

template<typename Real>
class TestSingle {
private:
  mutable ROL::VectorClone<Real> xclone_;
public:
 
  Real value( const ROL::Vector<Real>& x ) const {
    auto xc = xclone_(x); 
    xc->set(x);
    return xc->dot(x);
  }
};

template<typename Real> 
class TestMulti {
private:
  mutable ROL::VectorCloneMap<Real> clones_;
public:
  TestMulti() : clones_("x","y") {}

  Real value_x( const ROL::Vector<Real>& x ) const {
    auto xc = clones_(x,"x"); 
    xc->set(x);
    return xc->dot(x);
  }

  Real value_y( const ROL::Vector<Real>& y ) const {
    auto yc = clones_(y,"y"); 
    yc->set(y);
    return yc->dot(y);
  }

  Real value_z( const ROL::Vector<Real>& z ) const {
    auto zc = clones_(z,"z"); 
    zc->set(z);
    return zc->dot(z);
  }
};


int main( int argc, char* argv[] ) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;
  RealT errtol = ROL::ROL_THRESHOLD<RealT>();

  try {

    size_type N = 20;
    size_type M = 10;

    auto xp = ROL::makePtr<std::vector<RealT>>(N);
    auto x  = ROL::makePtr<ROL::StdVector<RealT>>(xp);
    auto yp = ROL::makePtr<std::vector<RealT>>(M);    // half the size
    auto y  = ROL::makePtr<ROL::StdVector<RealT>>(yp);
    auto z  = ROL::PartitionedVector<RealT>::create({y,y});

    ROL::ProfiledVector<size_type,RealT> xprofile(x);

    x->setScalar(1.0);
    y->setScalar(2.0);
 
    TestSingle<RealT> test_single_1;
    TestSingle<RealT> test_single_2;
    TestMulti<RealT>  test_multi;

    //------------------------------------------------------------------------
    // Test vector of same size and type 
    auto value1 = test_single_1.value(*x);
    RealT err1 = std::abs(N-value1);

    *outStream << "\nTesting single VectorClone of same size and type: ";
    if (err1<errtol) { *outStream << "Works!" << std::endl;   } 
    else {
      errorFlag += err1>errtol;
      *outStream << "Incorrect result!" << std::endl;
    }


    //------------------------------------------------------------------------
    // Test that exception is thrown if vector has wrong size
    bool passed_test_2 = false;
    *outStream << "Throw exception on mismatched dimension : ";

    try { test_single_1.value(*y); }
    catch( std::logic_error& size_mismatch ) { passed_test_2 = true;  }

    if( passed_test_2 ) { *outStream << "Works!" << std::endl; }
    else { 
      *outStream << "Failed to throw!" << std::endl;  
      errorFlag++;
    }

    //------------------------------------------------------------------------
    // Test that exception is thrown if vector has wrong type
    bool passed_test_3 = false;
    *outStream << "Throw exception on mismatched type : ";

    try { test_single_1.value(*y); }
    catch( std::logic_error& dim_mismatch ) { passed_test_3 = true;  }

    if( passed_test_3 ) { *outStream << "Works!" << std::endl; }
    else { 
      *outStream << "Failed to throw!" << std::endl;  
      errorFlag++;
    }

    //------------------------------------------------------------------------
    // Test that clone is only called once
    *outStream << "\n\nTesting with ProfiledVector. # calls to clone: ";
    for( int i=0; i<10; ++i ) {
      test_single_2.value(xprofile);
    }

    auto calls = getVectorFunctionCalls(xprofile);
    *outStream << calls.clone_ << std::endl;
    if( calls.clone_ > 1 ) { errorFlag++; }
      
    // Display number of function calls
    ROL::printVectorFunctionCalls(xprofile, *outStream);

    //------------------------------------------------------------------------
    // Test VectorCloneMap
    
    bool vcm_pass = true;

    *outStream << "Testing VectorCloneMap: ";

    auto x_value = test_multi.value_x(*x);
    auto y_value = test_multi.value_y(*y);
    auto z_value = test_multi.value_z(*z);
    
    auto errx = std::abs(x_value-20.0);
    auto erry = std::abs(y_value-40.0);
    auto errz = std::abs(z_value-80.0);

    if( errx>errtol ) { vcm_pass = false; errorFlag++; }
    if( erry>errtol ) { vcm_pass = false; errorFlag++; }
    if( errz>errtol ) { vcm_pass = false; errorFlag++; }
 
    if( vcm_pass ) { *outStream << "Works!" << std::endl; }
    else { 
      *outStream << "Error tolerance exceeded!" << std::endl; 
      *outStream << "x_value was " << x_value << ", should be 20." << std::endl;
      *outStream << "y_value was " << y_value << ", should be 40." << std::endl;
      *outStream << "z_value was " << z_value << ", should be 80." << std::endl;
    }

  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}

