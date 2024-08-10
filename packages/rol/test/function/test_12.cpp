// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_12.cpp
    \brief Validate that the Householder Reflector implmentation 
           works correctly.
*/


#include "ROL_HouseholderReflector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

template<class Real> 
void printVector( const ROL::Vector<Real> &x, std::ostream &outStream ) {

  ROL::Ptr<const std::vector<Real> > xp = 
    dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();

  outStream << "Standard Vector" << std::endl;
  for( size_t i=0; i<xp->size(); ++i ) {
    outStream << (*xp)[i] << std::endl;
  }
}



typedef double RealT;

int main(int argc, char *argv[]) {
  
  typedef ROL::Vector<RealT>    V;
  typedef ROL::StdVector<RealT> SV;

   

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag = 0;

  try {

    int dim = 10;

    RealT tol = std::sqrt(ROL::ROL_EPSILON<RealT>());  

    ROL::Ptr<V> v   = ROL::makePtr<SV>( ROL::makePtr<std::vector<RealT>>(dim) );
    ROL::Ptr<V> Hv  = v->clone();
    ROL::Ptr<V> HHv = v->clone();

    ROL::Ptr<V> e0 = v->basis(0);

    RandomizeVector(*v);

    ROL::HouseholderReflector<RealT> H(v,e0);

    // Reflect v about a vector to the x direction
    // Verify that the result is parallel to x

    printVector(*v,*outStream);

    H.apply(*Hv, *v, tol);
  
    printVector(*Hv,*outStream);

    H.apply(*HHv, *Hv, tol);
  
    printVector(*HHv,*outStream);

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


