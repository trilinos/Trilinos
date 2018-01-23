// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  test_12.cpp
    \brief Validate that the Householder Reflector implmentation 
           works correctly.
*/


#include "ROL_HouseholderReflector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_RandomVector.hpp"
#include "Teuchos_oblackholestream.hpp"
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
  Teuchos::oblackholestream bhs; // outputs nothing
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
  
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}


