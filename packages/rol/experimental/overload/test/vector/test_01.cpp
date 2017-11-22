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


/*! \file  test_01.cpp
    \brief Test std::vector interface using function overloading
*/
#include "XROL_StdVector.hpp"
#include "XROL.hpp"


int main( int argc, char *argv[] ) {

  using namespace std;
  using V        = vector<double>;
  using RealT    = XROL::magnitude_t<V>;
  using IndexT   = XROL::index_t<V>;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  default_random_engine generator;
  uniform_real_distribution<> dist(-1.0,1.0);
  
  ostream *os;

  auto errtol   = ROL::ROL_THRESHOLD<RealT>();
  IndexT dim    = 100;
  int errorFlag = 0;


  if( argc > 1 ) {
    os = &cout;
  }
  else {
    Teuchos::oblackholestream bhs;
    os = &bhs;
  }


  try {
   
    auto error_check = [&os,errtol,&errorFlag] ( auto val ) {
      if( abs(val) > errtol ) {
        *os << "---> POSSIBLE ERROR ABOVE!\n"; 
        ++errorFlag;
      }
      return 0;
    }; 
  
    // Create and randomize three vectors
    V x(dim);
    V y(dim);
    V w(dim);
  
    *os << "\n\nGenerating three vectors of length " << dim 
        << " with randomly generated element values from udf([-1,1])" << endl;

    XROL::randomize(generator,dist,x);
    XROL::randomize(generator,dist,y);
    XROL::randomize(generator,dist,w);

    *os << "\n ||w|| = " << XROL::norm(w) << endl;
    *os << "\n ||x|| = " << XROL::norm(x) << endl;
    *os << "\n ||y|| = " << XROL::norm(y) << endl;

     
    // Standard tests
    auto consistency = XROL::checkVector(x, y, w, *os);
    if( XROL::norm(consistency) > errtol ) {
      ++errorFlag;
    }
  
    auto z = *XROL::clone(x); XROL::basis(z,0);
    auto znorm = XROL::norm(z);
    *os << "Norm of ROL::Vector z (first basis vector):    " << znorm << endl;
    error_check( znorm-1 ); 
  
    XROL::basis(z,dim/2);
    znorm = XROL::norm(z);
    *os << "Norm of ROL::Vector z ('middle' basis vector): " << znorm << endl;
    error_check( znorm-1 ); 
  
    XROL::basis(z,dim-1);
    znorm = XROL::norm(z);
    *os << "Norm of ROL::Vector z (last basis vector):     " << znorm << endl;
    error_check( znorm-1 ); 
  
    // Repeat standard tests with a zero vector
    XROL::fill(x,0);
    consistency = XROL::checkVector(x,x,x,*os);
    if( XROL::norm(consistency) > 0 ) 
      ++errorFlag;
  }
  catch( logic_error err ) {
    *os << err.what() << endl;
    errorFlag = -1000;
  }; // end try

  cout << "End Result: TEST " << ( errorFlag ? "FAILED" : "PASSED" ) << endl;
  

  return 0;
}
