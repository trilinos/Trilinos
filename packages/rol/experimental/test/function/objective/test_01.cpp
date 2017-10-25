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
    \brief Test XROL::Objective interface - Compare values to those 
           using the ROL implementation of vectors
*/

#include "XROL.hpp"

#include "ZakharovObjective.hpp"

#include <iostream>



int main( int argc, char *argv[] ) {

  using namespace std;
  using V        = vector<double>;
  using IndexT   = XROL::index_t<V>;

  using XROL::Elementwise::evaluate_all;


  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  default_random_engine generator;
  uniform_real_distribution<> dist(-1.0,1.0);
  
  ostream *os;

  IndexT dim    = 10;
  int errorFlag = 0;

  if( argc > 1 ) {
    os = &cout;
  }
  else {
    Teuchos::oblackholestream bhs;
    os = &bhs;
  }

  try {
  
    Teuchos::ParameterList parlist;

    // Need to capture references to generator and distribution
    auto rnd = [&generator,&dist](auto &arg){
      XROL::randomize(generator,dist,arg); 
    };

    // Create vectors
    V x(dim);  V g(dim);  V v(dim);  V w(dim);
    V hv(dim);  
    V hw(dim);

    evaluate_all(rnd,x,g,v,w);

    auto k = make_unique<V>(dim);

    for( XROL::index_t<V> i=0; i<k->size(); ++i ) {
      (*k)[i] = 1.0*(1+i);
    }
   
    Zakharov::Objective<V> obj(move(k));
     
    XROL::checkGradient(obj,x,g,v,*os,parlist);
    XROL::checkHessVec(obj,x,g,v,*os,parlist);
    XROL::checkHessSym(obj,x,g,w,v,*os,parlist);

  }
  catch( logic_error err ) {
    *os << err.what() << endl;
    errorFlag = -1000;
  }; // end try

  cout << "End Result: TEST " << ( errorFlag ? "FAILED" : "PASSED" ) << endl;
  

  return 0;
}
