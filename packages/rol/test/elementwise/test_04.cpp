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

#include "ROL_StdVector.hpp"
#include <cstdlib>
#include <ctime>
#include <chrono>


/*! \file  test_04.cpp
    \brief Compare performance of elementwise functions with built-in
           ROL::Vector methods
*/

template<class Real> 
class Axpy : public ROL::Elementwise::BinaryFunction<Real> {
public:
  Axpy(const Real &scalar) : scalar_(scalar) { }
  Real apply( const Real &x, const Real &y ) const {
    return scalar_*x+y;  
  }
private:
  Real scalar_;
};


template<class Real> 
Real axpy_func( const Real &x, const Real &y ) {
  return 2.0*x+y;
}


typedef double RealT;

int main(int argc, char *argv[]) {

  using namespace Teuchos;

  typedef std::vector<RealT>     vec;
  typedef ROL::StdVector<RealT>  SV;

  // Number of elements in vector
  size_t dim = argc>1 ? atoi(argv[1]) : 1000; 

  size_t nTrials = 1000;

  RCP<vec> x_rcp = rcp( new vec(dim,1.0) );
  RCP<vec> y_rcp = rcp( new vec(dim,2.0) );

  SV x(x_rcp);
  SV y(y_rcp);

  double duration = 0;

  for(size_t i=0;i<nTrials;++i) { 
    std::clock_t startTime = std::clock();
    x.axpy(2.0,y);
    std::clock_t endTime = std::clock();
    duration += (static_cast<double>(endTime-startTime))/CLOCKS_PER_SEC;
  }

  double average_duration = duration/nTrials;
 
  std::cout << "Average time of Vector::axpy() for n = " << 
        dim <<  ": " << average_duration << " seconds" << std::endl;

  duration = 0;

  Axpy<RealT> axpy(2.0);

  for(size_t i=0;i<nTrials;++i) { 
    std::clock_t startTime = std::clock();
    x.applyBinary(axpy,y);
    std::clock_t endTime = std::clock();
    duration += (static_cast<double>(endTime-startTime))/CLOCKS_PER_SEC;
  }

  average_duration = duration/nTrials;
 
  std::cout << "Average time of Vector::applyBinary with inheritance for n = " << 
        dim <<  ": " << average_duration << " seconds" << std::endl;


  duration = 0;

  for(size_t i=0;i<nTrials;++i) { 
    std::clock_t startTime = std::clock();
    ROL::Elementwise::applyBinaryInPlace(x,y,axpy_func<RealT>);
    std::clock_t endTime = std::clock();
    duration += (static_cast<double>(endTime-startTime))/CLOCKS_PER_SEC;
  }

  average_duration = duration/nTrials;
 
  std::cout << "Average time of Vector::applyBinary with function wrapper for n = " << 
        dim <<  ": " << average_duration << " seconds" << std::endl;




}
