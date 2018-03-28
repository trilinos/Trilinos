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

#include <iostream>
#include <iomanip>
#include "ROL_Ptr.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "LowerBandedMatrix.hpp"

using RealT = double;
using size_type = std::vector<RealT>::size_type;

int main( int argc, char* argv[] ) {
  
  using std::vector;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);  

  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  auto print_vector = [outStream]( const vector<RealT>& x) { 
    *outStream << "\n";
    for( auto e : x ) *outStream << e << " ";
    *outStream << "\n";
  };

  try {

    size_type N = 8;
    vector<size_type> band_index{0,1};
    vector<vector<RealT>> A_band;

    vector<RealT> x{1, -1,  1, -1,   1, -1,  1, -1}; 
    vector<RealT> b{1, -1,  2, -7, -10,  5, -4,  2};
    vector<RealT> c{0,  0,  3, -6,  -9,  6, -3,  1};

    vector<RealT> y(N), z(N);

    vector<RealT> Ax(N);    
    vector<RealT> Atx(N);    

    vector<RealT> band0{ 1,  2,  4,  8, -8, -4, -2, -1 };
    vector<RealT> band1{ 1,  2,  1,  2,  1,  2,  1     };

    A_band.push_back(band0);
    A_band.push_back(band1);
 
    LowerBandedMatrix<RealT> A( band_index, A_band );

    A.apply(Ax,  x, 1.0, 0, N);
    A.applyTranspose(Atx, x, 1.0, 0, N);
    A.solve(y, Ax, 1.0, 0, N);
    A.solveTranspose(z, Atx, 1.0, 0, N);
    
//    RealT error2  = 0;
//    RealT errort2 = 0;

    for( size_type i=0; i<N; ++i ) {
      *outStream << std::setw(8) <<  x[i]    << " " 
                 << std::setw(8) <<  y[i]    << " "
                 << std::setw(8) <<  Ax[i]   << " "
                 << std::setw(8) <<  Atx[i]  << " "
                 << std::setw(8) <<  b[i]    << " "
                 << std::setw(8) <<  c[i]    << std::endl;
    }


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
