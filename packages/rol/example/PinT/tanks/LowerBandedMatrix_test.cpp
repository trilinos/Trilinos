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
#include <random>
#include <utility>
#include "ROL_Ptr.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "LowerBandedMatrix.hpp"

using RealT = double;
using size_type = std::vector<RealT>::size_type;

RealT& access( std::vector<RealT>& v, size_type i, size_type j, size_type cols ) { return v.at(i*cols+j); }

const RealT& access( const std::vector<RealT>& v, size_type i, size_type j, size_type cols ) { return v.at(i*cols+j); }

int main( int argc, char* argv[] ) {
  
  using namespace std;
  random_device r;
  default_random_engine gen(r());
  uniform_real_distribution<RealT> dist(1.0, 2.0);

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);  

  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
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

    size_type rows = 3;
    size_type cols = 3;
    size_type n    = rows*cols;
    size_type N    = n;
 
    auto zero = []( vector<RealT>& v ) { for( auto &e: v ) e = RealT{0}; };

    auto dot = [N]( const vector<RealT>& x1, const vector<RealT>& x2 ) {
                    RealT result{0}; 
                    for( size_type i=0; i<N; ++i )  
                      result += x1.at(i)*x2.at(i);
                    return result; };
    auto l1_diff = [N] ( const vector<RealT>& x1, const vector<RealT>& x2 ) {
                         RealT result{0}; 
                         for( size_type i=0; i<N; ++i )  
                           result += abs(x1.at(i)-x2.at(i));
                         return result; };

    auto alpha = dist(gen);
    auto beta  = dist(gen);
  
    RealT error_tol = sqrt(ROL::ROL_EPSILON<RealT>());

    LowerBandedMatrix<RealT> A( rows, cols, alpha, beta );
   
    vector<RealT> u(N,0.0), v(N,0.0), x(N,0.0), y(N,0.0), 
                  Ax(N,0.0), Aty(N,0.0), Aix(N,0.0), Aity(N,0.0);
    
    for( auto& e : x ) e = dist(gen);
    for( auto& e : y ) e = dist(gen);

    //--------------------------------------------------------------------------
    // Test of LowerBandedMatrix

    A.apply( Ax, x,  1.0, 0, 0 );
    A.applyTranspose( Aty, y, 1.0, 0, 0 );
    A.solve( u, Ax, 1.0, 0, 0 );
    A.solveTranspose( v, Aty, 1.0, 0, 0);

    A.solve( Aix, x, 1.0, 0, 0 );
    A.solveTranspose( Aity, y, 1.0, 0, 0 );

    RealT xnorm  = sqrt(dot(x,x));
    RealT ynorm  = sqrt(dot(y,y));

    RealT err_sol    = l1_diff(x,u);
    RealT err_sol_tr = l1_diff(y,v);

    RealT err_tr   = abs( dot(y,Ax)-dot(x,Aty) )/sqrt(ynorm*xnorm);
    RealT err_itr  = abs( dot(y,Aix)-dot(x,Aity) )/sqrt(ynorm*xnorm);

    *outStream << "\n\nTest of LowerBandedMatrix with alpha = " << alpha 
               << ", beta = " << beta << endl;
    *outStream << string(75,'-') << endl;
    *outStream << setw(60) << left 
                           << "Transpose error       | y'Ax-x'A'y |/sqrt(x'x y'x)" 
                           << " = " << err_tr << endl;
    *outStream << setw(60) << left 
                           << "Transpose-solve error | y'inv(A)x-x'inv(A)'y |/sqrt(x'x y'x)" 
                           << " = " << err_itr << endl;

    *outStream << setw(60) << left << "l1-norm of solve error" 
               << " = " << err_sol << endl;
    *outStream << setw(60) << left << "l1-norm of transpose solve error" 
               << " = " << err_sol_tr << endl;
    
    errorFlag += ( err_sol     > error_tol );
    errorFlag += ( err_sol_tr  > error_tol );
    errorFlag += ( err_itr     > error_tol );
    errorFlag += ( err_tr      > error_tol );

    //-------------------------------------------------------------------------
    // Test of SplitterMatrix
    
    SplitterMatrix<RealT> S(rows, cols);

    zero(Ax); zero(Aty);
    
    S.apply( Ax, x,  1.0, 0, 0 );
    S.applyTranspose( Aty, y, 1.0, 0, 0 );

    err_tr = abs( dot(y,Ax)-dot(x,Aty) )/(ynorm*xnorm);

    *outStream << "\nTest of SplitterMatrix" << endl;
    *outStream << string(75,'-') << endl;
    *outStream << setw(50) << left 
                           << "Transpose error | y'Sx-x'S'y |/sqrt(x'x y'x)" 
                           << " = " << err_tr << endl;
   
    //-------------------------------------------------------------------------
    // Test of TankLevelMatrix
    size_type Nhalf = N/2;
    vector<RealT> p(N,1.0); 
    p.at(Nhalf) = 0.0;

    TankLevelMatrix<RealT> T( rows, cols, alpha, p );
    zero(Ax); zero(Aty); zero(u); zero(v);

    T.apply( Ax, x,  1.0, 0, 0 );
    T.applyTranspose( Aty, y, 1.0, 0, 0 );
    T.solve( u, Ax, 1.0, 0, 0 );
    T.solveTranspose( v, Aty, 1.0, 0, 0);

    err_sol    = l1_diff(x,u);
    err_sol_tr = l1_diff(y,v);

    err_tr   = abs( dot(y,Ax)-dot(x,Aty) )/sqrt(ynorm*xnorm);

     *outStream << "\nTest of TankLevelMatrix with alpha = " << alpha 
                << " and pass-through in element " << Nhalf << endl;
     *outStream << string(75,'-') << endl;
     *outStream << setw(60) << left 
                            << "Transpose error       | y'Tx-x'T'y |/sqrt(x'x y'x)" 
                            << " = " << err_tr << endl;
    *outStream << setw(60) << left << "l1-norm of solve error" 
               << " = " << err_sol << endl;
    *outStream << setw(60) << left << "l1-norm of transpose solve error" 
               << " = " << err_sol_tr << endl;
    
    errorFlag += ( err_sol    > error_tol );
    errorFlag += ( err_sol_tr > error_tol );
    errorFlag += ( err_tr     > error_tol ); 
 

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
