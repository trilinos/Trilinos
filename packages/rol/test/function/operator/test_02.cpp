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

/*! \file  test_02.cpp
    \brief Test of StdTridiagonalOperator 
*/

#include "ROL_StdTridiagonalOperator.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

   

  using SV  = ROL::StdVector<RealT>;
  using MAT = ROL::StdLinearOperator<RealT>;
  using TRI = ROL::StdTridiagonalOperator<RealT>; 

  using vector = std::vector<RealT>;
  

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  int errorFlag  = 0;

  RealT tol = ROL::ROL_EPSILON<RealT>();
 
  // *** Test body.

  try {
   
    int dim = 3;    

    ROL::Ptr<vector> m_ptr = ROL::makePtr<vector>(
      std::initializer_list<RealT>{ 3.0,  1.0, 0.0, -2.0, 6.0, 2.0, 0.0, -1.0, 3.0 } );
    ROL::Ptr<vector> a_ptr = ROL::makePtr<vector>(
      std::initializer_list<RealT>{ 3.0,  6.0, 3.0} );
    ROL::Ptr<vector> b_ptr = ROL::makePtr<vector>(
      std::initializer_list<RealT>{ -2.0, -1.0 } );
    ROL::Ptr<vector> c_ptr = ROL::makePtr<vector>(
      std::initializer_list<RealT>{  1.0,  2.0 } );

    MAT M(m_ptr);
    TRI T(a_ptr,b_ptr,c_ptr);   

    SV xm( ROL::makePtr<vector>( std::initializer_list<RealT>{1.0, 2.0, -1.0} ) );
    SV ym( ROL::makePtr<vector>( dim ) );
    SV zm( ROL::makePtr<vector>( dim ) );
    
    SV xt( ROL::makePtr<vector>( dim ) );
    SV yt( ROL::makePtr<vector>( dim ) );
    SV zt( ROL::makePtr<vector>( dim ) );
  
    SV error( ROL::makePtr<vector>(dim) );
    RealT nerr = 0;
 
    xt.set(xm);
 
    M.apply(ym,xm,tol);
    M.applyInverse(zm,ym,tol);

    *outStream << "\nUsing StdLinearOperator - A is full matrix representation" << std::endl;
    *outStream << "x = "; xm.print(*outStream);
    *outStream << "y = Ax = "; ym.print(*outStream);
    *outStream << "z = inv(A)y = "; zm.print(*outStream);
     
    *outStream << "\nUsing StdTridiagonalOperator - T is tridiagonal representation" << std::endl;
    T.apply(yt,xt,tol);
    *outStream << "y = Tx = "; yt.print(*outStream);
    error.set(yt);
    error.axpy(-1.0,ym);
    nerr = error.norm();  
    errorFlag += static_cast<int>(nerr>tol);
    *outStream << "apply() error = " << nerr <<std::endl;

    T.applyInverse(zt,yt,tol);
    *outStream << "z = inv(T)y = "; yt.print(*outStream);
    error.set(zt);
    error.axpy(-1.0,zm);
    nerr = error.norm();
    errorFlag += static_cast<int>(nerr>tol);
    *outStream << "applyInverse() error = " << nerr <<std::endl;
         
    M.applyAdjoint(ym,xm,tol);
    M.applyAdjointInverse(zm,ym,tol);
    *outStream << "\nUsing StdLinearOperator - A is full matrix representation" << std::endl;
    *outStream << "x = "; xm.print(*outStream);
    *outStream << "y = A'x = "; ym.print(*outStream);
    *outStream << "z = inv(A')y = "; zm.print(*outStream);
        
    *outStream << "\nUsing StdTridiagonalOperator - T is tridiagonal representation" << std::endl;
    T.applyAdjoint(yt,xt,tol);
    *outStream << "y = T'x = "; yt.print(*outStream);
    error.set(yt);
    error.axpy(-1.0,ym);
    nerr = error.norm();  
    errorFlag += static_cast<int>(nerr>tol);
    *outStream << "applyAdjoint() error = " << nerr <<std::endl;

    T.applyAdjointInverse(zt,yt,tol);
    *outStream << "z = inv(T')y = "; yt.print(*outStream);
    error.set(zt);
    error.axpy(-1.0,zm);
    nerr = error.norm();
    errorFlag += static_cast<int>(nerr>tol);
    *outStream << "applyAdjointInverse() error = " << nerr <<std::endl;
 

    

  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return 0;

}

