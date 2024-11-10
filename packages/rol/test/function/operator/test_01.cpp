// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_01.cpp
    \brief Test of StdLinearOperator, its inverse and transpose

    \f$ A=\begin{pmatrix} 4 & 1 \\ 2 & 3 \end{pmatrix},\quad 
        A^{-1}=\frac{1}{10}\begin{pmatrix} 4 & -1 \\ -2 & 3 \end{pmatrix} \f$

    1) Compute \f$b\f$ in \f$Ax = b\f$, when \f$ x=\begin{pmatrix} 1 \\ -1 \end{pmatrix}\f$
 
    2) Solve for \f$x\f$ in the above when \f$b=\begin{pmatrix} 3 \\ -1 \end{pmatrix}\f$

    3) Compute \f$c\f$ in \f$A^\top y=c\f$ when \f$y=\begin{pmatrix} -2 \\ 1 \end{pmatrix}\f$ 

    4) Solve for \f$y\f$ in the above when \f$c=\begin{pmatrix} -6 \\ 1 \end{pmatrix}\f$

    Also ensure that the interface works with both ROL::Vector and std::vector arguments
*/

#include "ROL_StdLinearOperator.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

   

  typedef std::vector<RealT>            vector;

  typedef ROL::StdVector<RealT>         SV;

  typedef ROL::StdLinearOperator<RealT> StdLinearOperator;


  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  // Save the format state of the original std::cout.
  ROL::nullstream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  int errorFlag  = 0;

  // *** Test body.

  try {

    ROL::Ptr<vector> a_ptr  = ROL::makePtr<vector>(
      std::initializer_list<RealT>{4.0,2.0,1.0,3.0});
    ROL::Ptr<vector> ai_ptr = ROL::makePtr<vector>(
      std::initializer_list<RealT>{3.0/10.0, -2.0/10.0, -1.0/10.0, 4.0/10.0});

    ROL::Ptr<vector> x1_ptr  = ROL::makePtr<vector>(
      std::initializer_list<RealT>{1.0,-1.0});
    ROL::Ptr<vector> b1_ptr = ROL::makePtr<vector>(2);

    ROL::Ptr<vector> x2_ptr = ROL::makePtr<vector>(2);
    ROL::Ptr<vector> b2_ptr  = ROL::makePtr<vector>(
      std::initializer_list<RealT>{3.0,-1.0});
    
    ROL::Ptr<vector> y3_ptr = ROL::makePtr<vector>(
      std::initializer_list<RealT>{-2.0,1.0});
    ROL::Ptr<vector> c3_ptr = ROL::makePtr<vector>(2);

    ROL::Ptr<vector> y4_ptr = ROL::makePtr<vector>(2);
    ROL::Ptr<vector> c4_ptr = ROL::makePtr<vector>(
      std::initializer_list<RealT>{-6.0,1.0});

    StdLinearOperator A(a_ptr);
    StdLinearOperator Ai(ai_ptr);

    SV x1(x1_ptr); SV x2(x2_ptr); SV y3(y3_ptr); SV y4(y4_ptr);
    SV b1(b1_ptr); SV b2(b2_ptr); SV c3(c3_ptr); SV c4(c4_ptr);

    RealT tol = ROL::ROL_EPSILON<RealT>();

    // Test 1
    *outStream << "\nTest 1: Matrix multiplication" << std::endl;
    A.apply(b1,x1,tol);
    *outStream << "x = [" << (*x1_ptr)[0] << "," << (*x1_ptr)[1] << "]" << std::endl;
    *outStream << "b = [" << (*b1_ptr)[0] << "," << (*b1_ptr)[1] << "]" << std::endl;
    b1.axpy(-1.0,b2);

    RealT error1 = b1.norm();
    errorFlag += error1 > tol;
    *outStream << "Error = " << error1 << std::endl;

    // Test 2    
    *outStream << "\nTest 2: Linear solve" << std::endl;
    A.applyInverse(*x2_ptr,*b2_ptr,tol);
    *outStream << "x = [" << (*x2_ptr)[0] << "," << (*x2_ptr)[1] << "]" << std::endl;
    *outStream << "b = [" << (*b2_ptr)[0] << "," << (*b2_ptr)[1] << "]" << std::endl;
    x2.axpy(-1.0,x1);

    RealT error2 = x2.norm();
    errorFlag += error2 > tol;
    *outStream << "Error = " << error2 << std::endl;

    // Test 3
    *outStream << "\nTest 3: Transposed matrix multiplication" << std::endl;
    A.applyAdjoint(*c3_ptr,*y3_ptr,tol);
    *outStream << "y = [" << (*y3_ptr)[0] << "," << (*y3_ptr)[1] << "]" << std::endl;
    *outStream << "c = [" << (*c3_ptr)[0] << "," << (*c3_ptr)[1] << "]" << std::endl;
    c3.axpy(-1.0,c4);

    RealT error3 = c3.norm();
    errorFlag += error3 > tol;
    *outStream << "Error = " << error3 << std::endl;

    // Test 4
    *outStream << "\nTest 4: Linear solve with transpose" << std::endl;
    A.applyAdjointInverse(y4,c4,tol);
    *outStream << "y = [" << (*y4_ptr)[0] << "," << (*y4_ptr)[1] << "]" << std::endl;
    *outStream << "c = [" << (*c4_ptr)[0] << "," << (*c4_ptr)[1] << "]" << std::endl;
    y4.axpy(-1.0,y3);

    RealT error4 = y4.norm();
    errorFlag += error4 > tol;
    *outStream << "Error = " << error4 << std::endl;

    *outStream << "x1 = ";  x1.print(*outStream); 
    Ai.applyInverse(b1,x1,tol);
    *outStream << "b1 = ";  b1.print(*outStream);
    A.apply(b1,x1,tol);
    *outStream << "b1 = ";  b1.print(*outStream);
    A.applyInverse(x1,b1,tol);
    *outStream << "x1 = ";  x1.print(*outStream);
    Ai.apply(x1,b1,tol);
    *outStream << "x1 = ";  x1.print(*outStream);
    

  }
  catch (std::logic_error& err) {
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

