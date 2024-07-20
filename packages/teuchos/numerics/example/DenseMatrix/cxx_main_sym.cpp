// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_SerialSymDenseMatrix.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialSpdDenseSolver.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Version.hpp"

int main(int argc, char* argv[])
{
  std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  // Creating a double-precision matrix can be done in several ways:
  // Create an empty matrix with no dimension
  Teuchos::SerialSymDenseMatrix<int,double> Empty_Matrix;
  // Create an empty 4x4 matrix
  Teuchos::SerialSymDenseMatrix<int,double> My_Matrix( 4 );
  // Basic copy of My_Matrix
  Teuchos::SerialSymDenseMatrix<int,double> My_Copy1( My_Matrix ),
    // (Deep) Copy of principle 3x3 submatrix of My_Matrix
    My_Copy2( Teuchos::Copy, My_Matrix, 3 ),
    // (Shallow) Copy of 3x3 submatrix of My_Matrix
    My_Copy3( Teuchos::View, My_Matrix, 3, 1 );

  // The matrix dimensions and strided storage information can be obtained:
  int rows, cols, stride;
  rows = My_Copy3.numRows();  // number of rows
  cols = My_Copy3.numCols();  // number of columns
  stride = My_Copy3.stride(); // storage stride
  TEUCHOS_ASSERT_EQUALITY(rows, 3);
  TEUCHOS_ASSERT_EQUALITY(cols, 3);
  TEUCHOS_ASSERT_EQUALITY(stride, 4);

  // Matrices can change dimension:
  Empty_Matrix.shape( 3 );      // size non-dimensional matrices
  My_Matrix.reshape( 3 );       // resize matrices and save values

  // Filling matrices with numbers can be done in several ways:
  My_Matrix.random();             // random numbers
  My_Copy1.putScalar( 1.0 );      // every entry is 1.0
  My_Copy1 = 1.0;                 // every entry is 1.0 (still)
  My_Copy2(1,1) = 10.0;           // individual element access
  Empty_Matrix = My_Matrix;       // copy My_Matrix to Empty_Matrix

  // Basic matrix arithmetic can be performed:
  Teuchos::SerialDenseMatrix<int,double> My_Prod( 4, 3 ), My_GenMatrix( 4, 3 );
  My_GenMatrix = 1.0;
  // Matrix multiplication ( My_Prod = 1.0*My_GenMatrix*My_Matrix )
  My_Prod.multiply( Teuchos::RIGHT_SIDE, 1.0, My_Matrix, My_GenMatrix, 0.0 );
  My_Copy2 += My_Matrix;   // Matrix addition
  My_Copy2 *= 0.5;         // Matrix scaling

  // Matrices can be compared:
  // Check if the matrices are equal in dimension and values
  if (Empty_Matrix == My_Matrix) {
    std::cout<< "The matrices are the same!" <<std::endl;
  }
  // Check if the matrices are different in dimension or values
  if (My_Copy2 != My_Matrix) {
    std::cout<< "The matrices are different!" <<std::endl;
  }

  // The norm of a matrix can be computed:
  double norm_one, norm_inf, norm_fro;
  norm_one = My_Matrix.normOne();        // one norm
  norm_inf = My_Matrix.normInf();        // infinity norm
  norm_fro = My_Matrix.normFrobenius();  // frobenius norm

  std::cout << std::endl << "|| My_Matrix ||_1 = " << norm_one << std::endl;
  std::cout << "|| My_Matrix ||_Inf = " << norm_inf << std::endl;
  std::cout << "|| My_Matrix ||_F = " << norm_fro << std::endl << std::endl;

  // A matrix can be factored and solved using Teuchos::SerialDenseSolver.
  Teuchos::SerialSpdDenseSolver<int,double> My_Solver;
  Teuchos::SerialSymDenseMatrix<int,double> My_Matrix2( 3 );
  My_Matrix2.random();
  Teuchos::SerialDenseMatrix<int,double> X(3,1), B(3,1);
  X = 1.0;
  B.multiply( Teuchos::LEFT_SIDE, 1.0, My_Matrix2, X, 0.0 );
  X = 0.0;  // Make sure the computed answer is correct.

  int info = 0;
  My_Solver.setMatrix( Teuchos::rcp( &My_Matrix2, false ) );
  My_Solver.setVectors( Teuchos::rcp( &X, false ), Teuchos::rcp( &B, false ) );
  info = My_Solver.factor();
  if (info != 0)
    std::cout << "Teuchos::SerialSpdDenseSolver::factor() returned : " << info << std::endl;
  info = My_Solver.solve();
  if (info != 0)
    std::cout << "Teuchos::SerialSpdDenseSolver::solve() returned : " << info << std::endl;

  // A matrix triple-product can be computed:  C = alpha*W'*A*W
  double alpha=0.5;
  Teuchos::SerialDenseMatrix<int,double> W(3,2);
  Teuchos::SerialSymDenseMatrix<int,double> A1(2), A2(3);
  A1(0,0) = 1.0, A1(1,1) = 2.0;
  A2(0,0) = 1.0, A2(1,1) = 2.0, A2(2,2) = 3.00;
  W = 1.0;

  Teuchos::SerialSymDenseMatrix<int,double> C1(3), C2(2);

  Teuchos::symMatTripleProduct<int,double>( Teuchos::NO_TRANS, alpha, A1, W, C1);
  Teuchos::symMatTripleProduct<int,double>( Teuchos::TRANS, alpha, A2, W, C2 );

  // A matrix can be sent to the output stream:
  std::cout<< printMat(My_Matrix) << std::endl;
  std::cout<< printMat(X) << std::endl;

  return 0;
}
