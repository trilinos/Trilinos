// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_Version.hpp"

int main(int argc, char* argv[])
{
  std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  // Creating an instance of the LAPACK class for double-precision routines looks like:
  Teuchos::LAPACK<int, double> lapack;

  // This instance provides the access to all the LAPACK routines.
  Teuchos::SerialDenseMatrix<int, double> My_Matrix(4,4);
  Teuchos::SerialDenseVector<int, double> My_Vector(4);
  My_Matrix.random();
  My_Vector.random();

  // Perform an LU factorization of this matrix.
  int ipiv[4], info;
  char TRANS = 'N';
  lapack.GETRF( 4, 4, My_Matrix.values(), My_Matrix.stride(), ipiv, &info );

  // Solve the linear system.
  lapack.GETRS( TRANS, 4, 1, My_Matrix.values(), My_Matrix.stride(),
		ipiv, My_Vector.values(), My_Vector.stride(), &info );

  // Print out the solution.
  std::cout << printMat(My_Vector) << std::endl;

  return 0;
}
