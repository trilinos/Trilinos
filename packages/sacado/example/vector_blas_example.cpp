// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// vector_blas_example
//
//  usage: 
//     vector_blas_example
//
//  output:  
//     prints the results of differentiating a BLAS routine with forward
//     mode AD using the Sacado::Fad::DVFad class (uses dynamic memory
//     allocation for number of derivative components stored in a contiguous
//     array).

#include <iostream>
#include <iomanip>

#include "Sacado_No_Kokkos.hpp"
#include "Teuchos_BLAS.hpp"
#include "Sacado_Fad_BLAS.hpp"

typedef Sacado::Fad::DVFad<double> FadType;

int main(int argc, char **argv)
{
  const unsigned int n = 5;
  Sacado::Fad::Vector<unsigned int, FadType> A(n*n,0),B(n,n), C(n,n);
  for (unsigned int i=0; i<n; i++) {
    for (unsigned int j=0; j<n; j++)
      A[i+j*n] = FadType(Teuchos::ScalarTraits<double>::random());
    B[i] = FadType(n, Teuchos::ScalarTraits<double>::random());
    for (unsigned int j=0; j<n; j++)
      B[i].fastAccessDx(j) = Teuchos::ScalarTraits<double>::random();
    C[i] = 0.0;
  }

  double *a = A.vals();
  double *b = B.vals();
  double *bdx = B.dx();
  std::vector<double> c(n), cdx(n*n);

  Teuchos::BLAS<int,double> blas;
  blas.GEMV(Teuchos::NO_TRANS, n, n, 1.0, &a[0], n, &b[0], 1, 0.0, &c[0], 1);
  blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, n, n, n, 1.0, &a[0], n, &bdx[0], n, 0.0, &cdx[0], n);

  // Teuchos::BLAS<int,FadType> blas_fad;
  // blas_fad.GEMV(Teuchos::NO_TRANS, n, n, 1.0, &A[0], n, &B[0], 1, 0.0, &C[0], 1);

  Teuchos::BLAS<int,FadType> sacado_fad_blas(false);
  sacado_fad_blas.GEMV(Teuchos::NO_TRANS, n, n, 1.0, &A[0], n, &B[0], 1, 0.0, &C[0], 1);

  // Print the results
  int p = 4;
  int w = p+7;
  std::cout.setf(std::ios::scientific);
  std::cout.precision(p);

  std::cout << "BLAS GEMV calculation:" << std::endl;
  std::cout << "a = " << std::endl;
  for (unsigned int i=0; i<n; i++) {
    for (unsigned int j=0; j<n; j++)
      std::cout << " " << std::setw(w) << a[i+j*n];
    std::cout << std::endl;
  }
  std::cout << "b = " << std::endl;
  for (unsigned int i=0; i<n; i++) {
    std::cout << " " << std::setw(w) << b[i];
  }
  std::cout << std::endl;
  std::cout << "bdot = " << std::endl;
  for (unsigned int i=0; i<n; i++) {
    for (unsigned int j=0; j<n; j++)
      std::cout << " " << std::setw(w) << bdx[i+j*n];
    std::cout << std::endl;
  }
  std::cout << "c = " << std::endl;
  for (unsigned int i=0; i<n; i++) {
    std::cout << " " << std::setw(w) << c[i];
  }
  std::cout << std::endl;
  std::cout << "cdot = " << std::endl;
  for (unsigned int i=0; i<n; i++) {
    for (unsigned int j=0; j<n; j++)
      std::cout << " " << std::setw(w) << cdx[i+j*n];
    std::cout << std::endl;
  }
  std::cout << std::endl << std::endl;

  std::cout << "FAD BLAS GEMV calculation:" << std::endl;
  std::cout << "A.val() (should = a) = " << std::endl;
  for (unsigned int i=0; i<n; i++) {
    for (unsigned int j=0; j<n; j++)
      std::cout << " " << std::setw(w) << A[i+j*n].val();
    std::cout << std::endl;
  }
  std::cout << "B.val() (should = b) = " << std::endl;
  for (unsigned int i=0; i<n; i++) {
    std::cout << " " << std::setw(w) << B[i].val();
  }
  std::cout << std::endl;
  std::cout << "B.dx() (should = bdot) = " << std::endl;
  double *Bdx = B.dx();
  for (unsigned int i=0; i<n; i++) {
    for (unsigned int j=0; j<n; j++)
      std::cout << " " << std::setw(w) << Bdx[i+j*n];
    std::cout << std::endl;
  }
  std::cout << "C.val() (should = c) = " << std::endl;
  for (unsigned int i=0; i<n; i++) {
    std::cout << " " << std::setw(w) << C[i].val();
  }
  std::cout << std::endl;
  std::cout << "C.dx() (should = cdot) = " << std::endl;
  double *Cdx = C.dx();
  for (unsigned int i=0; i<n; i++) {
    for (unsigned int j=0; j<n; j++)
      std::cout << " " << std::setw(w) << Cdx[i+j*n];
    std::cout << std::endl;
  }

  double tol = 1.0e-14;
  bool failed = false;
  for (unsigned int i=0; i<n; i++) {
    if (std::fabs(C[i].val() - c[i]) > tol)
      failed = true;
    for (unsigned int j=0; j<n; j++) {
      if (std::fabs(C[i].dx(j) - cdx[i+j*n]) > tol) 
	failed = true;
    }
  }
  if (!failed) {
    std::cout << "\nExample passed!" << std::endl;
    return 0;
  }
  else {
    std::cout <<"\nSomething is wrong, example failed!" << std::endl;
    return 1;
  }
}
