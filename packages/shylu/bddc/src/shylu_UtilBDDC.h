
//@HEADER
// ************************************************************************
// 
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef UTILBDDC_H
#define UTILBDDC_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <algorithm>
#include <vector>

#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"

// Author: Clark R. Dohrmann
namespace bddc {
  
  template <class SX, class SM> 
  class UtilBDDC
{
public:
  UtilBDDC()
  {
  }

  static void calculateSchurComplement(std::vector<SX> & A, 
				       std::vector<int> & i1, 
				       std::vector<int> & i2, 
				       std::vector<SX> & Sc)
  {
    std::vector<SX> A21;
    calculateSchurComplement(A, i1, i2, Sc, A21);
  }

  static void calculateSchurComplement(std::vector<SX> & A, 
				       std::vector<int> & i1, 
				       std::vector<int> & i2, 
				       std::vector<SX> & Sc,
				       std::vector<SX> & A21)
  {
    // calculates Schur complement Sc = A11 - A12*inv(A22)*A21
    std::vector<SX> A12, A22;
    loadDenseMatrices(A, i1, i2, Sc, A12, A21, A22);
    Teuchos::BLAS<int, SX>  BLAS;
    Teuchos::LAPACK<int, SX> LAPACK;
    int numRows1 = i1.size();
    int numRows2 = i2.size();
    int INFO(0);
    if (numRows2 > 0) {
      // static condensation of "2" unknowns
      std::vector<int> IPIV(numRows2);
      LAPACK.GETRF(numRows2, numRows2, &A22[0], numRows2, &IPIV[0], &INFO);
      assert (INFO == 0);
      int NRHS = numRows1;
      LAPACK.GETRS('N', numRows2, NRHS, &A22[0], numRows2, &IPIV[0], 
		   &A21[0], numRows2, &INFO);
      assert (INFO == 0);
      SX ALPHA(-1), BETA(1);
      if (numRows1 > 0) {
	BLAS.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, numRows1, numRows1,
		  numRows2, ALPHA, &A12[0], numRows1, &A21[0], numRows2, 
		  BETA, &Sc[0], numRows1);
      }
    }
  }
  
  static void solveSymmetricGeneralizedEigenProblem(std::vector<SX> & A,
						    std::vector<SX> & B,
						    int N,
						    std::vector<SM> & W)
  {
    int ITYPE = 1; // A*x = lambda*B*x;
    char JOBZ('V'), UPLO('U');
    std::vector<SM> RWORK(std::max(1, 3*N-2));
    Teuchos::LAPACK<int, SX> LAPACK;
    int LWORK(-1), INFO(0);
    SX dumWORK[1];
    LAPACK.HEGV(ITYPE, JOBZ, UPLO, N, &A[0], N, &B[0], N, &W[0], 
		dumWORK, LWORK, &RWORK[0], &INFO);
    assert (INFO == 0);
    LWORK = int(real(dumWORK[0])+0.1);
    std::vector<SX> WORK(LWORK);
    LAPACK.HEGV(ITYPE, JOBZ, UPLO, N, &A[0], N, &B[0], N, &W[0], 
		&WORK[0], LWORK, &RWORK[0], &INFO);
    assert (INFO == 0);
  }

  static void solveRealEigenValueProblem(std::vector<SX> & A,
					 int N,
					 std::vector<SM> & W,
					 const char JOBZ)
  {
    if (N == 0) return;
    assert (int(A.size()) == N*N);
    char UPLO('U');
    W.resize(N);
    std::vector<SM> RWORK(std::max(1, 3*N));
    int LWORK = std::max(1, 3*N);
    std::vector<SX> WORK(LWORK);
    Teuchos::LAPACK<int, SX> LAPACK;
    int INFO(0);
    LAPACK.HEEV(JOBZ, UPLO, N, &A[0], N, &W[0], &WORK[0], LWORK, 
		&RWORK[0], &INFO);
    assert (INFO == 0);
  }

  static void calculateRealEigenValues(std::vector<SX> & A, 
				       int N,
				       std::vector<SM> & W)
  {
    char JOBZ('N');
    solveRealEigenValueProblem(A, N, W, JOBZ);
  }

  static void calculateRealEigenValuesAndVectors(std::vector<SX> & A, 
						 int N,
						 std::vector<SM> & W)
  {
    char JOBZ('V');
    solveRealEigenValueProblem(A, N, W, JOBZ);
  }

  static void printIndices(const int numIndices, 
			   const int* indices,
			   const char* fileName)
  {
    std::ofstream fout;
    fout.open(fileName);
    for (int j=0; j<numIndices; j++) {
      fout << indices[j]+1 << std::endl;;
    }
    fout.close();
  }

  static void printDenseMatrix(const int numRows,
			       const int numCols,
			       const SX* A,
			       const char* fileName)
  {
    std::ofstream fout;
    fout.open(fileName);
    for (int j=0; j<numCols; j++) {
      for (int i=0; i<numRows; i++) {
	fout << i+1 << " ";
	fout << j+1 << " ";
	fout << std::setw(22) << std::setprecision(15);
	SX value = A[i+numRows*j];
	fout << real(value);
	if (isComplex(value) == true) {
	  fout << " " << imag(value);
	}
	fout << std::endl;
      }
    }
    fout.close();
  }

  static void printSparseMatrix(const int numRows,
				const int* rowBegin,
				const int* columns,
				const SX* values,
				const char* fileName)
  {
    std::ofstream fout;
    fout.open(fileName);
    for (int i=0; i<numRows; i++) {
      for (int j=rowBegin[i]; j<rowBegin[i+1]; j++) {
	fout << i+1 << "  " << columns[j]+1 << " ";
	fout << std::setw(22) << std::setprecision(15);
	SX value = values[j];
	fout << real(value);
	if (isComplex(value) == true) {
	  fout << " " << imag(value);
	}
	fout << std::endl;
      }
    }
    fout.close();
  }

  static void printCoords(std::vector<SM> & xCoords,
			  std::vector<SM> & yCoords,
			  std::vector<SM> & zCoords,
			  const char* fileName)
  {
    std::ofstream fout;
    fout.open(fileName);
    int numRows = xCoords.size();
    for (int i=0; i<numRows; i++) {
      fout << std::setw(22) << std::setprecision(15);
      fout << xCoords[i] << " "
	   << yCoords[i] << " "
	   << zCoords[i] << std::endl;
    }
    fout.close();
  }

  static void printLocalDofs(std::vector<int> & localDofs,
			     const char* fileName)
  {
    std::ofstream fout;
    fout.open(fileName);
    int numRows = localDofs.size();
    for (int i=0; i<numRows; i++) {
      fout << localDofs[i] << std::endl;
    }
    fout.close();
  }

  static void calculateSchurComplement(int numRows1, 
				       int numRows2,
				       int ScOption,
				       std::vector<SX> & A,
				       std::vector<SX> & Sc)
  {
    std::vector<SX> A21;
    calculateSchurComplementHere(numRows1, numRows2, ScOption, A, Sc, A21);
  }

  static void calculateSchurComplement(int numRows1, 
				       int numRows2,
				       int ScOption,
				       std::vector<SX> & A,
				       std::vector<SX> & Sc,
				       std::vector<SX> & ExtensionMatrix)
  {
    calculateSchurComplementHere(numRows1, numRows2, ScOption, A, Sc, 
				 ExtensionMatrix);
  }

  static void calculateSchurComplementHere(int numRows1, 
					   int numRows2,
					   int ScOption,
					   std::vector<SX> & A,
					   std::vector<SX> & Sc,
					   std::vector<SX> & A22invA21)
  {
    assert ((ScOption == 1) || (ScOption == 2));
    int N = numRows1 + numRows2;
    if (N == 0) return;
    assert (int(A.size()) == N*N);
    std::vector<int> i1(numRows1), i2(numRows2);
    for (int i=0; i<numRows1; i++) i1[i] = i;
    for (int i=0; i<numRows2; i++) i2[i] = i + numRows1;
    if (ScOption == 1) {
      calculateSchurComplement(A, i1, i2, Sc, A22invA21);
    }
    else {
      calculateSchurComplement(A, i2, i1, Sc, A22invA21);
    }
  }

  static void loadDenseMatrices(std::vector<SX> & A, 
				std::vector<int> & i1,
				std::vector<int> & i2, 
				std::vector<SX> & A11, 
				std::vector<SX> & A12, 
				std::vector<SX> & A21, 
				std::vector<SX> & A22)
  {
    int numRows1 = i1.size();
    int numRows2 = i2.size();
    A11.resize(numRows1*numRows1);
    A12.resize(numRows1*numRows2);
    A21.resize(numRows2*numRows1);
    A22.resize(numRows2*numRows2);
    int LDA = numRows1 + numRows2;
    for (int i=0; i<numRows1; i++) {
      int row = i1[i];
      for (int j=0; j<numRows1; j++) {
	int col = i1[j];
	A11[i+j*numRows1] = A[row+col*LDA];
      }
      for (int j=0; j<numRows2; j++) {
	int col = i2[j];
	A12[i+j*numRows1] = A[row+col*LDA];
      }
    }
    for (int i=0; i<numRows2; i++) {
      int row = i2[i];
      for (int j=0; j<numRows1; j++) {
	int col = i1[j];
	A21[i+j*numRows2] = A[row+col*LDA];
      }
      for (int j=0; j<numRows2; j++) {
	int col = i2[j];
	A22[i+j*numRows2] = A[row+col*LDA];
      }
    }
  }

  static float imag(float a) {return 0;};
  static double imag(double a) {return 0;};
  static float imag(std::complex<float> a) {return std::imag(a);};
  static double imag(std::complex<double> a) {return std::imag(a);};

  static float real(float a) {return a;};
  static double real(double a) {return a;};
  static float real(std::complex<float> a) {return std::real(a);};
  static double real(std::complex<double> a) {return std::real(a);};

  static bool isComplex(float a) {return false;};
  static bool isComplex(double a) {return false;};
  static bool isComplex(std::complex<float> a) {return true;};
  static bool isComplex(std::complex<double> a) {return true;};

  static float conj(float a) {return a;};
  static double conj(double a) {return a;};
  static std::complex<float> conj(std::complex<float> a) 
  {
    return std::conj(a);
  };
  static std::complex<double> conj(std::complex<double> a) 
  {
    return std::conj(a);
  };

private:

};

} // namespace bddc

#endif // UTILBDDC_H
  
