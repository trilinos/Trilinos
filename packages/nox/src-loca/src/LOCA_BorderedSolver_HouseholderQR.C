// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "LOCA_BorderedSolver_HouseholderQR.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::BorderedSolver::HouseholderQR::HouseholderQR(
	 const Teuchos::RCP<LOCA::GlobalData>& global_data) : 
  globalData(global_data),
  dblas()
{
}

LOCA::BorderedSolver::HouseholderQR::~HouseholderQR()
{
}

void
LOCA::BorderedSolver::HouseholderQR::computeQR(
			    const NOX::Abstract::MultiVector::DenseMatrix& C,
			    const NOX::Abstract::MultiVector& B,
			    bool use_c_transpose,
			    NOX::Abstract::MultiVector::DenseMatrix& Y1,
			    NOX::Abstract::MultiVector& Y2,
			    NOX::Abstract::MultiVector::DenseMatrix& T,
			    NOX::Abstract::MultiVector::DenseMatrix& R)
{
  double beta;
  int m = B.numVectors();
  
  // Initialize
  Y1.putScalar(0.0);
  T.putScalar(0.0);
  Y2 = B;
  if (use_c_transpose) {
    for (int i=0; i<m; i++)
      for (int j=0; j<m; j++)
	R(i,j) = C(j,i);        // Copy transpose of C into R
  }
  else
    R.assign(C);

  // A temporary vector
  Teuchos::RCP<NOX::Abstract::MultiVector> v2 = Y2.clone(1);

  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> v1;
  Teuchos::RCP<NOX::Abstract::MultiVector> h2;
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> h1;
  Teuchos::RCP<NOX::Abstract::MultiVector> y2;
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> y1;
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> z;
  std::vector<int> h_idx;
  std::vector<int> y_idx;
  y_idx.reserve(m);

  for (int i=0; i<m; i++) {

    // Create view of column i of Y1 starting at row i
    v1 = 
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(Teuchos::View, 
							       Y1, 
							       m-i, 
							       1, i, i));

    // Create view of columns i through m-1 of Y2
    h_idx.resize(m-i);
    for (unsigned int j=0; j<h_idx.size(); j++)
      h_idx[j] = i+j;
    h2 = Y2.subView(h_idx);

    // Create view of columns i thru m-1 of R, starting at row i
    h1 = 
      Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(Teuchos::View, 
							       R,
							       m-i,
							       m-i,
							       i, i));

    if (i > 0) {

      // Create view of columns 0 through i-1 of Y2
      y_idx.push_back(i-1);
      y2 = Y2.subView(y_idx);
      
      // Create view of columns 0 through i-1 of Y1, starting at row i
      y1 = 
	Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(Teuchos::View,
								 Y1,
								 m-i,
								 i, i, 0));

      // Create view of column i, row 0 through i-1 of T
      z = 
	Teuchos::rcp(new NOX::Abstract::MultiVector::DenseMatrix(Teuchos::View,
								 T, 
								 i, 
								 1, 
								 0, i));
    }

    // Compute Householder Vector
    computeHouseholderVector(i, R, Y2, *v1, *v2, beta);

    // Apply Householder reflection
    applyHouseholderVector(*v1, *v2, beta, *h1, *h2);

    // Copy v2 into Y2
    Y2[i] = (*v2)[0];
    
    T(i,i) = -beta;

    if (i > 0) {

      // Compute z = y2^T * v2
      v2->multiply(1.0, *y2, *z);

      // Compute z = -beta * (z + y1^T * v1)
      z->multiply(Teuchos::TRANS, Teuchos::NO_TRANS, -beta, *y1, *v1, -beta);

      // Compute z = T * z
      dblas.TRMV(Teuchos::UPPER_TRI, Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG,
		 i, T.values(), m, z->values(), 1);

    }
  }

}

void
LOCA::BorderedSolver::HouseholderQR::computeHouseholderVector(
			  int col,
			  const NOX::Abstract::MultiVector::DenseMatrix& A1,
			  const NOX::Abstract::MultiVector& A2,
			  NOX::Abstract::MultiVector::DenseMatrix& V1,
			  NOX::Abstract::MultiVector& V2,
			  double& beta)
{
  double houseP = A1(col,col);
  
  V1(0,0) = 1.0;
  V2[0] = A2[col];

  double sigma = A2[col].innerProduct(A2[col]);
  for (int i=col+1; i<A1.numRows(); i++)    
    sigma += A1(i,col)*A1(i,col);

  if (sigma == 0.0)
    beta = 0.0;
  else {
    double mu = sqrt(houseP*houseP + sigma);
    if (houseP <= 0.0)
      houseP = houseP - mu;
    else
      houseP = -sigma / (houseP + mu);
    beta = 2.0*houseP*houseP/(sigma + houseP*houseP);
    
    V2.scale(1.0/houseP);
    for (int i=1; i<V1.numRows(); i++)
      V1(i,0) = A1(i+col,col) / houseP;
  }


  return;
}

void
LOCA::BorderedSolver::HouseholderQR::applyHouseholderVector(
			   const NOX::Abstract::MultiVector::DenseMatrix& V1,
			   const NOX::Abstract::MultiVector& V2,
			   double beta,
			   NOX::Abstract::MultiVector::DenseMatrix& A1,
			   NOX::Abstract::MultiVector& A2)
{
  int nColsA = A2.numVectors();

  // Compute u = V2^T * A2
  NOX::Abstract::MultiVector::DenseMatrix u(1, nColsA);
  A2.multiply(1.0, V2, u);

  // Compute u = u + V1^T * A_P
  u.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, V1, A1, 1.0);

  // Compute A1 = A1 - b*V1*u
  A1.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, -beta, V1, u, 1.0);

  // Compute A2 = A2 - b*V2*u
  A2.update(Teuchos::NO_TRANS, -beta, V2, u, 1.0);
}

void
LOCA::BorderedSolver::HouseholderQR::applyCompactWY(
			    const NOX::Abstract::MultiVector::DenseMatrix& Y1,
			    const NOX::Abstract::MultiVector& Y2,
			    const NOX::Abstract::MultiVector::DenseMatrix& T,
			    NOX::Abstract::MultiVector::DenseMatrix& X1,
			    NOX::Abstract::MultiVector& X2,
			    bool isZeroX1, bool isZeroX2,
			    bool useTranspose) const
{
  if (isZeroX1 && isZeroX2) {
    X1.putScalar(0.0);
    X2.init(0.0);
    return;
  }

  int m = Y2.numVectors();

  Teuchos::ETransp T_flag;
  if (useTranspose)
    T_flag = Teuchos::TRANS;
  else
    T_flag = Teuchos::NO_TRANS;

  NOX::Abstract::MultiVector::DenseMatrix tmp(m, X2.numVectors());

  // Compute Y1^T*X1 + Y2^T*X2
  if (!isZeroX2)
    X2.multiply(1.0, Y2, tmp);

  // Opportunity for optimization here since Y1 is a lower-triangular
  // matrix with unit diagonal
  if (!isZeroX2 && !isZeroX1)
    tmp.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Y1, X1, 1.0);
  else if (!isZeroX1)
    tmp.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Y1, X1, 0.0);

  // Compute op(T)*(Y1^T*X1 + Y2^T*X2)
  dblas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, T_flag, 
	     Teuchos::NON_UNIT_DIAG, tmp.numRows(), tmp.numCols(), 1.0, 
	     T.values(), T.numRows(), tmp.values(), tmp.numRows());

  // Compute X1 = X1 + Y1*op(T)*(Y1^T*X1 + Y2^T*X2)
  // Opportunity for optimization here since Y1 is a lower-triangular
  // matrix with unit diagonal
  if (isZeroX1)
    X1.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Y1, tmp, 0.0);
  else
    X1.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Y1, tmp, 1.0);

  // Compute X2 = X2 + Y1*op(T)*(Y1^T*X1 + Y2^T*X2)
  if (isZeroX2)
    X2.update(Teuchos::NO_TRANS, 1.0, Y2, tmp, 0.0);
  else
    X2.update(Teuchos::NO_TRANS, 1.0, Y2, tmp, 1.0); 
}

void
LOCA::BorderedSolver::HouseholderQR::applyCompactWY(
		     const NOX::Abstract::MultiVector::DenseMatrix& Y1,
		     const NOX::Abstract::MultiVector& Y2,
		     const NOX::Abstract::MultiVector::DenseMatrix& T,
		     const NOX::Abstract::MultiVector::DenseMatrix* input1,
		     const NOX::Abstract::MultiVector* input2,
		     NOX::Abstract::MultiVector::DenseMatrix& result1,
		     NOX::Abstract::MultiVector& result2,
		     bool useTranspose) const
{
  bool isZeroX1 = (input1 == NULL);
  bool isZeroX2 = (input2 == NULL);

  if (!isZeroX1)
    result1.assign(*input1);
  if (!isZeroX2)
    result2 = *input2;

  applyCompactWY(Y1, Y2, T, result1, result2, isZeroX1, isZeroX2, 
		 useTranspose);
}

