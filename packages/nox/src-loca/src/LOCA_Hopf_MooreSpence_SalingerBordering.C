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

#include "LOCA_Hopf_MooreSpence_SalingerBordering.H"
#include "LOCA_Hopf_MooreSpence_ExtendedGroup.H"
#include "LOCA_Hopf_MooreSpence_AbstractGroup.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::Hopf::MooreSpence::SalingerBordering::SalingerBordering(
	 const Teuchos::RCP<LOCA::GlobalData>& global_data,
	 const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& slvrParams) : 
  globalData(global_data),
  solverParams(slvrParams),
  group(),
  hopfGroup(),
  yVector(),
  zVector(),
  CeRealVector(),
  CeImagVector(),
  dfdp(),
  dCedpReal(),
  dCedpImag(),
  ByVector(),
  minusBzVector()
{
}

LOCA::Hopf::MooreSpence::SalingerBordering::~SalingerBordering()
{
}

void
LOCA::Hopf::MooreSpence::SalingerBordering::setBlocks(
	 const Teuchos::RCP<LOCA::Hopf::MooreSpence::AbstractGroup>& group_,
	 const Teuchos::RCP<LOCA::Hopf::MooreSpence::ExtendedGroup>& hopfGroup_,
	 const Teuchos::RCP<const NOX::Abstract::Vector>& yVector_,
	 const Teuchos::RCP<const NOX::Abstract::Vector>& zVector_,
	 const Teuchos::RCP<const NOX::Abstract::Vector>& CeRealVector_,
	 const Teuchos::RCP<const NOX::Abstract::Vector>& CeImagVector_,
	 const Teuchos::RCP<const NOX::Abstract::Vector>& dfdp_,
	 const Teuchos::RCP<const NOX::Abstract::Vector>& dCedpReal_,
	 const Teuchos::RCP<const NOX::Abstract::Vector>& dCedpImag_,
	 const Teuchos::RCP<const NOX::Abstract::Vector>& ByVector_,
	 const Teuchos::RCP<const NOX::Abstract::Vector>& mBzVector_,
	 double w_)
{
  group = group_;
  hopfGroup = hopfGroup_;
  yVector = yVector_;
  zVector = zVector_;
  CeRealVector = CeRealVector_;
  CeImagVector = CeImagVector_;
  dfdp = dfdp_;
  dCedpReal = dCedpReal_;
  dCedpImag = dCedpImag_;
  ByVector = ByVector_;
  minusBzVector = mBzVector_;
  w = w_;
}

NOX::Abstract::Group::ReturnType 
LOCA::Hopf::MooreSpence::SalingerBordering::solve(
	   Teuchos::ParameterList& params,
	   const LOCA::Hopf::MooreSpence::ExtendedMultiVector& input,
           LOCA::Hopf::MooreSpence::ExtendedMultiVector& result) const
{
  std::string callingFunction = 
    "LOCA::Hopf::MooreSpence::SalingerBordering::solve()";
  NOX::Abstract::Group::ReturnType status;
  
  // Get components of input
  Teuchos::RCP<const NOX::Abstract::MultiVector> input_x = 
    input.getXMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector> input_y = 
    input.getRealEigenMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector> input_z = 
    input.getImagEigenMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> input_w 
    = input.getFrequencies();
  Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> input_p 
    = input.getBifParams();

  // Get components of result
  Teuchos::RCP<NOX::Abstract::MultiVector> result_x = 
    result.getXMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector> result_y = 
    result.getRealEigenMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector> result_z = 
    result.getImagEigenMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> result_w = 
    result.getFrequencies();
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> result_p = 
    result.getBifParams();

  int m = input.numVectors();

  std::vector<int> index_input(m);
  for (int i=0; i<m; i++)
    index_input[i] = i;
  
  // Create new multivectors with m+2 columns
  // First m columns store input_x, input_y, input_z, result_x, result_y, 
  // and result_z respectively, next column stores df/dp, (Jy-wBz)_p, 
  // (Jz+wBy)_p, J^-1 dfdp, C^-1 (Jy-wBz)_p, C^-1 (Jz+wBy)_p
  // respectively.  Last column stores -Bz, By C^-1 (-Bz) and C^-1 (By)
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_input_x = 
    input_x->clone(m+1);
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_input_y = 
    input_y->clone(m+2);
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_input_z = 
    input_z->clone(m+2);
  
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_result_x = 
    result_x->clone(m+1);
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_result_y = 
    result_y->clone(m+2);
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_result_z = 
    result_y->clone(m+2);
  
  // Set first m columns to input_x
  cont_input_x->setBlock(*input_x, index_input);

  // Set last column to dfdp
  (*cont_input_x)[m] = *dfdp;
  
  // Set first m columns to input_y
  cont_input_y->setBlock(*input_y, index_input);
  
  // Set next column to (Jy-wBz)_p
  (*cont_input_y)[m] = *dCedpReal;

  // Set last column to -Bz
  (*cont_input_y)[m+1] = *minusBzVector;

  // Set first m columns to input_z
  cont_input_z->setBlock(*input_z, index_input);
  
  // Set next column to (Jz+wBy)_p
  (*cont_input_z)[m] = *dCedpImag;

  // Set last column to By
  (*cont_input_z)[m+1] = *ByVector;
  
  // Initialize result multivectors to 0
  cont_result_x->init(0.0);
  cont_result_y->init(0.0);
  cont_result_z->init(0.0);
    
  // Solve
  status = solveContiguous(params, *cont_input_x, *cont_input_y, *cont_input_z,
			   *input_w, *input_p, *cont_result_x, *cont_result_y,
			   *cont_result_z, *result_w, *result_p);
  
  // Create views of first m columns for result_x, result_y, result_z
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_result_x_view = 
    cont_result_x->subView(index_input);
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_result_y_view = 
    cont_result_y->subView(index_input);
  Teuchos::RCP<NOX::Abstract::MultiVector> cont_result_z_view = 
    cont_result_z->subView(index_input);
  
  // Copy first m columns back into result_x, result_null
  *result_x = *cont_result_x_view;
  *result_y = *cont_result_y_view;
  *result_z = *cont_result_z_view;

   return status;
}

// Solves Hopf equations via classic Salinger bordering
// The first m columns of input_x, input_y, input_z store the RHS, the
// next column stores df/dp, (Jy-wBz)_p and (Jz+wBy)_p respectively, the
// last column of input_y and input_z store Bz and -By respectively.  Note 
// input_x has m+1 columns, input_y and input_z have m+2, and input_w and
// input_p have m columns.  result_x, result_y, result_z, result_w and 
// result_param have the same dimensions as their input counterparts
NOX::Abstract::Group::ReturnType 
LOCA::Hopf::MooreSpence::SalingerBordering::solveContiguous(
		      Teuchos::ParameterList& params,
		      const NOX::Abstract::MultiVector& input_x,
		      const NOX::Abstract::MultiVector& input_y,
		      const NOX::Abstract::MultiVector& input_z,
		      const NOX::Abstract::MultiVector::DenseMatrix& input_w,
		      const NOX::Abstract::MultiVector::DenseMatrix& input_p,
		      NOX::Abstract::MultiVector& result_x,
		      NOX::Abstract::MultiVector& result_y,
		      NOX::Abstract::MultiVector& result_z,
		      NOX::Abstract::MultiVector::DenseMatrix& result_w,
	              NOX::Abstract::MultiVector::DenseMatrix& result_p) const
{
  std::string callingFunction = 
    "LOCA::Hopf::MooreSpence::SalingerBordering::solveContiguous()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  int m = input_x.numVectors()-1;
  std::vector<int> index_input(m);
  std::vector<int> index_dp(1);
  std::vector<int> index_B(1);
  std::vector<int> index_ip(m+1);
  for (int i=0; i<m; i++) {
    index_input[i] = i;
    index_ip[i] = i;
  }
  index_ip[m] = m;
  index_dp[0] = m;
  index_B[0] = m+1;

  // verify underlying Jacobian is valid
  if (!group->isJacobian()) {
    status = group->computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }
  
  // compute [A b] = J^-1 [F df/dp]
  status = group->applyJacobianInverseMultiVector(params, input_x, result_x);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  Teuchos::RCP<NOX::Abstract::MultiVector> A = 
    result_x.subView(index_input);
  Teuchos::RCP<NOX::Abstract::MultiVector> b = 
    result_x.subView(index_dp);

  // verify underlying complex matrix is valid
   if (!group->isComplex()) {
    status = group->computeComplex(w);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  // compute (J+iwB)(y+iz)_x [A b]
  Teuchos::RCP<NOX::Abstract::MultiVector> tmp_real = 
    result_y.clone(NOX::ShapeCopy);
  Teuchos::RCP<NOX::Abstract::MultiVector> tmp_real_sub =
    tmp_real->subView(index_ip);
  Teuchos::RCP<NOX::Abstract::MultiVector> tmp_imag = 
    result_y.clone(NOX::ShapeCopy);
  Teuchos::RCP<NOX::Abstract::MultiVector> tmp_imag_sub =
    tmp_imag->subView(index_ip);
  tmp_real->init(0.0);
  tmp_imag->init(0.0);
  status = group->computeDCeDxa(*yVector, *zVector, w, result_x,
				*CeRealVector, *CeImagVector, *tmp_real_sub,
				*tmp_imag_sub);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  // compute [G+iH d(J+iwB)(y+iz)/dp iB(y+iz)] - [(J+iwB)_x[A b] 0+i0]
  tmp_real->update(1.0, input_y, -1.0);
  tmp_imag->update(1.0, input_z, -1.0);

  // verify underlying complex matrix is valid
  if (!group->isComplex()) {
    status = group->computeComplex(w);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  // compute [C+iD e+if g+ih] = (J+iwB)^-1 (tmp_real + i tmp_imag)
  status = group->applyComplexInverseMultiVector(params, *tmp_real, *tmp_imag,
						 result_y, result_z);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  Teuchos::RCP<NOX::Abstract::MultiVector> C = 
    result_y.subView(index_input);
  Teuchos::RCP<NOX::Abstract::MultiVector> D = 
    result_z.subView(index_input);
  Teuchos::RCP<NOX::Abstract::MultiVector> e = 
    result_y.subView(index_dp);
  Teuchos::RCP<NOX::Abstract::MultiVector> f = 
    result_z.subView(index_dp);
  Teuchos::RCP<NOX::Abstract::MultiVector> g = 
    result_y.subView(index_B);
  Teuchos::RCP<NOX::Abstract::MultiVector> h = 
    result_z.subView(index_B);

  // compute lambda = ((phi^T h)(phi^T C-u) - (phi^T g)(phi^T D-v)) /
  //                  ((phi^T h)(phi^T e)-(phi^T g)(phi^T f))
  NOX::Abstract::MultiVector::DenseMatrix ltC(1,m);
  NOX::Abstract::MultiVector::DenseMatrix ltD(1,m);
  double lte = hopfGroup->lTransNorm((*e)[0]);
  double ltf = hopfGroup->lTransNorm((*f)[0]);
  double ltg = hopfGroup->lTransNorm((*g)[0]);
  double lth = hopfGroup->lTransNorm((*h)[0]);
  double denom = lth*lte - ltg*ltf;
  hopfGroup->lTransNorm(*C, ltC); 
  ltC -= input_w; 
  ltC.scale(lth);
  hopfGroup->lTransNorm(*D, ltD); 
  ltD -= input_p; 
  result_p.assign(ltD);
  result_p.scale(-ltg);
  result_p += ltC;
  result_p.scale(1.0/denom);

  // compute omega = (phi^T D-v - (phi^T f)lambda)/(phi^T h)
  result_w.assign(result_p);
  result_w.scale(-ltf);
  result_w += ltD;
  result_w.scale(1.0/lth);

  // compute A = A - b*lambda (remember A is a sub-view of result_x)
  A->update(Teuchos::NO_TRANS, -1.0, *b, result_p, 1.0);

  // compute C = C - e*lambda - g*omega (remember C is a sub-view of result_y)
  C->update(Teuchos::NO_TRANS, -1.0, *e, result_p, 1.0);
  C->update(Teuchos::NO_TRANS, -1.0, *g, result_w, 1.0);

  // compute D = D - f*lambda - h*omega (remember D is a sub-view of result_z)
  D->update(Teuchos::NO_TRANS, -1.0, *f, result_p, 1.0);
  D->update(Teuchos::NO_TRANS, -1.0, *h, result_w, 1.0);

  return finalStatus;
}
