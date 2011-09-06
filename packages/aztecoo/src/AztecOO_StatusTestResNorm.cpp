
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#include "AztecOO_StatusTestResNorm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"

AztecOO_StatusTestResNorm::AztecOO_StatusTestResNorm(const Epetra_Operator & Operator, 
						     const Epetra_Vector & LHS, const Epetra_Vector & RHS,
						     double Tolerance)
  : AztecOO_StatusTest(),
    operator_(Operator),
    lhs_(LHS),
    rhs_(RHS),
    tolerance_(Tolerance),
    maxNumExtraIterations_(0),
    numExtraIterations_(0),
    restype_(Implicit),
    resnormtype_(TwoNorm),
    scaletype_(NormOfInitRes),
    scalenormtype_(TwoNorm),
    resweights_(0),
    scaleweights_(0),
    scalevalue_(1.0),
    resvalue_(0.0),
    convergedOnce_(false),
    status_(Unchecked),
    resvecrequired_(false),
    firstcallCheckStatus_(true),
    firstcallDefineResForm_(true),
    firstcallDefineScaleForm_(true),
    localresvector_(0),
    curresvecexplicit_(false)
{
  // At this point we are prepared to use ||r||/||r0|| <= tolerance using the 2-norm of the 
  // residual vector (or its estimate received from the iterative method).
}
AztecOO_StatusTestResNorm::~AztecOO_StatusTestResNorm() {

  if (localresvector_!=0) { delete localresvector_; localresvector_ = 0;}

}
int AztecOO_StatusTestResNorm::DefineResForm( ResType TypeOfResidual, NormType TypeOfNorm, 
					      Epetra_Vector * Weights)
{

  if (!firstcallDefineResForm_) EPETRA_CHK_ERR(-1); // We can only have this routine called once.
  firstcallDefineResForm_ = false;

  restype_ = TypeOfResidual;
  resnormtype_ = TypeOfNorm;
  resweights_ = Weights;

  // These conditions force the residual vector to be computed
  if (restype_==Explicit ||
      resnormtype_!=TwoNorm) resvecrequired_ = true;
  
  return(0);
}

int AztecOO_StatusTestResNorm::DefineScaleForm(ScaleType TypeOfScaling, NormType TypeOfNorm, 
					       Epetra_Vector * Weights, 
					       double ScaleValue )
{

  if (!firstcallDefineScaleForm_) EPETRA_CHK_ERR(-1); // We can only have this routine called once.
  firstcallDefineScaleForm_ = false;

  scaletype_ = TypeOfScaling;
  scalenormtype_ = TypeOfNorm;
  scaleweights_ = Weights;
  scalevalue_ = ScaleValue;

  // These conditions force the residual vector to be computed
  if (scaletype_==NormOfInitRes && scalenormtype_!=TwoNorm) resvecrequired_ = true;
  
  return(0);
}

bool AztecOO_StatusTestResNorm::ResidualVectorRequired() const
{
  return(resvecrequired_);
}

AztecOO_StatusType
AztecOO_StatusTestResNorm::CheckStatus(int CurrentIter, 
                                       Epetra_MultiVector * CurrentResVector, 
                                       double CurrentResNormEst,  
                                       bool SolutionUpdated)
{
  (void)CurrentIter;
  Epetra_Vector * crv = dynamic_cast<Epetra_Vector *>(CurrentResVector);

  // This section computes the norm of the residual vector
  if (restype_==Implicit && resnormtype_==TwoNorm && CurrentResNormEst!=-1.0) 
    resvalue_ = CurrentResNormEst;
  else if (crv==0) { // Cannot proceed because there is no norm est or res vector
    status_ = Failed;
    return(status_);
  }
  else if (restype_==Explicit && SolutionUpdated) {
    curresvecexplicit_ = true;
    if (localresvector_==0) localresvector_ = new Epetra_Vector(crv->Map());
    // Compute explicit residual
    operator_.Apply(lhs_, *localresvector_);
    localresvector_->Update(1.0, rhs_, -1.0); // localresvector_ = rhs_ - operator_* lhs_
    if (resweights_!=0) { // Check if we should scale the vector
      // localresvector_ = resweights_ * localresvector_
      localresvector_->Multiply(1.0, *resweights_, *localresvector_, 0.0);
    }
    resvalue_ = ComputeNorm(*localresvector_, resnormtype_);
  }
  else {
    curresvecexplicit_ = false;
    if (resweights_!=0) { // Check if we should scale the vector
      if (localresvector_==0) localresvector_ = new Epetra_Vector(crv->Map());
      // localresvector_ = resweights_ * localresvector_
      localresvector_->Multiply(1.0, *resweights_, *crv, 0.0);
      resvalue_ = ComputeNorm(*localresvector_, resnormtype_);
    }
    else
      resvalue_ = ComputeNorm(*crv, resnormtype_);
  }


  // Compute scaling term (done once)
  if (firstcallCheckStatus_) {
    if (scaletype_==NormOfRHS) {
      if (scaleweights_!=0) { // Check if we should scale the vector
	if (localresvector_==0) localresvector_ = new Epetra_Vector(rhs_.Map());
	// localresvector = scaleweights_ * rhs_
	localresvector_->Multiply(1.0, *scaleweights_, rhs_, 0.0);
	scalevalue_ = ComputeNorm(*localresvector_, resnormtype_);
      }
      else {
	scalevalue_ = ComputeNorm(rhs_, scalenormtype_);
      }
    }
    else if (scaletype_==NormOfInitRes) {
      if (restype_==Implicit && scalenormtype_==TwoNorm && CurrentResNormEst!=-1.0) 
	scalevalue_ = CurrentResNormEst;
      else {
	if (scaleweights_!=0) { // Check if we should scale the vector
	  if (localresvector_==0) localresvector_ = new Epetra_Vector(crv->Map());
	  // weightedrhs = scaleweights_ * initial residual
	  localresvector_->Multiply(1.0, *scaleweights_, *crv, 0.0);
	  scalevalue_ = ComputeNorm(*localresvector_, resnormtype_);
	}
	else {
	  scalevalue_ = ComputeNorm(rhs_, scalenormtype_);
	}
      }
    }
    if (scalevalue_==0.0) {
      status_ = Failed;
      return(status_);
    } 
  }

  testvalue_ = resvalue_/scalevalue_;
  if (testvalue_>tolerance_) 
    if (convergedOnce_) {
      if (numExtraIterations_<maxNumExtraIterations_) {
	numExtraIterations_++;
	status_ = Unconverged;
      }
      else {
	status_ = PartialFailed;
      }
    }
    else {
      status_ = Unconverged;
    }
  else if (testvalue_<=tolerance_) {
    convergedOnce_ = true;	  
    status_ = Converged;
  }
  else
    status_ = NaN;

  firstcallCheckStatus_ = false;
  return status_;
}

 double AztecOO_StatusTestResNorm::ComputeNorm(const Epetra_Vector & vec, NormType typeofnorm) {

   double result = 0.0;
   if (typeofnorm==TwoNorm) vec.Norm2(&result);
   else if (typeofnorm==OneNorm) vec.Norm1(&result);
   else vec.NormInf(&result);
   return(result);
 }

ostream& AztecOO_StatusTestResNorm::Print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  PrintStatus(stream, status_);
  stream << "(";
  if (resweights_!=0) stream << "Weighted ";
  stream << ((resnormtype_==OneNorm) ? "1-Norm" : (resnormtype_==TwoNorm) ? "2-Norm" : "Inf-Norm");
  stream << ((curresvecexplicit_) ? " Exp" : " Imp");
  stream << " Res Vec) ";
  if (scaletype_!=None)
    stream << "/";
  if (scaletype_==UserProvided)
    stream << " (User Scale)";
  else {
    stream << "(";
    if (scaleweights_!=0) stream << "Weighted ";
    stream << ((scalenormtype_==OneNorm) ? "1-Norm" : (resnormtype_==TwoNorm) ? "2-Norm" : "Inf-Norm");
    if (scaletype_==NormOfInitRes)
      stream << " Res0";
    else
      stream << " RHS ";
    stream << ")";
  }
  if (status_==Unchecked)
    stream << " Unchecked << ";
  else {
    stream << " = " << testvalue_;
    stream << ((testvalue_<tolerance_) ? " < " : (testvalue_==tolerance_) ? " = " : (testvalue_>tolerance_) ? " > " : " <> ");
  }
  stream << tolerance_;
  stream << endl;
    
  return stream;
}
