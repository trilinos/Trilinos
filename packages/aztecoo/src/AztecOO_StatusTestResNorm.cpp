
/* Copyright (2001) Sandia Corportation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 *
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

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
    restype_(Implicit),
    resnormtype_(TwoNorm),
    scaletype_(NormOfInitRes),
    scalenormtype_(TwoNorm),
    weights_(0),
    scalevalue_(1.0),
    resvalue_(0.0),
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
  weights_ = Weights;

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
  weights_ = Weights;
  scalevalue_ = ScaleValue;

  // These conditions force the residual vector to be computed
  if (scaletype_==NormOfInitRes && scalenormtype_!=TwoNorm) resvecrequired_ = true;
  
  return(0);
}

bool AztecOO_StatusTestResNorm::ResidualVectorRequired() const
{
  return(resvecrequired_);
}

AztecOO_StatusType AztecOO_StatusTestResNorm::CheckStatus(int CurrentIter, 
							      Epetra_MultiVector * CurrentResVector, 
							      double CurrentResNormEst,  
							      bool SolutionUpdated)
{

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
    resvalue_ = ComputeNorm(*localresvector_, resnormtype_);
  }
  else {
    curresvecexplicit_ = false;
    resvalue_ = ComputeNorm(*crv, resnormtype_);
  }

  if (firstcallCheckStatus_) {
    if (scaletype_==NormOfRHS) scalevalue_ = ComputeNorm(rhs_, scalenormtype_);
    else if (scaletype_==NormOfInitRes) 
      if (restype_==Implicit && scalenormtype_==TwoNorm && CurrentResNormEst!=-1.0) 
	scalevalue_ = CurrentResNormEst;
      else
	scalevalue_ = ComputeNorm(rhs_, scalenormtype_);
    if (scalevalue_==0.0) {
      status_ = Failed;
      return(status_);
    } 
  }

  testvalue_ = resvalue_/scalevalue_;
  if (testvalue_>tolerance_)
    status_ = Unconverged;
  else if (testvalue_<=tolerance_)
    status_ = Converged;
  else
    status_ = NaN;

  firstcallCheckStatus_ = false;
  return status_;
}

AztecOO_StatusType AztecOO_StatusTestResNorm::GetStatus() const
{
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
  stream << ((resnormtype_==OneNorm) ? "1-Norm" : (resnormtype_==TwoNorm) ? "2-Norm" : "Inf-Norm");
  stream << ((curresvecexplicit_) ? " Exp" : " Imp");
  stream << " Res Vec) ";
  if (scaletype_!=None)
    stream << "/";
  if (scaletype_==UserProvided)
    stream << " (User Scale)";
  else {
    stream << "(";
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
