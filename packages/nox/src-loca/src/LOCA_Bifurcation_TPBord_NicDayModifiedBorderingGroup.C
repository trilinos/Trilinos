// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "LOCA_Bifurcation_TPBord_NicDayModifiedBorderingGroup.H"
#include "LOCA_Bifurcation_TPBord_AbstractGroup.H"
#include "LOCA_Parameter_Vector.H"
#include "NOX_Parameter_List.H"
#include "NOX_LAPACK_Wrappers.H"
#include "LOCA_ErrorCheck.H"

LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup::NicDayModifiedBorderingGroup(
			       LOCA::Bifurcation::TPBord::AbstractGroup& g,
			       NOX::Parameter::List& bifParamList)   
  : LOCA::Bifurcation::TPBord::ExtendedGroup(g, bifParamList)
{
}

LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup::NicDayModifiedBorderingGroup(
			      LOCA::Bifurcation::TPBord::AbstractGroup& g,
			      const NOX::Abstract::Vector& lenVec,
			      const NOX::Abstract::Vector& nullVec,
			      int paramId)
  : LOCA::Bifurcation::TPBord::ExtendedGroup(g, lenVec, nullVec, paramId)
{
}

LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup::NicDayModifiedBorderingGroup(
			  const LOCA::Bifurcation::TPBord::AbstractGroup& g,
			  const NOX::Abstract::Vector& lenVec,
			  const NOX::Abstract::Vector& nullVec,
			  int paramId)
  : LOCA::Bifurcation::TPBord::ExtendedGroup(g, lenVec, nullVec, paramId)
{
}

LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup::NicDayModifiedBorderingGroup(
	     const LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup& source, 
	     NOX::CopyType type)
  : LOCA::Bifurcation::TPBord::ExtendedGroup(source, type)
{
}


LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup::~NicDayModifiedBorderingGroup() 
{
}

NOX::Abstract::Group&
LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup::operator=(
		       const NOX::Abstract::Group& source)
{
  return *this = 
    dynamic_cast<const LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup&>(source);
}

LOCA::Continuation::AbstractGroup&
LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup::operator=(
		       const LOCA::Continuation::AbstractGroup& source)
{
  return *this = 
    dynamic_cast<const LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup&>(source);
}

LOCA::Extended::AbstractGroup&
LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup::operator=(
		       const LOCA::Extended::AbstractGroup& source)
{
  return *this = 
    dynamic_cast<const LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup&>(source);
}

LOCA::Bifurcation::TPBord::ExtendedGroup&
LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup::operator=(
		       const LOCA::Bifurcation::TPBord::ExtendedGroup& source)
{
  return *this = 
    dynamic_cast<const LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup&>(source);
}

LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup&
LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup::operator=(
	     const LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup& source) 
{
  LOCA::Bifurcation::TPBord::ExtendedGroup::operator=(source);

  return *this;
}

NOX::Abstract::Group*
LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup::clone(NOX::CopyType type) const 
{
  return new LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup(*this, type);
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup::applyJacobianInverse(
					NOX::Parameter::List& params,
					const NOX::Abstract::Vector& input,
					NOX::Abstract::Vector& result) const 
{
  string callingFunction = 
    "LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup::applyJacobianInverse()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  if (!isJacobian()) {
    LOCA::ErrorCheck::throwError(callingFunction,
				 "Called with invalid Jacobian!");
  }

  // cast vectors to turning point vectors
  const LOCA::Bifurcation::TPBord::ExtendedVector& tp_input = 
    dynamic_cast<const LOCA::Bifurcation::TPBord::ExtendedVector&>(input);
  LOCA::Bifurcation::TPBord::ExtendedVector& tp_result = 
    dynamic_cast<LOCA::Bifurcation::TPBord::ExtendedVector&>(result);

  // Get componenets of input vector
  const NOX::Abstract::Vector& input_x = tp_input.getXVec();
  const NOX::Abstract::Vector& input_y = tp_input.getNullVec();
  double input_p = tp_input.getBifParam();

  // Get components of result vector.  Note these are references.
  NOX::Abstract::Vector& result_x = tp_result.getXVec();
  NOX::Abstract::Vector& result_y = tp_result.getNullVec();
  double& result_p = tp_result.getBifParam();

  // Temporary vectors
  NOX::Abstract::Vector* a = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector* b = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector* c = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector* d = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector* e = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector* tmp1 = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector* tmp2 = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector* tmp3 = input_x.clone(NOX::ShapeCopy);
  double aa, bb, cc, dd, ee;

  // Get reference to null vector, Jn vectors
  const NOX::Abstract::Vector& v = tpXVec.getNullVec();
  const NOX::Abstract::Vector& Jv = tpFVec.getNullVec();

  // Compute s = ||Jv||_2, u = v/s
  double s = Jv.norm(NOX::Abstract::Vector::TwoNorm);
  NOX::Abstract::Vector *u = Jv.clone(NOX::DeepCopy);
  u->scale(1.0/s);

  // verify underlying Jacobian is valid
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Compute a, aa
  status = deflatedJacobianInverse(params, *u, input_x, *a, aa);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute b, bb
  status = deflatedJacobianInverse(params, *u, *derivResidualParamPtr, *b, bb);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute input_y - (dJv/dx)*a
  status = grpPtr->computeDJnDxa(v, *a, Jv, *tmp1);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  tmp1->update(1.0, input_y, -1.0);

  // Compute d(Jv)/dp - (dJv/dx)*b
  status = grpPtr->computeDJnDxa(v, *b, Jv, *tmp2);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  tmp2->update(1.0, *derivNullResidualParamPtr, -1.0);

  // Compute (dJv/dx)*v
  status = grpPtr->computeDJnDxa(v, v, Jv, *tmp3);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // verify underlying Jacobian is valid
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Compute c, cc
  status = deflatedJacobianInverse(params, *u, *tmp1, *c, cc);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute d, dd
  status = deflatedJacobianInverse(params, *u, *tmp2, *d, dd);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute e, ee
  status = deflatedJacobianInverse(params, *u, *tmp3, *e, ee);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Fill coefficient arrays
  double A[9], B[3];
  A[0] = s;   A[1] = ee;  A[2] = -lTransNorm(*e);
  A[3] = 0.0; A[4] = s;   A[5] =  lTransNorm(v);
  A[6] = bb;  A[7] = dd;  A[8] = -lTransNorm(*d);

  B[0] = aa;  B[1] = cc;  B[2] = input_p - lTransNorm(*c);

  // Solve A*C = B
  int one = 1;
  int three = 3;
  int piv[3];
  int info;
  DGESV_F77(&three, &one, A, &three, piv, B, &three, &info);
  if (info != 0) {
    LOCA::ErrorCheck::throwError(callingFunction,
				 "Solve of 3x3 coefficient matrix failed!");
    return NOX::Abstract::Group::Failed;
  }

  double alpha = B[0];
  double beta = B[1];
  result_p = B[2];

  // result_x = a - result_p*b + alpha*v
  result_x = *a;
  result_x.update(-result_p, *b, alpha, v, 1.0);

  // result_y = c - result_p*d - alpha*e + beta*v
  result_y.update(1.0, *c, -result_p, *d, 0.0);
  result_y.update(-alpha, *e, beta, v, 1.0);

  delete a;
  delete b;
  delete c;
  delete d;
  delete e;
  delete tmp1;
  delete tmp2;
  delete tmp3;
  delete u;
 
  return NOX::Abstract::Group::Ok;
}

// NOX::Abstract::Group::ReturnType
// LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup::applyJacobianInverseMulti(
// 			    NOX::Parameter::List& params,
// 			    const NOX::Abstract::Vector* const* inputs,
// 			    NOX::Abstract::Vector** results, int nVecs) const
// {
//   string callingFunction = 
//     "LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup::applyJacobianInverseMulti()";
//   NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
//   NOX::Abstract::Group::ReturnType status;

//   if (!isJacobian()) {
//     LOCA::ErrorCheck::throwError(callingFunction,
// 				 "Called with invalid Jacobian!");
//   }

//   // Number of input vectors
//   int m = nVecs; 
  
//   // Build arrays of solution, null vector and parameter components
//   const NOX::Abstract::Vector** inputs_x = 
//     new const NOX::Abstract::Vector*[m+1];
//   const NOX::Abstract::Vector** inputs_null =
//     new const NOX::Abstract::Vector*[m+2];
//   double *inputs_params = new double[m];
//   double *zeros = new double[m+2];

//   NOX::Abstract::Vector** tmp1 = new NOX::Abstract::Vector*[m+2];
//   NOX::Abstract::Vector** tmp2 = new NOX::Abstract::Vector*[m+2];
//   NOX::Abstract::Vector** tmp3 = new NOX::Abstract::Vector*[m+2];
//   double *scalars_x = new double[m+1];
//   double *scalars_null = new double[m+2];

//   // Get reference to null vector, Jn vectors
//   const NOX::Abstract::Vector& v = tpXVec.getNullVec();
//   const NOX::Abstract::Vector& Jv = tpFVec.getNullVec();

//   // Compute s = ||Jv||_2, u = v/s
//   double s = Jv.norm(NOX::Abstract::Vector::TwoNorm);
//   NOX::Abstract::Vector *u = Jv.clone(NOX::DeepCopy);
//   u->scale(1.0/s);

//   const LOCA::Bifurcation::TPBord::ExtendedVector* constTPVecPtr;

//   for (int i=0; i<m; i++) {
//     constTPVecPtr = 
//       dynamic_cast<const LOCA::Bifurcation::TPBord::ExtendedVector*>(inputs[i]);
//     inputs_x[i] = &(constTPVecPtr->getXVec());
//     inputs_null[i] = &(constTPVecPtr->getNullVec());
//     inputs_params[i] = constTPVecPtr->getBifParam();

//     tmp1[i] = inputs_x[i]->clone(NOX::ShapeCopy); tmp1[i]->init(0.0);
//     tmp2[i] = inputs_x[i]->clone(NOX::ShapeCopy); tmp2[i]->init(0.0);
//     tmp3[i] = inputs_x[i]->clone(NOX::ShapeCopy); tmp3[i]->init(0.0);
//     zeros[i] = 0.0;
//     scalars_x[i] = 0.0;
//     scalars_null[i] = 0.0;
//   }

//   // Set last components to deriv. w.r.t. parameter
//   inputs_x[m] = derivResidualParamPtr;
//   inputs_null[m] = derivNullResidualParamPtr;
//   tmp1[m] = inputs_x[m]->clone(NOX::ShapeCopy); tmp1[m]->init(0.0);
//   tmp2[m] = inputs_x[m]->clone(NOX::ShapeCopy); tmp2[m]->init(0.0);
//   tmp3[m] = inputs_x[m]->clone(NOX::ShapeCopy); tmp3[m]->init(0.0);
//   zeros[m] = 0.0;
//   scalars_x[m] = 0.0;
//   scalars_null[m] = 0.0;

//   tmp1[m+1] = v.clone(NOX::DeepCopy);
//   tmp2[m+1] = inputs_x[m]->clone(NOX::ShapeCopy); tmp2[m+1]->init(0.0);
//   tmp3[m+1] = inputs_x[m]->clone(NOX::ShapeCopy); tmp3[m+1]->init(0.0);
//   zeros[m+1] = 0.0;
//   scalars_null[m+1] = 0.0;

//   // verify underlying Jacobian is valid
//   if (!grpPtr->isJacobian()) {
//     status = grpPtr->computeJacobian();
//     finalStatus = 
//       LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
// 						   callingFunction);
//   }

//   // Solve J*tmp1 = inputs_x
//   finalStatus = grpPtr->applyBorderedJacobianInverseMulti(false, params, *u, v,
// 							  inputs_x, zeros, 
// 							  tmp1,
// 							  scalars_x, m+1);
//   LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

//   // Compute tmp2 = inputs_null - (dJv/dx)*tmp1 
//   for (int i=0; i<m+2; i++) {
//     status = grpPtr->computeDJnDxa(v, *tmp1[i], Jv, *tmp2[i]);
//     finalStatus = 
//       LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
// 						   callingFunction);

//     if (i < m+1)
//       tmp2[i]->update(1.0, *inputs_null[i], -1.0);
//   } 

//   // verify underlying Jacobian is valid
//   if (!grpPtr->isJacobian()) {
//     status = grpPtr->computeJacobian();
//     finalStatus = 
//       LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
// 						   callingFunction);
//   }

//   // Solve J*tmp3 = tmp2
//   status = grpPtr->applyBorderedJacobianInverseMulti(false, params, *u, v,
// 						     tmp2, zeros, tmp3,
// 						     scalars_null, m+2);
//   finalStatus = 
//     LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
// 						 callingFunction);

//   // Fill coefficient arrays
//   double A[9];
//   double *B = new double[3*m];
//   NOX::Abstract::Vector& b = *tmp1[m];   double bb = scalars_x[m];
//   NOX::Abstract::Vector& d = *tmp3[m];   double dd = scalars_null[m];
//   NOX::Abstract::Vector& e = *tmp3[m+1]; double ee = scalars_null[m+1];
//   A[0] = s;   A[1] = ee;  A[2] = -lTransNorm(e);
//   A[3] = 0.0; A[4] = s;   A[5] =  lTransNorm(v);
//   A[6] = bb;  A[7] = dd;  A[8] = -lTransNorm(d);

//   for (int i=0; i<m; i++) {
//     B[0 + 3*i] = scalars_x[i];  
//     B[1 + 3*i] = scalars_null[i];  
//     B[2 + 3*i] = inputs_params[i] - lTransNorm(*tmp3[i]);
//   }
//   // Solve A*C = B
//   int three = 3;
//   int piv[3];
//   int info;
//   DGESV_F77(&three, &m, A, &three, piv, B, &three, &info);
//   if (info != 0)
//     return NOX::Abstract::Group::Failed;
  
//   // Compute and set results
//   double alpha, beta, w;
//   LOCA::Bifurcation::TPBord::ExtendedVector* tpVecPtr;
//   for (int i=0; i<m; i++) {
//     tpVecPtr = 
//       dynamic_cast<LOCA::Bifurcation::TPBord::ExtendedVector*>(results[i]);
//     alpha = B[0 + 3*i];
//     beta = B[1 + 3*i];
//     w = B[2 + 3*i];

//     tpVecPtr->getXVec() = *tmp1[i];
//     (tpVecPtr->getXVec()).update(-w, b, alpha, v, 1.0);

//     (tpVecPtr->getNullVec()).update(1.0, *tmp3[i], -w, d, 0.0);
//     (tpVecPtr->getNullVec()).update(-alpha, e, beta, v, 1.0);
//     tpVecPtr->getBifParam() = w;

//     delete tmp1[i];
//     delete tmp2[i];
//     delete tmp3[i];
//   }

//   delete tmp1[m];
//   delete tmp2[m];
//   delete tmp3[m];

//   delete tmp1[m+1];
//   delete tmp2[m+1];
//   delete tmp3[m+1];

//   delete [] tmp1;
//   delete [] tmp2;
//   delete [] tmp3;
//   delete [] inputs_x;
//   delete [] inputs_null;
//   delete [] inputs_params;
//   delete [] zeros;
//   delete [] scalars_x;
//   delete [] scalars_null;
//   delete [] B;

//   return finalStatus;
// }

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup::deflatedJacobianInverse(
					NOX::Parameter::List& params,
					const NOX::Abstract::Vector& u, 
					const NOX::Abstract::Vector& input,
					NOX::Abstract::Vector& result,
					double& r) const 
{
  string callingFunction = 
    "LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup::deflatedJacobianInverse()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // verify underlying Jacobian is valid
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Compute deflatedInput = input - (u^T*input)*u
  NOX::Abstract::Vector *deflatedInput = input.clone(NOX::DeepCopy);
  r = u.dot(input);
  deflatedInput->update(-r, u, 1.0);

  // Solve J*result = deflatedInput
  status = grpPtr->applyJacobianInverse(params, *deflatedInput, result);
  finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);

  delete deflatedInput;

  return finalStatus;
}

  
