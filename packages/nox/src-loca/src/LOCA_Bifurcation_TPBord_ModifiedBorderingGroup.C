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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "LOCA_Bifurcation_TPBord_ModifiedBorderingGroup.H"
#include "LOCA_Bifurcation_TPBord_AbstractGroup.H"
#include "LOCA_Parameter_Vector.H"
#include "NOX_Parameter_List.H"
#include "NOX_LAPACK_Wrappers.H"

LOCA::Bifurcation::TPBord::ModifiedBorderingGroup::ModifiedBorderingGroup(
			      LOCA::Bifurcation::TPBord::AbstractGroup& g,
			      const NOX::Abstract::Vector& lenVec,
			      const NOX::Abstract::Vector& nullVec,
			      int paramId)
  : LOCA::Bifurcation::TPBord::ExtendedGroup(g, lenVec, nullVec, paramId)
{
}

LOCA::Bifurcation::TPBord::ModifiedBorderingGroup::ModifiedBorderingGroup(
			  const LOCA::Bifurcation::TPBord::AbstractGroup& g,
			  const NOX::Abstract::Vector& lenVec,
			  const NOX::Abstract::Vector& nullVec,
			  int paramId)
  : LOCA::Bifurcation::TPBord::ExtendedGroup(g, lenVec, nullVec, paramId)
{
}

LOCA::Bifurcation::TPBord::ModifiedBorderingGroup::ModifiedBorderingGroup(
	     const LOCA::Bifurcation::TPBord::ModifiedBorderingGroup& source, 
	     NOX::CopyType type)
  : LOCA::Bifurcation::TPBord::ExtendedGroup(source, type)
{
}


LOCA::Bifurcation::TPBord::ModifiedBorderingGroup::~ModifiedBorderingGroup() 
{
}

NOX::Abstract::Group&
LOCA::Bifurcation::TPBord::ModifiedBorderingGroup::operator=(
		       const NOX::Abstract::Group& source)
{
  return *this = 
    dynamic_cast<const LOCA::Bifurcation::TPBord::ModifiedBorderingGroup&>(source);
}

LOCA::Continuation::AbstractGroup&
LOCA::Bifurcation::TPBord::ModifiedBorderingGroup::operator=(
		       const LOCA::Continuation::AbstractGroup& source)
{
  return *this = 
    dynamic_cast<const LOCA::Bifurcation::TPBord::ModifiedBorderingGroup&>(source);
}

LOCA::Extended::AbstractGroup&
LOCA::Bifurcation::TPBord::ModifiedBorderingGroup::operator=(
		       const LOCA::Extended::AbstractGroup& source)
{
  return *this = 
    dynamic_cast<const LOCA::Bifurcation::TPBord::ModifiedBorderingGroup&>(source);
}

LOCA::Bifurcation::TPBord::ExtendedGroup&
LOCA::Bifurcation::TPBord::ModifiedBorderingGroup::operator=(
		       const LOCA::Bifurcation::TPBord::ExtendedGroup& source)
{
  return *this = 
    dynamic_cast<const LOCA::Bifurcation::TPBord::ModifiedBorderingGroup&>(source);
}

LOCA::Bifurcation::TPBord::ModifiedBorderingGroup&
LOCA::Bifurcation::TPBord::ModifiedBorderingGroup::operator=(
	     const LOCA::Bifurcation::TPBord::ModifiedBorderingGroup& source) 
{
  LOCA::Bifurcation::TPBord::ExtendedGroup::operator=(source);

  return *this;
}

NOX::Abstract::Group*
LOCA::Bifurcation::TPBord::ModifiedBorderingGroup::clone(NOX::CopyType type) const 
{
  return new LOCA::Bifurcation::TPBord::ModifiedBorderingGroup(*this, type);
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBord::ModifiedBorderingGroup::applyJacobianInverse(
					NOX::Parameter::List& params,
					const NOX::Abstract::Vector& input,
					NOX::Abstract::Vector& result) const 
{
  NOX::Abstract::Group::ReturnType res;

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
  NOX::Abstract::Vector* tmp = input_x.clone(NOX::ShapeCopy);
  double aa, bb, cc, dd, ee;

  // Get reference to null vector, Jn vectors
  const NOX::Abstract::Vector& v = tpXVec.getNullVec();
  const NOX::Abstract::Vector& Jv = tpFVec.getNullVec();

  // Compute s = ||Jv||_2, u = v/s
  double s = Jv.norm(NOX::Abstract::Vector::TwoNorm);
  NOX::Abstract::Vector *u = Jv.clone(NOX::DeepCopy);
  u->scale(1.0/s);

  // Compute a
  res = grpPtr->applyBorderedJacobianInverse(false, params, *u, v, input_x, 
					     0.0, *a, aa);

  // Compute b
  res = grpPtr->applyBorderedJacobianInverse(false, params, *u, v, 
					     *derivResidualParamPtr, 0.0,
					     *b, bb);

  // Compute input_y - (dJv/dx)*a
  res = grpPtr->computeDJnDxa(v, *a, Jv, *tmp);
  tmp->update(1.0, input_y, -1.0);

  // Compute c
  res = grpPtr->applyBorderedJacobianInverse(false, params, *u, v, *tmp, 0.0,
					     *c, cc);

  // Compute d(Jv)/dp - (dJv/dx)*b
  res = grpPtr->computeDJnDxa(v, *b, Jv, *tmp);
  tmp->update(1.0, *derivNullResidualParamPtr, -1.0);

  // Compute d
  res = grpPtr->applyBorderedJacobianInverse(false, params, *u, v, *tmp, 0.0,
					     *d, dd);

  // Compute (dJv/dx)*v
  res = grpPtr->computeDJnDxa(v, v, Jv, *tmp);

  // Compute e
  res = grpPtr->applyBorderedJacobianInverse(false, params, *u, v, *tmp, 0.0,
					     *e, ee);

  // Fill coefficient arrays
  double A[9], B[3];
  A[0] = s;   A[1] = ee;  A[2] = -lengthVecPtr->dot(*e)/ lengthVecPtr->length();
  A[3] = 0.0; A[4] = s;   A[5] =  lengthVecPtr->dot(v)/ lengthVecPtr->length();
  A[6] = bb;  A[7] = dd;  A[8] = -lengthVecPtr->dot(*d)/ lengthVecPtr->length();

  B[0] = aa;  B[1] = cc;  B[2] = input_p - lengthVecPtr->dot(*c)/ lengthVecPtr->length();

  // Solve A*C = B
  int one = 1;
  int three = 3;
  int piv[3];
  int info;
  DGESV_F77(&three, &one, A, &three, piv, B, &three, &info);
  if (info != 0)
    return NOX::Abstract::Group::Failed;

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
  delete tmp;
  delete u;

//   NOX::Abstract::Vector *residual = input.clone(NOX::ShapeCopy);
//   applyJacobian(result, *residual);
//   residual->update(-1.0, input, 1.0);
//   residual->scale(1.0/input.norm(NOX::Abstract::Vector::TwoNorm));
//   cout << "Scaled L2 Norm of linear solve residual = " 
//        << residual->norm(NOX::Abstract::Vector::TwoNorm) << endl;

//   delete residual;
 
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBord::ModifiedBorderingGroup::applyJacobianInverseMulti(
			    NOX::Parameter::List& params,
			    const NOX::Abstract::Vector* const* inputs,
			    NOX::Abstract::Vector** results, int nVecs) const 
{
  NOX::Abstract::Group::ReturnType res;

  for (int i=0; i<nVecs; i++) {
    res = applyJacobianInverse(params, *(inputs[i]), *(results[i]));
    if (res != NOX::Abstract::Group::Ok)
      return res;
  }

  return res;
}

