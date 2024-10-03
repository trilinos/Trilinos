// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCALinearConstraint.H"
#include "LOCA_Parameter_Vector.H"

LinearConstraint::LinearConstraint(int NumConstraints,
                   const LOCA::ParameterVector& pVec,
                   NOX::Abstract::Vector& cloneVec) :
  m(NumConstraints),
  constraints(m,1),
  isValidConstraints(false),
  x(),
  dgdx(),
  isZeroDgDx(false),
  p(pVec),
  pvec(Teuchos::View, p.getDoubleArrayPointer(), p.length(), p.length(), 1),
  dgdp(NumConstraints, pVec.length())
{
  constraints.putScalar(0.0);
  x = cloneVec.createMultiVector(1);
  dgdx = cloneVec.createMultiVector(m);
  dgdp.putScalar(0.0);
}

LinearConstraint::LinearConstraint(const LinearConstraint& source,
                   NOX::CopyType type) :
  m(source.m),
  constraints(source.constraints),
  isValidConstraints(false),
  x(source.x->clone(type)),
  dgdx(source.dgdx->clone(type)),
  isZeroDgDx(source.isZeroDgDx),
  p(source.p),
  pvec(source.pvec),
  dgdp(source.dgdp)
{
  if (source.isValidConstraints && type == NOX::DeepCopy)
    isValidConstraints = true;
}

LinearConstraint::~LinearConstraint()
{
}

void
LinearConstraint::copy(const LOCA::MultiContinuation::ConstraintInterface& src)
{
  const LinearConstraint& source = dynamic_cast<const LinearConstraint&>(src);

  if (this != &source) {
    m = source.m;
    constraints = source.constraints;
    isValidConstraints = source.isValidConstraints;
    *x = *source.x;
    *dgdx = *source.dgdx;
    isZeroDgDx = source.isZeroDgDx;
    p = source.p;
    pvec.assign(source.pvec);
    dgdp = source.dgdp;
  }
}

Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
LinearConstraint::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new LinearConstraint(*this, type));
}

int
LinearConstraint::numConstraints() const
{
  return m;
}

void
LinearConstraint::setX(const NOX::Abstract::Vector& y)
{
  (*x)[0] = y;
  isValidConstraints = false;
}

void
LinearConstraint::setParam(int paramID, double val)
{
  p[paramID] = val;
}

void
LinearConstraint::setParams(
             const std::vector<int>& paramIDs,
             const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  for (unsigned int i=0; i<paramIDs.size(); i++)
    p[paramIDs[i]] = vals(i,0);
}

NOX::Abstract::Group::ReturnType
LinearConstraint::computeConstraints()
{
  // compute dg/dx*x + dg/dp*p
  x->multiply(1.0, *dgdx, constraints);
  if (pvec.numRows() > 0)
    constraints.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, dgdp,
             pvec, 1.0);
  isValidConstraints = true;

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LinearConstraint::computeDX()
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LinearConstraint::computeDP(const std::vector<int>& paramIDs,
              NOX::Abstract::MultiVector::DenseMatrix& dp,
              bool isValidG)
{
  if (!isValidG) {
    for (int i=0; i<m; i++)
      dp(i,0) = constraints(i,0);
  }

  for (unsigned int i=0; i<paramIDs.size(); i++) {
    for (int j=0; j<m; j++)
      dp(j,i+1) = dgdp(j,paramIDs[i]);
  }

  return NOX::Abstract::Group::Ok;
}

bool
LinearConstraint::isConstraints() const
{
  return isValidConstraints;
}

bool
LinearConstraint::isDX() const
{
  return true;
}

const NOX::Abstract::MultiVector::DenseMatrix&
LinearConstraint::getConstraints() const
{
  return constraints;
}

const NOX::Abstract::MultiVector*
LinearConstraint::getDX() const
{
  return dgdx.get();
}

bool
LinearConstraint::isDXZero() const
{
  return isZeroDgDx;
}

void
LinearConstraint::setDgDx(const NOX::Abstract::MultiVector& A)
{
  *dgdx = A;
  isValidConstraints = false;
}

void
LinearConstraint::setDgDp(const NOX::Abstract::MultiVector::DenseMatrix& A)
{
  dgdp.assign(A);
  isValidConstraints = false;
}

void
LinearConstraint::setIsZeroDX(bool flag)
{
  isZeroDgDx = flag;
}
