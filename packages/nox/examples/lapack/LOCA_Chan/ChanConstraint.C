// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ChanConstraint.H"
#include "LOCA_Parameter_Vector.H"

ChanConstraint::ChanConstraint(int N, const LOCA::ParameterVector& pVec) :
  n(N),
  constraints(1,1),
  isValidConstraints(false),
  p(pVec),
  x()
{
  constraints.putScalar(0.0);
  NOX::LAPACK::Vector xx(n);
  x = xx.createMultiVector(1);
}

ChanConstraint::ChanConstraint(const ChanConstraint& source,
                   NOX::CopyType type) :
  n(source.n),
  constraints(source.constraints),
  isValidConstraints(false),
  p(source.p),
  x(source.x->clone(type))
{
  if (source.isValidConstraints && type == NOX::DeepCopy)
    isValidConstraints = true;
}

ChanConstraint::~ChanConstraint()
{
}

void
ChanConstraint::copy(const LOCA::MultiContinuation::ConstraintInterface& src)
{
  const ChanConstraint& source = dynamic_cast<const ChanConstraint&>(src);

  if (this != &source) {
    n = source.n;
    constraints = source.constraints;
    isValidConstraints = source.isValidConstraints;
    p = source.p;
    *x = *source.x;
  }
}

Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
ChanConstraint::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new ChanConstraint(*this, type));
}

int
ChanConstraint::numConstraints() const
{
  return 1;
}

void
ChanConstraint::setX(const NOX::Abstract::Vector& y)
{
  (*x)[0] = y;
  isValidConstraints = false;
}

void
ChanConstraint::setParam(int paramID, double val)
{
  p[paramID] = val;
  isValidConstraints = false;
}

void
ChanConstraint::setParams(
             const std::vector<int>& paramIDs,
             const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  for (unsigned int i=0; i<paramIDs.size(); i++)
    p[paramIDs[i]] = vals(i,0);
  isValidConstraints = false;
}

NOX::Abstract::Group::ReturnType
ChanConstraint::computeConstraints()
{
  constraints(0,0) = 0.5 * (*x)[0].innerProduct((*x)[0]) / n - p.getValue("gamma");
  isValidConstraints = true;

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
ChanConstraint::computeDX()
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
ChanConstraint::computeDP(const std::vector<int>& paramIDs,
              NOX::Abstract::MultiVector::DenseMatrix& dgdp,
              bool isValidG)
{
  if (!isValidG) {
    dgdp(0,0) = constraints(0,0);
  }

  for (unsigned int i=0; i<paramIDs.size(); i++) {
    if (p.getLabel(paramIDs[i]) == "gamma") {
      dgdp(0,i+1) = -1.0;
    }
    else {
      dgdp(0,i+1) = 0.0;
    }
  }

  return NOX::Abstract::Group::Ok;
}

bool
ChanConstraint::isConstraints() const
{
  return isValidConstraints;
}

bool
ChanConstraint::isDX() const
{
  return true;
}

const NOX::Abstract::MultiVector::DenseMatrix&
ChanConstraint::getConstraints() const
{
  return constraints;
}

NOX::Abstract::Group::ReturnType
ChanConstraint::multiplyDX(
            double alpha,
            const NOX::Abstract::MultiVector& input_x,
            NOX::Abstract::MultiVector::DenseMatrix& result_p) const
{
  input_x.multiply(alpha/n, *x, result_p);
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
ChanConstraint::addDX(Teuchos::ETransp transb,
              double alpha,
              const NOX::Abstract::MultiVector::DenseMatrix& b,
              double beta,
              NOX::Abstract::MultiVector& result_x) const
{
  result_x.update(transb, alpha/n, *x, b, beta);
  return NOX::Abstract::Group::Ok;
}

bool
ChanConstraint::isDXZero() const
{
  return false;
}
