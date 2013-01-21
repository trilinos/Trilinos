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
