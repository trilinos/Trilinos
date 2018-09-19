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

#include "NOX_Abstract_Group.H"

#include "NOX_Abstract_MultiVector.H"
#include "Teuchos_ParameterList.hpp"

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::computeJacobian()
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::computeGradient()
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::computeNewton(Teuchos::ParameterList& /* params */)
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyJacobian(const NOX::Abstract::Vector& /* input */,
                    NOX::Abstract::Vector& /* result */) const
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyJacobianTranspose(const NOX::Abstract::Vector& /* input */,
                         NOX::Abstract::Vector& /* result */) const
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyJacobianInverse(Teuchos::ParameterList& /* params */,
                       const NOX::Abstract::Vector& /* input */,
                       NOX::Abstract::Vector& /* result */) const
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyRightPreconditioning(bool /* useTranspose */,
                        Teuchos::ParameterList& /* params */,
                        const NOX::Abstract::Vector& /* input */,
                        NOX::Abstract::Vector& /* result */
                        ) const
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyJacobianMultiVector(
                    const NOX::Abstract::MultiVector& input,
                    NOX::Abstract::MultiVector& result) const
{
  NOX::Abstract::Group::ReturnType status, finalStatus;
  finalStatus = NOX::Abstract::Group::Ok;

  for (int i=0; i<input.numVectors(); i++) {
    status = applyJacobian(input[i], result[i]);

    if (status == NotDefined || status == BadDependency)
      return status;
    else if (status == Failed)
      finalStatus = Failed;
    else if (status == NotConverged && finalStatus != Failed)
      finalStatus = NotConverged;
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyJacobianTransposeMultiVector(
                    const NOX::Abstract::MultiVector& input,
                    NOX::Abstract::MultiVector& result) const
{
  NOX::Abstract::Group::ReturnType status, finalStatus;
  finalStatus = NOX::Abstract::Group::Ok;

  for (int i=0; i<input.numVectors(); i++) {
    status = applyJacobianTranspose(input[i], result[i]);

    if (status == NotDefined || status == BadDependency)
      return status;
    else if (status == Failed)
      finalStatus = Failed;
    else if (status == NotConverged && finalStatus != Failed)
      finalStatus = NotConverged;
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyJacobianInverseMultiVector(
                                    Teuchos::ParameterList& params,
                    const NOX::Abstract::MultiVector& input,
                    NOX::Abstract::MultiVector& result) const
{
  NOX::Abstract::Group::ReturnType status, finalStatus;
  finalStatus = NOX::Abstract::Group::Ok;

  for (int i=0; i<input.numVectors(); i++) {
    status = applyJacobianInverse(params, input[i], result[i]);

    if (status == NotDefined || status == BadDependency)
      return status;
    else if (status == Failed)
      finalStatus = Failed;
    else if (status == NotConverged && finalStatus != Failed)
      finalStatus = NotConverged;
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyRightPreconditioningMultiVector(
                   bool useTranspose,
                   Teuchos::ParameterList& params,
                   const NOX::Abstract::MultiVector& input,
                   NOX::Abstract::MultiVector& result) const
{
  NOX::Abstract::Group::ReturnType status, finalStatus;
  finalStatus = NOX::Abstract::Group::Ok;

  for (int i=0; i<input.numVectors(); i++) {
    status = applyRightPreconditioning(useTranspose, params, input[i],
                       result[i]);

    if (status == NotDefined || status == BadDependency)
      return status;
    else if (status == Failed)
      finalStatus = Failed;
    else if (status == NotConverged && finalStatus != Failed)
      finalStatus = NotConverged;
  }

  return finalStatus;
}

bool NOX::Abstract::Group::isJacobian() const
{
  return false;
}

bool NOX::Abstract::Group::isGradient() const
{
  return false;
}

bool NOX::Abstract::Group::isNewton() const
{
  return false;
}

void NOX::Abstract::Group::logLastLinearSolveStats(NOX::SolverStats& ) const {}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::getNormLastLinearSolveResidual(double& /* residual */) const
{
  return NOX::Abstract::Group::NotDefined;
}

