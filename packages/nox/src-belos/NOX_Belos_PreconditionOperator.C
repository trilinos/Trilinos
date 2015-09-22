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

#include "NOX_Belos_PreconditionOperator.H"
#include "NOX_Belos_MultiVector.H"
#include "NOX_Abstract_MultiVector.H"
#include "NOX_Parameter_List.H"

NOX::Belos::PreconditionOperator::PreconditionOperator(
                   NOX::Abstract::Group& g,
                   NOX::Parameter::List& preconditionerParameters)
  : grp(g),
    precondParams(preconditionerParameters)
{
}


NOX::Belos::PreconditionOperator::~PreconditionOperator()
{
}

::Belos::ReturnType
NOX::Belos::PreconditionOperator::Apply(
                 const ::Belos::MultiVec<double>& x,
                 ::Belos::MultiVec<double>& y,
                 ::Belos::ETrans trans) const
{
  // Cast x and y to NOX::Belos::MultiVec's
  const NOX::Belos::MultiVector& nox_belos_x =
    dynamic_cast<const NOX::Belos::MultiVector&>(x);
  NOX::Belos::MultiVector& nox_belos_y =
    dynamic_cast<NOX::Belos::MultiVector&>(y);

  // Get underlying NOX::Abstract::MultiVector's
  const NOX::Abstract::MultiVector& nox_x = nox_belos_x.getNoxMultiVector();
  NOX::Abstract::MultiVector& nox_y = nox_belos_y.getNoxMultiVector();

  // NOX return type
  NOX::Abstract::Group::ReturnType nox_status;

  bool useTranspose = false;
  if (trans == ::Belos::TRANS)
    useTranspose = true;

  nox_status = grp.applyRightPreconditioningMultiVector(useTranspose,
                            precondParams,
                            nox_x,
                            nox_y);

  return noxReturnTypeToBelos(nox_status);
}

::Belos::ReturnType
NOX::Belos::PreconditionOperator::ApplyInverse(
                 const ::Belos::MultiVec<double>& x,
                 ::Belos::MultiVec<double>& y,
                 ::Belos::ETrans trans) const
{
  return ::Belos::Undefined;
}

::Belos::ReturnType
NOX::Belos::PreconditionOperator::noxReturnTypeToBelos(
                NOX::Abstract::Group::ReturnType noxStatus) const
{
  if (noxStatus == NOX::Abstract::Group::Ok ||
      noxStatus == NOX::Abstract::Group::NotConverged)
    return ::Belos::Ok;
  else if (noxStatus == NOX::Abstract::Group::NotDefined)
    return ::Belos::Undefined;
  else
    return ::Belos::Error;
}


