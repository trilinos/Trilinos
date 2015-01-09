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

#include "LOCA_Continuation_AbstractGroup.H"
#include "LOCA_SingularJacobianSolve_Default.H"

LOCA::SingularJacobianSolve::Default::Default(Teuchos::ParameterList& params)
{
  reset(params);
}

LOCA::SingularJacobianSolve::Default::Default(
              const LOCA::SingularJacobianSolve::Default& source)
{
}

LOCA::SingularJacobianSolve::Default::~Default()
{
}

LOCA::SingularJacobianSolve::Generic*
LOCA::SingularJacobianSolve::Default::clone() const
{
  return new Default(*this);
}

LOCA::SingularJacobianSolve::Generic&
LOCA::SingularJacobianSolve::Default::operator=(
              const LOCA::SingularJacobianSolve::Generic& source)
{
  return operator=(dynamic_cast<const LOCA::SingularJacobianSolve::Default&>(source));
}

LOCA::SingularJacobianSolve::Default&
LOCA::SingularJacobianSolve::Default::operator=(
              const LOCA::SingularJacobianSolve::Default& source)
{
  return *this;
}

NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::Default::reset(Teuchos::ParameterList& params)
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::Default::compute(
                Teuchos::ParameterList& params,
                LOCA::Continuation::AbstractGroup& grp,
                const NOX::Abstract::Vector& input,
                    const NOX::Abstract::Vector& approxNullVec,
                const NOX::Abstract::Vector& jacApproxNullVec,
                NOX::Abstract::Vector& result)
{
  return grp.applyJacobianInverse(params, input, result);
}

NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::Default::computeMulti(
                Teuchos::ParameterList& params,
                LOCA::Continuation::AbstractGroup& grp,
                const NOX::Abstract::Vector*const* inputs,
                const NOX::Abstract::Vector& approxNullVec,
                const NOX::Abstract::Vector& jacApproxNullVec,
                NOX::Abstract::Vector** results,
                int nVecs)
{
  return grp.applyJacobianInverseMulti(params, inputs, results, nVecs);
}
