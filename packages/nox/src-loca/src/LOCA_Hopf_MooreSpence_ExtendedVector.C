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

#include "LOCA_Hopf_MooreSpence_ExtendedVector.H"  // Class definition
#include "LOCA_Hopf_MooreSpence_ExtendedMultiVector.H"

LOCA::Hopf::MooreSpence::ExtendedVector::ExtendedVector(
            const Teuchos::RCP<LOCA::GlobalData>& global_data,
            const NOX::Abstract::Vector& xVec,
            const NOX::Abstract::Vector& realEigenVec,
            const NOX::Abstract::Vector& imagEigenVec,
            double frequency,
            double bifParam) :
  LOCA::Extended::Vector(global_data,3,2)
{
  setVector(0, xVec);
  setVector(1, realEigenVec);
  setVector(2, imagEigenVec);
  setScalar(0, frequency);
  setScalar(1, bifParam);
}

LOCA::Hopf::MooreSpence::ExtendedVector::ExtendedVector(
                    const LOCA::Hopf::MooreSpence::ExtendedVector& source,
            NOX::CopyType type) :
  LOCA::Extended::Vector(source, type)
{
}


LOCA::Hopf::MooreSpence::ExtendedVector::~ExtendedVector()
{
}

NOX::Abstract::Vector&
LOCA::Hopf::MooreSpence::ExtendedVector::operator=(
                          const NOX::Abstract::Vector& y)
{
  return operator=(dynamic_cast<const LOCA::Hopf::MooreSpence::ExtendedVector&>(y));
}

LOCA::Extended::Vector&
LOCA::Hopf::MooreSpence::ExtendedVector::operator=(
                          const LOCA::Extended::Vector& y)
{
  return operator=(dynamic_cast<const LOCA::Hopf::MooreSpence::ExtendedVector&>(y));
}

LOCA::Hopf::MooreSpence::ExtendedVector&
LOCA::Hopf::MooreSpence::ExtendedVector::operator=(
                         const LOCA::Hopf::MooreSpence::ExtendedVector& y)
{
  LOCA::Extended::Vector::operator=(y);
  return *this;
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::Hopf::MooreSpence::ExtendedVector::clone(NOX::CopyType type) const
{
  return
    Teuchos::rcp(new LOCA::Hopf::MooreSpence::ExtendedVector(*this, type));
}

void
LOCA::Hopf::MooreSpence::ExtendedVector::setVec(
                   const NOX::Abstract::Vector& xVec,
                   const NOX::Abstract::Vector& realEigenVec,
                   const NOX::Abstract::Vector& imagEigenVec,
                   double frequency,
                   double bifPar)
{
  setVector(0, xVec);
  setVector(1, realEigenVec);
  setVector(2, imagEigenVec);
  setScalar(0, frequency);
  setScalar(1, bifPar);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::Hopf::MooreSpence::ExtendedVector::getXVec() const
{
  return getVector(0);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::Hopf::MooreSpence::ExtendedVector::getRealEigenVec() const
{
  return getVector(1);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::Hopf::MooreSpence::ExtendedVector::getImagEigenVec() const
{
  return getVector(2);
}

double
LOCA::Hopf::MooreSpence::ExtendedVector::getFrequency() const
{
  return getScalar(0);
}

double
LOCA::Hopf::MooreSpence::ExtendedVector::getBifParam() const
{
  return getScalar(1);
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::Hopf::MooreSpence::ExtendedVector::getXVec()
{
  return getVector(0);
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::Hopf::MooreSpence::ExtendedVector::getRealEigenVec()
{
  return getVector(1);
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::Hopf::MooreSpence::ExtendedVector::getImagEigenVec()
{
  return getVector(2);
}

double&
LOCA::Hopf::MooreSpence::ExtendedVector::getFrequency()
{
  return getScalar(0);
}

double&
LOCA::Hopf::MooreSpence::ExtendedVector::getBifParam()
{
  return getScalar(1);
}

LOCA::Hopf::MooreSpence::ExtendedVector::ExtendedVector(
          const Teuchos::RCP<LOCA::GlobalData>& global_data) :
  LOCA::Extended::Vector(global_data,3,2)
{
}

Teuchos::RCP<LOCA::Extended::MultiVector>
LOCA::Hopf::MooreSpence::ExtendedVector::generateMultiVector(
                            int nColumns,
                            int /* nVectorRows */,
                            int /* nScalarRows */) const
{
  return
    Teuchos::rcp(new LOCA::Hopf::MooreSpence::ExtendedMultiVector(globalData,
                                  nColumns));
}
