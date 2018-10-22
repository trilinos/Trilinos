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

#include "NOX_GlobalData.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_MeritFunction_SumOfSquares.H"
#include "NOX_SolverStats.hpp"

NOX::GlobalData::
GlobalData(const Teuchos::RCP<Teuchos::ParameterList>& noxParams)
{ this->initialize(noxParams); }

NOX::GlobalData::
GlobalData(const Teuchos::RCP<NOX::Utils>& utils,
           const Teuchos::RCP<NOX::MeritFunction::Generic>& mf) :
  utilsPtr(utils),
  meritFunctionPtr(mf)
{
  solverStatsPtr = Teuchos::rcp(new NOX::SolverStats);
}

NOX::GlobalData::~GlobalData() {}

void NOX::GlobalData::
initialize(const Teuchos::RCP<Teuchos::ParameterList>& noxParams)
{
  paramListPtr = noxParams;
  utilsPtr = Teuchos::rcp(new NOX::Utils(noxParams->sublist("Printing")));
  if (is_null(solverStatsPtr))
    solverStatsPtr = Teuchos::rcp(new NOX::SolverStats);

  Teuchos::ParameterList& so = noxParams->sublist("Solver Options");

  if (so.isType< Teuchos::RCP<NOX::MeritFunction::Generic> >("User Defined Merit Function")) {
    meritFunctionPtr = so.get< Teuchos::RCP<NOX::MeritFunction::Generic> >("User Defined Merit Function");
  }

  // PL validator sets a default null RCP. If it is null, allocate
  // a concrete default.
  if (is_null(meritFunctionPtr))
    meritFunctionPtr = Teuchos::rcp(new NOX::MeritFunction::SumOfSquares(utilsPtr)); 
}

Teuchos::RCP<NOX::Utils> NOX::GlobalData::getUtils() const
{ return utilsPtr; }

Teuchos::RCP<NOX::MeritFunction::Generic>
NOX::GlobalData::getMeritFunction() const
{ return meritFunctionPtr; }

Teuchos::RCP<Teuchos::ParameterList>
NOX::GlobalData::getNoxParameterList() const
{ return paramListPtr; }

Teuchos::RCP<const NOX::SolverStats>
NOX::GlobalData::getSolverStatistics() const
{ return solverStatsPtr; }

Teuchos::RCP<NOX::SolverStats>
NOX::GlobalData::getNonConstSolverStatistics()
{ return solverStatsPtr; }
