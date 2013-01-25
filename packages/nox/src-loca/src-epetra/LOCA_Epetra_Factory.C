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

#include "Teuchos_ParameterList.hpp"

#include "LOCA_Epetra_Factory.H"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_BorderedSolver_EpetraHouseholder.H"
#include "LOCA_BorderedSolver_EpetraAugmented.H"
#ifdef HAVE_NOX_EPETRAEXT
#ifdef HAVE_MPI
#include "LOCA_Epetra_AnasaziOperator_Floquet.H"
#endif
#endif

LOCA::Epetra::Factory::Factory() :
  globalData()
{
}

LOCA::Epetra::Factory::~Factory()
{
}

void
LOCA::Epetra::Factory::init(
		   const Teuchos::RCP<LOCA::GlobalData>& global_data)
{
  globalData = global_data;
}

bool
LOCA::Epetra::Factory::createBorderedSolverStrategy(
       const std::string& strategyName,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
       Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy>& strategy)
{
  // Instantiate Householder strategy if requested
  if (strategyName == "Householder") {
    strategy = 
      Teuchos::rcp(new LOCA::BorderedSolver::EpetraHouseholder(globalData,
							       topParams,
							       solverParams));
    return true;
  }
  // Instantiate augmented strategy if requested
  else if (strategyName == "Augmented") {
    strategy = 
      Teuchos::rcp(new LOCA::BorderedSolver::EpetraAugmented(globalData,
							     topParams,
							     solverParams));
    return true;
  }
  else
    return false;
}

bool
LOCA::Epetra::Factory::createAnasaziOperatorStrategy(
       const std::string& strategyName,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& eigenParams,
       const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
       const Teuchos::RCP<NOX::Abstract::Group>& grp,
       Teuchos::RCP<LOCA::AnasaziOperator::AbstractStrategy>& strategy)
{
#ifdef HAVE_NOX_EPETRAEXT
#ifdef HAVE_MPI
if (strategyName == "Floquet") {

    strategy = 
      Teuchos::rcp(new LOCA::Epetra::AnasaziOperator::Floquet(globalData,
                                topParams, eigenParams, solverParams, grp));
    return true;
  }
 else
#endif
#endif
    return false;
}


