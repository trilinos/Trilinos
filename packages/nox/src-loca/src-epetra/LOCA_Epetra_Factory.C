// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
       const string& strategyName,
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
       const string& strategyName,
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


