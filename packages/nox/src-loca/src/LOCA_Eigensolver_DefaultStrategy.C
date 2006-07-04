// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov) or Eric Phipps
// (etphipp@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#include "LOCA_Eigensolver_DefaultStrategy.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"


LOCA::Eigensolver::DefaultStrategy::DefaultStrategy(
	const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RefCountPtr<Teuchos::ParameterList>& eigenParams) :
  globalData(global_data)
{
}

LOCA::Eigensolver::DefaultStrategy::~DefaultStrategy()
{
}

NOX::Abstract::Group::ReturnType
LOCA::Eigensolver::DefaultStrategy::computeEigenvalues(
		 NOX::Abstract::Group& group,
		 Teuchos::RefCountPtr< std::vector<double> >& evals_r,
		 Teuchos::RefCountPtr< std::vector<double> >& evals_i,
		 Teuchos::RefCountPtr< NOX::Abstract::MultiVector >& evecs_r,
	         Teuchos::RefCountPtr< NOX::Abstract::MultiVector >& evecs_i)
{
  // Print a warning that this eigensolver strategy doesn't do anything
  globalData->locaErrorCheck->printWarning(
   "LOCA::Eigensolver::DefaultStrategy::computeEigenvalues()",
   "\nThe default Eigensolver strategy does not compute eigenvalues.\nSet the \"Method\" parameter of the \"Eigensolver\" sublist to chose an \neigensolver method.");
  return NOX::Abstract::Group::Ok;
}
