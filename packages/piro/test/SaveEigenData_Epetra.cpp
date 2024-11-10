// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "SaveEigenData_Epetra.hpp"
#include "NOX_Abstract_MultiVector.H"
#include "NOX_Epetra_MultiVector.H"
#include "Epetra_Vector.h"

#include <algorithm>
#include <iostream>
#include <iomanip>

SaveEigenData_Epetra::
SaveEigenData_Epetra(Teuchos::ParameterList& params)
   :  nsave(params.get("Save Eigenvectors", 0))
{
  std::cout
    << "\nSaveEigenData_Epetra: Will save up to "
    << nsave << " eigenvectors." << std::endl;
}

NOX::Abstract::Group::ReturnType
SaveEigenData_Epetra::save(
		 Teuchos::RCP< std::vector<double> >& evals_r,
		 Teuchos::RCP< std::vector<double> >& evals_i,
		 Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_r,
	         Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_i)
{
  if (nsave==0) return NOX::Abstract::Group::Ok;

  Teuchos::RCP<NOX::Epetra::MultiVector> ne_r =
    Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(evecs_r);
  Teuchos::RCP<NOX::Epetra::MultiVector> ne_i =
    Teuchos::rcp_dynamic_cast<NOX::Epetra::MultiVector>(evecs_i);
  Epetra_MultiVector& e_r = ne_r->getEpetraMultiVector();
  Epetra_MultiVector& e_i = ne_i->getEpetraMultiVector();

  const int ns = std::min(nsave, evecs_r->numVectors());

  for (int i=0; i<ns; i++) {
    if ((*evals_i)[i]==0) {
      std::cout << std::setprecision(8)
           << "Eigenvalue " << i << " with value: " << (*evals_r)[i]
           << "\n   Has Eigenvector: " << *(e_r(i)) << "\n" << std::endl;
    }
    else {
      std::cout << std::setprecision(8)
           << "Eigenvalue " << i << " with value: " << (*evals_r)[i]
           << " +  " << (*evals_i)[i] << " i \nHas Eigenvector Re, Im"
           << *(e_r(i)) << "\n" << *(e_i(i)) << std::endl;
    }
  }

  return NOX::Abstract::Group::Ok;
}
