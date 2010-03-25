// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include "SaveEigenData_Epetra.hpp"
#include "NOX_Abstract_MultiVector.H"
#include "NOX_Epetra_MultiVector.H"
#include "Epetra_Vector.h"
#include <string>

SaveEigenData_Epetra::
SaveEigenData_Epetra(Teuchos::ParameterList& locaParams)
   :  nsave(0)
{
  bool doEig = locaParams.sublist("Stepper").get("Compute Eigenvalues", false);
  if (doEig) nsave = locaParams.sublist("Stepper").
              sublist("Eigensolver").get("Save Eigenvectors",0);

  cout << "\nSaveEigenData_Epetra: Will save up to " 
       << nsave << " eigenvectors." << endl;
}

SaveEigenData_Epetra::~SaveEigenData_Epetra()
{
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

  int ns = nsave;
  if (ns > evecs_r->numVectors())
    ns = evecs_r->numVectors();

  for (int i=0; i<ns; i++) {
    if ((*evals_i)[i]==0) {
      cout << setprecision(8) 
           << "Eigenvalue " << i << " with value: " << (*evals_r)[i] 
           << "\n   Has Eigenvector: " << *(e_r(i)) << "\n" << endl;
    }
    else {
      cout << setprecision(8) 
           << "Eigenvalue " << i << " with value: " << (*evals_r)[i] 
           << " +  " << (*evals_i)[i] << " i \nHas Eigenvector Re, Im" 
           << *(e_r(i)) << "\n" << *(e_i(i)) << endl;
    }
  }

  return NOX::Abstract::Group::Ok;
}
