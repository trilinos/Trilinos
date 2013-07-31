// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
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
