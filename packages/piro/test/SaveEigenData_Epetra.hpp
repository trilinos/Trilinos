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

#ifndef SAVEEIGENDATA_EPETRA_H
#define SAVEEIGENDATA_EPETRA_H

#include "NOX_Common.H" // <string> and more
#include "LOCA_SaveEigenData_AbstractStrategy.H" // base class
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

//! Strategy for saving eigenvector/value data
/*!
 * Saves eigenvectors and corresponding eigenvalues
 * using a custom strategy.
 */
class SaveEigenData_Epetra : public LOCA::SaveEigenData::AbstractStrategy {

public:

  //! Constructor
  SaveEigenData_Epetra(Teuchos::ParameterList& locaParams);
    
  //! Destructor
  virtual ~SaveEigenData_Epetra();

  //! Save eigenvalues/eigenvectors
  virtual NOX::Abstract::Group::ReturnType 
  save(Teuchos::RCP< std::vector<double> >& evals_r,
	 Teuchos::RCP< std::vector<double> >& evals_i,
	 Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_r,
	 Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_i);

private:

  //! Private to prohibit copying
  SaveEigenData_Epetra(const SaveEigenData_Epetra&);

  //! Private to prohibit copying
  SaveEigenData_Epetra& operator = (const SaveEigenData_Epetra&);

protected:

  //! number of eigenvalues/vectors to save
  int nsave;

}; // Class SaveEigenData_Epetra
#endif
