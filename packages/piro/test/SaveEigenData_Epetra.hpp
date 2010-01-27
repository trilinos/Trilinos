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
