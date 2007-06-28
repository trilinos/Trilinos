#ifndef ML_VALIDATEPARAMETERS_H
#define ML_VALIDATEPARAMETERS_H

#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) 
#include "Teuchos_ParameterList.hpp"

namespace ML_Epetra
{
  //! Builds a list of "valid" parameters for parameter validation
  Teuchos::ParameterList * GetValidMLPParameters();
  
  //! Validates the parameters of inList (warning: level-specific parameters will not be validated).
  bool ValidateMLPParameters(const Teuchos::ParameterList &inList);  
}
#endif
#endif
