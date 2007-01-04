#ifndef IFPACK_VALIDPARAMETERS_H
#define IFPACK_VALIDPARAMETERS_H

#include "Ifpack_ConfigDefs.h"
#include "Teuchos_ParameterList.hpp"

//! Returns a list which contains all the parameters possibly used by IFPACK.
Teuchos::ParameterList Ifpack_GetValidParameters();

#endif
