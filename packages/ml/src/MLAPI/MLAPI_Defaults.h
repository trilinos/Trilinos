#ifndef MLAPI_DEFAULTS_H
#define MLAPI_DEFAULTS_H

#include "MLAPI_Error.h"
#include "Teuchos_ParameterList.hpp"

namespace MLAPI 
{

// ====================================================================== 
//! Sets default values in input \c List.
// ====================================================================== 

void SetDefaults(Teuchos::ParameterList& List);

} // namespace MLAPI

#endif // MLAPI_DEFAULTS_H
