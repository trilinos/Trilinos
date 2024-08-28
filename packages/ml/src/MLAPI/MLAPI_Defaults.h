#ifndef MLAPI_DEFAULTS_H
#define MLAPI_DEFAULTS_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*!
\file MLAPI_Defaults.h

\brief Function to set default values in a  parameter list.

\author Marzio Sala, D-INFK/ETHZ.

\date Last updated on Mar-06.
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

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
