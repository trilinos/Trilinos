// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_MAPS_H
#define GALERI_MAPS_H

#include "Galeri_ConfigDefs.h"
#include "Galeri_Exception.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Teuchos_ParameterList.hpp"

namespace Galeri {

// Map creation functions, contained in galeri/src/Maps

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
Epetra_Map* CreateMap(std::string MapType, Epetra_Comm& Comm,
                        Teuchos::ParameterList& List);
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
Epetra_Map* CreateMap64(std::string MapType, Epetra_Comm& Comm,
                        Teuchos::ParameterList& List);
#endif

} // namespace Galeri

#endif
