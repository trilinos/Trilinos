// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_MATRICES_H
#define GALERI_MATRICES_H

#include "Galeri_ConfigDefs.h"
#include "Galeri_Exception.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_ParameterList.hpp"

namespace Galeri {

// Matrix creation function, contained in galeri/src/Maps
Epetra_CrsMatrix* CreateCrsMatrix(std::string MatrixType,
                                  const Epetra_Map* Map,
                                  Teuchos::ParameterList& List);

} // namespace Galeri

#endif
