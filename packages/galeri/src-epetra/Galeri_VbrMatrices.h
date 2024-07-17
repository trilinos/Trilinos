// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_VBRMATRICES_H
#define GALERI_VBRMATRICES_H

#include "Galeri_ConfigDefs.h"
#include "Galeri_Exception.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"

namespace Galeri {

Epetra_VbrMatrix* CreateVbrMatrix(const Epetra_CrsMatrix* CrsMatrix,
                                  const int NumPDEs);

} // namespace Galeri
#endif
