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
