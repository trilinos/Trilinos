#ifndef GALERI_VBRMATRICES_H
#define GALERI_VBRMATRICES_H

#include "Galeri_ConfigDefs.h"
#include "Galeri_Exception.h"
class Epetra_VbrMatrix;
class Epetra_CrsMatrix;

namespace Galeri {

Epetra_VbrMatrix* CreateVbrMatrix(const Epetra_CrsMatrix* CrsMatrix,
                                  const int NumPDEs);

} // namespace Galeri
#endif
