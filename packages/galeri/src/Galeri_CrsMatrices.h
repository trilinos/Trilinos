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
Epetra_CrsMatrix* CreateCrsMatrix(string MatrixType,
                                  const Epetra_Map* Map,
                                  Teuchos::ParameterList& List);

}; // namespace Galeri

#endif
