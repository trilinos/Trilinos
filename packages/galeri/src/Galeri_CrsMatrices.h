#ifndef GALERI_MATRICES_H
#define GALERI_MATRICES_H

#include "Galeri_ConfigDefs.h"
#include "Galeri_Exception.h"

class Epetra_Map;
class Epetra_CrsMatrix;
namespace Teuchos {
  class ParameterList;
}

namespace Galeri {

// Matrix creation function, contained in galeri/src/Maps
Epetra_CrsMatrix* CreateCrsMatrix(string MatrixType,
                                  const Epetra_Map* Map,
                                  Teuchos::ParameterList& List);

}; // namespace Galeri

#endif
