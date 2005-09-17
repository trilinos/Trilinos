#ifndef GALERI_MAPS_H
#define GALERI_MAPS_H

#include "Galeri_ConfigDefs.h"
#include "Galeri_Exception.h"

class Epetra_Comm;
class Epetra_Map;
namespace Teuchos {
  class ParameterList;
}

namespace Galeri {

// Map creation function, contained in galeri/src/Maps
Epetra_Map* CreateMap(string MapType, Epetra_Comm& Comm,
                        Teuchos::ParameterList& List);

}; // namespace Galeri

#endif
