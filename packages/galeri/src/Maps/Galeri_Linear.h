#ifndef GALERI_LINEAR_H
#define GALERI_LINEAR_H

#include "Galeri_Exception.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"

namespace Galeri {
namespace Maps {

inline
Epetra_Map* 
Linear(Epetra_Comm& Comm, int NumGlobalElements)
{
  if (NumGlobalElements <= 0)
    throw(Exception(__FILE__, __LINE__,
                    "Incorrect input parameter to Maps::Linear()",
                    "n = " + toString(NumGlobalElements)));

  return(new Epetra_Map (NumGlobalElements, 0, Comm));

} // Linear()

} // namespace Linear
} // namespace Galeri
#endif
