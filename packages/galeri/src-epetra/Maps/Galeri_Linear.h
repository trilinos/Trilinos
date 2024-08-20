// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_LINEAR_H
#define GALERI_LINEAR_H

#include "Galeri_Exception.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"

namespace Galeri {
namespace Maps {

template<typename int_type>
inline
Epetra_Map* 
TLinear(Epetra_Comm& Comm, int_type NumGlobalElements)
{
  if (NumGlobalElements <= 0)
    throw(Exception(__FILE__, __LINE__,
                    "Incorrect input parameter to Maps::Linear()",
                    "n = " + toString(NumGlobalElements)));

  return(new Epetra_Map (NumGlobalElements, (int_type) 0, Comm));

} // TLinear()

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
inline
Epetra_Map* 
Linear(Epetra_Comm& Comm, int NumGlobalElements) {
  return TLinear<int>(Comm, NumGlobalElements);
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
inline
Epetra_Map* 
Linear64(Epetra_Comm& Comm, long long NumGlobalElements) {
  return TLinear<long long>(Comm, NumGlobalElements);
}
#endif

} // namespace Linear
} // namespace Galeri
#endif
