// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Galeri_ConfigDefs.h"

class Epetra_Comm;
class Epetra_Map;
class Epetra_CrsMatrix;
class Epetra_Vector;

namespace Galeri 
{
void ReadHB(const char* data_file, const Epetra_Comm& comm, 
            Epetra_Map*& map,  Epetra_CrsMatrix*& A, 
            Epetra_Vector*& x, Epetra_Vector*& b,
            Epetra_Vector*& xexact);
}
