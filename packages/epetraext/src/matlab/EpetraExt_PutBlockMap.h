//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#ifndef PUT_ROWBLOCKMAP_H
#define PUT_ROWBLOCKMAP_H
#include <Epetra_ConfigDefs.h>

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Comm.h"

// the following deal with matlab provided headers:
#include "engine.h"
#include "mex.h"
#undef printf

class Epetra_BlockMap;

namespace Matlab {
 
int CopyBlockMap(mxArray* matlabA, const Epetra_BlockMap& map)


  // Internal function
  int DoCopyBlockMap(mxArray* matlabA, int& valueCount, int length, const int * v1, const int * v2, bool doSizes);

} // namespace Matlab

#endif /* PUT_ROWBLOCKMAP_H */
