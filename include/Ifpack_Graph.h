/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK_GRAPH_H
#define IFPACK_GRAPH_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif
class Epetra_Comm;

//! Ifpack_Graph: a pure virtual class that defines graphs for IFPACK.
/*!
Class Ifpack_Graph defines the abstract interface to use graphs in
IFPACK. This class contains all the functions that are required by
IFPACK classes.

\author Marzio Sala, SNL 9214.

\date Last modified on Nov-04.

*/

#include "Epetra_ConfigDefs.h"

class Ifpack_Graph {

public:

  //! Destructor.
  virtual ~Ifpack_Graph() {};

  //! Returns the number of local rows.
  virtual int NumMyRows() const = 0;

  //! Returns the number of local columns.
  virtual int NumMyCols() const = 0;

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  //! Returns the number of global rows.
  virtual int NumGlobalRows() const = 0;

  //! Returns the number of global columns.
  virtual int NumGlobalCols() const = 0;
#endif

  virtual long long NumGlobalRows64() const = 0;

  virtual long long NumGlobalCols64() const = 0;

  //! Returns the maximun number of entries for row.
  virtual int MaxMyNumEntries() const = 0;

  //! Returns the number of local nonzero entries.
  virtual int NumMyNonzeros() const = 0;

  //! Returns \c true is graph is filled.
  virtual bool Filled() const = 0;

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  //! Returns the global row ID of input local row.
  virtual int GRID(int) const = 0;

  //! Returns the global column ID of input local column.
  virtual int GCID(int) const = 0;
#endif

  virtual long long GRID64(int) const = 0;

  //! Returns the global column ID of input local column.
  virtual long long GCID64(int) const = 0;

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  //! Returns the local row ID of input global row.
  virtual int LRID(int) const = 0;

  //! Returns the local column ID of input global column.
  virtual int LCID(int) const = 0;
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  //! Returns the local row ID of input global row.
  virtual int LRID(long long) const = 0;

  //! Returns the local column ID of input global column.
  virtual int LCID(long long) const = 0;
#endif

  //! Extracts a copy of input local row.
  virtual int ExtractMyRowCopy(int MyRow, int LenOfIndices,
                               int &NumIndices, int *Indices) const = 0;

  //! Returns the communicator object of the graph.
  virtual const Epetra_Comm& Comm() const = 0;

  //! Prints basic information about the graph object.
  virtual std::ostream& Print(std::ostream& os) const = 0;

};

inline std::ostream& operator<<(std::ostream& os, const Ifpack_Graph& obj)
{
  return(obj.Print(os));
}

#endif // iFPACK_GRAPH_H
