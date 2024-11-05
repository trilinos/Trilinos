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

#ifndef IFPACK_EPETRA_CRSGRAPH_H
#define IFPACK_EPETRA_CRSGRAPH_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif
#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Graph.h"
#include "Epetra_CrsGraph.h"
#include "Teuchos_RefCountPtr.hpp"

class Epetra_Comm;
class Epetra_CrsGraph;

//! Ifpack_Graph_Epetra_CrsGraph: a class to define Ifpack_Graph as a light-weight conversion of Epetra_CrsGraph's.

/*!
Class Ifpack_Graph_Epetra_CrsGraph enables the construction of an
Ifpack_Graph based on the input Epetra_CrsGraph. Note that data are
not copied to \e this object; instead, wrappers are furnished.

\date Set-04.
*/
class Ifpack_Graph_Epetra_CrsGraph : public Ifpack_Graph {

public:

  //! Constructor.
  Ifpack_Graph_Epetra_CrsGraph(const Teuchos::RefCountPtr<const Epetra_CrsGraph>& CrsGraph);

  //! Destructor.
  virtual ~Ifpack_Graph_Epetra_CrsGraph() {};

  //! Returns the number of local rows.
  int NumMyRows() const
  {
    return(NumMyRows_);
  }

  //! Returns the number of local columns.
  int NumMyCols() const
  {
    return(NumMyCols_);
  }

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  //! Returns the number of global rows.
  int NumGlobalRows() const
  {
    if(CrsGraph_->RowMap().GlobalIndicesInt())
      return (int) (NumGlobalRows_);
        else
      throw "Ifpack_Graph_Epetra_CrsGraph::NumGlobalRows: GlobalIndices not int.";
  }
#endif
  long long NumGlobalRows64() const
  {
    return(NumGlobalRows_);
  }

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  //! Returns the number of global columns.
  int NumGlobalCols() const
  {
    if(CrsGraph_->ColMap().GlobalIndicesInt())
      return (int) (NumGlobalCols_);
        else
      throw "Ifpack_Graph_Epetra_CrsGraph::NumGlobalCols: GlobalIndices not int.";
  }
#endif
  long long NumGlobalCols64() const
  {
    return(NumGlobalCols_);
  }

  //! Returns the maximun number of entries for row.
  int MaxMyNumEntries() const
  {
    return(MaxNumIndices_);
  }

  //! Returns the number of local nonzero entries.
  int NumMyNonzeros() const;

  //! Returns \c true is graph is filled.
  bool Filled() const;

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  //! Returns the global row ID of input local row.
  int GRID(int) const;
#endif
  long long GRID64(int) const;

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  //! Returns the global column ID of input local column.
  int GCID(int) const;
#endif
  long long GCID64(int) const;

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  //! Returns the local row ID of input global row.
  int LRID(int) const;

  //! Returns the local column ID of input global column.
  int LCID(int) const;
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  //! Returns the local row ID of input global row.
  int LRID(long long) const;

  //! Returns the local column ID of input global column.
  int LCID(long long) const;
#endif

  //! Extracts a copy of input local row.
  int ExtractMyRowCopy(int GlobalRow, int LenOfIndices,
                       int &NumIndices, int *Indices) const;

  //! Returns the communicator object of the graph.
  const Epetra_Comm& Comm() const;

  //! Prints basic information about the graph object.
  virtual std::ostream& Print(std::ostream& os) const;

private:

  //! Number of local rows.
  int NumMyRows_;
  //! Number of local columns.
  int NumMyCols_;
  //! Number of global rows.
  long long NumGlobalRows_;
  //! Number of global columns.
  long long NumGlobalCols_;
  //! Maximum number of indices per row.
  int MaxNumIndices_;
  //! Pointer to the wrapped Epetra_CrsGraph.
  Teuchos::RefCountPtr<const Epetra_CrsGraph> CrsGraph_;
};

#endif
