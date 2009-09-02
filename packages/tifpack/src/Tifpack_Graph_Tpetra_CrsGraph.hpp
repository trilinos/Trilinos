/*@HEADER
// ***********************************************************************
//
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
*/

#ifndef TIFPACK_EPETRA_CRSGRAPH_HPP
#define TIFPACK_EPETRA_CRSGRAPH_HPP
#include "Tifpack_ConfigDefs.hpp"
#include "Tifpack_Graph.hpp"
#include "Teuchos_RefCountPtr.hpp"

class Tpetra_Comm;
class Tpetra_CrsGraph;

//! Tifpack_Graph_Tpetra_CrsGraph: a class to define Tifpack_Graph as a light-weight conversion of Tpetra_CrsGraph's.

/*! 
Class Tifpack_Graph_Tpetra_CrsGraph enables the construction of an
Tifpack_Graph based on the input Tpetra_CrsGraph. Note that data are
not copied to \e this object; instead, wrappers are furnished.

\date Set-04.
*/
class Tifpack_Graph_Tpetra_CrsGraph : public Tifpack_Graph {

public:
    
  //! Constructor.
  Tifpack_Graph_Tpetra_CrsGraph(const Teuchos::RefCountPtr<const Tpetra_CrsGraph>& CrsGraph);

  //! Destructor.
  virtual ~Tifpack_Graph_Tpetra_CrsGraph() {};

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

  //! Returns the number of global rows.
  int NumGlobalRows() const
  {
    return(NumGlobalRows_);
  }

  //! Returns the number of global columns.
  int NumGlobalCols() const
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

  //! Returns the global row ID of input local row.
  int GRID(int) const;

  //! Returns the global column ID of input local column.
  int GCID(int) const;
  
  //! Returns the local row ID of input global row.
  int LRID(int) const;

  //! Returns the local column ID of input global column.
  int LCID(int) const;

  //! Extracts a copy of input local row.
  int ExtractMyRowCopy(int GlobalRow, int LenOfIndices, 
		       int &NumIndices, int *Indices) const;

  //! Returns the communicator object of the graph.
  const Tpetra_Comm& Comm() const;
  
  //! Prints basic information about the graph object.
  virtual ostream& Print(std::ostream& os) const;

private:

  //! Number of local rows.
  int NumMyRows_;
  //! Number of local columns.
  int NumMyCols_;
  //! Number of global rows.
  int NumGlobalRows_;
  //! Number of global columns.
  int NumGlobalCols_;
  //! Maximum number of indices per row.
  int MaxNumIndices_;
  //! Pointer to the wrapped Tpetra_CrsGraph.
  Teuchos::RefCountPtr<const Tpetra_CrsGraph> CrsGraph_;
};

#endif
