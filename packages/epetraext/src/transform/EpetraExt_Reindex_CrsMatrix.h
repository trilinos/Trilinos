//@HEADER
// ************************************************************************
//
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
                                                                                                    
#ifndef EpetraExt_CRSMATRIX_REINDEX_H
#define EpetraExt_CRSMATRIX_REINDEX_H

#include <EpetraExt_Transform.h>

class Epetra_CrsMatrix;
class Epetra_Map;

namespace EpetraExt {

///
/** Given an Epetra_CrsMatrix, a "reindexed" version is returned based on
 *  the new row map.  The row map must be conformal to the original.  The
 *  Matrix data will be shared by the new Matrix using the new indexing
 */
class CrsMatrix_Reindex : public ViewTransform<Epetra_CrsMatrix> {

  const Epetra_Map & NewRowMap_;
  Epetra_Map * NewColMap_;

 public:

  ///
  /** Destructor
   */
  ~CrsMatrix_Reindex();

  ///
  /** Constructor
   */
  CrsMatrix_Reindex( const Epetra_Map & new_row_map )
  : NewRowMap_(new_row_map),
    NewColMap_(0)
  {}

  ///
  /** Constructs "reindexed" Matrix
   */
  NewTypeRef operator()( OriginalTypeRef orig );

};

} //namespace EpetraExt

#endif //EpetraExt_CRSMATRIX_REINDEX_H
