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

#ifndef EDT_CRSGRAPH_VIEW_H
#define EDT_CRSGRAPH_VIEW_H

#include <Epetra_Transform.h>

class Epetra_CrsGraph;
class Epetra_BlockMap;

namespace EpetraExt {

///
/** Generates a "view" object of a contiguous subset of a original Epetra_CrsGraph
 */

class CrsGraph_View : public ViewTransform<Epetra_CrsGraph> {

  const Epetra_BlockMap * NewRowMap_;
  const Epetra_BlockMap * NewColMap_;

 public:

  ///
  /** Destructor
   */
  ~CrsGraph_View();

  ///
  /** Constructor
   * input param new_row_map - Row Map for the new "view" object
   * input param new_col_map - Col Map for the new "view" object
   */
  CrsGraph_View( const Epetra_BlockMap * new_row_map,
                 const Epetra_BlockMap * new_col_map = 0 )
  : NewRowMap_(new_row_map),
    NewColMap_(new_col_map)
  {}

  ///
  /** Generates the contiguous subset "view" of the input Epetra_CrsGraph
   */
  NewTypeRef operator()( OriginalTypeRef orig );

};

} //namespace EpetraExt

#endif //EDT_CRSGRAPH_VIEW_H
