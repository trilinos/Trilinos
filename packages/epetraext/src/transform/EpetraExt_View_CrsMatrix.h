// @HEADER
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
// @HEADER

#ifndef EpetraExt_CRSMATRIX_VIEW_H
#define EpetraExt_CRSMATRIX_VIEW_H

#include <EpetraExt_Transform.h>

class Epetra_CrsGraph;
class Epetra_CrsMatrix;

namespace EpetraExt {

//! Generates a sub-block view of a Epetra_CrsMatrix
class CrsMatrix_View : public ViewTransform<Epetra_CrsMatrix> {

  const Epetra_CrsGraph & OrigGraph_;
  const Epetra_CrsGraph & NewGraph_;

 public:

  //! Destructor
  ~CrsMatrix_View();

  //! Constructor
  CrsMatrix_View( const Epetra_CrsGraph & orig_graph,
                  const Epetra_CrsGraph & new_graph )
  : OrigGraph_(orig_graph),
    NewGraph_(new_graph)
  { /*Should test graphs for requirements*/ }

  //! Transformation Operator
  NewTypeRef operator()( OriginalTypeRef orig );

};

} //namespace EpetraExt

#endif //EpetraExt_CRSMATRIX_VIEW_H
