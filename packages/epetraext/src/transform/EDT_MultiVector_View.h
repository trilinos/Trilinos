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

#ifndef EDT_MULTIVECTOR_VIEW_H
#define EDT_MULTIVECTOR_VIEW_H

#include <Epetra_Transform.h>

class Epetra_MultiVector;
class Epetra_BlockMap;

namespace EpetraExt {

class MultiVector_View : public ViewTransform<Epetra_MultiVector> {

  const Epetra_BlockMap & OrigMap_;
  const Epetra_BlockMap & NewMap_;

  const int NumVec_;

 public:

  ~MultiVector_View();

  MultiVector_View( const Epetra_BlockMap & orig_map,
                    const Epetra_BlockMap & new_map,
                    const int num_vec = -1 )
  : OrigMap_(orig_map),
    NewMap_(new_map),
    NumVec_(num_vec)
  {}

  NewTypeRef operator()( OriginalTypeRef orig );

};

} //namespace EpetraExt

#endif //EDT_MULTIVECTOR_VIEW_H
