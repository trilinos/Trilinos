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

#ifndef EpetraExt_CRSMATRIX_BTF_H
#define EpetraExt_CRSMATRIX_BTF_H

#include <EpetraExt_Transform.h>

class Epetra_CrsMatrix;
class Epetra_CrsGraph;
class Epetra_Map;
class Epetra_Import;

namespace EpetraExt {

class CrsMatrix_BTF : public SameTypeTransform<Epetra_CrsMatrix> {

 public:

  ~CrsMatrix_BTF();

  CrsMatrix_BTF( double thres = 0.0,
                 bool verbose = false )
  : NewRowMap_(0),
    NewColMap_(0),
    NewMatrix_(0),
    NewGraph_(0),
    Importer_(0),
    threshold_(thres),
    verbose_(verbose)
  {}

  NewTypeRef operator()( OriginalTypeRef orig );

  bool fwd();
  bool rvs();

 private:

  Epetra_Map * NewRowMap_;
  Epetra_Map * NewColMap_;
  
  Epetra_CrsMatrix * NewMatrix_;
  Epetra_CrsGraph * NewGraph_;

  Epetra_Import * Importer_;

  const double threshold_;

  const bool verbose_;
};

} //namespace EpetraExt

#endif //EpetraExt_CRSMATRIX_BTF_H
