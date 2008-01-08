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

#ifndef EpetraExt_AMESOS_BTF_CRSMATRIX_H
#define EpetraExt_AMESOS_BTF_CRSMATRIX_H

#include <EpetraExt_Transform.h>
#include <Teuchos_RCP.hpp>

class Epetra_CrsMatrix;
class Epetra_CrsGraph;
class Epetra_Map;
class Epetra_Import;

namespace EpetraExt {

///
/** Block Triangular Factorization (Reordering) of Epetra_CrsMatrix
 *
 * Uses Tim Davis' BTF algorithm to find a block upper triangular
 * ordering form a Epetra_CrsMatrix.
 */

class AmesosBTF_CrsMatrix : public SameTypeTransform<Epetra_CrsMatrix> {

 public:

  ~AmesosBTF_CrsMatrix();

  AmesosBTF_CrsMatrix( double thres = 0.0,
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

  Teuchos::RCP<Epetra_Map> NewRowMap_;
  Teuchos::RCP<Epetra_Map> NewColMap_;
  
  Teuchos::RCP<Epetra_CrsMatrix> NewMatrix_;
  Teuchos::RCP<Epetra_CrsGraph> NewGraph_;

  Teuchos::RCP<Epetra_Import> Importer_;

  const double threshold_;

  const bool verbose_;
};

} //namespace EpetraExt

#endif //EpetraExt_AMESOS_BTF_CRSMATRIX_H
