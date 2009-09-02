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
#include <vector>

class Epetra_CrsMatrix;
class Epetra_CrsGraph;
class Epetra_Map;
class Epetra_Import;

namespace EpetraExt {

///
/** Block Triangular Factorization (Reordering) of Epetra_CrsMatrix
 *
 * Uses Tim Davis' BTF algorithm to find a block lower or upper triangular
 * ordering form a Epetra_CrsMatrix.
 */

class AmesosBTF_CrsMatrix : public SameTypeTransform<Epetra_CrsMatrix> {

 public:

  ~AmesosBTF_CrsMatrix();

  AmesosBTF_CrsMatrix( double thres = 0.0, bool upperTri = false,
                 bool verbose = false, bool debug = false )
  : numBlocks_(0), 
    threshold_(thres),
    upperTri_(upperTri),
    verbose_(verbose),
    debug_(debug)
  {}

  NewTypeRef operator()( OriginalTypeRef orig );

  bool fwd();
  bool rvs();

  Teuchos::RCP<Epetra_Import> Importer() { return Importer_; }
  std::vector<int> RowPerm() { return rowPerm_; }
  std::vector<int> ColPerm() { return colPerm_; }
  std::vector<int> BlockPtr() { return blockPtr_; }
  int NumBlocks() { return numBlocks_; }

 private:

  Teuchos::RCP<Epetra_Map> NewRowMap_;
  Teuchos::RCP<Epetra_Map> NewColMap_;
  
  Teuchos::RCP<Epetra_CrsMatrix> NewMatrix_;
  Teuchos::RCP<Epetra_CrsGraph> NewGraph_;

  Teuchos::RCP<Epetra_Import> Importer_;

  std::vector<int> rowPerm_, colPerm_, blockPtr_;

  int numBlocks_;

  const double threshold_;

  const bool upperTri_, verbose_, debug_;

};

} //namespace EpetraExt

#endif //EpetraExt_AMESOS_BTF_CRSMATRIX_H
