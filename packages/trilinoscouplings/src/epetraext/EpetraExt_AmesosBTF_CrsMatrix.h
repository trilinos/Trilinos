// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
