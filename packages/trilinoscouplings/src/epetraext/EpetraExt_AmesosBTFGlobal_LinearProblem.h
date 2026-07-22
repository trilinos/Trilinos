// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef EpetraExt_AMESOS_BTF_GLOBAL_LINEARPROBLEM_H
#define EpetraExt_AMESOS_BTF_GLOBAL_LINEARPROBLEM_H

#include <EpetraExt_Transform.h>
#include <Teuchos_RCP.hpp>
#include <vector>
#include <string>

class Epetra_CrsMatrix;
class Epetra_MultiVector;
class Epetra_CrsGraph;
class Epetra_Map;
class Epetra_Import;
class Epetra_LinearProblem;

namespace EpetraExt {

///
/** Block Triangular Factorization (Reordering) of an Epetra_CrsMatrix
 *
 * Uses Tim Davis' BTF algorithm to find a block lower or upper triangular
 * ordering form a Epetra_CrsMatrix.
 */

class AmesosBTFGlobal_LinearProblem: public SameTypeTransform<Epetra_LinearProblem> {

 public:

  ~AmesosBTFGlobal_LinearProblem();

  AmesosBTFGlobal_LinearProblem( double thresh = 0.0, const std::string& balanceType="linear",
			         bool upperTri = false, bool verbose = false, bool debug = false )
  : numBlocks_(0),
    analysisDone_(false),
    threshold_(thresh),
    upperTri_(upperTri),
    verbose_(verbose),
    debug_(debug),
    balance_(balanceType)
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
  
  // Old linear system (pre-transform)
  Teuchos::RCP<Epetra_CrsMatrix> OldMatrix_;
  Teuchos::RCP<Epetra_MultiVector> OldRHS_;
  Teuchos::RCP<Epetra_MultiVector> OldLHS_;

  Teuchos::RCP<Epetra_CrsMatrix> NewMatrix_;
  Teuchos::RCP<Epetra_MultiVector> NewRHS_;
  Teuchos::RCP<Epetra_MultiVector> NewLHS_;
  
  Teuchos::RCP<Epetra_CrsGraph> NewGraph_;

  Teuchos::RCP<Epetra_Import> Importer_, Importer2_;

  std::vector<int> rowPerm_, colPerm_, blockPtr_;

  int numBlocks_;
  bool analysisDone_; 

  const double threshold_;

  const bool upperTri_, verbose_, debug_;

  std::string balance_;

};

} //namespace EpetraExt

#endif //EpetraExt_AMESOS_BTF_GLOBAL_LINEARPROBLEM_H
