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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

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
