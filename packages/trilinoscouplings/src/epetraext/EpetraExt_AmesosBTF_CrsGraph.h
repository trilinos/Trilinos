// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef EpetraExt_AMESOS_BTF_CRSGRAPH_H
#define EpetraExt_AMESOS_BTF_CRSGRAPH_H

#include <EpetraExt_Transform.h>
#include <Teuchos_RCP.hpp>
#include <vector>

class Epetra_CrsGraph;
class Epetra_Map;

namespace EpetraExt {

///
/** Block Triangular Factorization (Reordering) of Epetra_CrsGraph
 *
 * Uses Tim Davis' BTF algorithm to find a block lower or upper triangular
 * ordering form a Epetra_CrsGraph.
 */

class AmesosBTF_CrsGraph : public StructuralSameTypeTransform<Epetra_CrsGraph> {

 public:

  ///
  /** Destructor
   */
  ~AmesosBTF_CrsGraph();

  ///
  /** Default Constructor
   */
  AmesosBTF_CrsGraph( bool upperTri = false, bool verbose = false )
    : numBlocks_(0),
      upperTri_(upperTri),
      verbose_(verbose)
    {}
    
  ///
  /** Construction of BTF ordered Epetra_CrsGraph from <tt>orig</tt> object.
   *
   * Preconditions:<ul>
   * </ul>
   *
   * Invariants:<ul>
   * </ul>
   *
   * Postconditions:<ul>
   * </ul>
   *
   */
  NewTypeRef operator()( OriginalTypeRef orig );

  std::vector<int> RowPerm() { return rowPerm_; }
  std::vector<int> ColPerm() { return colPerm_; }
  std::vector<int> BlockPtr() { return blkPtr_; }
  int NumBlocks() { return numBlocks_; }

 private:

  Teuchos::RCP<Epetra_Map> NewRowMap_;
  Teuchos::RCP<Epetra_Map> NewColMap_;
  
  Teuchos::RCP<Epetra_CrsGraph> NewGraph_;

  std::vector<int> rowPerm_, colPerm_, blkPtr_;

  int numBlocks_;

  const bool upperTri_, verbose_;
};

} //namespace EpetraExt

#endif //EpetraExt_AMESOS_BTF_CRSGRAPH_H
