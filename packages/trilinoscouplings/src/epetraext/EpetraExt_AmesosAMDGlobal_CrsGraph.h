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
//                                                                                                    
#ifndef EpetraExt_AMESOS_AMD_CRSGRAPH_H
#define EpetraExt_AMESOS_AMD_CRSGRAPH_H

#include <EpetraExt_Transform.h>
#include <Teuchos_RCP.hpp>
#include <vector>

class Epetra_CrsGraph;
class Epetra_Map;

namespace EpetraExt {

///
/** Block Triangular Factorization (Reordering) of Epetra_CrsGraph
 *
 * Uses Tim Davis' AMDGlobal algorithm to find a block lower or upper triangular
 * ordering form a Epetra_CrsGraph.
 */

class AmesosAMDGlobal_CrsGraph : public StructuralSameTypeTransform<Epetra_CrsGraph> {

 public:

  ///
  /** Destructor
   */
  ~AmesosAMDGlobal_CrsGraph();

  ///
  /** Default Constructor
   */
  AmesosAMDGlobal_CrsGraph( bool verbose = false, bool debug = false )
  : verbose_(verbose),
    debug_(debug)
  {}
    
  ///
  /** Construction of AMDGlobal ordered Epetra_CrsGraph from <tt>orig</tt> object.
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
  std::vector<int> Perm() { return perm_; }
  std::vector<int> BlockPtr() { return blkPtr_; }
  int NumBlocks() { return numBlocks_; }

 private:

  Teuchos::RCP<Epetra_Map> NewRowMap_;
  Teuchos::RCP<Epetra_Map> NewColMap_;
  
  Teuchos::RCP<Epetra_CrsGraph> NewGraph_;

  std::vector<int> perm_, rowPerm_, colPerm_, blkPtr_;

  int numBlocks_;

  const bool verbose_, debug_;
};

} //namespace EpetraExt

#endif //EpetraExt_AMESOS_AMD_CRSGRAPH_H
