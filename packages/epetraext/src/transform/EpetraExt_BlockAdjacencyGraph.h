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
                                                                                                    
#ifndef EpetraExt_BLOCK_ADJACENCY_GRAPH_H
#define EpetraExt_BLOCK_ADJACENCY_GRAPH_H

class Epetra_CrsGraph;

#include <Teuchos_RCP.hpp>
#include <vector>

namespace EpetraExt {

///
/** Given an Epetra_CrsGraph that has block structure, an adjacency graph
 *  is constructed representing the block connectivity of the original graph.
 *
 *  \authors David Day, Heidi Thornquist
 *
 */

int compare_ints(const void *a, const void *b);

class BlockAdjacencyGraph {

public:

  ///
  /** Destructor
   */
  ~BlockAdjacencyGraph() {}

  ///
  /** Constructor
   */
  BlockAdjacencyGraph() {}

  ///
  /** Constructs an adjacency graph representing the block connectivity of 
      the input graph, where \c nbrr is the number of block rows in \c B and \c r
      contains the row index where each block begins.  A reference-counted
      pointer to an Epetra_CrsGraph that has \c nbrr rows is returned as well as the 
      vector of \c weights.  This vector is of length \c nbrr returns some weighting 
      on the block adjacency graph that can be used to balance the original graph \c B.
      Right now, that weighting is just the number of rows in each block.
   */
   Teuchos::RCP<Epetra_CrsGraph> compute( Epetra_CrsGraph& B, int nbrr, std::vector<int>&r, std::vector<double>& weights, bool verbose = false);

private:

  // Some binary search tree helper functions.
  int* csr_bst( int n );
  int csr_bstrootindex( int n );

};

} //namespace EpetraExt

#endif //EpetraExt_BLOCK_ADJACENCY_GRAPH_H
