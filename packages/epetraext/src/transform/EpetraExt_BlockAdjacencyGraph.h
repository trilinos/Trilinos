//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
