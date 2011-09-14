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
  
#include <EpetraExt_BlockAdjacencyGraph.h>

#include <Epetra_CrsMatrix.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>

#include <math.h>
#include <cstdlib>

namespace EpetraExt {

  int compare_ints(const void *a, const void *b)
  {
    return(*((int *) a) - *((int *) b));
  }

  int ceil31log2(int n)
  { // Given 1 <= n < 2^31, find l such that 2^(l-1) < n <= 2^(l)
   int l=0, m = 1;
   while( n > m & l < 31 )
   {
      m = 2*m;
      ++l;
   }
   return(l);
  }
//  Purpose: Compute the block connectivity graph of a matrix.
//  An nrr by nrr sparse matrix admits a (Dulmage-Mendelsohn)
//  permutation to block triangular form with nbrr blocks.
//  Abstractly, the block triangular form corresponds to a partition
//  of the set of integers {0,...,n-1} into nbrr disjoint sets.
//  The graph of the sparse matrix, with nrr vertices, may be compressed
//  into the graph of the blocks, a graph with nbrr vertices, that is
//  called here the block connectivity graph.
//     The partition of the rows and columns of B is represented by
//  r(0:nbrr),  0 = r(0) < r(1) < .. < r(nbrr) = nrr, 
//  The graph (Mp,Mj) of the nbrr x nbrr matrix is represened by
//  a sparse matrix in sparse coordinate format.
//  Mp: row indices, dimension determined here (nzM).
//  Mj: column indices, dimension determined here (nzM).
//  The integer vector, weights, of block sizes  (dimension nbrr) is also
//  computed here.  
//     The case of nbrr proportional to nrr is critical.  One must
//  efficiently look up the column indices of B in the partition.
//  This is done here using a binary search tree, so that the
//  look up cost is nzB*log2(nbrr).  

  Teuchos::RCP<Epetra_CrsGraph> BlockAdjacencyGraph::compute( Epetra_CrsGraph& B, int nbrr, std::vector<int>&r, std::vector<double>& weights, bool verbose)
  {
    // Check if the graph is on one processor.
    int myMatProc = -1, matProc = -1;
    int myPID = B.Comm().MyPID();
    for (int proc=0; proc<B.Comm().NumProc(); proc++)
      {
	if (B.NumGlobalEntries() == B.NumMyEntries())
	  myMatProc = myPID;
      }
    B.Comm().MaxAll( &myMatProc, &matProc, 1 );
    
    if( matProc == -1)
      { cout << "FAIL for Global!  All CrsGraph entries must be on one processor!\n"; abort(); }
    
    int i= 0, j = 0, k, l = 0, p, pm, q = -1, ns;
    int tree_height;
    int error = -1;    /* error detected, possibly a problem with the input */
    int nrr;           /* number of rows in B */
    int nzM = 0;       /* number of edges in graph */
    int m = 0;         /* maximum number of nonzeros in any block row of B */
    int* colstack = 0; /* stack used to process each block row */
    int* bstree = 0;   /* binary search tree */
    std::vector<int> Mi, Mj, Mnum(nbrr+1,0);
    nrr = B.NumMyRows();
    if ( matProc == myPID && verbose )
      std::printf(" Matrix Size = %d      Number of Blocks = %d\n",nrr, nbrr);
    else
      nrr = -1;     /* Prevent processor from doing any computations */
    bstree = csr_bst(nbrr);  /* 0 : nbrr-1 */
    tree_height = ceil31log2(nbrr) + 1;
    error = -1;

    l = 0; j = 0; m = 0;
    for( i = 0; i < nrr; i++ ){
      if( i >= r[l+1] ){
	++l;                 /* new block row */
	m = EPETRA_MAX(m,j) ;   /* nonzeros in block row */
	j = B.NumGlobalIndices(i);
      }else{
	j += B.NumGlobalIndices(i);
      }
    }
    /* one more time for the final block */
     m = EPETRA_MAX(m,j) ;   /* nonzeros in block row */

    colstack = (int*) malloc( EPETRA_MAX(m,1) * sizeof(int) );
    // The compressed graph is actually computed twice,
    // due to concerns about memory limitations.  First, 
    // without memory allocation, just nzM is computed.  
    // Next Mj is allocated. Then, the second time, the
    // arrays are actually populated.
    nzM = 0; q = -1; l = 0;
    int * indices;
    int numEntries;
    for( i = 0; i <= nrr; i++ ){
      if( i >= r[l+1] ){
	if( q > 0 ) std::qsort(colstack,q+1,sizeof(int),compare_ints); /* sort stack */
	if( q >= 0 ) ns = 1; /* l, colstack[0] M */
	for( j=1; j<=q ; j++ ){ /* delete copies */
	  if( colstack[j] > colstack[j-1] ) ++ns;
	}
	nzM += ns; /*M->p[l+1] = M->p[l] + ns;*/
	++l;
	q = -1;
      }
      if( i < nrr ){
	B.ExtractMyRowView( i, numEntries, indices );
	for( k = 0; k < numEntries; k++){
	  j = indices[k];  ns = 0; p = 0;
	  while( (r[bstree[p]] > j)  ||  (j >= r[bstree[p]+1])  ){
	    if( r[bstree[p]] > j){
	      p = 2*p+1;
	    }else{
	      if( r[bstree[p]+1] <= j) p = 2*p+2;
	    }
	    ++ns;
	    if( p > nbrr || ns > tree_height ) {
	      error = j;
	      std::printf("error: p %d  nbrr %d  ns %d %d\n",p,nbrr,ns,j); break;
	    }
	  }
	  colstack[++q] = bstree[p];
	}
	//if( error >-1 ){ std::printf("%d\n",error); break; }
        // p > nbrr is a fatal error that is ignored
      }
    }
    
    if ( matProc == myPID && verbose )
      std::printf("nzM =  %d \n", nzM );
    Mi.resize( nzM );
    Mj.resize( nzM );
    q = -1; l = 0; pm = -1;
    for( i = 0; i <= nrr; i++ ){
      if( i >= r[l+1] ){
	if( q > 0 ) std::qsort(colstack,q+1,sizeof(colstack[0]),compare_ints); /* sort stack */
	if( q >= 0 ){
	  Mi[++pm] = l;
	  Mj[pm] = colstack[0];
	}
	for( j=1; j<=q ; j++ ){ /* delete copies */
	  if( colstack[j] > colstack[j-1] ){ /* l, colstack[j] */
	    Mi[++pm] = l;
	    Mj[pm] = colstack[j];
	  }
	}
	++l;
	Mnum[l] = pm + 1;
	
	/* sparse row format: M->p[l+1] = M->p[l] + ns; */
	q = -1;
      }
      if( i < nrr ){
	B.ExtractMyRowView( i, numEntries, indices );
	for( k = 0; k < numEntries; k++){
	  j = indices[k]; ns = 0; p = 0;
	  while( (r[bstree[p]] > j)  ||  (j >= r[bstree[p]+1])  ){
	    if( r[bstree[p]] > j){
	      p = 2*p+1;
	    }else{
	      if( r[bstree[p]+1] <= j) p = 2*p+2;
	    }
	    ++ns;
	  }
	  colstack[++q] = bstree[p];
	}
      }
    }
    if ( bstree ) free ( bstree );
    if ( colstack ) free( colstack );
    
    // Compute weights as number of rows in each block.
    weights.resize( nbrr );
    for( l=0; l<nbrr; l++) weights[l] = r[l+1] - r[l];
    
    // Compute Epetra_CrsGraph and return
    Teuchos::RCP<Epetra_Map> newMap;
    if ( matProc == myPID )
      newMap = Teuchos::rcp( new Epetra_Map(nbrr, nbrr, 0, B.Comm() ) );
    else
      newMap = Teuchos::rcp( new Epetra_Map( nbrr, 0, 0, B.Comm() ) );
    Teuchos::RCP<Epetra_CrsGraph> newGraph = Teuchos::rcp( new Epetra_CrsGraph( Copy, *newMap, 0 ) );
    for( l=0; l<newGraph->NumMyRows(); l++) {
      newGraph->InsertGlobalIndices( l, Mnum[l+1]-Mnum[l], &Mj[Mnum[l]] );
    }
    newGraph->FillComplete();
    
    return (newGraph);  
  }
  
  /*
   * bst(n) returns the complete binary tree, stored in an integer array of dimension n.
   * index i has children 2i+1 and 2i+2, root index 0.
   * A binary search of a sorted n-array: find array(k)
   * i=0; k = tree(i);
   * if < array( k ), i = 2i+1;
   * elif > array( k ), i = 2i+2;
   * otherwise exit
   */
  int* BlockAdjacencyGraph::csr_bst( int n )
  {
    int i, l=1, m, nstack = 0, nexp=0, os, *array, *stack;
    int max_nstack = 0;
    if( n == 0 ) return(NULL);
    while( l <= n ){
      l = 2*l;
      ++nexp;
    } /* n < l = 2^nexp */
    array = (int *) malloc( n * sizeof(int) );
    stack = (int *) malloc(3*nexp * sizeof(int) );
    stack[3*nstack] = 0; stack[3*nstack+1] = 0; stack[3*nstack+2] = n;
    ++nstack;
    /*if( debug ) std::printf("stack : %d %d %d\n", stack[0] , stack[1], stack[2] );*/
    while( nstack > 0 ){
      --nstack;
      i = stack[3*nstack]; os = stack[3*nstack+1]; m = stack[3*nstack+2];
      array[i] = csr_bstrootindex(m) + os; /* 5 */
      if( 2*i+2 < n){   /* right child                4     3      1                      3 - 5 - 1     */
	stack[3*nstack] = 2*i+2; stack[3*nstack+1] = array[i] + 1 ; stack[3*nstack+2] = m-array[i]-1+os; /* 6 10 -3 */
	++nstack;
      }
      if( 2*i+1 < n){   /* left  child */
	stack[3*nstack] = 2*i+1; stack[3*nstack+1] = os; stack[3*nstack+2] = array[i] - os; /* 5 4 5 */
	++nstack;
      }
      if( nstack > max_nstack ) max_nstack =  nstack;
    }
    free( stack );
    return(array);
  }

  /*
    Given a complete binary search tree with n nodes, return
    the index of the root node.  A nodeless tree, n=0, has no
    root node (-1).  Nodes are indexed from 0 to n-1.
  */
  int BlockAdjacencyGraph::csr_bstrootindex( int n )
  {
    int l = 1, nexp = 0, i;
    if( n == 0) return(-1);
    while( l <= n ){
      l = 2*l;
      ++nexp;
    } /* n < l = 2^nexp */
    i = l/2 - 1;
    if(n<4) return(i);
    if (n-i < l/4)
      return( n - l/4  );
    else
      return(i);  
  }
  
} //namespace EpetraExt

