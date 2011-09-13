//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006,2011) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA

************************************************************************
*/
//@HEADER

#ifndef _Isorropia_EpetraMatcher_hpp_
#define _Isorropia_EpetraMatcher_hpp_

#include <Isorropia_Epetra.hpp>

#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Import.h>

#ifdef ISORROPIA_HAVE_OMP
#include <omp.h>
#endif

#include <Isorropia_ConfigDefs.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#ifdef MIN
#undef MIN
#endif

#define MIN(a,b) ((a)<(b)?(a):(b))

namespace Isorropia {
namespace Epetra {

/** An implementation of the Matcher interface that operates on
    Epetra matrices and Graphs.
    
    One of the use of Maximum Cardinality Matching is that it provides nonzero
    diagonal for the sparse matrices. This multithreaded version of the
    matching algorithms provides an interface to solve the Bipartite Matching
    problem as well as it also provides permutation to get zero free diagonal
    for input sparse matrices.
*/
class Matcher {
private:
    int* mateU_;
    int* mateV_;
    int* Queue_;
    int* LV_;
    int* LU_;
    int* unmatchedU_;
    int *parent_;
    int *lookahead_;
    int *CRS_pointers_;
    int *CRS_indices_;
    double *CRS_vals_;
    bool finish_;
    int U_,V_,E_,avgDegU_,k_star_,icm_,BFSInd_,numThread_,Qst_,Qend_,matched_,choice_;
   const Epetra_CrsMatrix *A_;

#ifdef ISORROPIA_HAVE_OMP
    omp_lock_t *scannedV_;
#endif

    void delete_matched_v();
    void complete_nonperfect_permutation();
    int SGM();  
    int match_dfs();
    int match_hk();
    int construct_layered_graph();
    int find_set_del_M();
    int recursive_path_finder(int, int);
    int dfs_path_finder(int);
    int dfs_augment();
    int augment_matching(int);
    int DW_phase();

public:
    /** Constructor
    \param[in] Epetra_CRSMatrix* compressed row matrix
    \param[in] paramlist this parameter list is used to select the algorithm.
     The input is a string. "PHK","PHKDW", "PDFS" and "PPF" for parallel
     Hopcroft-Karp, parallel Hopcroft-Karp with Duff-Wiberg , parallel DFS and
     parallel Pothen-Fan respectively. The default is the "PHKDW".
    */
    Matcher(const Epetra_CrsMatrix*, const Teuchos::ParameterList&
    paramlist=Teuchos::ParameterList("EmptyParameterList"));
    
    /** Constructor
    \param[in] RCP of the CrsMatrix
    \param[in] paramlist this parameter list is used to select the algorithm.
     The input is a string. "PHK","PHKDW", "PDFS" and "PPF" for parallel
     Hopcroft-Karp, parallel Hopcroft-Karp with Duff-Wiberg , parallel DFS and
     parallel Pothen-Fan respectively. The default is the "PHKDW".
    */
    Matcher(Teuchos::RCP<const Epetra_CrsMatrix>,const Teuchos::ParameterList&
    paramlist=Teuchos::ParameterList("EmptyParameterList"));

    /** Constructor
    \param[in] Pointer to the Epetra_CrsGraph
    \param[in] paramlist this parameter list is used to select the algorithm.
     The input is a string. "PHK","PHKDW", "PDFS" and "PPF" for parallel
     Hopcroft-Karp, parallel Hopcroft-Karp with Duff-Wiberg , parallel DFS and
     parallel Pothen-Fan respectively. The default is the "PHKDW".
    */
    Matcher(const Epetra_CrsGraph *,const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"));
    
    /** Constructor
    \param[in] RCP of the Epetra_CrsGraph
    \param[in] paramlist this parameter list is used to select the algorithm.
     The input is a string. "PHK","PHKDW", "PDFS" and "PPF" for parallel
     Hopcroft-Karp, parallel Hopcroft-Karp with Duff-Wiberg , parallel DFS and
     parallel Pothen-Fan respectively. The default is the "PHKDW".
    */
    Matcher(Teuchos::RCP<const Epetra_CrsGraph>,const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"));
    
    /** Destructor
    */
     virtual ~Matcher();

    /** Copies the matched columns corresponding to each row in a user given
    array.
   
    \param[in] The length of the array provided by the user
    \param[out] The amount of data copied into the array
    \param[in/out] The user array
    */
    void getMatchedColumnsForRowsCopy(int, int&, int* ) const;

    /** Copies the matched rows corresponding to each column in a user given
    array.

    \param[in] The length of the array provided by the user
    \param[out] The amount of data copied into the array
    \param[in/out] The user array
    */
    void getMatchedRowsForColumnsCopy(int, int&, int*) const;

    /* Returns the number of matched vertices from Maximum Cardinality
     Matching set
    */
    int getNumberOfMatchedVertices();

    /** Applies the row permutation from matching and returns a new matrix.

    */
    Teuchos::RCP<Epetra_CrsMatrix> applyRowPermutation();

    /** Applies the column permutation from matching and returns a new matrix.

    */
    Teuchos::RCP<Epetra_CrsMatrix> applyColumnPermutation();

    //Epetra_Map* getPermutedRowMap();

    /* Provides the row map which is actually the complete row permutation
    */
    //Epetra_Map* getPermutedColumnMap();

    /** Provides the column map which is actually the complete column
     * permutation
    */

    /** Computes the maximum cardinality matching, user should call this
     * function after calling the constructor.
    */
    int match();
};
}
}
#endif
