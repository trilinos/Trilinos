//@HEADER
//************************************************************************
//
//              Isorropia: Partitioning and Load Balancing Package
//                Copyright (2006) Sandia Corporation
//
//Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
//license for use of this work by or on behalf of the U.S. Government.
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
//************************************************************************
//@HEADER

#ifndef _ispatest_lbeval_utils_hpp_
#define _ispatest_lbeval_utils_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <vector>

#ifdef HAVE_EPETRA

#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>
#include <Isorropia_EpetraCostDescriber.hpp>

/*  Load balance evaluation (As calculated in Zoltan_LB_Eval)

  Given an Epetra graph and possible vertex and hyperedge weights,
  calculate the hypergraph balance, cutn and cutl.  We think of
  the rows of the graph as the objects (vertices) to be partitioned
  and the columns of the graph as the hyperedges.

  ====================================

  Definition of hypergraph balance:

  Suppose Size_p is target size of partition p.  (Sum of all Size_p is 1.0)
  (For now Size_p is always 1/p in isorropia - we don't have an interface
  to set different desired partition sizes.)

  Let W_p be the total row weight on process p and WT be the sum of the
  row weights on all processes.  Then

  imbalance on process p =   |W_p - Size_p*WT| / Size_p*WT

  balance = (1 + maximum imbalance over all processes p)

  ====================================

  Definition of hypergraph cutn and cutl:

  Suppose Cut_k is the number of cuts in column k.  (So column k is in 
  Cut_k + 1 different partitions.)  Suppose W_k is the weight of column k.

  cutn = Sum over all k of W_K * ((Cut_k > 0) ? 1 : 0)

  cutl = Sum over all k of Cut_k * W_k

  ====================================

  TODO explain metrics computed in compute_graph_metrics
*/
namespace ispatest {

/** Given an Epetra_Vector of weights, compute the weight imbalance on each process.  Set 
    the minimum imbalance across all processes, the maximum, and the average.  
    myShare is between 0 and 1, and usually would be 1.0/numProc.
  */

int compute_balance(const Epetra_Vector &wgts, double myGoalWeight,
                              double &min, double &max, double &avg);

/** Compute Zoltan-style hypergraph metrics given a partitioned
    CrsGraph and a CostDescriber (weight) object.
    If the CostDescriber has no weights in it, reasonable defaults
    will be used. 
 */
int compute_hypergraph_metrics(const Epetra_CrsGraph &graph,
            Isorropia::Epetra::CostDescriber &costs,
            double &myGoalWeight,
            double &balance, double &cutn, double &cutl);

/** Compute Zoltan-style hypergraph metrics given a partitioned
    RowMatrix and a CostDescriber (weight) object.
    If the CostDescriber has no weights in it, reasonable defaults
    will be used. 
 */
int compute_hypergraph_metrics(const Epetra_RowMatrix &matrix,
            Isorropia::Epetra::CostDescriber &costs,
            double &myGoalWeight,
            double &balance, double &cutn, double &cutl);

/** Compute graph metrics given an Epetra_RowMatrix.
    A CostDescriber object may provide vertex (row) and/or edge (non-zeroes)
    weights, or it may be an initialized object with no weights.  If no vertex
    weights are provided, each vertex is assumed to be weight 1.  If no
    edge weights are provided, each edge is assumed to be weight 1.

    The goal weight for a process is the proportion of the total vertex (row)
    weights that were to be assigned to this process under repartitioning.  If
    all processes are to get an equal proportion of the weight, set this 
    value to (1.0 / # processes).  This value is needed in order to compute
    how close the repartitioning is to being perfectly balanced.

    If the CostDescriber has no weights in it, reasonable defaults
    will be used. 
  */
int compute_graph_metrics(const Epetra_RowMatrix &matrix,
            Isorropia::Epetra::CostDescriber &costs,
            double &myGoalWeight,
            double &balance, int &numCuts, double &cutWgt, double &cutn, double &cutl);

/** Compute graph metrics given an Epetra_CrsGraph.
    A CostDescriber object may provide vertex (row) and/or edge (non-zeroes)
    weights, or it may be an initialized object with no weights.  If no vertex
    weights are provided, each vertex is assumed to be weight 1.  If no
    edge weights are provided, each edge is assumed to be weight 1.

    The goal weight for a process is the proportion of the total vertex (row)
    weights that were to be assigned to this process under repartitioning.  If
    all processes are to get an equal proportion of the weight, set this 
    value to (1.0 / # processes).  This value is needed in order to compute
    how close the repartitioning is to being perfectly balanced.

    If the CostDescriber has no weights in it, reasonable defaults
    will be used. 
  */
int compute_graph_metrics(const Epetra_CrsGraph &graph,
            Isorropia::Epetra::CostDescriber &costs,
            double &myGoalWeight,
            double &balance, int &numCuts, double &cutWgt, double &cutn, double &cutl);

/** Print out a distributed RowMatrix.  This only works for small test
    matrices of 1s and 0s, and 10 or fewer processes.
  */

void show_matrix(const char *txt, const Epetra_RowMatrix &matrix, const Epetra_Comm &comm);

/** Print out a distributed CrsGraph.  This only works for small test
    matrices of 1s and 0s and 10 or fewer processes.
  */

void show_matrix(const char *txt, const Epetra_CrsGraph &graph, const Epetra_Comm &comm);

/** Print out a distributed LinearProblem.  This only works for small test
    matrices of 1s and 0s and 10 or fewer processes.
  */

void show_matrix(const char *txt, const Epetra_LinearProblem &problem, const Epetra_Comm &comm);



}//namespace ispatest

#endif //HAVE_EPTERA

#endif

