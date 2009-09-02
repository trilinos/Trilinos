//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

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

void compute_balance(const Epetra_Vector &wgts, double myGoalWeight,
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
    value to (1.0 / #processes).  This value is needed in order to compute
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
    value to (1.0 / #processes).  This value is needed in order to compute
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

