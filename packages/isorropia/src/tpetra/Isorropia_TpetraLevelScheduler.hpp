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

#ifndef _Isorropia_TpetraLevelScheduler_hpp_
#define _Isorropia_TpetraLevelScheduler_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <Isorropia_TpetraOperator.hpp>
#include <Isorropia_LevelScheduler.hpp>
#include <Teuchos_RCP.hpp>

#ifdef HAVE_ISORROPIA_TPETRA

#include <Tpetra_CrsGraph_decl.hpp>
#include <Kokkos_DefaultNode.hpp>

namespace Isorropia {

namespace Tpetra {

/** An implementation of the LevelScheduler interface that operates on
    and Tpetra::CrsGraph, representing the non-zeros in a matrix.  The
    elements to be partitioned into levels are matrix rows.  Assumption
    is that matrix is lower triangular or upper triangular.
*/

template <class Node=Kokkos::DefaultNode::DefaultNodeType>
class LevelScheduler :  public Isorropia::LevelScheduler, public Isorropia::Tpetra::Operator<Node> {

public:

    /** Constructor

    \param[in] input_graph the graph representing the non-zeros of the matrix
    \param[in] paramlist list of parameters
    \param[in] compute_now  if @c true, the scheduling is computed in the constructor, otherwise call Isorropia::Tpetra::LevelScheduler::schedule when you want to compute the scheduling, defaults to @c false
    */

  LevelScheduler(Teuchos::RCP<const ::Tpetra::CrsGraph<int,int,Node> > input_graph, 
                 const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
	         bool compute_now=true);

  /** Destructor */
  ~LevelScheduler(); // {} ;

  /** Compute the scheduling if it has not already been computed, same effect as
       Isorropia::Tpetra::LevelScheduler::compute

    \param[in] force_scheduling if @c true recompute the scheduling even if it has already been computed, defaults to @c false
    */

  void schedule(bool force_scheduling=false);

  /** Compute the scheduling if it has not already been computed, same effect as
       Isorropia::Tpetra::LevelScheduler::schedule

    \param[in] force_compute if @c true recompute the scheduling even if it has already been computed, defaults to @c false
    */

  void compute(bool force_compute=false) {
    schedule(force_compute);
  }

private:
  Teuchos::RCP<const ::Tpetra::CrsGraph<int,int,Node> > input_graph_;

};//class LevelScheduler

}//namespace Tpetra
}//namespace Isorropia

#endif //HAVE_ISORROPIA_TPETRA

#endif

