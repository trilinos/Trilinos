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

template <class Node = ::Tpetra::Map<int, int>::node_type >
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

