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

#ifndef _Isorropia_TpetraPartitioner_hpp_
#define _Isorropia_TpetraPartitioner_hpp_


/*! \file Isorropia_TpetraPartitioner.hpp
 *
 * \brief Contains functionality used for partitioning Tpetra
 * matrices, maps, multivectors, etc.
 */


#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>



#ifdef HAVE_ISORROPIA_TPETRA
#include <Isorropia_TpetraOperator.hpp>
#include <Isorropia_Partitioner.hpp>
#include <Isorropia_TpetraZoltanLib.hpp>


#include <Kokkos_DefaultNode.hpp>
#include <Tpetra_CrsGraph_decl.hpp>
#include <Tpetra_RowMatrix.hpp>
#include <Tpetra_MultiVector_decl.hpp>

namespace Isorropia {

namespace Tpetra {




template <class Node = ::Tpetra::Map<int, int>::node_type>
class Partitioner : public Isorropia::Partitioner, public Isorropia::Tpetra::Operator<Node>  {
public:


  Partitioner(Teuchos::RCP<const ::Tpetra::CrsGraph<int,int,Node> > inputGraph,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  Partitioner(const ::Tpetra::CrsGraph<int,int,Node> *inputGraph,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  Partitioner(Teuchos::RCP<const ::Tpetra::CrsGraph<int,int,Node> > inputGraph,
              Teuchos::RCP<CostDescriber<Node> > costs,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  Partitioner(const ::Tpetra::CrsGraph<int,int,Node> *inputGraph,
              CostDescriber<Node> * costs,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  Partitioner(Teuchos::RCP<const ::Tpetra::RowMatrix<double,int,int,Node> > inputMatrix,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  Partitioner(const ::Tpetra::RowMatrix<double,int,int,Node> *inputMatrix,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  Partitioner(Teuchos::RCP<const ::Tpetra::RowMatrix<double,int,int,Node> > inputMatrix,
              Teuchos::RCP<CostDescriber<Node> > costs,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  Partitioner(const ::Tpetra::RowMatrix<double,int,int,Node> *inputMatrix,
              CostDescriber<Node>  *costs,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  Partitioner(Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > coords,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  Partitioner(const ::Tpetra::MultiVector<double,int,int,Node> *coords,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  Partitioner(Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > coords,
              Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > weights,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  Partitioner(const ::Tpetra::MultiVector<double,int,int,Node> *coords,
              const ::Tpetra::MultiVector<double,int,int,Node> *weights,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  Partitioner(Teuchos::RCP<const ::Tpetra::Map<int,int,Node> > inputMap,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  Partitioner(const ::Tpetra::Map<int,int,Node> *inputMap,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);


  // MMW: currently not supporting hiearchical partitioning





  virtual ~Partitioner();

   /* @ingroup partitioning_grp
    * Set the relative number of objects in each part.  The default is to
    * evenly divide objects across parts.  The numbers can be fractions of
    * one, or whole numbers.  Zoltan adds the values supplied and takes the sizes
    * as proportional to that whole.
    *
    * We make a copy of id and size lists.
    *
    * Caller should supply either global part IDs or local part IDs.
    * Part IDs are integers beginning at zero for the first part.
    *
    * No communication is done during this call.  One process can make the call
    * for all parts, or many processes can make the call.  Zoltan checks the
    * consistency of the information provided.
    */

  void setPartSizes(int len, int *global_part_id, float *part_size);

  /* @ingroup partitioning_grp
   * Free the memory allocated to store part sizes.
   */
  void clearPartSizes();

  /** @ingroup partitioning_grp
      partition is the method that computes
       a rebalanced partitioning for the data in the object
      that this class was constructed with.

      \param force_repartitioning Optional argument defaults to false. By
         default, compute_partitioning() only does anything the first time
         it is called, and subsequent repeated calls are no-ops. If the user's
         intent is to re-compute the partitioning (e.g., if parameters
         or other inputs have been changed), then setting this flag to
         true will force a new partitioning to be computed.
   */
  void partition(bool force_repartitioning=false);

  /** @ingroup partitioning_grp
   */
  virtual void compute(bool forceRecomputing=false);

  /** @ingroup partitioning_grp
   */
  int numElemsInPart(int part) const {
    return (numElemsWithProperty(part));
  }

  /** @ingroup partitioning_grp
   */
  void elemsInPart(int part, int* elementList, int len) const {
    elemsWithProperty(part, elementList, len);
  }

  /**
      Create a new @c ::Tpetra::Map corresponding to the new partition.

      This method is essentially used by the
      Isorropia::Tpetra::Redistributor object.

      \return @c ::Tpetra::Map that contains the new distribution of elements.

      \pre The number of parts might be the same or lower than the
      number of processors.
  */
  Teuchos::RCP< ::Tpetra::Map<int,int,Node> > createNewMap();

  /**
      Create a new @c ::Tpetra::Map corresponding to the new partition.

      This method is essentially used by the
      Isorropia::Tpetra::Redistributor object.

      \param[out] outputMap @c ::Tpetra::Map that contains the new distribution of elements.

      \pre The number of parts might be the same or lower than the
      number of processors.
  */
  void createNewMap(::Tpetra::Map<int,int,Node> *&outputMap);

private:
  int *partGIDs;
  float *partSizes;
  int numPartSizes;

};//class Partitioner

}//namespace Tpetra
}//namespace Isorropia

#endif //HAVE_ISORROPIA_TPETRA

#endif

