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

#ifndef _Isorropia_EpetraPartitioner_hpp_
#define _Isorropia_EpetraPartitioner_hpp_


/*! \file Isorropia_EpetraPartitioner.hpp
 *
 * \brief Contains functionality used for partitioning Epetra
 * matrices, maps, multivectors, etc.
 */

/** \example matrix_1.cpp
 *
 * This example demonstrates creating a balanced copy of Epetra_CrsGraph and 
 * Epetra_CrsMatrix objects using Isorropia::Epetra::createBalancedCopy functions.
 *
 */

/** \example geometric/example_rcb.cpp
 *
 * This example shows how to use geometric partitioning methods with multivectors.
 *
 */

/** \example vert_weights.cpp
 *
 * This example shows how to use graph partitioning methods with vertex weights.
 *
 */

/** \example graphedge_weights.cpp
 *
 * This example shows how to use graph partitioning methods with edge weights.
 *
 */

/** \example hgedge_weights.cpp
 *
 * This example shows how to use hypergraph partitioning methods with hyperedge weights.
 *
 */

/** \example part_redist.cpp
 *
 * This example demonstrates repartitioning and redistributing the
 * contents of an Epetra_LinearProblem object, using the Isorropia::Partitioner
 * and Isorropia::Redistributor classes. This program does not use
 * user-specified weights/costs.
 *
 */



#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Isorropia_EpetraCostDescriber.hpp>
#include <Isorropia_EpetraOperator.hpp>
#include <Isorropia_Partitioner.hpp>

#ifdef HAVE_EPETRA
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_RowMatrix;
class Epetra_LinearProblem;

namespace Isorropia {

namespace Epetra {
  class CostDescriber;

/** An implementation of the Partitioner interface that operates on
    Epetra matrices and linear systems.

\ingroup partitioning_grp partitioning_rcp_grp partitioning_ptr_grp
*/

class Partitioner : public Isorropia::Partitioner, public Isorropia::Epetra::Operator  {
public:
  
  /** @ingroup partitioning_rcp_grp 
  */
  Partitioner(Teuchos::RCP<const Epetra_CrsGraph> inputGraph,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),  
              bool compute_partitioning_now=true);

  /** 
      \ingroup partitioning_ptr_grp 
  */
  Partitioner(const Epetra_CrsGraph *inputGraph,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),  
              bool compute_partitioning_now=true);

  /** 
      \ingroup partitioning_rcp_grp 
  */
  Partitioner(Teuchos::RCP<const Epetra_CrsGraph> inputGraph,
              Teuchos::RCP<CostDescriber> costs,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  /** 
      \ingroup partitioning_ptr_grp 
  */
  Partitioner(const Epetra_CrsGraph *inputGraph,
              CostDescriber* costs,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  /** 
      \ingroup partitioning_rcp_grp 
  */
  Partitioner(Teuchos::RCP<const Epetra_RowMatrix> inputMatrix,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  /** 
      \ingroup partitioning_ptr_grp 
  */
  Partitioner(const Epetra_RowMatrix *inputMatrix,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  /** 
      \ingroup partitioning_rcp_grp 
  */
  Partitioner(Teuchos::RCP<const Epetra_RowMatrix> inputMatrix,
              Teuchos::RCP<CostDescriber> costs,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  /** 
      \ingroup partitioning_ptr_grp 
  */
  Partitioner(const Epetra_RowMatrix *inputMatrix,
              CostDescriber *costs,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  /** 
      \ingroup partitioning_rcp_grp 
  */
  Partitioner(Teuchos::RCP<const Epetra_MultiVector> coords,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  /** 
      \ingroup partitioning_ptr_grp 
  */
  Partitioner(const Epetra_MultiVector *coords,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  /** 
      \ingroup partitioning_rcp_grp 
  */
  Partitioner(Teuchos::RCP<const Epetra_MultiVector> coords,
              Teuchos::RCP<const Epetra_MultiVector> weights,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  /** @ingroup partitioning_ptr_grp 
  */
  Partitioner(const Epetra_MultiVector *coords,
              const Epetra_MultiVector *weights,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  /** @ingroup partitioning_rcp_grp 
  */
  Partitioner(Teuchos::RCP<const Epetra_BlockMap> inputMap,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  /** @ingroup partitioning_ptr_grp 
  */
  Partitioner(const Epetra_BlockMap *inputMap,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  /** @ingroup partitioning_rcp_grp 
  */
  Partitioner(Teuchos::RCP<const Epetra_CrsGraph> inputGraph,
	      Teuchos::RCP<const Epetra_MultiVector> coords,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  /** @ingroup partitioning_ptr_grp 
  */
  Partitioner(const Epetra_CrsGraph *inputGraph,
	      const Epetra_MultiVector *coords,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);
  

  /** @ingroup partitioning_rcp_grp 
  */
  Partitioner(Teuchos::RCP<const Epetra_CrsGraph> inputGraph,
              Teuchos::RCP<CostDescriber> costs,
	      Teuchos::RCP<const Epetra_MultiVector> coords,
              Teuchos::RCP<const Epetra_MultiVector> weights,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  /** @ingroup partitioning_ptr_grp 
  */
  Partitioner(const Epetra_CrsGraph *inputGraph,
              CostDescriber *costs,
	      const Epetra_MultiVector *coords,
              const Epetra_MultiVector *weights,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  /** @ingroup partitioning_rcp_grp 
  */
  Partitioner(Teuchos::RCP<const Epetra_RowMatrix> inputMatrix,
	      Teuchos::RCP<const Epetra_MultiVector> coords,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  /** @ingroup partitioning_ptr_grp 
  */
  Partitioner(const Epetra_RowMatrix *inputMatrix,
	      const Epetra_MultiVector *coords,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  /** @ingroup partitioning_rcp_grp 
  */
  Partitioner(Teuchos::RCP<const Epetra_RowMatrix> inputMatrix,
              Teuchos::RCP<CostDescriber> costs,
	      Teuchos::RCP<const Epetra_MultiVector> coords,
              Teuchos::RCP<const Epetra_MultiVector> weights,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);

  /** @ingroup partitioning_ptr_grp 
  */
  Partitioner(const Epetra_RowMatrix *inputMatrix,
              CostDescriber *costs,
	      const Epetra_MultiVector *coords,
              const Epetra_MultiVector *weights,
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
              bool compute_partitioning_now=true);


  /** @ingroup partitioning_grp
      Destructor 
  */
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

  /** @ingroup partitioning_rcp_grp
      Create a new @c Epetra_Map corresponding to the new partition.

      This method is essentially used by the
      Isorropia::Epetra::Redistributor object.

      \return @c Epetra_Map that contains the new distribution of elements.

      \pre The number of parts might be the same or lower than the
      number of processors.
  */
  Teuchos::RCP<Epetra_Map> createNewMap();

  /** @ingroup partitioning_ptr_grp
      Create a new @c Epetra_Map corresponding to the new partition.

      This method is essentially used by the
      Isorropia::Epetra::Redistributor object.

      \param[out] outputMap @c Epetra_Map that contains the new distribution of elements.

      \pre The number of parts might be the same or lower than the
      number of processors.
  */
  void createNewMap(Epetra_Map *&outputMap);

  int printZoltanMetrics() { return printMetrics; }

private:
  int *partGIDs;
  float *partSizes;
  int numPartSizes;
  int printMetrics;

};//class Partitioner

}//namespace Epetra
}//namespace Isorropia

#endif //HAVE_EPETRA

#endif


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

