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

//#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

//#include <EEP_Isorropia_TpetraCostDescriber.hpp>
#include <EEP_Isorropia_TpetraOperator.hpp>
#include <Isorropia_Partitioner.hpp>

#ifdef HAVE_ISORROPIA_ZOLTAN
#include <EEP_Isorropia_TpetraZoltanLib.hpp>
#endif
#include <EEP_Isorropia_Tpetra.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

//#include <Epetra_Comm.h>
//#include <Epetra_Map.h>
//#include <Epetra_Import.h>
//#include <Epetra_Vector.h>
//#include <Epetra_MultiVector.h>
#include <Tpetra_CrsGraph_decl.hpp>
//#include <Epetra_CrsMatrix.h>
//#include <Epetra_LinearProblem.h>

//#include <cstring>
//#include <iostream>
//#include <sstream>
//#include <string>
//#include <ctype.h>

#if 1 // def HAVE_TPETRA // AquiToDo

namespace Isorropia {

namespace Tpetra {
  template <class LocalOrdinal,
            class GlobalOrdinal,
            class Node> class CostDescriber;

/** An implementation of the Partitioner interface that operates on
    Epetra matrices and linear systems.

\ingroup partitioning_grp partitioning_rcp_grp partitioning_ptr_grp
*/

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class Partitioner : public Isorropia::Partitioner, public Operator<LocalOrdinal, GlobalOrdinal, Node>  {
public:
  
  /** @ingroup partitioning_rcp_grp 
  */
  Partitioner(Teuchos::RCP< const ::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > inputGraph, // EEP__
              const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),  
              bool compute_partitioning_now=true);
#if 0 // EEP
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
#endif // EEP

  /** @ingroup partitioning_grp
      Destructor 
  */
  virtual ~Partitioner();

#if 0 // EEP
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
#endif // EEP
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
  void partition(bool force_repartitioning=false); // EEP__
#if 0 // EEP
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
#endif // EEP

private:
  int *partGIDs;
  float *partSizes;
  int numPartSizes;
  int printMetrics;

};//class Partitioner

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Partitioner<LocalOrdinal, GlobalOrdinal, Node>::Partitioner(Teuchos::RCP< const ::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > input_graph, // EEP__
                                                            const Teuchos::ParameterList& paramlist,
                                                            bool compute_partitioning_now):
  ::Isorropia::Tpetra::Operator<LocalOrdinal, GlobalOrdinal, Node> (input_graph, paramlist, 0),
  partGIDs(NULL), partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(1)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(1)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

#if 0 // EEP
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_CrsGraph *input_graph,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_CrsGraph>(input_graph,false), paramlist, 0),
  partGIDs(NULL), partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(2)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(2)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
			 Teuchos::RCP<CostDescriber> costs,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_graph, costs, paramlist, 0) ,
  partGIDs(NULL), partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(3)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(3)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_CrsGraph *input_graph,
			 CostDescriber *costs,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_CrsGraph>(input_graph,false), 
            Teuchos::RCP<CostDescriber>(costs,false), paramlist, 0),
  partGIDs(NULL), partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(4)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(4)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_matrix, paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(5)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(5)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_RowMatrix *input_matrix,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_RowMatrix>(input_matrix,false), paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(6)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(6)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
			 Teuchos::RCP<CostDescriber> costs,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_matrix, costs, paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(7)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(7)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_RowMatrix *input_matrix,
			 CostDescriber *costs,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_RowMatrix>(input_matrix,false), 
            Teuchos::RCP<CostDescriber>(costs,false), paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(8)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(8)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(Teuchos::RCP<const Epetra_MultiVector> coords,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (coords, paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(9)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(9)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_MultiVector *coords,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_MultiVector>(coords,false), paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(10)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(10)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(Teuchos::RCP<const Epetra_MultiVector> coords,
                         Teuchos::RCP<const Epetra_MultiVector> weights,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (coords, weights, paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(11)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(11)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_MultiVector *coords,
                         const Epetra_MultiVector *weights,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_MultiVector>(coords,false), 
            Teuchos::RCP<const Epetra_MultiVector>(weights,false), paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(12)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(12)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(Teuchos::RCP<const Epetra_BlockMap> input_map,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_map, paramlist, 0),
  partGIDs(NULL), partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(13)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(13)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_BlockMap *input_map,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_BlockMap>(input_map,false), paramlist, 0),
  partGIDs(NULL), partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(14)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(14)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
                         Teuchos::RCP<const Epetra_MultiVector> coords,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_graph, coords, paramlist, 0),
  partGIDs(NULL), partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(15)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(15)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_CrsGraph *input_graph,
                         const Epetra_MultiVector *coords,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_CrsGraph>(input_graph,false), 
            Teuchos::RCP<const Epetra_MultiVector>(coords,false), paramlist, 0),
  partGIDs(NULL), partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(16)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(16)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
			 Teuchos::RCP<CostDescriber> costs,
                         Teuchos::RCP<const Epetra_MultiVector> coords,
                         Teuchos::RCP<const Epetra_MultiVector> weights,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_graph, costs, coords, weights, paramlist, 0) ,
  partGIDs(NULL), partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(17)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(17)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_CrsGraph *input_graph,
			 CostDescriber *costs,
                         const Epetra_MultiVector *coords,
                         const Epetra_MultiVector *weights,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_CrsGraph>(input_graph,false), 
            Teuchos::RCP<CostDescriber>(costs,false), 
            Teuchos::RCP<const Epetra_MultiVector>(coords,false), 
            Teuchos::RCP<const Epetra_MultiVector>(weights,false), paramlist, 0) ,
  partGIDs(NULL), partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(18)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(18)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
                         Teuchos::RCP<const Epetra_MultiVector> coords,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_matrix, coords, paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(19)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(19)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_RowMatrix *input_matrix,
                         const Epetra_MultiVector *coords,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_RowMatrix>(input_matrix,false), 
            Teuchos::RCP<const Epetra_MultiVector>(coords,false), paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(20)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(20)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
			 Teuchos::RCP<CostDescriber> costs,
                         Teuchos::RCP<const Epetra_MultiVector> coords,
                         Teuchos::RCP<const Epetra_MultiVector> weights,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_matrix, costs, coords, weights, paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(21)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(21)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Partitioner::Partitioner(const Epetra_RowMatrix *input_matrix,
			 CostDescriber *costs,
                         const Epetra_MultiVector *coords,
                         const Epetra_MultiVector *weights,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (Teuchos::RCP<const Epetra_RowMatrix>(input_matrix,false), 
            Teuchos::RCP<CostDescriber>(costs,false), 
            Teuchos::RCP<const Epetra_MultiVector>(coords,false), 
            Teuchos::RCP<const Epetra_MultiVector>(weights,false), paramlist, 0) ,
  partGIDs(NULL),  partSizes(NULL), numPartSizes(0), printMetrics(0)
{
  std::cout << "EEP Entering isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(22)..." << std::endl;
  if (compute_partitioning_now)
    partition(true);
  std::cout << "EEP Leaving isorropia/src/tpetra/EEP_Isorropia_TpetraPartitioner.hpp Partitioner::constructor(22)" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////

Partitioner::~Partitioner(){}

void Partitioner::
clearPartSizes()
{
  if (partGIDs){
    delete [] partGIDs;
    partGIDs = NULL;
  }
  if (partSizes){
    delete [] partSizes;
    partSizes = NULL;
  }
  numPartSizes = 0;
}

void Partitioner::
setPartSizes(int len, int *global_part_id, float *part_size)
{
  clearPartSizes();

  if (len < 1) return;

  numPartSizes = len;

  partGIDs = new int [len];
  memcpy(partGIDs, global_part_id, sizeof(int) * len);

  partSizes = new float [len];
  memcpy(partSizes, part_size, sizeof(float) * len);
}
#endif // EEP

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
Partitioner<LocalOrdinal, GlobalOrdinal, Node>::partition(bool force_repartitioning) // EEP__
{
  std::cout << "EEP Entering isorropia/src/tpetra/Isorropia_TpetraPartitioner.hpp Partitioner::partition()..." << std::endl;
  int input_type = Library<LocalOrdinal, GlobalOrdinal, Node>::unspecified_input_; // EEP__

  std::string partitioning_method_str("PARTITIONING METHOD"); // Aqui
  std::string partitioning_method =
    this->paramlist_.get(partitioning_method_str, "UNSPECIFIED");

  std::string zoltan("ZOLTAN");

  if (this->alreadyComputed() && !force_repartitioning)
    return;

#ifdef HAVE_ISORROPIA_ZOLTAN
  std::cout << "EEP In Partitioner::partition(): pos 001" << std::endl;
  if (partitioning_method == "SIMPLE_LINEAR") {
    throw std::runtime_error("Partitioner::partition - Only Zoltan Partitionner is now supported.");
  }

#if 0 // EEP
  if (input_graph_.get() != 0 && input_coords_.get() != 0)
  {
    std::cout << "EEP In Partitioner::partition(): pos 001.1" << std::endl;
    if (weights_.get())
    {
      lib_ = Teuchos::rcp(new ZoltanLibClass(input_graph_,costs_,input_coords_, weights_));
    }
    else
    {
      lib_ = Teuchos::rcp(new ZoltanLibClass(input_graph_,input_coords_));
    }
  }
  else if (input_matrix_.get() != 0 && input_coords_.get() != 0)
  {
    std::cout << "EEP In Partitioner::partition(): pos 001.2" << std::endl;
    if (weights_.get())
    {
      lib_ = Teuchos::rcp(new ZoltanLibClass(input_matrix_,costs_,input_coords_, weights_));
    }
    else
    {
      lib_ = Teuchos::rcp(new ZoltanLibClass(input_matrix_,input_coords_));
    }
  }
  else if (input_graph_.get() != 0) {
#endif // EEP
    std::cout << "EEP In Partitioner::partition(): pos 001.3" << std::endl;
    this->lib_ = Teuchos::rcp(new ZoltanLibClass<LocalOrdinal, GlobalOrdinal, Node>(this->input_graph_, this->costs_)); // EEP__
#if 0 // EEP
  }
  else if (input_matrix_.get() != 0) {
    std::cout << "EEP In Partitioner::partition(): pos 001.4" << std::endl;
    lib_ = Teuchos::rcp(new ZoltanLibClass(input_matrix_, costs_));
  }
  else if (input_coords_.get() != 0)
  {
    std::cout << "EEP In Partitioner::partition(): pos 001.5" << std::endl;
    if (weights_.get())
    {
      lib_ = Teuchos::rcp(new ZoltanLibClass(input_coords_, weights_));
    }
    else
    {
      lib_ = Teuchos::rcp(new ZoltanLibClass(input_coords_));
    }
  }
  else if (input_map_.get() != 0)
  {
    std::cout << "EEP In Partitioner::partition(): pos 001.6" << std::endl;
    lib_ = Teuchos::rcp(new ZoltanLibClass(input_map_));
  }
  else
  {
    throw Isorropia::Exception("Partitioner::partition - no input object.");
  }
#endif // EEP
  this->lib_->numPartSizes = numPartSizes; // EEP__
  this->lib_->partGIDs = partGIDs;
  this->lib_->partSizes = partSizes;

#endif /* HAVE_ISORROPIA_ZOLTAN */
  std::cout << "EEP In Partitioner::partition(): pos 002" << std::endl;
  Teuchos::ParameterList sublist = this->paramlist_.sublist(zoltan);
  // TODO: Add "block" and "random" partitioning.

  if (partitioning_method == "UNSPECIFIED" && sublist.isParameter("LB_METHOD")) // AquiExceptioThrown
  {
    std::cout << "EEP In Partitioner::partition(): pos 003 temporarylly comenting out an exception" << std::endl;
    //throw Isorropia::Exception("Isorropia \"PARTITIONING METHOD\" as to be set\n"
    //                           "ZOLTAN/LB_METHOD is no longer supported.\n"
    //                           "See readme and release notes for details.");
  }

#if 0 // EEP
  if (input_coords_.get() != 0)
  {
    std::cout << "EEP In Partitioner::partition(): pos 004.1" << std::endl;
    if (partitioning_method == "UNSPECIFIED")
    {
      sublist.set("LB_METHOD", "RCB");
      input_type = Library::geometric_input_;
    }
    else if (partitioning_method == "BLOCK")
    {
      input_type = Library::simple_input_;
      sublist.set("LB_METHOD", "BLOCK");
    }
    else if (partitioning_method == "CYCLIC")
    {
      input_type = Library::simple_input_;
      sublist.set("LB_METHOD", "CYCLIC");
    }
    else if (partitioning_method == "RANDOM")
    {
      input_type = Library::simple_input_;
      sublist.set("LB_METHOD", "RANDOM");
    }
    else if (partitioning_method == "HIER_GRAPH_GEOM") // Can perhaps simply this partitioning method name by using another parameter
    {
      sublist.set("LB_METHOD", "HIER");
      input_type = Library::graph_geometric_input_;
    }
    else if (partitioning_method == "HIER_HGRAPH_GEOM") // Can perhaps simply this partitioning method name by using another parameter
    {
      sublist.set("LB_METHOD", "HIER");
      input_type = Library::hgraph_geometric_input_;
    }
    else if (partitioning_method == "HIER_HGRAPH_GRAPH_GEOM") // Can perhaps simply this partitioning method name by using another parameter
    {
      sublist.set("LB_METHOD", "HIER");
      input_type = Library::hgraph_graph_geometric_input_;
    }
    else // Default case: copy partitioning_method (e.g., RCB, RIB, HSFC)
    {
      sublist.set("LB_METHOD", partitioning_method);
      input_type = Library::geometric_input_;
    }
  }
  else if (input_graph_.get() != 0 || input_matrix_.get() != 0) // graph or matrix input
  {
    std::cout << "EEP In Partitioner::partition(): pos 004.2" << std::endl;
    if (partitioning_method == "GRAPH")
    {
      input_type = Library::graph_input_;
      sublist.set("LB_METHOD", "GRAPH");
    }
    else if (partitioning_method == "BLOCK")
    {
      input_type = Library::simple_input_;
      sublist.set("LB_METHOD", "BLOCK");
    }
    else if (partitioning_method == "CYCLIC")
    {
      input_type = Library::simple_input_;
      sublist.set("LB_METHOD", "CYCLIC");
    }
    else if (partitioning_method == "RANDOM")
    {
      input_type = Library::simple_input_;
      sublist.set("LB_METHOD", "RANDOM");
    }
    else if (partitioning_method == "HIER_GRAPH")
    {
      input_type = Library::graph_input_;
      sublist.set("LB_METHOD", "HIER");
    }
    else if (partitioning_method == "HIER_HGRAPH_GRAPH") // Can perhaps simplify this partitioning method name by using another parameter
    {
      sublist.set("LB_METHOD", "HIER");
      input_type = Library::hgraph_graph_input_;
    }
    else //Hypergraph by default
    {
#endif // EEP
      std::cout << "EEP In Partitioner::partition(): pos 004.3" << std::endl;
      input_type = Library<LocalOrdinal, GlobalOrdinal, Node>::hgraph_input_; // EEP__
      sublist.set("LB_METHOD", "HYPERGRAPH");
#if 0 // EEP
    }
  }
  else  // BlockMap input
  {
    std::cout << "EEP In Partitioner::partition(): pos 004.4" << std::endl;
    if (partitioning_method == "CYCLIC")
    {
      input_type = Library::simple_input_;
      sublist.set("LB_METHOD", "CYCLIC");
    }
    else if (partitioning_method == "RANDOM")
    {
      input_type = Library::simple_input_;
      sublist.set("LB_METHOD", "RANDOM");
    }
    else // BLOCK by default
    {
      input_type = Library::simple_input_;
      sublist.set("LB_METHOD", "BLOCK");
    }
  }
#endif
  std::cout << "EEP In Partitioner::partition(): pos 005" << std::endl;

  if (this->paramlist_.isParameter("NUM PARTS")) {
    sublist.set("NUM_GLOBAL_PARTS", this->paramlist_.template get<std::string>("NUM PARTS"));
  }
  if (this->paramlist_.isParameter("IMBALANCE TOL")) {
    sublist.set("IMBALANCE_TOL", this->paramlist_.template get<std::string>("IMBALANCE TOL"));
  }
  if (this->paramlist_.isParameter("BALANCE OBJECTIVE")
      && this->paramlist_.template get<std::string>("BALANCE OBJECTIVE") == "NONZEROS") {
    sublist.set("ADD_OBJ_WEIGHT", "NONZEROS");
  }

  std::string print_metrics_str("PRINT ZOLTAN METRICS");
  std::string print_metrics =
    this->paramlist_.get(print_metrics_str, "UNSPECIFIED");

  if (print_metrics == ("UNSPECIFIED")){
    printMetrics = 0;
  }
  else if (print_metrics == ("1")){
    printMetrics = 1;
  }
  else if (print_metrics == ("2")){
    printMetrics = 2;
  }
  else {
    printMetrics = 1;
  }

  this->lib_->input_type_ = input_type; // EEP__
  this->lib_->repartition(sublist, this->properties_, this->exportsSize_, this->imports_);
  std::cout << "EEP In Partitioner::partition(): pos 006" << std::endl;
  this->computeNumberOfProperties();
  this->operation_already_computed_ = true;
  std::cout << "EEP Leaving isorropia/src/tpetra/Isorropia_TpetraPartitioner.hpp Partitioner::partition()" << std::endl;
}

#if 0 // EEP
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void Partitioner::
compute(bool force_repartitioning)
{
  partition(force_repartitioning);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Create a new RowMap 
////////////////////////////////////////////////////////////////////////////////
Teuchos::RCP<Epetra_Map>
Partitioner::createNewMap()
{
  Epetra_Map *outputMap;

  createNewMap(outputMap);

  return( Teuchos::RCP<Epetra_Map>(outputMap) );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Create a new RowMap 
////////////////////////////////////////////////////////////////////////////////
void
Partitioner::createNewMap(Epetra_Map * &outputMap)
{
  if (!alreadyComputed()) {
    partition();
  }

  //Generate New Element List
  int myPID = input_map_->Comm().MyPID();
  int numMyElements = input_map_->NumMyElements();
  std::vector<int> elementList( numMyElements );
  if (numMyElements > 0)
    input_map_->MyGlobalElements( &elementList[0] );
  else
    input_map_->MyGlobalElements((int*)NULL); // disambiguate int/long long

  int newGIDSize = numMyElements - exportsSize_;

  std::vector<int> myNewGID;

  if (newGIDSize > 0){
    myNewGID.resize(newGIDSize);
    std::vector<int>::iterator newElemsIter;
    std::vector<int>::const_iterator elemsIter;

    for (elemsIter = properties_.begin(), newElemsIter= myNewGID.begin() ;
         elemsIter != properties_.end() ; elemsIter ++) {
      if ((*elemsIter) == myPID) {
        (*newElemsIter) = elementList[elemsIter - properties_.begin()];
        newElemsIter ++;
      }
    }
  }
  //Add imports to end of list
  myNewGID.insert(myNewGID.end(), imports_.begin(), imports_.end());

  int *gidptr;
  if (myNewGID.size() > 0)
    gidptr = &myNewGID[0];
  else
    gidptr = NULL;

  outputMap = new Epetra_Map(-1, myNewGID.size(), gidptr, 0, input_map_->Comm());

  return;
}
////////////////////////////////////////////////////////////////////////////////

#endif // EEP

}//namespace Tpetra
}//namespace Isorropia

#endif //HAVE_TPETRA

#endif


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

