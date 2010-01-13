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

#ifndef _Isorropia_EpetraPartitioner_hpp_
#define _Isorropia_EpetraPartitioner_hpp_

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

*/

class Partitioner : public Isorropia::Partitioner, public Isorropia::Epetra::Operator  {
public:
  
  Partitioner(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
              bool compute_partitioning_now=true);

  Partitioner(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
              const Teuchos::ParameterList& paramlist,
              bool compute_partitioning_now=true);
  
  Partitioner(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
              Teuchos::RCP<CostDescriber> costs,
              const Teuchos::ParameterList& paramlist,
              bool compute_partitioning_now=true);

  Partitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
              bool compute_partitioning_now=true);

  Partitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
              const Teuchos::ParameterList& paramlist,
              bool compute_partitioning_now=true);

  Partitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
              Teuchos::RCP<CostDescriber> costs,
              const Teuchos::ParameterList& paramlist,
              bool compute_partitioning_now=true);

  Partitioner(Teuchos::RCP<const Epetra_MultiVector> coords,
              bool compute_partitioning_now=true);

  Partitioner(Teuchos::RCP<const Epetra_MultiVector> coords,
              const Teuchos::ParameterList& paramlist,
              bool compute_partitioning_now=true);

  Partitioner(Teuchos::RCP<const Epetra_MultiVector> coords,
              Teuchos::RCP<const Epetra_MultiVector> weights,
              const Teuchos::ParameterList& paramlist,
              bool compute_partitioning_now=true);

  Partitioner(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
	      Teuchos::RCP<const Epetra_MultiVector> coords,
              bool compute_partitioning_now=true);

  Partitioner(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
	      Teuchos::RCP<const Epetra_MultiVector> coords,
              const Teuchos::ParameterList& paramlist,
              bool compute_partitioning_now=true);
  
  Partitioner(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
              Teuchos::RCP<CostDescriber> costs,
	      Teuchos::RCP<const Epetra_MultiVector> coords,
              Teuchos::RCP<const Epetra_MultiVector> weights,
              const Teuchos::ParameterList& paramlist,
              bool compute_partitioning_now=true);

  Partitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
	      Teuchos::RCP<const Epetra_MultiVector> coords,
              bool compute_partitioning_now=true);

  Partitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
	      Teuchos::RCP<const Epetra_MultiVector> coords,
              const Teuchos::ParameterList& paramlist,
              bool compute_partitioning_now=true);

  Partitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
              Teuchos::RCP<CostDescriber> costs,
	      Teuchos::RCP<const Epetra_MultiVector> coords,
              Teuchos::RCP<const Epetra_MultiVector> weights,
              const Teuchos::ParameterList& paramlist,
              bool compute_partitioning_now=true);


  /** Destructor */
  virtual ~Partitioner();

   /* Set the relative number of objects in each part.  The default is to
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

  /*
   * Free the memory allocated to store part sizes.
   */
  void clearPartSizes();

  /**  partition is the method that computes 
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

  virtual void compute(bool forceRecomputing=false);

  int numElemsInPart(int part) const {
    return (numElemsWithProperty(part));
  }

  void elemsInPart(int part, int* elementList, int len) const {
    elemsWithProperty(part, elementList, len);
  }

  /** Create a new @c Epetra_Map corresponding to the new partition.

      This method is essentially used by the
      Isorropia::Epetra::Redistributor object.

      \return @c Epetra_Map that contains the new distribution of elements.

      \pre The number of parts might be the same or lower than the
      number of processors.
  */
  Teuchos::RCP<Epetra_Map> createNewMap();

private:
  int *partGIDs;
  float *partSizes;
  int numPartSizes;

};//class Partitioner

}//namespace Epetra
}//namespace Isorropia

#endif //HAVE_EPETRA

#endif

