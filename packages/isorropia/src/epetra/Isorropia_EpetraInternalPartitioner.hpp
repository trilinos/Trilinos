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

#ifndef _Isorropia_EpetraInternal_hpp_
#define _Isorropia_EpetraInternal_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Isorropia_EpetraCostDescriber.hpp>
#include <Isorropia_EpetraLibrary.hpp>

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

class InternalPartitioner : public Library {
public:

  /** Constructor
   */
  InternalPartitioner(Teuchos::RCP<const Epetra_CrsGraph> input_graph);
  /** Constructor
   */
  InternalPartitioner(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
		      Teuchos::RCP<CostDescriber> costs);
  /** Constructor
   */
  InternalPartitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix);
  /** Constructor
   */
  InternalPartitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
		      Teuchos::RCP<CostDescriber> costs);
  /** Constructor
   */
  InternalPartitioner(Teuchos::RCP<const Epetra_MultiVector> coords);
  /** Constructor
   */
  InternalPartitioner(Teuchos::RCP<const Epetra_MultiVector> coords,
                      Teuchos::RCP<const Epetra_MultiVector> weights);


  /** Destructor
   */
  virtual ~InternalPartitioner();

  /** Method to partition the object that the InternalPartitioner was contructed with.
    
      \param[in] paramlist  Parameters to govern partitioning.  At this point in
                            time the parameter list is ignored.

      \param[out]  myNewElements  The new partition for each of my objects, in
                                   local ID order.  The objects may be rows (for
                               CrsGraph and RowMatrix input) or coordinates (for
                               MultiVector input).  Partition numbers can range from
                               zero to numProcs-1.
      \param[out]  exportsSize  The number of my objects that will be exported to
                              another process under the new partitioning.  This is
                             also the number of elements in myNewElements that are
                             not equal to my process rank.
      \param[out]  imports   A list of the global IDs of the objects that will be
                            imported to my process under the new partitioning
   */
  virtual int
  repartition(Teuchos::ParameterList& paramlist,
	      std::vector<int>& myNewElements,
	      int& exportsSize,
	      std::vector<int>& imports);

  /** Coloring is not implemented in InternalPartitioner
    */
  virtual int
  color(Teuchos::ParameterList& paramlist,
	std::vector<int>& myNewElements);

  /** Ordering is not implemented in InternalPartitioner
    */
  virtual int
  order(Teuchos::ParameterList& paramlist,
	std::vector<int>& myNewElements);

protected:

  virtual int precompute();
  virtual int postcompute();

};//class InternalPartitioner

}//namespace Epetra
}//namespace Isorropia

#endif //HAVE_EPETRA

#endif

