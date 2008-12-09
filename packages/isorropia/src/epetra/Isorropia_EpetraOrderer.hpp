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

#ifndef _Isorropia_EpetraOrderer_hpp_
#define _Isorropia_EpetraOrderer_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Isorropia_EpetraCostDescriber.hpp>
#include <Isorropia_EpetraOperator.hpp>
#include <Isorropia_Orderer.hpp>


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

/** An implementation of the Partitioner interface that operates on
    Epetra matrices and linear systems.

*/

class Orderer : virtual public Isorropia::Orderer, public Isorropia::Epetra::Operator {
public:

  Orderer(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
	  const Teuchos::ParameterList& paramlist,
	  bool compute_now=true);

  Orderer(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
	  const Teuchos::ParameterList& paramlist,
	  bool compute_now=true);

  /** Destructor */
  ~Orderer() {} ;

  /** Method which computes a new ordering.

     \param force_ordering Optional argument defaults to false.
        Depending on the implementation, order() should
        only perform a reordering the first time it is called, and
        subsequent repeated calls are no-ops. If the user's intent is
        to re-compute the ordering (e.g., if parameters or other
        inputs have been changed), then setting this flag to true
        will force a new ordering to be computed.
   */
  void order(bool force_ordering=false);


  void compute(bool forceOrdering=false) {
    return (order(forceOrdering));
  }

};//class Orderer

}//namespace Epetra
}//namespace Isorropia

#endif //HAVE_EPETRA

#endif

