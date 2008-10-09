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

#include <Isorropia_EpetraOrderer.hpp>
#ifdef HAVE_ISORROPIA_ZOLTAN
#include <Isorropia_EpetraZoltanLib.hpp>
#endif
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>

#endif

#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <ctype.h>

namespace Isorropia {

#ifdef HAVE_EPETRA

namespace Epetra {


Orderer::Orderer(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
		 const Teuchos::ParameterList& paramlist,
		 bool compute_now):
  Operator (input_graph, paramlist) {
#ifdef HAVE_ISORROPIA_ZOLTAN
  lib_ = Teuchos::rcp(new ZoltanLibClass(input_graph, Library::graph_input_));
#else /* HAVE_ISORROPIA_ZOLTAN */
  throw Isorropia::Exception("Ordering only available in Zoltan");
  return ;
#endif /* HAVE_ISORROPIA_ZOLTAN */

  if (compute_now)
    order(true);
}

Orderer::Orderer(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
		 const Teuchos::ParameterList& paramlist,
		 bool compute_now):
  Operator (input_matrix, paramlist) {
#ifdef HAVE_ISORROPIA_ZOLTAN
  lib_ = Teuchos::rcp(new ZoltanLibClass(input_matrix, Library::graph_input_));
#else /* HAVE_ISORROPIA_ZOLTAN */
  throw Isorropia::Exception("Ordering only available in Zoltan");
  return ;
#endif /* HAVE_ISORROPIA_ZOLTAN */

  if (compute_now)
    order(true);
}


void
Orderer::order(bool force_ordering)
{
  if (alreadyComputed() && !force_ordering)
    return;
  lib_->order(paramlist_, properties_);
  operation_already_computed_ = true;
}



} // namespace EPETRA

#endif //HAVE_EPETRA

}//namespace Isorropia

