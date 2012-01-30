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


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Orderer::Orderer(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
		 const Teuchos::ParameterList& paramlist,
		 bool compute_now):
  Operator (input_graph, paramlist, 0) {
#ifdef HAVE_ISORROPIA_ZOLTAN
  lib_ = Teuchos::rcp(new ZoltanLibClass(input_graph, Library::graph_input_));
#else /* HAVE_ISORROPIA_ZOLTAN */
  throw Isorropia::Exception("Ordering only available in Zoltan");
  return ;
#endif /* HAVE_ISORROPIA_ZOLTAN */

  if (compute_now)
    order(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Orderer::Orderer(const Epetra_CrsGraph *input_graph,
		 const Teuchos::ParameterList& paramlist,
		 bool compute_now):
  Operator (Teuchos::RCP<const Epetra_CrsGraph>(input_graph,false), paramlist, 0) 
{
#ifdef HAVE_ISORROPIA_ZOLTAN
  lib_ = Teuchos::rcp(new ZoltanLibClass(input_graph_, Library::graph_input_));
#else /* HAVE_ISORROPIA_ZOLTAN */
  throw Isorropia::Exception("Ordering only available in Zoltan");
  return ;
#endif /* HAVE_ISORROPIA_ZOLTAN */

  if (compute_now)
    order(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Orderer::Orderer(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
		 const Teuchos::ParameterList& paramlist,
		 bool compute_now):
  Operator (input_matrix, paramlist, 0) {
#ifdef HAVE_ISORROPIA_ZOLTAN
  lib_ = Teuchos::rcp(new ZoltanLibClass(input_matrix, Library::graph_input_));  
#else /* HAVE_ISORROPIA_ZOLTAN */
  throw Isorropia::Exception("Ordering only available in Zoltan");
  return ;
#endif /* HAVE_ISORROPIA_ZOLTAN */

  if (compute_now)
    order(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Orderer::Orderer(const Epetra_RowMatrix *input_matrix,
		 const Teuchos::ParameterList& paramlist,
		 bool compute_now):
  Operator (Teuchos::RCP<const Epetra_RowMatrix>(input_matrix,false), paramlist, 0) 
{
#ifdef HAVE_ISORROPIA_ZOLTAN
  lib_ = Teuchos::rcp(new ZoltanLibClass(input_matrix_, Library::graph_input_));
#else /* HAVE_ISORROPIA_ZOLTAN */
  throw Isorropia::Exception("Ordering only available in Zoltan");
  return ;
#endif /* HAVE_ISORROPIA_ZOLTAN */

  if (compute_now)
    order(true);
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void Orderer::order(bool force_ordering)
{
  if (alreadyComputed() && !force_ordering)
    return;
  lib_->order(paramlist_.sublist("ZOLTAN"), properties_);
  operation_already_computed_ = true;
}
////////////////////////////////////////////////////////////////////////////////



} // namespace EPETRA

#endif //HAVE_EPETRA

}//namespace Isorropia

