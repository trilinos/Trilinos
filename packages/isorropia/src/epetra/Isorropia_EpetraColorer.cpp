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

#include <Isorropia_EpetraColorer.hpp>
#ifdef HAVE_ISORROPIA_ZOLTAN
#include <Isorropia_EpetraZoltanLib.hpp>
#endif /* HAVE_ISORROPIA_ZOLTAN */
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
Colorer::Colorer(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
		 const Teuchos::ParameterList& paramlist, bool compute_now):
  Operator(input_graph, paramlist, 1)
{
#ifdef HAVE_EPETRAEXT
  colmap_ = Teuchos::rcp(&(input_graph->ColMap()),false);
#endif /* HAVE_EPETRAEXT */

#ifdef HAVE_ISORROPIA_ZOLTAN
  lib_ = Teuchos::rcp(new ZoltanLibClass(input_graph, Library::graph_input_));
#else /* HAVE_ISORROPIA_ZOLTAN */
  throw Isorropia::Exception("Coloring only available in Zoltan");
  return ;
#endif /* HAVE_ISORROPIA_ZOLTAN */
  if (compute_now)
    color(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Colorer::Colorer(const Epetra_CrsGraph * input_graph,
		 const Teuchos::ParameterList& paramlist, bool compute_now):
  Operator(Teuchos::RCP<const Epetra_CrsGraph>(input_graph,false), paramlist, 1)
{
#ifdef HAVE_EPETRAEXT
  colmap_ = Teuchos::rcp(&(input_graph_->ColMap()),false);
#endif /* HAVE_EPETRAEXT */

#ifdef HAVE_ISORROPIA_ZOLTAN
  lib_ = Teuchos::rcp(new ZoltanLibClass(input_graph_, Library::graph_input_));
#else /* HAVE_ISORROPIA_ZOLTAN */
  throw Isorropia::Exception("Coloring only available in Zoltan");
  return ;
#endif /* HAVE_ISORROPIA_ZOLTAN */
  if (compute_now)
    color(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Colorer::Colorer(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
		 const Teuchos::ParameterList& paramlist,
		 bool compute_now):
  Operator (input_matrix, paramlist, 1) {

#ifdef HAVE_EPETRAEXT
  colmap_ = Teuchos::rcp(&(input_matrix->RowMatrixColMap()),false);
#endif /* HAVE_EPETRAEXT */

#ifdef HAVE_ISORROPIA_ZOLTAN
  lib_ = Teuchos::rcp(new ZoltanLibClass(input_matrix, Library::graph_input_));
#else /* HAVE_ISORROPIA_ZOLTAN */
  throw Isorropia::Exception("Coloring only available in Zoltan");
  return ;
#endif /* HAVE_ISORROPIA_ZOLTAN */

  if (compute_now)
    color(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Colorer::Colorer(const Epetra_RowMatrix *input_matrix,
		 const Teuchos::ParameterList& paramlist, bool compute_now):
  Operator (Teuchos::RCP<const Epetra_RowMatrix>(input_matrix,false), paramlist, 1) 
{
#ifdef HAVE_EPETRAEXT
  colmap_ = Teuchos::rcp(&(input_matrix_->RowMatrixColMap()),false);
#endif /* HAVE_EPETRAEXT */

#ifdef HAVE_ISORROPIA_ZOLTAN
  lib_ = Teuchos::rcp(new ZoltanLibClass(input_matrix_, Library::graph_input_));
#else /* HAVE_ISORROPIA_ZOLTAN */
  throw Isorropia::Exception("Coloring only available in Zoltan");
  return ;
#endif /* HAVE_ISORROPIA_ZOLTAN */

  if (compute_now)
    color(true);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void
Colorer::color(bool force_coloring)
{
  if (alreadyComputed() && !force_coloring)
    return;

  std::string zoltan("ZOLTAN");
  Teuchos::ParameterList sublist = paramlist_.sublist(zoltan);

  if (paramlist_.isParameter("DISTANCE")) {
    sublist.set("COLORING_PROBLEM", "DISTANCE-"+paramlist_.get<std::string>("DISTANCE"));
  }

  lib_->color(sublist, properties_);
  operation_already_computed_ = true;
  computeNumberOfProperties();
}
////////////////////////////////////////////////////////////////////////////////

#ifdef HAVE_EPETRAEXT

Teuchos::RCP<Epetra_MapColoring>
Colorer::generateRowMapColoring()
{
  Teuchos::RCP<Epetra_MapColoring> colorMap;

  color(false);
  colorMap = Teuchos::rcp(new Epetra_MapColoring(*input_map_, 1));

  for(unsigned int i = 0; i < properties_.size(); i++ ) {
    (*colorMap)[i] = properties_[i];
  }
  return (colorMap);
}


Teuchos::RCP<Epetra_MapColoring>
Colorer::generateColMapColoring()
{
  Teuchos::RCP<Epetra_MapColoring> rowColorMap = generateRowMapColoring();

  // Color map has colored rows -- need colored columns
  Epetra_Import importer(*colmap_, *input_map_);

  Teuchos::RCP<Epetra_MapColoring> colorMap =
    Teuchos::rcp(new Epetra_MapColoring(*colmap_));

  colorMap->Import(*rowColorMap, importer, Insert);
  return (colorMap);
}


#endif /* HAVE_EPETRAEXT */

} // namespace EPETRA

#endif //HAVE_EPETRA

}//namespace Isorropia

