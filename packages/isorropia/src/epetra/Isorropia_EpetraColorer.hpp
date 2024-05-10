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

#ifndef _Isorropia_EpetraColorer_hpp_
#define _Isorropia_EpetraColorer_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Isorropia_EpetraCostDescriber.hpp>
#include <Isorropia_EpetraOperator.hpp>
#include <Isorropia_Colorer.hpp>

#ifdef HAVE_EPETRAEXT
#include <Epetra_MapColoring.h>
#endif /* HAVE_EPETRAEXT */

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

/** An implementation of the Colorer interface that operates on
    Epetra matrices and linear systems.

\ingroup coloring_grp
*/

class Colorer : public Isorropia::Colorer, public Isorropia::Epetra::Operator {
public:

    /** Constructor

    \param input_graph (in) the graph which is to have colors assigned to its rows
    \param paramlist (in) this parameter list may be used to pass parameters to Zoltan
    \param compute_now (in) if @c true, the coloring is computed in the constructor, otherwise call Isorropia::Epetra::Colorer::color when you want to compute the coloring, defaults to @c true

\ingroup coloring_grp
    */

  Colorer(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
	  const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
	  bool compute_now=true);

    /** Constructor

    \param input_graph (in) the graph which is to have colors assigned to its rows
    \param paramlist (in) this parameter list may be used to pass parameters to Zoltan
    \param compute_now (in) if @c true, the coloring is computed in the constructor, otherwise call Isorropia::Epetra::Colorer::color when you want to compute the coloring, defaults to @c true

\ingroup coloring_grp
    */

  Colorer(const Epetra_CrsGraph * input_graph,
	  const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
	  bool compute_now=true);

    /** Constructor

    \param input_matrix (in) the matrix which is to have colors assigned to its rows
    \param paramlist (in) this parameter list may be used to pass parameters to Zoltan
    \param compute_now (in) if @c true, the coloring is computed in the constructor, otherwise call Isorropia::Epetra::Colorer::color when you want to compute the coloring, defaults to @c true

\ingroup coloring_grp
  */

  Colorer(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
	  const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
	  bool compute_now=true);


    /** Constructor

    \param input_matrix (in) the matrix which is to have colors assigned to its rows
    \param paramlist (in) this parameter list may be used to pass parameters to Zoltan
    \param compute_now (in) if @c true, the coloring is computed in the constructor, otherwise call Isorropia::Epetra::Colorer::color when you want to compute the coloring, defaults to @c true

\ingroup coloring_grp
  */

  Colorer(const Epetra_RowMatrix * input_matrix,
	  const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
	  bool compute_now=true);


  /** Destructor */
  ~Colorer() {} ;

  /** Compute the coloring if it has not already been computed, same effect as
       Isorropia::Epetra::Colorer::compute

    \param force_coloring (in) if @c true recompute the coloring even if it has already been computed, defaults to @c false

\ingroup coloring_grp
    */

  void color(bool force_coloring=false);

  /** Compute the coloring if it has not already been computed, same effect as
       Isorropia::Epetra::Colorer::color

    \param force_compute (in) if @c true recompute the coloring even if it has already been computed, defaults to @c false

\ingroup coloring_grp
    */

  void compute(bool force_compute=false) {
    color(force_compute);
  }

#ifdef HAVE_EPETRAEXT
  /** Generate an @c Epetra_MapColoring object.

  Provide access on the coloring thru the EpetraEXT color class @c Epetra_MapColoring.
  This methods requires EpetraEXT support.

  \deprecated It's recommended to use @see generateRowMapColoring() to this.
  */
  __deprecated Teuchos::RCP<Epetra_MapColoring> generateMapColoring() {
    return generateRowMapColoring();
  }


  /** Generate an @c Epetra_MapColoring object corresponding of rows color.

  Provide access on the coloring thru the EpetraEXT color class @c Epetra_MapColoring.
  This methods requires EpetraEXT support.

\ingroup coloring_grp
  */
  Teuchos::RCP<Epetra_MapColoring> generateRowMapColoring() ;

  /** Generate an @c Epetra_MapColoring object corresponding of columns color.

  Provide access on the coloring thru the EpetraEXT color class @c Epetra_MapColoring.
  This methods requires EpetraEXT support.

\ingroup coloring_grp
  */
  Teuchos::RCP<Epetra_MapColoring> generateColMapColoring() ;

private:
  Teuchos::RCP<const Epetra_BlockMap> colmap_;

#endif /* HAVE_EPETRAEXT */

};//class Colorer

}//namespace Epetra
}//namespace Isorropia

#endif //HAVE_EPETRA

#endif


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

