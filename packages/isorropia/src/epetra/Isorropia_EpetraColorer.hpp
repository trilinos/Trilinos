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

*/

class Colorer : public Isorropia::Colorer, public Isorropia::Epetra::Operator {
public:

    /** Constructor

    \param[in] input_graph the graph which is to have colors assigned to its rows
    \param[in] paramlist this parameter list may be used to pass parameters to Zoltan
    \param[in] compute_now  if @c true, the coloring is computed in the constructor, otherwise call Isorropia::Epetra::Colorer::color when you want to compute the coloring, defaults to @c false
    */

  Colorer(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
	  const Teuchos::ParameterList& paramlist,
	  bool compute_now=true);

    /** Constructor

    \param[in] input_matrix the matrix which is to have colors assigned to its rows
    \param[in] paramlist this parameter list may be used to pass parameters to Zoltan
    \param[in] compute_now  if @c true, the coloring is computed in the constructor, otherwise call Isorropia::Epetra::Colorer::color when you want to compute the coloring, defaults to @c true
  */

  Colorer(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
	  const Teuchos::ParameterList& paramlist,
	  bool compute_now=true);

  /** Destructor */
  ~Colorer() {} ;

  /** Compute the coloring if it has not already been computed, same effect as
       Isorropia::Epetra::Colorer::compute

    \param[in] force_coloring if @c true recompute the coloring even if it has already been computed, defaults to @c false
    */

  void color(bool force_coloring=false);

  /** Compute the coloring if it has not already been computed, same effect as
       Isorropia::Epetra::Colorer::color

    \param[in] force_compute if @c true recompute the coloring even if it has already been computed, defaults to @c false
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
  */
  Teuchos::RCP<Epetra_MapColoring> generateRowMapColoring() ;

  /** Generate an @c Epetra_MapColoring object corresponding of columns color.

  Provide access on the coloring thru the EpetraEXT color class @c Epetra_MapColoring.
  This methods requires EpetraEXT support.
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

