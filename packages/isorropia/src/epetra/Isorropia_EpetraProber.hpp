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

#ifndef _Isorropia_EpetraProber_hpp_
#define _Isorropia_EpetraProber_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Isorropia_Colorer.hpp>
#include <Isorropia_EpetraColorer.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_CrsGraph.h>
#include <Teuchos_RCP.hpp>
class Epetra_MultiVector;
class Epetra_CrsMatrix;
class Epetra_Operator;

namespace Isorropia {

namespace Epetra {

/** An implementation of the Prober interface that operates on
    Epetra matrices and linear systems.  The Prober currently works only on structurally
    symmetric problems.  Support for structually non-symmetric problems is under development.

\ingroup probing_grp
*/

class Prober {
public:

    /** Constructor

    \param[in] input_graph the graph whose sparsity pattern is to guide the probing.
    \param[in] paramlist this parameter list may be used to pass parameters to the colorer.
    \param[in] compute_now  if @c true, the coloring is computed in the constructor, otherwise call Isorropia::Epetra::Prober::color when you want to compute the coloring, defaults to @c false

\ingroup probing_grp
    */

  Prober(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
         const Teuchos::ParameterList& paramlist,
         bool compute_now=true);


    /** Constructor

    \param[in] input_graph the graph whose sparsity pattern is to guide the probing.
    \param[in] paramlist this parameter list may be used to pass parameters to the colorer.
    \param[in] compute_now  if @c true, the coloring is computed in the constructor, otherwise call Isorropia::Epetra::Prober::color when you want to compute the coloring, defaults to @c false

\ingroup probing_grp
    */

  Prober(const Epetra_CrsGraph * input_graph,
         const Teuchos::ParameterList& paramlist,
         bool compute_now=true);

    /** Constructor

    \param[in] input_matrix the matrix whose sparsity pattern is to guide the probing.
    \param[in] paramlist this parameter list may be used to pass parameters to the colorer.
    \param[in] compute_now  if @c true, the coloring is computed in the constructor, otherwise call Isorropia::Epetra::Prober::color when you want to compute the coloring, defaults to @c true

\ingroup probing_grp
  */

  Prober(Teuchos::RCP<const Epetra_CrsMatrix> input_matrix,
         const Teuchos::ParameterList & paramlist,
         bool compute_now=true);


    /** Constructor

    \param[in] input_matrix the matrix whose sparsity pattern is to guide the probing.
    \param[in] paramlist this parameter list may be used to pass parameters to the colorer.
    \param[in] compute_now  if @c true, the coloring is computed in the constructor, otherwise call Isorropia::Epetra::Prober::color when you want to compute the coloring, defaults to @c true

\ingroup probing_grp
  */

  Prober(const Epetra_CrsMatrix * input_matrix,
         const Teuchos::ParameterList & paramlist,
         bool compute_now=true);

  

  /** Default Constructor
  */
  
  Prober();

  
  /** Destructor */
  ~Prober(){delete colorer_;}


  /** Sets the parameter list */
  void setList(const Teuchos::ParameterList& paramlist);

  /** Sets the graph */
  void setGraph(Teuchos::RCP<const Epetra_CrsGraph> input_graph){input_graph_=input_graph; has_colored=false;}  

  /** Compute the coloring.
\ingroup probing_grp
    */
  void color();

  /** Get the number of orthogonal vectors (or the number of colors from
      coloring)
   \return number of orthogonal vectors
\ingroup probing_grp
    */
  int getNumOrthogonalVectors();
  
  /** @ingroup probing_grp  
      Perform the actual probing.
   \param[in] op is the operator we are probing
   \param[in,out] out_matrix is the matrix 
    */
  int probe(const Epetra_Operator & op, Epetra_CrsMatrix & out_matrix);

 /** @ingroup probing_grp
     Perform the actual probing.
   \param[in] op is the operator we are probing
   \return RCP to the matrix
    */
  Teuchos::RCP<Epetra_CrsMatrix> probe(const Epetra_Operator & op);
  
private:
  Teuchos::RCP<const Epetra_CrsGraph> input_graph_;
  Colorer *colorer_;
  Teuchos::ParameterList List_;
  bool has_colored;
};//class Prober

}//namespace Epetra
}//namespace Isorropia

#endif //HAVE_EPETRA

#endif

