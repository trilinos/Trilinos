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


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

