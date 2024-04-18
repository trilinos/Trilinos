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

/** An implementation of the Orderer interface that operates on
    Epetra matrices and linear systems.
\ingroup ordering_grp
*/

class Orderer : public Isorropia::Orderer, public Isorropia::Epetra::Operator {
public:

  Orderer(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
	  const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
	  bool compute_now=true);

  Orderer(const Epetra_CrsGraph *input_graph,
	  const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
	  bool compute_now=true);

  Orderer(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
	  const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
	  bool compute_now=true);

  Orderer(const Epetra_RowMatrix * input_matrix,
	  const Teuchos::ParameterList& paramlist=Teuchos::ParameterList("EmptyParameterList"),
	  bool compute_now=true);

  /** Destructor */
  ~Orderer() {} ;

  /** Method which does the work of computing a new ordering.

     \param force_ordering Optional argument defaults to false.

\ingroup ordering_grp
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


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

