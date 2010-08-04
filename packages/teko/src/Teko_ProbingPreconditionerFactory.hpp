/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
//  
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software. 
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
// 
// ***********************************************************************
// 
// @HEADER
*/

#ifndef __Teko_ProbingPreconditionerFactory_hpp__
#define __Teko_ProbingPreconditionerFactory_hpp__

#include "Teko_Config.h"

#ifdef Teko_ENABLE_Isorropia

// Teko includes
#include "Teko_PreconditionerState.hpp"
#include "Teko_PreconditionerFactory.hpp"

// Isorropia includes
#include "Isorropia_EpetraProber.hpp"

namespace Teko {

/** \brief Preconditioner factory that for (block) diagonals of explicit operators.
  *
  * Preconditioner factory that for (block) diagonals of explicit operators.
  * These operators need to be Epetra_CrsMatrices under the hood or this will bomb.
  * Uses EpetraExt_PointToBlockDiagPermute.
  */
class ProbingPreconditionerFactory 
   : public virtual Teko::PreconditionerFactory {
public:

  //! @name Constructors.
  //@{
  
  /** Build an empty probing preconditioner factory
   */
  ProbingPreconditionerFactory();

  //@}
  
  /** Create the probed preconditioner operator.
   */
  LinearOp buildPreconditionerOperator(LinearOp & lo,PreconditionerState & state) const;

  //! Initialize from a parameter list
  virtual void initializeFromParameterList(const Teuchos::ParameterList & pl);

  void setGraphOperator(const Teko::LinearOp & graphOp);
  void setGraph(const Teuchos::RCP<const Epetra_CrsGraph> & graph);

  void setProberList(const Teuchos::ParameterList & list);

  void setInverseFactory(const Teuchos::RCP<Teko::InverseFactory> & invFactory)
  { invFactory_ = invFactory; }

protected: 
  //! some members
  Teuchos::RCP<Isorropia::Epetra::Prober> prober;
  Teuchos::RCP<Teko::InverseFactory> invFactory_;
};

} // end namespace Teko

#endif
#endif
