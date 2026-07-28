// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _MiniEM_SOLVERS_OperatorRequestCallback_hpp_
#define _MiniEM_SOLVERS_OperatorRequestCallback_hpp_

#include <string>

#include "Teko_Utilities.hpp"
#include "Teko_RequestHandler.hpp"
#include "Teko_RequestCallback.hpp"

#include "Panzer_GlobalEvaluationDataContainer.hpp"

namespace mini_em {

/** \brief Implements the Teko request handler interface to look up
  * assembled operators (matrices) by name from a
  * panzer::GlobalEvaluationDataContainer, for use as preconditioner
  * building blocks.
  */
class OperatorRequestCallback : public Teko::RequestCallback<Teko::LinearOp> {
public:
   /// \brief Construct from the container to look operators up from, and whether to write out the matrices requested.
   OperatorRequestCallback(const Teuchos::RCP<const panzer::GlobalEvaluationDataContainer> & gedc, const bool & matrix_out);

   // RequestCallback member functions
   ///////////////////////////////////////////

   /// \brief Returns the requested operator, looked up by name from the global evaluation data container.
   Teko::LinearOp request(const Teko::RequestMesg & rm);

   /// \brief Returns true if the requested operator's name is present in the global evaluation data container.
   bool handlesRequest(const Teko::RequestMesg & rm);

   /// \brief Asserts that the requested operator's name is present in the global evaluation data container.
   void preRequest(const Teko::RequestMesg & rm);

private:
   Teuchos::RCP<const panzer::GlobalEvaluationDataContainer> gedc_;

   bool matrix_output;
};

}

#endif
