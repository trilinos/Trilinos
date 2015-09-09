// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_EVALUATORS_REGISTRAR_HPP
#define PANZER_EVALUATORS_REGISTRAR_HPP

#include "Phalanx_FieldManager.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {

/** Classes that call PHX::FieldManager::registerEvaluator on
  * panzer::EvaluatorWithBaseImpl objects inherit from this class to wrap the
  * registerEvaluator call. This class injects the WorksetDetails index into the
  * evaluator.
 */
class EvaluatorsRegistrar {
public:
  //! Set the WorksetDetails index in all evaluators registered through
  //! EquationSetBase::registerEvaluator. The details index can be set multiple
  //! times. The current value applies at registration.
  void setDetailsIndex(const int details_index) { details_index_ = details_index; }
  //! Get the WorksetDetails index.
  int getDetailsIndex() const { return details_index_; }

protected:
  //! Default ctor initializes WorksetDetails index to 0.
  EvaluatorsRegistrar() : details_index_(0) {}
  virtual ~EvaluatorsRegistrar() {}

  //! Register the evaluator and initialize it with any information in
  //! panzer::EvaluatorWithBaseImpl.
  template <typename EvalT>
  void registerEvaluator(PHX::FieldManager<panzer::Traits>& fm,
                         const Teuchos::RCP< PHX::Evaluator<panzer::Traits> >& op) const;

private:
  int details_index_;
};

template<typename EvalT>
void EvaluatorsRegistrar::
registerEvaluator(PHX::FieldManager<panzer::Traits>& fm,
                  const Teuchos::RCP< PHX::Evaluator<panzer::Traits> >& op) const
{
  Teuchos::RCP< panzer::EvaluatorWithBaseImpl<panzer::Traits> >
    pop = Teuchos::rcp_dynamic_cast< panzer::EvaluatorWithBaseImpl<panzer::Traits> >(op);
  // Temporarily allow casting failure so that Charon Charon continues to work.
#if 0
  TEUCHOS_TEST_FOR_EXCEPTION(pop.is_null(), std::runtime_error,
                             op->getName() + " does not inherit from panzer::EvaluatorWithBaseImpl.");
  pop->setDetailsIndex(details_index_);
#else
  if (Teuchos::nonnull(pop))
    pop->setDetailsIndex(details_index_);
#endif
  fm.template registerEvaluator<EvalT>(op);
}

}

#endif
