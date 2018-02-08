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

#ifndef PANZER_CLOSURE_MODEL_FACTORY_HPP
#define PANZER_CLOSURE_MODEL_FACTORY_HPP

#include "PanzerDiscFE_config.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Phalanx_Evaluator.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_ClosureModel_Factory_Base.hpp"

#include <string>
#include <vector>

namespace panzer {

  class FieldLayoutLibrary;
  class IntegrationRule;
  struct GlobalData;

  template<typename EvalT>
  class ClosureModelFactory : public panzer::ClosureModelFactoryBase {

  protected:
    bool m_throw_if_model_not_found;
  public:

    ClosureModelFactory(bool throw_if_model_not_found=true) : m_throw_if_model_not_found(throw_if_model_not_found) {}

    virtual ~ClosureModelFactory() {}

    Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
    virtual  buildClosureModels(const std::string& model_id,
                                const Teuchos::ParameterList& models,
                                const panzer::FieldLayoutLibrary& fl,
                                const Teuchos::RCP<panzer::IntegrationRule>& ir,
                                const Teuchos::ParameterList& equation_set_params,
                                const Teuchos::ParameterList& user_data,
                                const Teuchos::RCP<panzer::GlobalData>& global_data,
                                PHX::FieldManager<panzer::Traits>& fm) const = 0;

    /** This a convenience function for registering the evaluators. Essentially this
      * facilitates better usage of the ClosureModel TM and allows an easy registration
      * process externally without knowning the compile-time evaluation type.
      *
      * \param[in] evaluators Evaluators to register
      * \param[in] fm Field manager where the evaluators will be registered on completion.
      */
    virtual void registerEvaluators(const std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > & evaluators,
                                    PHX::FieldManager<panzer::Traits>& fm) const
    {
      for (std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > >::size_type i=0; i < evaluators.size(); ++i)
        this->template registerEvaluator<EvalT>(fm, evaluators[i]);
    }

    virtual void setThrowOnModelNotFound(bool do_throw) {
      m_throw_if_model_not_found=do_throw;
    }

  };

}

#endif
