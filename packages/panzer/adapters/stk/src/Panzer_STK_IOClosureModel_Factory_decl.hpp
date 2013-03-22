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

#ifndef PANZER_STK_IOCLOSURE_MODEL_FACTORY_DECL_HPP
#define PANZER_STK_IOCLOSURE_MODEL_FACTORY_DECL_HPP

#include "Panzer_config.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_ClosureModel_Factory.hpp"

#include "Panzer_STK_Interface.hpp"

#include <vector>
#include <string>

namespace panzer {
  class InputEquationSet;
}

namespace panzer_stk {

  template<typename EvalT>
  class IOClosureModelFactory : public panzer::ClosureModelFactory<EvalT> {
  public:

    IOClosureModelFactory(const Teuchos::RCP<const panzer::ClosureModelFactory<EvalT> > userCMF_,
                          const Teuchos::RCP<STK_Interface> & mesh,
                          const Teuchos::ParameterList & outputList);

    Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
    buildClosureModels(const std::string& model_id,
		       const Teuchos::ParameterList& models,
		       const panzer::FieldLayoutLibrary& fl,
		       const Teuchos::RCP<panzer::IntegrationRule>& ir,
		       const Teuchos::ParameterList& default_params,
		       const Teuchos::ParameterList& user_data,
		       const Teuchos::RCP<panzer::GlobalData>& global_data,
		       PHX::FieldManager<panzer::Traits>& fm) const;

  private:
    void parseOutputList(const Teuchos::ParameterList & pl,
                         std::map<std::string,std::vector<std::string> > & blockIdToFields) const;

    //! Mesh pointer, will be passed around
    Teuchos::RCP<STK_Interface> mesh_;
 
    //! Map showing which cell averaged fields need to be written out for each element block
    std::map<std::string,std::vector<std::string> > blockIdToCellAvgFields_;

    //! Map showing which cell fields need to be written out for each element block
    std::map<std::string,std::vector<std::string> > blockIdToCellFields_;

    /** Map stating if an evaluator for a particular block ID has been included.
      *
      * This is a bit of hack that is done to gurantee only one evaluator is
      * added to each field manager. However, if an instantiation of this closure model factory is
      * used in multiple places then the appropriate evaluator will be added and
      * required only once. So its likely that not every field manager will have (and require)
      * the scatter cell evaluators.
      */    
    mutable std::map<std::string,bool> blockIdEvaluated_;

    //! we will reuse the drekar closure model factory
    Teuchos::RCP<const panzer::ClosureModelFactory<EvalT> > userCMF_;
  };

  template < >
  Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
  panzer_stk::IOClosureModelFactory<panzer::Traits::Residual>::buildClosureModels(const std::string& model_id,
		       const Teuchos::ParameterList& models,
		       const panzer::FieldLayoutLibrary& fl,
		       const Teuchos::RCP<panzer::IntegrationRule>& ir,
		       const Teuchos::ParameterList& default_params,
		       const Teuchos::ParameterList& user_data,
		       const Teuchos::RCP<panzer::GlobalData>& global_data,
		       PHX::FieldManager<panzer::Traits>& fm) const;

}

#endif
