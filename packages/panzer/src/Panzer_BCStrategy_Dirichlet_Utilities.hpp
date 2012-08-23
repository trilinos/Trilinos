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

#ifndef PANZER_BCSTRATEGY_DIRCHLET_UTILITIES_HPP
#define PANZER_BCSTRATEGY_DIRCHLET_UTILITIES_HPP

#include "Phalanx_FieldManager.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Teuchos_ParameterList.hpp"
#include <vector>
#include <list>
#include <string>
#include <algorithm>

namespace panzer {


/** \brief Builds the closure models for a particular physics block for a dirichlet bc.

    In many cases, we want the material models for a particular
    physics block for a boundary condition.  We could make a direct
    call to the physics block function
    buildAndRegisterClosureModelEvaluatorsForType<EvalT>().  However,
    in many cases this registers certain unneeded evaluators that
    cause problems.  Some evaluators may instert required fields at
    integration points that can not be met for dirichlet conditions
    where we work at basis points.  Use this class to build a subset
    of the closure models for a dirichlet bc.

    \param fm [in/out] FieldManager - on exit will contain registerd evaluators
    \param pb [in] Physics block
    \param factory [in] closure model factory
    \param comma_sep_closure_model_list [in] a string containing a comma separated list of the models to build
    \param models [in] list of models avaialble to the simulation
    \param user_data [in] User suppplied data that optionally can be used to build models
*/
template<typename EvalT>
void 
buildAndRegisterSubsetOfClosureModelEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
						      const panzer::PhysicsBlock& pb,
						      const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
						      const std::string comma_sep_closure_model_list,
						      const Teuchos::ParameterList& models,
						      const Teuchos::ParameterList& user_data)
{
  
  std::vector<std::string> closure_model_vector;
  panzer::StringTokenizer(closure_model_vector,comma_sep_closure_model_list,",");
  
  // copy into list
  std::list<std::string> closure_model_list;
  for (std::vector<std::string>::iterator i=closure_model_vector.begin(); i != closure_model_vector.end(); ++i)
    closure_model_list.push_back(*i);
  
  Teuchos::ParameterList models_to_build;
  
  for (Teuchos::ParameterList::ConstIterator model = models.begin(); model != models.end(); ++model) {
    
    std::list<std::string>::iterator search = 
      std::find(closure_model_list.begin(), closure_model_list.end(), model->first);
    
    if (search != closure_model_list.end()) {
      closure_model_list.erase(search);
      models_to_build.sublist(model->first) = models.sublist(model->first);
    }
    else
      models_to_build.sublist(model->first);
  }
  
  TEUCHOS_TEST_FOR_EXCEPTION(closure_model_list.size() != 0, std::logic_error,
			     "Error - the list of closure models \"" << comma_sep_closure_model_list << "\" contains an invalid model.");

  pb.buildAndRegisterClosureModelEvaluatorsForType<EvalT>(fm,factory,models_to_build,user_data);

}
  
} // namespace panzer

#endif
