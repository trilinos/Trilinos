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

#ifndef   __myClosureModelFactoryImpl_hpp__
#define   __myClosureModelFactoryImpl_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <iostream>
#include <sstream>
#include <typeinfo>

// Panzer
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Constant.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Integrator_Scalar.hpp"

// Phalanx
#include "Phalanx_FieldTag_Tag.hpp"

// Teuchos
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_TypeNameTraits.hpp"

// Files for this specific example.
#include "mySourceTerm.hpp"

///////////////////////////////////////////////////////////////////////////////
//
//  buildClosureModels()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT>
Teuchos::RCP<std::vector<Teuchos::RCP<PHX::Evaluator<panzer::Traits>>>>
MyClosureModelFactory<EvalT>::
buildClosureModels(
  const std::string&                           modelId,
  const Teuchos::ParameterList&                models,
  const panzer::FieldLayoutLibrary&            /* fl            */,
  const Teuchos::RCP<panzer::IntegrationRule>& ir,
  const Teuchos::ParameterList&                /* defaultParams */,
  const Teuchos::ParameterList&                /* userData      */,
  const Teuchos::RCP<panzer::GlobalData>&      /* globalData    */,
  PHX::FieldManager<panzer::Traits>&           /* fm            */) const
{
  using   panzer::BasisIRLayout;
  using   panzer::Constant;
  using   panzer::PureBasis;
  using   panzer::Traits;
  using   PHX::Evaluator;
  using   std::cout;
  using   std::endl;
  using   std::logic_error;
  using   std::string;
  using   std::vector;
  using   Teuchos::getValue;
  using   Teuchos::ParameterEntry;
  using   Teuchos::ParameterList;
  using   Teuchos::RCP;
  using   Teuchos::rcp;
  typedef Teuchos::ParameterList::ConstIterator ParamIter;

  // Create a vector of evaluators to return.
  RCP<vector<RCP<Evaluator<Traits>>>> evaluators =
    rcp(new vector<RCP<Evaluator<Traits>>>);

  // Ensure that the requested closure model is actually specified in the list
  // of all closure models.
  if (not models.isSublist(modelId))
  {
    models.print(cout);
    TEUCHOS_TEST_FOR_EXCEPTION(true, logic_error, "Failed to find requested " \
      "model:  \"" << modelId << "\"." << endl)
  } // end if (not models.isSublist(modelId)

  // Loop over all the closure models that correspond to the given model ID.
  const ParameterList& myModels = models.sublist(modelId);
  for (ParamIter itr = myModels.begin(); itr != myModels.end(); ++itr)
  {
    bool found(false);

    // Extract the current closure model information, which consists of the
    // "key" (its name) and the ParameterList that defines the closure model
    // itself.
    const string key(itr->first);
    ParameterList input;
    const ParameterEntry& entry = itr->second;
    const ParameterList&  plist = getValue<ParameterList>(entry);

    // Create an Evaluator for MySourceTerm.
    if ((plist.isType<string>("Type")             )  and
        (plist.get<string>("Type") == "MySourceTerm"))
    {
      RCP<Evaluator<Traits>> e =
        rcp(new MySourceTerm<EvalT, Traits>(key, *ir));
      evaluators->push_back(e);
      found = true;
    } // end check that the input is correct

    // Fail if we didn't create an Evaluator corresponding to a closure model
    // specified in the input ParameterList.
    if (not found)
      TEUCHOS_TEST_FOR_EXCEPTION(true, logic_error, "ClosureModelFactory "    \
        "failed to build evaluator for key \"" << key << "\"\nin model \""
        << modelId << "\".  Please correct the type or add support to the "   \
        "factory." << endl)
  } // end loop over models in parameter list
  return evaluators;
} // end of buildClosureModels()

#endif // __myClosureModelFactoryImpl_hpp__
