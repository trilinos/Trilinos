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

#ifndef   __myClosureModelFactory_hpp__
#define   __myClosureModelFactory_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Panzer
#include "Panzer_ClosureModel_Factory.hpp"

/**
 *  \brief Our closure model factory.
 *
 *  This class is used to build our closure models.
 */
template<typename EvalT>
class MyClosureModelFactory
  :
  public panzer::ClosureModelFactory<EvalT>
{
  public:
    
    /**
     *  \brief Build the closure models.
     *
     *  This routine builds all the evaluators for all the various closure
     *  models we'll support.  In our case, there is only one, corresponding to
     *  our source term.
     *
     *  \param[in] modelId       The closure model ID, which is the "name" of a
     *                           `ParameterList` in the "Closure Models"
     *                           `ParameterList` in the input XML file.
     *  \param[in] models        The "Closure Models" `ParameterList` from the
     *                           input XML file.
     *  \param[in] fl            This is unused in this routine, though it is
     *                           part of the `ClosureModelFactory` interface.
     *  \param[in] ir            The integration rule that is used in creating
     *                           our closure model objects.
     *  \param[in] defaultParams This is unused in this routine, though it is
     *                           part of the `ClosureModelFactory` interface.
     *  \param[in] userData      This is unused in this routine, though it is
     *                           part of the `ClosureModelFactory` interface.
     *  \param[in] globalData    This is unused in this routine, though it is
     *                           part of the `ClosureModelFactory` interface.
     *  \param[in] fm            This is unused in this routine, though it is
     *                           part of the `ClosureModelFactory` interface.
     *
     *  \returns A list of evaluators corresponding to all the various closure
     *           models we support.
     */
    Teuchos::RCP<std::vector<Teuchos::RCP<PHX::Evaluator<panzer::Traits>>>>
    buildClosureModels(
      const std::string&                           modelId,
      const Teuchos::ParameterList&                models,
      const panzer::FieldLayoutLibrary&            fl,
      const Teuchos::RCP<panzer::IntegrationRule>& ir,
      const Teuchos::ParameterList&                defaultParams,
      const Teuchos::ParameterList&                userData,
      const Teuchos::RCP<panzer::GlobalData>&      globalData,
      PHX::FieldManager<panzer::Traits>&           fm) const;

}; // end of class MyClosureModelFactory

#include "myClosureModelFactoryImpl.hpp"

#endif // __myClosureModelFactory_hpp__
