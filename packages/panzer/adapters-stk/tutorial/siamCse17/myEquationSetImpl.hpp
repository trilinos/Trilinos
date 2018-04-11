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

#ifndef   __myEquationSetImpl_hpp__
#define   __myEquationSetImpl_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Panzer
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Integrator_BasisTimesScalar.hpp"
#include "Panzer_Integrator_GradBasisDotVector.hpp"

// Phalanx
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"

// Teuchos
#include "Teuchos_Assert.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

///////////////////////////////////////////////////////////////////////////////
//
//  Default Constructor
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
MyEquationSet<EvalT>::
MyEquationSet(
  const Teuchos::RCP<Teuchos::ParameterList>& params,
  const int&                                  defaultIntegrationOrder,
  const panzer::CellData&                     cellData,
  const Teuchos::RCP<panzer::GlobalData>&     globalData,
  const bool                                  buildTransientSupport)
  :
  panzer::EquationSet_DefaultImpl<EvalT>(params, defaultIntegrationOrder,
    cellData, globalData, buildTransientSupport)
{
  using std::string;
  using Teuchos::ParameterList;

  // Validate and parse the unnamed ParameterList from the input.xml file that
  // corresponds to the "domain" Physics Block.
  {
    ParameterList validParameters;
    this->setDefaultValidParameters(validParameters);
    validParameters.set("Model ID", "", "The closure model ID associated "    \
      "with this equaiton set.");
    validParameters.set("Prefix", "", "A prefix that can be tacked on if "    \
      "you'd like to use multiple instatiations of this equation set.");
    validParameters.set("Basis Type", "HGrad", "The type of basis to use.");
    validParameters.set("Basis Order", 1, "The order of the basis.");
    validParameters.set("Integration Order", -1, "The order of the "          \
      "integration rule.");
    params->validateParametersAndSetDefaults(validParameters);
  }

  // Extract some information from that ParameterList.
  string modelId(params->get<string>("Model ID")),
    prefix(params->get<string>("Prefix")),
    basisType(params->get<string>("Basis Type"));
  int basisOrder(params->get<int>("Basis Order")),
    integrationOrder(params->get<int>("Integration Order"));

  // Setup the degrees of freedom (along with any necessary gradients) and the
  // closure models.
  {
    dofName_ = prefix + "U";
    this->addDOF(dofName_, basisType, basisOrder, integrationOrder);
    this->addDOFGrad(dofName_);
  }
  this->addClosureModel(modelId);
  this->setupDOFs();
} // end of Default Constructor

///////////////////////////////////////////////////////////////////////////////
//
//  buildAndRegisterEquationSetEvaluators()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
void
MyEquationSet<EvalT>::
buildAndRegisterEquationSetEvaluators(
  PHX::FieldManager<panzer::Traits>& fm,
  const panzer::FieldLibrary&        /* fl       */,
  const Teuchos::ParameterList&      /* userData */) const
{
  using panzer::BasisIRLayout;
  using panzer::IntegrationRule;
  using panzer::Integrator_BasisTimesScalar;
  using panzer::Integrator_GradBasisDotVector;
  using panzer::Traits;
  using PHX::Evaluator;
  using std::string;
  using std::vector;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Get the integration rule and basis layout corresponding to the degree of
  // freedom.
  RCP<IntegrationRule> ir    = this->getIntRuleForDOF(dofName_);
  RCP<BasisIRLayout>   basis = this->getBasisIRLayoutForDOF(dofName_);

  // Create the Laplacian term:  \int_\Omega \nabla u \cdot \nabla v\,d\Omega.
  const string laplacianName("RESIDUAL_" + dofName_ + "_LAPLACIAN");
  {
    ParameterList p;
    p.set("Residual Name", laplacianName     );
    p.set("Flux Name",     "GRAD_" + dofName_);
    p.set("IR",            ir                );
    p.set("Basis",         basis             );
    p.set("Multiplier",    1.0               );
    RCP<Evaluator<Traits>> op =
      rcp(new Integrator_GradBasisDotVector<EvalT, Traits>(p));
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // Create the Helmholtz term:  (1 - 8\pi^2) \int_\Omega uv\,d\Omega.
  const string helmholtzName("RESIDUAL_" + dofName_ + "_HELMHOLTZ");
  {
    ParameterList p;
    p.set("Residual Name", helmholtzName            );
    p.set("Value Name",    dofName_                 );
    p.set("IR",            ir                       );
    p.set("Basis",         basis                    );
    p.set("Multiplier",    (1.0 - 8.0 * M_PI * M_PI));
    RCP<Evaluator<Traits>> op =
      rcp(new Integrator_BasisTimesScalar<EvalT, Traits>(p));
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // Create the source term:  -\int_\Omega \sin(2\pi x)\sin(2\pi y)\,d\Omega.
  const string sourceName("RESIDUAL_" + dofName_ + "_SOURCE");
  {
    ParameterList p;
    p.set("Residual Name", sourceName          );
    p.set("Value Name",    dofName_ + "_SOURCE");
    p.set("IR",            ir                  );
    p.set("Basis",         basis               );
    p.set("Multiplier",    -1.0                );
    RCP<Evaluator<Traits>> op =
      rcp(new Integrator_BasisTimesScalar<EvalT, Traits>(p));
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // Use a sum operator to form the overall residual for the equation.
  {
    vector<string> residualOperatorNames{laplacianName, helmholtzName,
      sourceName};
    this->buildAndRegisterResidualSummationEvaluator(fm, dofName_,
      residualOperatorNames);
  }
} // end of buildAndRegisterEquationSetEvaluators()

#endif // __myEquationSetImpl_hpp__
