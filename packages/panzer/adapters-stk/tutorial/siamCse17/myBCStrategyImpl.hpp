// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Panzer
#include "Panzer_Constant.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_PureBasis.hpp"

// Phalanx
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Phalanx_MDField.hpp"

// Teuchos
#include "Teuchos_Assert.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

///////////////////////////////////////////////////////////////////////////////
//
//  Default Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT>
MyBCStrategy<EvalT>::
MyBCStrategy(
  const panzer::BC&                       bc,
  const Teuchos::RCP<panzer::GlobalData>& globalData)
  :
  panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>(bc, globalData)
{
  // Ensure that the "Strategy" from the input XML file matches this boundary
  // condition.
  TEUCHOS_ASSERT(this->m_bc.strategy() == "MyBCStrategy");
} // end of Default Constructor

///////////////////////////////////////////////////////////////////////////////
//
//  setup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT>
void
MyBCStrategy<EvalT>::
setup(
  const panzer::PhysicsBlock&   sidePB,
  const Teuchos::ParameterList& /* userData */)
{
  using   std::runtime_error;
  using   Teuchos::is_null;
  typedef std::vector<std::pair<std::string, Teuchos::RCP<panzer::PureBasis>>>
    DofVec;
  typedef DofVec::const_iterator DofVecIter;

  // Alert the Panzer library of a degree of freedom that is required by this
  // boundary condition.
  this->addDOF(this->m_bc.equationSetName());
  this->addTarget("MyBCStrategy_" + this->m_bc.equationSetName(),
    this->m_bc.equationSetName(), "Residual_" + this->m_bc.identifier());

  // Find the basis for this degree of freedom.
  const DofVec& dofs = sidePB.getProvidedDOFs();
  for (DofVecIter dof = dofs.begin(); dof != dofs.end(); ++dof)
  {
    if (dof->first == this->m_bc.equationSetName())
      this->basis = dof->second;
  } // end loop over the degrees of freedom

  // Ensure that we did in fact find a basis.
  TEUCHOS_TEST_FOR_EXCEPTION(is_null(this->basis), runtime_error, "Error:  "  \
    "The name \"" << this->m_bc.equationSetName() << "\" is not a valid DOF " \
    "for the boundary condition:\n" << this->m_bc << "\n");
} // end of setup()

///////////////////////////////////////////////////////////////////////////////
//
//  buildAndRegisterEvaluators()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT>
void
MyBCStrategy<EvalT>::
buildAndRegisterEvaluators(
  PHX::FieldManager<panzer::Traits>&                                 fm,
  const panzer::PhysicsBlock&                                        /* pb  */,
  const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& /* fac */,
  const Teuchos::ParameterList&                                      /* mod */,
  const Teuchos::ParameterList&                                      /* ud  */)
  const
{
  using panzer::Constant;
  using panzer::Traits;
  using PHX::Evaluator;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Create a Constant evaluator with a value of zero for the boundary
  // condition.
  ParameterList p("BC Constant Dirichlet");
  p.set("Name",        "MyBCStrategy_" + this->m_bc.equationSetName());
  p.set("Data Layout", basis->functional                             );
  p.set("Value",       0.0                                           );
  RCP<Evaluator<Traits>> op = rcp(new Constant<EvalT, Traits>(p));
  this->template registerEvaluator<EvalT>(fm, op);
} // end of buildAndRegisterEvaluators()
