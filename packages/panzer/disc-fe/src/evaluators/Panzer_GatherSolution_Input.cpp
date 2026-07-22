// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_GatherSolution_Input.hpp"

namespace panzer {

GatherSolution_Input::
GatherSolution_Input()
{
  setParameterList(Teuchos::rcp(new Teuchos::ParameterList(*getValidParameters())));
}

void 
GatherSolution_Input::
setParameterList(const Teuchos::ParameterList & p)
{
  setParameterList(Teuchos::rcp(new Teuchos::ParameterList(p)));
}

void 
GatherSolution_Input::
setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & p)
{
  using Teuchos::RCP;

  // correction for non const basis
  if(p->isType<RCP<PureBasis> >("Basis")) {
    RCP<const PureBasis> basis = p->get<RCP<PureBasis> >("Basis");
    p->remove("Basis");
    p->set("Basis",basis);
  }

  // set complete state
  p->validateParametersAndSetDefaults(*getValidParameters());

  // for posterity save the modified list
  setMyParamList(p);

  dofNames_            = *p->get<RCP< std::vector<std::string> > >("DOF Names");
  indexerNames_        = *p->get<RCP< std::vector<std::string> > >("Indexer Names");
  useTimeDerivSolnVec_ = p->get<bool>("Use Time Derivative Solution Vector");
  globalDataKey_       = p->get<std::string>("Global Data Key");
  basis_               = p->get<RCP<const panzer::PureBasis> >("Basis");

  // required by Tangent types
  tangentNames_ = *p->get<RCP<std::vector<std::vector<std::string> > > >("Tangent Names"); 

  // required by Jacobian types
  sensName_        = p->get<std::string>("Sensitivities Name");
  gatherSeedIndex_ = p->get<int>("Gather Seed Index");
  firstSensAvail_  = p->get<bool>("First Sensitivities Available");

  // required by Hessian types
  secondSensAvail_         = p->get<bool>("Second Sensitivities Available");    
  secondSensDataKeyPrefix_ = p->get<std::string>("Second Sensitivities Data Key Prefix");
}
  
Teuchos::RCP<const Teuchos::ParameterList> 
GatherSolution_Input::
getValidParameters() const
{
  using Teuchos::RCP;

  RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  RCP<std::vector<std::string> > emptyList = Teuchos::rcp(new std::vector<std::string>);
  RCP<std::vector<std::vector<std::string> > > emptyDblList = Teuchos::rcp(new std::vector<std::vector<std::string> >);

  // required by Residual types (and all others)
  p->set<RCP< std::vector<std::string> > >("DOF Names",emptyList);
  p->set<RCP< std::vector<std::string> > >("Indexer Names",emptyList);
  p->set<RCP<const panzer::PureBasis> >("Basis",Teuchos::null);
  p->set<bool>("Use Time Derivative Solution Vector",false);
  p->get<std::string>("Global Data Key","Solution Gather Container");

  // required by Tangent types
  p->set<RCP< std::vector<std::vector<std::string> > > >("Tangent Names",emptyDblList); // only Tangent

  // required by Jacobian types
  p->set<std::string>("Sensitivities Name","");        // Hessian and Jacobian
  p->set<bool>("First Sensitivities Available",true);  // Hessian and Jacobian
  p->set<int>("Gather Seed Index",-1);                 // Hessian and Jacobian

  // required by Hessian types
  p->set<bool>("Second Sensitivities Available",true);          // Hessian only
  p->set<std::string>("Second Sensitivities Data Key Prefix","DELTA_"); // Hessian only

  return p;
}

}
