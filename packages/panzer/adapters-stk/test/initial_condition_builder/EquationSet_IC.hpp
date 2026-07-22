// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __EquationSet_IC_hpp__
#define __EquationSet_IC_hpp__

#include <vector>
#include <string>

#include "Panzer_EquationSet_DefaultImpl.hpp"
#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_Factory_Defines.hpp"

template <typename EvalT>
class EquationSet_IC : public panzer::EquationSet_DefaultImpl<EvalT> {
public:    

   /** In the constructor you set all the fields provided by this
     * equation set. 
     */
   EquationSet_IC(const Teuchos::RCP<Teuchos::ParameterList>& params,
                  const int& default_integration_order,
                  const panzer::CellData& cell_data,
                  const Teuchos::RCP<panzer::GlobalData>& global_data,
                  const bool build_transient_support);
    
   /** The specific evaluators are registered with the field manager argument.
     */
   void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                              const panzer::FieldLibrary& field_library,
                                              const Teuchos::ParameterList& user_data) const {}

};

// ***********************************************************************
template <typename EvalT>
EquationSet_IC<EvalT>::
EquationSet_IC(const Teuchos::RCP<Teuchos::ParameterList>& params,
               const int& default_integration_order,
               const panzer::CellData& cell_data,
               const Teuchos::RCP<panzer::GlobalData>& global_data,
               const bool build_transient_support) :
  panzer::EquationSet_DefaultImpl<EvalT>(params, default_integration_order, cell_data, global_data, build_transient_support )
{
  // ********************
  // Validate and parse parameter list
  // ********************
  Teuchos::ParameterList valid_parameters_sublist;
  valid_parameters_sublist.set("Basis Type","HGrad","Type of Basis to use");
  valid_parameters_sublist.set("Basis Order",1,"Order of the basis");

  for(auto itr=params->begin();itr!=params->end();++itr) {
     
    const std::string field = params->name(itr);
    const Teuchos::ParameterEntry & entry = params->entry(itr);

    // only process lists
    if(!entry.isList()) continue;

    Teuchos::ParameterList & basisPL = entry.getValue((Teuchos::ParameterList *) 0);
    basisPL.validateParametersAndSetDefaults(valid_parameters_sublist);

    std::string basis_type = basisPL.get<std::string>("Basis Type");
    int basis_order = basisPL.get<int>("Basis Order");

    this->addDOF(field,basis_type,basis_order,default_integration_order);
  }
  
  this->addClosureModel("");

  this->setupDOFs();
}

PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(EquationSet_IC, EquationSet_IC)

// A user written factory that creates each equation set.  The key member here
// is buildEquationSet
class IC_EquationSetFactory : public panzer::EquationSetFactory {
public:

   Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> >
   buildEquationSet(const Teuchos::RCP<Teuchos::ParameterList>& params,
		    const int& default_integration_order,
                    const panzer::CellData& cell_data,
		    const Teuchos::RCP<panzer::GlobalData>& global_data,
                    const bool build_transient_support) const
   {
      Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > eq_set= 
         Teuchos::rcp(new panzer::EquationSet_TemplateManager<panzer::Traits>);
         
      bool found = false; // this is used by PANZER_BUILD_EQSET_OBJECTS
         
      // macro checks if(ies.name=="Poisson") then an EquationSet_Energy object is constructed
      PANZER_BUILD_EQSET_OBJECTS("IC", EquationSet_IC)
         
      // make sure your equation set has been found
      if(!found) {
	std::string msg = "Error - the \"Equation Set\" called \"" + params->get<std::string>("Type") +
                           "\" is not a valid equation set identifier. Please supply the correct factory.\n";
         TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
      }
         
      return eq_set;
   }
    
};

#endif
