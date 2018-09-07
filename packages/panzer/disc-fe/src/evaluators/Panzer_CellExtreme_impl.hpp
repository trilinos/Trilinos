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

#ifndef __Panzer_CellExtreme_impl_hpp__
#define __Panzer_CellExtreme_impl_hpp__

#include <limits>

#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
CellExtreme<EvalT, Traits>::
CellExtreme(
  const Teuchos::ParameterList& p) : quad_index(-1)
{
  Teuchos::RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  // setup default value
  use_max = true;
  if(p.isType<bool>("Use Max")) 
    use_max = p.get<bool>("Use Max");

  Teuchos::RCP<panzer::IntegrationRule> ir = p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR");
  quad_order = ir->cubature_degree;

  Teuchos::RCP<PHX::DataLayout> dl_cell = Teuchos::rcp(new PHX::MDALayout<Cell>(ir->dl_scalar->extent(0)));
  extreme = PHX::MDField<ScalarT>( p.get<std::string>("Extreme Name"), dl_cell);
  scalar = PHX::MDField<const ScalarT,Cell,IP>( p.get<std::string>("Field Name"), ir->dl_scalar);

  this->addEvaluatedField(extreme);
  this->addDependentField(scalar);
    
  multiplier = 1.0;
  if(p.isType<double>("Multiplier"))
     multiplier = p.get<double>("Multiplier");

  if (p.isType<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers")) {
    const std::vector<std::string>& field_multiplier_names = 
      *(p.get<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers"));

    for (std::vector<std::string>::const_iterator name = 
	   field_multiplier_names.begin(); 
	 name != field_multiplier_names.end(); ++name) {
      PHX::MDField<const ScalarT,Cell,IP> tmp_field(*name, p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
      field_multipliers.push_back(tmp_field);
    }
  }

  for (typename std::vector<PHX::MDField<const ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->addDependentField(*field);

  std::string n = "CellExtreme: " + extreme.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
CellExtreme<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  num_qp = scalar.extent(1);
  quad_index =  panzer::getIntegrationRuleIndex(quad_order,(*sd.worksets_)[0], this->wda);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
CellExtreme<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 
  for (index_t cell = 0; cell < workset.num_cells; ++cell) {
    
    for (std::size_t qp = 0; qp < num_qp; ++qp) {
      ScalarT current= multiplier * scalar(cell,qp);
      for (typename std::vector<PHX::MDField<const ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
	   field != field_multipliers.end(); ++field)
        current *= (*field)(cell,qp);  

      // take first quad point value
      if(qp==0)
        extreme(cell) = current;

      // take largest value in the cell
      if(use_max)
        extreme(cell) = extreme(cell)<current ? current : extreme(cell);
      else // use_min
        extreme(cell) = extreme(cell)>current ? current : extreme(cell);
    }
  }
}

//**********************************************************************
template<typename EvalT, typename TRAITS>
Teuchos::RCP<Teuchos::ParameterList> 
CellExtreme<EvalT, TRAITS>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Extreme Name", "?");
  p->set<std::string>("Field Name", "?");

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);
  p->set<bool>("Use Max", true);
  p->set<double>("Multiplier", 1.0);

  Teuchos::RCP<const std::vector<std::string> > fms;
  p->set("Field Multipliers", fms);
  return p;
}

//**********************************************************************

}

#endif
