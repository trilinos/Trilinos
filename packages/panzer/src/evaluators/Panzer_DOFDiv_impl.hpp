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

#ifndef PANZER_DOF_DIV_IMPL_HPP
#define PANZER_DOF_DIV_IMPL_HPP

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"

namespace panzer {

namespace {

//**********************************************************************
template<typename ScalarT,typename ArrayT>                   
void evaluateDiv_withSens(int numCells,
                          int basis_dimension,
                          PHX::MDField<ScalarT,Cell,IP> & dof_div, 
                          PHX::MDField<ScalarT,Cell,Point> & dof_value,
                          const ArrayT & div_basis)
{ 
  if(numCells>0) {
    // evaluate at quadrature points

    int numFields = div_basis.dimension(1);
    int numPoints = div_basis.dimension(2);

    for (int cell=0; cell<numCells; cell++) {
      for (int pt=0; pt<numPoints; pt++) {
        // first initialize to the right thing (prevents over writing with 0)
        // then loop over one less basis function
        // ScalarT & div = dof_div(cell,pt);
        dof_div(cell,pt) = dof_value(cell, 0) * div_basis(cell, 0, pt);
        for (int bf=1; bf<numFields; bf++)
          dof_div(cell,pt) += dof_value(cell, bf) * div_basis(cell, bf, pt);
      }
    }
  }
}

}

//**********************************************************************
// MOST EVALUATION TYPES
//**********************************************************************

//**********************************************************************
template<typename EvalT, typename TRAITS>                   
DOFDiv<EvalT, TRAITS>::
DOFDiv(const Teuchos::ParameterList & p) :
  dof_value( p.get<std::string>("Name"), 
	     p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->functional),
  basis_name(p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->name())
{
  Teuchos::RCP<const PureBasis> basis 
     = p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->getBasis();

  // Verify that this basis supports the div operation
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->supportsDiv(),std::logic_error,
                             "DOFDiv: Basis of type \"" << basis->name() << "\" does not support DIV");
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->requiresOrientations(),std::logic_error,
                             "DOFDiv: Basis of type \"" << basis->name() << "\" in DOF Div should require orientations. So we are throwing.");

  // build dof_div
  basis_dimension = basis->dimension();
  dof_div = PHX::MDField<ScalarT,Cell,IP>(p.get<std::string>("Div Name"), 
      	                                  p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar );

  // add to evaluation graph
  this->addEvaluatedField(dof_div);
  this->addDependentField(dof_value);
  
  std::string n = "DOFDiv: " + dof_div.fieldTag().name() + " ("+PHX::typeAsString<EvalT>()+")";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename TRAITS>                   
void DOFDiv<EvalT, TRAITS>::
postRegistrationSetup(typename TRAITS::SetupData sd,
                      PHX::FieldManager<TRAITS>& fm)
{
  this->utils.setFieldData(dof_value,fm);
  this->utils.setFieldData(dof_div,fm);

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0]);
}

//**********************************************************************
template<typename EvalT, typename TRAITS>                   
void DOFDiv<EvalT, TRAITS>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
  evaluateDiv_withSens<ScalarT>(workset.num_cells,basis_dimension,dof_div,dof_value,workset.bases[basis_index]->div_basis);
}

//**********************************************************************

//**********************************************************************
// JACOBIAN EVALUATION TYPES
//**********************************************************************

//**********************************************************************
template<typename TRAITS>                   
DOFDiv<panzer::Traits::Jacobian, TRAITS>::
DOFDiv(const Teuchos::ParameterList & p) :
  dof_value( p.get<std::string>("Name"), 
	     p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->functional),
  basis_name(p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->name())
{
  Teuchos::RCP<const PureBasis> basis 
     = p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->getBasis();

  // do you specialize because you know where the basis functions are and can
  // skip a large number of AD calculations?
  if(p.isType<Teuchos::RCP<const std::vector<int> > >("Jacobian Offsets Vector")) {
    offsets = *p.get<Teuchos::RCP<const std::vector<int> > >("Jacobian Offsets Vector");
    accelerate_jacobian = true;  // short cut for identity matrix
  }
  else
    accelerate_jacobian = false; // don't short cut for identity matrix
  accelerate_jacobian = false; // don't short cut for identity matrix

  // Verify that this basis supports the div operation
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->supportsDiv(),std::logic_error,
                             "DOFDiv: Basis of type \"" << basis->name() << "\" does not support DIV");
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->requiresOrientations(),std::logic_error,
                             "DOFDiv: Basis of type \"" << basis->name() << "\" in DOF Div should require orientations. So we are throwing.");

  // build dof_div
  basis_dimension = basis->dimension();
  dof_div = PHX::MDField<ScalarT,Cell,IP>(p.get<std::string>("Div Name"), 
      	                                  p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar );

  // add to evaluation graph
  this->addEvaluatedField(dof_div);
  this->addDependentField(dof_value);
  
  std::string n = "DOFDiv: " + dof_div.fieldTag().name() + " ("+PHX::typeAsString<panzer::Traits::Jacobian>()+")";
  this->setName(n);
}

//**********************************************************************
template<typename TRAITS>                   
void DOFDiv<panzer::Traits::Jacobian, TRAITS>::
postRegistrationSetup(typename TRAITS::SetupData sd,
                      PHX::FieldManager<TRAITS>& fm)
{
  this->utils.setFieldData(dof_value,fm);
  this->utils.setFieldData(dof_div,fm);

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0]);
}

template<typename TRAITS>                   
void DOFDiv<panzer::Traits::Jacobian,TRAITS>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
  if(!accelerate_jacobian) {
    // do the case where we use the AD types to determine the derivatives
    evaluateDiv_withSens(workset.num_cells,basis_dimension,dof_div,dof_value,workset.bases[basis_index]->div_basis);
    return;
  }

  TEUCHOS_ASSERT(false);
}

}

#endif
