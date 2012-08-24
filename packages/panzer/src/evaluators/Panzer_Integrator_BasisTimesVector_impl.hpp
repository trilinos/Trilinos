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

#ifndef PANZER_EVALUATOR_BASISTIMESSCALAR_IMPL_HPP
#define PANZER_EVALUATOR_BASISTIMESSCALAR_IMPL_HPP

#include "Intrepid_FunctionSpaceTools.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset_Utilities.hpp"

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(Integrator_BasisTimesVector,p) :
  residual( p.get<std::string>("Residual Name"), 
	    p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->functional),
  vectorField( p.get<std::string>("Value Name"), 
 	       p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_vector),
  basis_name(p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->name())
{
  Teuchos::RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  Teuchos::RCP<const PureBasis> basis 
     = p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->getBasis();

  // Verify that this basis supports the gradient operation
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->isVectorBasis(),std::logic_error,
                             "Integrator_BasisTimesVector: Basis of type \"" << basis->name() << "\" is not "
                             "a vector basis.");
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->requiresOrientations(),std::logic_error,
                             "Integrator_BasisTimesVector: Basis of type \"" << basis->name() << "\" does not "
                             "require orientation. This seems very strange, so I'm failing.");

  dof_orientation = PHX::MDField<ScalarT,Cell,BASIS>(p.get<std::string>("Test Field Name")+" Orientation", 
                                                     basis->functional);

  this->addEvaluatedField(residual);
  this->addDependentField(vectorField);
  this->addDependentField(dof_orientation);
    
  multiplier = p.get<double>("Multiplier");


  if (p.isType<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers")) {
    const std::vector<std::string>& field_multiplier_names = 
      *(p.get<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers"));

    for (std::vector<std::string>::const_iterator name = 
	   field_multiplier_names.begin(); 
	 name != field_multiplier_names.end(); ++name) {
      PHX::MDField<ScalarT,Cell,IP> tmp_field(*name, p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
      field_multipliers.push_back(tmp_field);
    }
  }

  for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->addDependentField(*field);

  std::string n = "Integrator_BasisTimesVector: " + residual.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(Integrator_BasisTimesVector,sd,fm)
{
  this->utils.setFieldData(residual,fm);
  this->utils.setFieldData(vectorField,fm);
  this->utils.setFieldData(dof_orientation,fm);
  
  for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->utils.setFieldData(*field,fm);

  basis_card = residual.dimension(1); // basis cardinality
  num_qp = vectorField.dimension(1); 
  num_dim = vectorField.dimension(2); // dimension of a vector

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0]);

  tmp = Intrepid::FieldContainer<ScalarT>(vectorField.dimension(0), num_qp, num_dim);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(Integrator_BasisTimesVector,workset)
{ 
  for (int i=0; i < residual.size(); ++i)
    residual[i] = 0.0;

  for (std::size_t cell = 0; cell < workset.num_cells; ++cell) {
    for (std::size_t qp = 0; qp < num_qp; ++qp) {
      for (std::size_t d = 0; d < num_dim; ++d) {
        tmp(cell,qp,d) = multiplier * vectorField(cell,qp,d);
        for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
 	     field != field_multipliers.end(); ++field)
	  tmp(cell,qp,d) *= (*field)(cell,qp);  
      }
    }
  }

  if(workset.num_cells>0) {
    Intrepid::FieldContainer<double> weighted_basis = (workset.bases[basis_index])->weighted_basis;

    // assign ScalarT "dof_orientation" to double "orientation"
    Intrepid::FieldContainer<double> orientation(dof_orientation.dimension(0),
                                                 dof_orientation.dimension(1));
    for(int i=0;i<dof_orientation.dimension(0);i++)
       for(int j=0;j<dof_orientation.dimension(1);j++)
          orientation(i,j) = Sacado::ScalarValue<ScalarT>::eval(dof_orientation(i,j));

    // make sure things are orientated correctly
    Intrepid::FunctionSpaceTools::
       applyFieldSigns<ScalarT>(weighted_basis,orientation);

    // evaluate at quadrature points
    Intrepid::FunctionSpaceTools::
      integrate<ScalarT>(residual, tmp, 
  		            weighted_basis,
		            Intrepid::COMP_BLAS);
  }
}

//**********************************************************************
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList> 
Integrator_BasisTimesVector<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Residual Name", "?");
  p->set<std::string>("Value Name", "?");
  p->set<std::string>("Test Field Name", "?");
  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);
  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);
  p->set<double>("Multiplier", 1.0);
  Teuchos::RCP<const std::vector<std::string> > fms;
  p->set("Field Multipliers", fms);
  return p;
}

//**********************************************************************

}

#endif

