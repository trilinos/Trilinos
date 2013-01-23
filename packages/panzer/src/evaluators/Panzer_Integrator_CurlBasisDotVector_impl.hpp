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

#ifndef PANZER_EVALUATOR_CURLBASISDOTVECTOR_IMPL_HPP
#define PANZER_EVALUATOR_CURLBASISDOTVECTOR_IMPL_HPP

#include "Intrepid_FunctionSpaceTools.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset_Utilities.hpp"

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(Integrator_CurlBasisDotVector,p) :
  residual( p.get<std::string>("Residual Name"), 
	    p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->functional),
  basis_name(p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->name())
{
  Teuchos::RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  Teuchos::RCP<const PureBasis> basis 
     = p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->getBasis();

  // Verify that this basis supports the curl operation
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->supportsCurl(),std::logic_error,
                             "Integrator_CurlBasisDotVector: Basis of type \"" << basis->name() << "\" does not support CURL.");
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->requiresOrientations(),std::logic_error,
                             "Integration_CurlBasisDotVector: Basis of type \"" << basis->name() << "\" should require orientations. So we are throwing.");
  TEUCHOS_TEST_FOR_EXCEPTION(!(basis->dimension()==2 || basis->dimension()==3),std::logic_error,
                             "Integrator_CurlBasisDotVector: Evaluator requires 2D or 3D basis types, the basis \"" << basis->name() << "\" is neither.");

  // use a scalar field only if dimension is 2D
  useScalarField = (basis->dimension()==2);
  
  // determine if using scalar field for curl or a vector field (2D versus 3D)
  if(!useScalarField)
     flux = PHX::MDField<ScalarT>( p.get<std::string>("Value Name"), 
	                                  p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_vector );
  else 
     flux = PHX::MDField<ScalarT>( p.get<std::string>("Value Name"), 
   	                                  p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar );

  // build dof_orientation
  dof_orientation = PHX::MDField<ScalarT,Cell,BASIS>(p.get<std::string>("Test Field Name")+" Orientation", 
                                                     basis->functional);


  this->addEvaluatedField(residual);
  this->addDependentField(flux);
  this->addDependentField(dof_orientation);
  
  multiplier = p.get<double>("Multiplier");
  if (p.isType<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers")) 
  {
    const std::vector<std::string>& field_multiplier_names = 
      *(p.get<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers"));

    for (std::vector<std::string>::const_iterator name = field_multiplier_names.begin(); 
      name != field_multiplier_names.end(); ++name) 
    {
      PHX::MDField<ScalarT,Cell,IP> tmp_field(*name, p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
      field_multipliers.push_back(tmp_field);
    }
  }

  for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->addDependentField(*field);

  std::string n = 
    "Integrator_CurlBasisDotVector: " + residual.fieldTag().name();

  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(Integrator_CurlBasisDotVector,sd,fm)
{
  this->utils.setFieldData(residual,fm);
  this->utils.setFieldData(flux,fm);
  this->utils.setFieldData(dof_orientation,fm);

  for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->utils.setFieldData(*field,fm);

  num_nodes = residual.dimension(1);
  num_qp = flux.dimension(1);
  num_dim = (useScalarField ? 2 : 3);  // this only works in 2D or 3D

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0]);

  if(!useScalarField)
     tmp = Intrepid::FieldContainer<ScalarT>(flux.dimension(0), num_qp, num_dim); 
  else
     tmp = Intrepid::FieldContainer<ScalarT>(flux.dimension(0), num_qp); 
}

//**********************************************************************
PHX_EVALUATE_FIELDS(Integrator_CurlBasisDotVector,workset)
{ 
  for (int i=0; i < residual.size(); ++i)
    residual[i] = 0.0;
  
  for (std::size_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (std::size_t qp = 0; qp < num_qp; ++qp)
    {
      ScalarT tmpVar = 1.0;
      for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
           field != field_multipliers.end(); ++field)
        tmpVar = tmpVar * (*field)(cell,qp);  

      if(!useScalarField) {
        // for vector fields loop over dimension
        for (std::size_t dim = 0; dim < num_dim; ++dim)
          tmp(cell,qp,dim) = multiplier * tmpVar * flux(cell,qp,dim);
      }
      else {
        // no dimension to loop over for scalar fields
        tmp(cell,qp) = multiplier * tmpVar * flux(cell,qp);
      }
    }
  }
  
  if(workset.num_cells>0) {
    Intrepid::FieldContainer<double> weighted_curl_basis = (workset.bases[basis_index])->weighted_curl_basis;

    // assign ScalarT "dof_orientation" to double "orientation"
    Intrepid::FieldContainer<double> orientation(dof_orientation.dimension(0),
                                                 dof_orientation.dimension(1));
    for(int i=0;i<dof_orientation.dimension(0);i++)
       for(int j=0;j<dof_orientation.dimension(1);j++)
          orientation(i,j) = Sacado::ScalarValue<ScalarT>::eval(dof_orientation(i,j));

    // make sure things are orientated correctly
    Intrepid::FunctionSpaceTools::
       applyFieldSigns<ScalarT>(weighted_curl_basis,orientation);

    // evaluate at quadrature points
     Intrepid::FunctionSpaceTools::
       integrate<ScalarT>(residual, tmp, 
                       weighted_curl_basis, 
		       Intrepid::COMP_BLAS);
  }
}

//**********************************************************************

template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList> 
Integrator_CurlBasisDotVector<EvalT, Traits>::getValidParameters() const
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
