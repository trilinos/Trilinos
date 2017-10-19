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

#ifndef PANZER_BasisValues_Evaluator_IMPL_HPP
#define PANZER_BasisValues_Evaluator_IMPL_HPP

#include <algorithm>
#include "Panzer_PointRule.hpp"
#include "Panzer_Workset_Utilities.hpp"

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(BasisValues_Evaluator,p)
  : derivativesRequired_(true)
{
  Teuchos::RCP<const panzer::PointRule> pointRule 
     = p.get< Teuchos::RCP<const panzer::PointRule> >("Point Rule");
  Teuchos::RCP<const panzer::PureBasis> inBasis
     = p.get<Teuchos::RCP<const panzer::PureBasis> >("Basis");

  bool derivativesRequired = true;
  if(p.isType<bool>("Derivatives Required"))
    derivativesRequired = p.get<bool>("Derivatives Required");

  initialize(pointRule,inBasis,derivativesRequired);
}

//**********************************************************************
template <typename EvalT, typename TRAITST>
BasisValues_Evaluator<EvalT,TRAITST>::BasisValues_Evaluator(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                                                            const Teuchos::RCP<const panzer::PureBasis> & inBasis)
  : derivativesRequired_(true)
{
  bool derivativesRequired = true;
  initialize(pointRule,inBasis,derivativesRequired);
}

//**********************************************************************
template <typename EvalT, typename TRAITST>
BasisValues_Evaluator<EvalT,TRAITST>::BasisValues_Evaluator(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                                                            const Teuchos::RCP<const panzer::PureBasis> & inBasis,
                                                            bool derivativesRequired)
  : derivativesRequired_(true)
{
  initialize(pointRule,inBasis,derivativesRequired);
}

//**********************************************************************
template <typename EvalT, typename TRAITST>
void BasisValues_Evaluator<EvalT,TRAITST>::initialize(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                                                      const Teuchos::RCP<const panzer::PureBasis> & inBasis,
                                                      bool derivativesRequired)
{
  basis = inBasis;
  derivativesRequired_ = derivativesRequired;

  int space_dim = basis->dimension();

  // setup all fields to be evaluated and constructed
  pointValues = PointValues2<ScalarT>(pointRule->getName()+"_",false);
  pointValues.setupArrays(pointRule);

  // the field manager will allocate all of these field
  {
    constPointValues = pointValues;
    this->addDependentField(constPointValues.coords_ref);
    this->addDependentField(constPointValues.jac);
    this->addDependentField(constPointValues.jac_inv);
    this->addDependentField(constPointValues.jac_det);
  }

  // setup all fields to be evaluated and constructed
  Teuchos::RCP<panzer::BasisIRLayout> layout = Teuchos::rcp(new panzer::BasisIRLayout(basis,*pointRule));
  basisValues = Teuchos::rcp(new BasisValues2<ScalarT>(basis->name()+"_"+pointRule->getName()+"_",false));
  basisValues->setupArrays(layout,derivativesRequired_);

  // the field manager will allocate all of these field

  if(basis->getElementSpace()==panzer::PureBasis::HGRAD) {
    this->addEvaluatedField(basisValues->basis_ref_scalar);      
    this->addEvaluatedField(basisValues->basis_scalar);           

    if(derivativesRequired) {
      this->addEvaluatedField(basisValues->grad_basis_ref);   
      this->addEvaluatedField(basisValues->grad_basis);        
    }
  }

  if(basis->getElementSpace()==panzer::PureBasis::HCURL) {
    this->addEvaluatedField(basisValues->basis_ref_vector);      
    this->addEvaluatedField(basisValues->basis_vector);           

    if(derivativesRequired && space_dim==2) {
      this->addEvaluatedField(basisValues->curl_basis_ref_scalar);     
      this->addEvaluatedField(basisValues->curl_basis_scalar);          
    }
    else if(derivativesRequired && space_dim==3) {
      this->addEvaluatedField(basisValues->curl_basis_ref_vector);     
      this->addEvaluatedField(basisValues->curl_basis_vector);          
    }
  }

  if(basis->getElementSpace()==panzer::PureBasis::HDIV) {
    this->addEvaluatedField(basisValues->basis_ref_vector);      
    this->addEvaluatedField(basisValues->basis_vector);           

    if(derivativesRequired) {
      this->addEvaluatedField(basisValues->div_basis_ref);     
      this->addEvaluatedField(basisValues->div_basis);          
    }
  }

  std::string n = "BasisValues_Evaluator: " +basis->name() + "_" + pointRule->getName();
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(BasisValues_Evaluator, sd, fm)
{
  int space_dim = basis->dimension();

  orientations = sd.orientations_;

  basisValues->setExtendedDimensions(fm.template getKokkosExtendedDataTypeDimensions<EvalT>());

  // setup the pointers for the point values data structure
  this->utils.setFieldData(pointValues.coords_ref,fm);
  this->utils.setFieldData(pointValues.jac,fm);
  this->utils.setFieldData(pointValues.jac_inv,fm);
  this->utils.setFieldData(pointValues.jac_det,fm);

  // setup the pointers for the basis values data structure

  if(basis->getElementSpace()==panzer::PureBasis::HGRAD) {
    this->utils.setFieldData(basisValues->basis_ref_scalar,fm);      
    this->utils.setFieldData(basisValues->basis_scalar,fm);           

    if(derivativesRequired_) {
      this->utils.setFieldData(basisValues->grad_basis_ref,fm);   
      this->utils.setFieldData(basisValues->grad_basis,fm);        
    }
  }

  if(basis->getElementSpace()==panzer::PureBasis::HCURL) {
    this->utils.setFieldData(basisValues->basis_ref_vector,fm);      
    this->utils.setFieldData(basisValues->basis_vector,fm);           

    if(derivativesRequired_ && space_dim==2) {
      this->utils.setFieldData(basisValues->curl_basis_ref_scalar,fm);     
      this->utils.setFieldData(basisValues->curl_basis_scalar,fm);          
    }
    else if(derivativesRequired_ && space_dim==3) {
      this->utils.setFieldData(basisValues->curl_basis_ref_vector,fm);     
      this->utils.setFieldData(basisValues->curl_basis_vector,fm);          
    }
  }

  if(basis->getElementSpace()==panzer::PureBasis::HDIV) {
    this->utils.setFieldData(basisValues->basis_ref_vector,fm);      
    this->utils.setFieldData(basisValues->basis_vector,fm);           

    if(derivativesRequired_) {
      this->utils.setFieldData(basisValues->div_basis_ref,fm);     
      this->utils.setFieldData(basisValues->div_basis,fm);          
    }
  }

}

//**********************************************************************
PHX_EVALUATE_FIELDS(BasisValues_Evaluator, workset)
{ 
  // evaluate the point values (construct jacobians etc...)
  basisValues->evaluateValues(pointValues.coords_ref,
                              pointValues.jac,
                              pointValues.jac_det,
                              pointValues.jac_inv);

  // this can be over-ridden in basisValues e.g., DG element setting
  if(basis->requiresOrientations()) {
    const WorksetDetails & details = workset;
    
    std::vector<Intrepid2::Orientation> ortPerWorkset;
    for (index_t c=0;c<workset.num_cells;++c)
      ortPerWorkset.push_back((*orientations)[details.cell_local_ids[c]]);
    
    basisValues->applyOrientations(ortPerWorkset);
  }
}

//**********************************************************************

}

#endif
