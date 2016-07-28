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

#ifndef __Panzer_Integerator_BasisTimesVector_impl_hpp__
#define __Panzer_Integerator_BasisTimesVector_impl_hpp__

#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_CommonArrayFactories.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"

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

  this->addEvaluatedField(residual);
  this->addDependentField(vectorField);
    
  multiplier = p.get<double>("Multiplier");


  if (p.isType<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers")) {
    const std::vector<std::string>& field_multiplier_names = 
      *(p.get<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers"));

    int i=0;
    field_multipliers.resize(field_multiplier_names.size());
    kokkos_field_multipliers = Kokkos::View<Kokkos::View<ScalarT**>* >("BasisTimesVector::KokkosFieldMultipliers", field_multiplier_names.size());
    Kokkos::fence();
    for (std::vector<std::string>::const_iterator name = 
           field_multiplier_names.begin(); 
         name != field_multiplier_names.end(); ++name) {
      PHX::MDField<ScalarT,Cell,IP> tmp_field(*name, p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
      Kokkos::fence();
      field_multipliers[i++] = tmp_field;
      this->addDependentField(tmp_field);
    }
  }

  std::string n = "Integrator_BasisTimesVector: " + residual.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(Integrator_BasisTimesVector,sd,fm)
{
  this->utils.setFieldData(residual,fm);
  this->utils.setFieldData(vectorField,fm);
  // this->utils.setFieldData(dof_orientation,fm);
  
  for (std::size_t i=0; i<field_multipliers.size(); ++i) {
    this->utils.setFieldData(field_multipliers[i],fm);
    kokkos_field_multipliers(i) = field_multipliers[i].get_static_view();
  }

  basis_card = residual.dimension(1); // basis cardinality
  num_qp = vectorField.dimension(1); 
  num_dim = vectorField.dimension(2); // dimension of a vector

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0], this->wda);

  // tmp = Intrepid2:FieldContainer<ScalarT>(vectorField.dimension(0), num_qp, num_dim);
  MDFieldArrayFactory af("",fm.template getKokkosExtendedDataTypeDimensions<EvalT>(),true);
}


//**********************************************************************
template<typename EvalT, typename TRAITS>
template<int NUM_FIELD_MULT>
KOKKOS_INLINE_FUNCTION
void Integrator_BasisTimesVector<EvalT, TRAITS>::operator()(const FieldMultTag<NUM_FIELD_MULT> &, const size_t &cell) const {
  const int nqp = vectorField.extent_int(1), ndim = vectorField.extent_int(2);
  const int nfm = kokkos_field_multipliers.extent_int(0);
  const int nbf = weighted_basis_vector.extent_int(1);

  for (int lbf = 0; lbf < nbf; lbf++)
    residual(cell,lbf) = 0.0;

  ScalarT tmp, fmm=1;
  if ( NUM_FIELD_MULT == 0 ){
    for (int qp = 0; qp < nqp; ++qp) {
      for (int d = 0; d < ndim; ++d) {
        tmp = multiplier * vectorField(cell,qp,d);
        for (int lbf = 0; lbf < nbf; lbf++)
          residual(cell,lbf) += weighted_basis_vector(cell, lbf, qp, d)*tmp;
      }
    }
  } else if ( NUM_FIELD_MULT == 1 ){
    for (int qp = 0; qp < nqp; ++qp) {
      for (int d = 0; d < ndim; ++d) {
        tmp = multiplier * vectorField(cell,qp,d)*kokkos_field_multipliers(0)(cell,qp);
        for (int lbf = 0; lbf < nbf; lbf++)
          residual(cell,lbf) += weighted_basis_vector(cell, lbf, qp, d)*tmp;
      }
    }
  } else {
    for (int qp = 0; qp < nqp; ++qp) {
      for (int i = 0; i < nfm; ++i)
        fmm *= kokkos_field_multipliers(i)(cell,qp);
      for (int d = 0; d < ndim; ++d) {
        tmp = multiplier * vectorField(cell,qp,d)*fmm;
        for (int lbf = 0; lbf < nbf; lbf++)
          residual(cell,lbf) += weighted_basis_vector(cell, lbf, qp, d)*tmp;
      }
    }
  }
}

//**********************************************************************
PHX_EVALUATE_FIELDS(Integrator_BasisTimesVector,workset)
{ 
  // residual.deep_copy(ScalarT(0.0));
  weighted_basis_vector = this->wda(workset).bases[basis_index]->weighted_basis_vector;
  if ( field_multipliers.size() == 0)
    Kokkos::parallel_for(Kokkos::RangePolicy<FieldMultTag<0> >(0, workset.num_cells),*this);
  else if ( field_multipliers.size() == 1)
    Kokkos::parallel_for(Kokkos::RangePolicy<FieldMultTag<1> >(0, workset.num_cells),*this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<FieldMultTag<-1> >(0, workset.num_cells),*this);
}

//**********************************************************************
template<typename EvalT, typename TRAITS>
Teuchos::RCP<Teuchos::ParameterList> 
Integrator_BasisTimesVector<EvalT, TRAITS>::getValidParameters() const
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

