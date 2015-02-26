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

#include "Intrepid_FunctionSpaceTools.hpp"
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
  // this->utils.setFieldData(dof_orientation,fm);
  
  for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->utils.setFieldData(*field,fm);

  basis_card = residual.dimension(1); // basis cardinality
  num_qp = vectorField.dimension(1); 
  num_dim = vectorField.dimension(2); // dimension of a vector

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0]);

  // tmp = Intrepid:FieldContainer<ScalarT>(vectorField.dimension(0), num_qp, num_dim);
  MDFieldArrayFactory af("",fm.template getKokkosExtendedDataTypeDimensions<EvalT>(),true);
  scratch = af.buildStaticArray<ScalarT,Cell,IP,Dim>("btv_scratch",vectorField.dimension(0),num_qp,num_dim);
}


//**********************************************************************

namespace {

template <typename ScalarT>
class FillScratch_Initialize {
public:
  typedef typename PHX::Device execution_space;
  double multiplier;
  PHX::MDField<ScalarT,Cell,IP,Dim> vectorField;
  PHX::MDField<ScalarT,Cell,IP,Dim> scratch;
  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned cell) const
  {
    for (std::size_t qp = 0; qp < scratch.dimension_1(); ++qp) 
      for (std::size_t d = 0; d < scratch.dimension_2(); ++d)
        scratch(cell,qp,d) = multiplier * vectorField(cell,qp,d);
  }
};

template <typename ScalarT>
class FillScratch_FieldMultipliers {
public:
  typedef typename PHX::Device execution_space;
  PHX::MDField<ScalarT,Cell,IP,Dim> scratch;
  PHX::MDField<ScalarT,Cell,IP> field;
  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned cell) const
  {
    for (std::size_t qp = 0; qp < scratch.dimension_1(); ++qp) 
      for (std::size_t d = 0; d < scratch.dimension_2(); ++d)
        scratch(cell,qp,d) *= field(cell,qp);  
  }
};

template <typename ScalarT>
class Integrate_Values {
public:
  typedef typename PHX::Device execution_space;
  PHX::MDField<ScalarT,Cell,IP,Dim> scratch;
  PHX::MDField<ScalarT,Cell,BASIS> residual;
  PHX::MDField<double,Cell,BASIS,IP,Dim> weighted_basis_vector;
  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned cell) const
  {
    for (int lbf = 0; lbf < weighted_basis_vector.dimension_1(); lbf++) {
      residual(cell,lbf) = 0.0;
      for (int qp = 0; qp < weighted_basis_vector.dimension_2(); qp++) {
        for (int iVec = 0; iVec < weighted_basis_vector.dimension_3(); iVec++) {
          residual(cell,lbf) += weighted_basis_vector(cell, lbf, qp, iVec)*scratch(cell, qp, iVec);
        } // D-loop
      } // P-loop
    } // F-loop
  }
};

} // end internal namespace

//**********************************************************************
PHX_EVALUATE_FIELDS(Integrator_BasisTimesVector,workset)
{ 
  // residual.deep_copy(ScalarT(0.0));

  // initialize the scratch field with vectorField times the 
  // multiplier
  {
    FillScratch_Initialize<ScalarT> fillScratch_I;
    fillScratch_I.scratch = scratch;
    fillScratch_I.multiplier = multiplier;
    fillScratch_I.vectorField = vectorField;
    Kokkos::parallel_for(workset.num_cells,fillScratch_I);
  }

  // this multipliers each entry by the field
  {
    FillScratch_FieldMultipliers<ScalarT> fillScratch_FM;
    fillScratch_FM.scratch = scratch;
    for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
         field != field_multipliers.end(); ++field) {
      fillScratch_FM.field = *field;
      Kokkos::parallel_for(workset.num_cells,fillScratch_FM);
    }
  }

  if(workset.num_cells>0) {

    Integrate_Values<ScalarT> integrate_V;
    integrate_V.weighted_basis_vector = workset.bases[basis_index]->weighted_basis_vector;
    integrate_V.residual = residual;
    integrate_V.scratch = scratch;

    Kokkos::parallel_for(workset.num_cells,integrate_V);
  }
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

