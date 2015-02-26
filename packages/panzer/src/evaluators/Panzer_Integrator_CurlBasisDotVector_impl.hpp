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
#include "Panzer_CommonArrayFactories.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"

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
  if(!useScalarField) {
     flux_vector = PHX::MDField<const ScalarT,Cell,IP,Dim>( p.get<std::string>("Value Name"), 
	                                  p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_vector );
     this->addDependentField(flux_vector);
  }
  else {
     flux_scalar = PHX::MDField<const ScalarT,Cell,IP>( p.get<std::string>("Value Name"), 
   	                                  p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar );
     this->addDependentField(flux_scalar);
  }

  this->addEvaluatedField(residual);
  
  multiplier = p.get<double>("Multiplier");
  if (p.isType<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers")) 
  {
    const std::vector<std::string>& field_multiplier_names = 
      *(p.get<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers"));

    for (std::vector<std::string>::const_iterator name = field_multiplier_names.begin(); 
      name != field_multiplier_names.end(); ++name) 
    {
      PHX::MDField<const ScalarT,Cell,IP> tmp_field(*name, p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
      field_multipliers.push_back(tmp_field);
    }
  }

  for (typename std::vector<PHX::MDField<const ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->addDependentField(*field);

  std::string n = 
    "Integrator_CurlBasisDotVector: " + residual.fieldTag().name();

  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(Integrator_CurlBasisDotVector,sd,fm)
{
  // setup the output field
  this->utils.setFieldData(residual,fm);

  // initialize the input field multiplier fields
  for (typename std::vector<PHX::MDField<const ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->utils.setFieldData(*field,fm);

  // setup the input fields, mainly differentiate between 2D and 3D
  MDFieldArrayFactory af("",fm.template getKokkosExtendedDataTypeDimensions<EvalT>(),true);
  if(!useScalarField) {
    this->utils.setFieldData(flux_vector,fm);

    num_nodes = residual.dimension(1);
    num_qp = flux_vector.dimension(1);
    num_dim = 3;

    basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0]);

    scratch_vector = af.buildStaticArray<ScalarT,Cell,IP,Dim>("btv_scratch",flux_vector.dimension(0),num_qp,num_dim);
    // tmp = Intrepid::FieldContainer<ScalarT>(flux.dimension(0), num_qp, num_dim); 
  }
  else {
    this->utils.setFieldData(flux_scalar,fm);

    num_nodes = residual.dimension(1);
    num_qp = flux_scalar.dimension(1);
    num_dim = 2;

    basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0]);

    scratch_scalar = af.buildStaticArray<ScalarT,Cell,IP>("btv_scratch",flux_scalar.dimension(0),num_qp);
    // tmp = Intrepid::FieldContainer<ScalarT>(flux.dimension(0), num_qp); 
  }
}

namespace {

template <typename ScalarT>
class FillScratchVector {
public:
  typedef typename PHX::Device execution_space;

  // Required for all functors
  PHX::MDField<ScalarT,Cell,IP,Dim> scratch;

  // Required for "Initialize" functor
  double multiplier;
  PHX::MDField<const ScalarT,Cell,IP,Dim> vectorField;

  // Required for "FieldMultipliers" functor
  PHX::MDField<const ScalarT,Cell,IP> field;

  struct Initialize {};
  struct FieldMultipliers {};

  KOKKOS_INLINE_FUNCTION
  void operator()(const Initialize,const unsigned cell) const
  {
    for (std::size_t qp = 0; qp < scratch.dimension_1(); ++qp) 
      for (std::size_t d = 0; d < scratch.dimension_2(); ++d)
        scratch(cell,qp,d) = multiplier * vectorField(cell,qp,d);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const FieldMultipliers,const unsigned cell) const
  {
    for (std::size_t qp = 0; qp < scratch.dimension_1(); ++qp) 
      for (std::size_t d = 0; d < scratch.dimension_2(); ++d)
        scratch(cell,qp,d) *= field(cell,qp)*scratch(cell,qp,d);  
  }
};

template <typename ScalarT>
class FillScratchScalar {
public:
  typedef typename PHX::Device execution_space;

  // Required for all functors
  PHX::MDField<ScalarT,Cell,IP> scratch;

  // Required for "Initialize" functor
  double multiplier;
  PHX::MDField<const ScalarT,Cell,IP> vectorField;

  // Required for "FieldMultipliers" functor
  PHX::MDField<const ScalarT,Cell,IP> field;

  struct Initialize {};
  struct FieldMultipliers {};

  KOKKOS_INLINE_FUNCTION
  void operator()(const Initialize,const unsigned cell) const
  {
    for (std::size_t qp = 0; qp < scratch.dimension_1(); ++qp) 
      scratch(cell,qp) *= vectorField(cell,qp);  
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const FieldMultipliers,const unsigned cell) const
  {
    for (std::size_t qp = 0; qp < scratch.dimension_1(); ++qp) 
      scratch(cell,qp) *= field(cell,qp)*scratch(cell,qp);  
  }
};

template <typename ScalarT,int spaceDim>
class IntegrateValuesVector {
public:
  typedef typename PHX::Device execution_space;
  PHX::MDField<ScalarT,Cell,IP,Dim> scratch;
  PHX::MDField<ScalarT,Cell,BASIS> residual;
  PHX::MDField<double,Cell,BASIS,IP,Dim> weighted_curl_basis;
  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned cell) const
  {
    for (int lbf = 0; lbf < weighted_curl_basis.dimension_1(); lbf++) {
      residual(cell,lbf) = 0.0;
      for (int qp = 0; qp < weighted_curl_basis.dimension_2(); qp++) {
        for (int d = 0; d < spaceDim; d++) {
          residual(cell,lbf) += scratch(cell, qp, d)*weighted_curl_basis(cell, lbf, qp, d);
        } // D-loop
      } // P-loop
    } // F-loop
  }
};

template <typename ScalarT>
class IntegrateValuesScalar {
public:
  typedef typename PHX::Device execution_space;
  PHX::MDField<ScalarT,Cell,IP> scratch;
  PHX::MDField<ScalarT,Cell,BASIS> residual;
  PHX::MDField<double,Cell,BASIS,IP> weighted_curl_basis;
  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned cell) const
  {
    for (int lbf = 0; lbf < weighted_curl_basis.dimension_1(); lbf++) {
      residual(cell,lbf) = 0.0;
      for (int qp = 0; qp < weighted_curl_basis.dimension_2(); qp++) {
          residual(cell,lbf) += scratch(cell,qp)*weighted_curl_basis(cell,lbf,qp);
      } // P-loop
    } // F-loop
  }
};

} // end internal namespace

//**********************************************************************
PHX_EVALUATE_FIELDS(Integrator_CurlBasisDotVector,workset)
{ 
  // residual.deep_copy(ScalarT(0.0));

  const BasisValues2<double> & bv = *workset.bases[basis_index];

  if(!useScalarField) {
    typedef FillScratchVector<ScalarT> FillScratch;
    FillScratch fillScratch;
    fillScratch.scratch     = scratch_vector;
    fillScratch.multiplier  = multiplier;
    fillScratch.vectorField = flux_vector;

    // initialize in first pass (basically a scaled copy)
    Kokkos::parallel_for(Kokkos::RangePolicy<PHX::Device,typename FillScratch::Initialize>(0,workset.num_cells), fillScratch);

    // multiply agains all the fields in the next one
    for (typename std::vector<PHX::MDField<const ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
         field != field_multipliers.end(); ++field) {
      fillScratch.field = *field;

      Kokkos::parallel_for(Kokkos::RangePolicy<PHX::Device,typename FillScratch::FieldMultipliers>(0,workset.num_cells), fillScratch);
    }

    // Build integrate values functor (hard code spatial dimensions)
    IntegrateValuesVector<ScalarT,3> intValues;
    intValues.scratch     = scratch_vector;
    intValues.residual    = residual;
    intValues.weighted_curl_basis = bv.weighted_curl_basis_vector;
    
    Kokkos::parallel_for(workset.num_cells, intValues);
  }
  else {
    typedef FillScratchScalar<ScalarT> FillScratch;
    FillScratch fillScratch;
    fillScratch.scratch     = scratch_scalar;
    fillScratch.multiplier  = multiplier;
    fillScratch.vectorField = flux_scalar;

    // initialize in first pass (basically a scaled copy)
    Kokkos::parallel_for(Kokkos::RangePolicy<PHX::Device,typename FillScratch::Initialize>(0,workset.num_cells), fillScratch);

    // multiply agains all the fields in the next one
    for (typename std::vector<PHX::MDField<const ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
         field != field_multipliers.end(); ++field) {
      fillScratch.field = *field;

      Kokkos::parallel_for(Kokkos::RangePolicy<PHX::Device,typename FillScratch::FieldMultipliers>(0,workset.num_cells), fillScratch);
    }

    // Build integrate values functor
    IntegrateValuesScalar<ScalarT> intValues;
    intValues.scratch     = scratch_scalar;
    intValues.residual    = residual;
    intValues.weighted_curl_basis = bv.weighted_curl_basis_scalar;

    Kokkos::parallel_for(workset.num_cells, intValues);
  }

/*
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
          scratch_vector(cell,qp,dim) = multiplier * tmpVar * flux_vector(cell,qp,dim);
      }
      else {
        // no dimension to loop over for scalar fields
        scratch_scalar(cell,qp) = multiplier * tmpVar * flux_scalar(cell,qp);
      }
    }
  }
*/

/*
  if(!useScalarField) {
    auto weighted_curl_basis_vector = bv.weighted_curl_basis_vector;

    for (std::size_t cell = 0; cell < workset.num_cells; ++cell)
      for (std::size_t basis = 0; basis < num_nodes; ++basis)
        for (std::size_t qp = 0; qp < num_qp; ++qp)
          for (std::size_t dim = 0; dim < num_dim; ++dim)
            residual(cell,basis) += scratch_vector(cell,qp,dim)*weighted_curl_basis_vector(cell,basis,qp,dim);
  }
  else { // useScalarField
    auto weighted_curl_basis_scalar = bv.weighted_curl_basis_scalar;

    for (std::size_t cell = 0; cell < workset.num_cells; ++cell)
      for (std::size_t basis = 0; basis < num_nodes; ++basis)
        for (std::size_t qp = 0; qp < num_qp; ++qp)
          residual(cell,basis) += scratch_scalar(cell,qp)*weighted_curl_basis_scalar(cell,basis,qp);
  }
*/
}

//**********************************************************************

template<typename EvalT, typename TRAITS>
Teuchos::RCP<Teuchos::ParameterList> 
Integrator_CurlBasisDotVector<EvalT, TRAITS>::getValidParameters() const
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
