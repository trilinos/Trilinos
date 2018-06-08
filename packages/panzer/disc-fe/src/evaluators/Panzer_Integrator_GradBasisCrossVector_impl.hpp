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

#ifndef PANZER_EVALUATOR_GRADBASISCROSSVECTOR_IMPL_HPP
#define PANZER_EVALUATOR_GRADBASISCROSSVECTOR_IMPL_HPP

#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Kokkos_ViewFactory.hpp"

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
Integrator_GradBasisCrossVector<EvalT, Traits>::
Integrator_GradBasisCrossVector(
  const Teuchos::ParameterList& p):
  _num_basis_nodes(0),
  _num_quadrature_points(0),
  _basis_index(-1)
{

  // TODO: Get all of the panzer modules to read in const integration rules and basisIRlayouts
  Teuchos::RCP<const BasisIRLayout> basis_layout = p.get< Teuchos::RCP<BasisIRLayout> >("Basis").getConst();
  Teuchos::RCP<const panzer::IntegrationRule> ir = p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR").getConst();
  Teuchos::RCP<PHX::DataLayout> vector_dl = p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout Vector");

  Teuchos::RCP<const PureBasis> basis = basis_layout->getBasis();
  _basis_name = basis_layout->name();

  // Verify that this basis supports the gradient operation
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->supportsGrad(),std::logic_error,"Integrator_GradBasisCrossVector: Basis of type \"" << basis->name() << "\" does not support GRAD");

  // Setup residuals
  _residuals.clear();
  if (p.isType<Teuchos::RCP<const std::vector<std::string> > >("Residual Names")){
    const std::vector<std::string> & names = *(p.get<Teuchos::RCP<const std::vector<std::string> > >("Residual Names"));
    for(const auto & name : names){
      PHX::MDField<ScalarT,Cell,BASIS> res(name, basis_layout->functional);
      _residuals.push_back(res);
    }
  }
  _num_dims = _residuals.size();

  // Currently we only allow for vectors of length 3
  TEUCHOS_TEST_FOR_EXCEPTION(_num_dims != vector_dl->dimension(2),std::logic_error,"Vector must be same length as number of residuals.");

  // Setup vector
  _vector = PHX::MDField<const ScalarT,Cell,IP,Dim>(p.get<std::string>("Vector Name"), vector_dl);

  // TODO: I'm assuming that the gradient is held in the same number of dimensions as ir->dl_vector (dl_vector is initialized from mesh dimensions)
  _num_grad_dims = ir->dl_vector->dimension(2);

  // Vector must be at least the same dimension as the GRAD operation
  TEUCHOS_TEST_FOR_EXCEPTION(_num_grad_dims > _num_dims,std::logic_error,"Vector size must have at least as many components as there are dimensions in the mesh.");

  // Read the scalar multiplier
  _multiplier = ScalarT(p.get<double>("Multiplier"));

  // Setup field multipliers
  if (p.isType<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers")){
    const std::vector<std::string> & field_multiplier_names = *(p.get<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers"));
    for (const std::string & name : field_multiplier_names){
      _field_multipliers.push_back(PHX::MDField<const ScalarT,Cell,IP>(name, ir->dl_scalar));
    }
  }

  // Register evaluated fields
  for (auto & residual : _residuals){
    this->addEvaluatedField(residual);
  }

  // Register dependent fields
  this->addDependentField(_vector);
  for (auto & field : _field_multipliers){
    this->addDependentField(field);
  }

  // Name the module
  this->setName("Integrator_GradBasisCrossVector: " + _vector.fieldTag().name());
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
Integrator_GradBasisCrossVector<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& fm)
{

  for (auto & residual : _residuals){
    this->utils.setFieldData(residual,fm);
  }

  this->utils.setFieldData(_vector,fm);

  for (auto & field : _field_multipliers){
    this->utils.setFieldData(field,fm);
  }

  _num_basis_nodes = _residuals[0].dimension(1);
  _num_quadrature_points = _vector.dimension(1);

  _basis_index = panzer::getBasisIndex(_basis_name, (*sd.worksets_)[0], this->wda);

  // TODO: figure out a clean way of cloning _vector
  _tmp = Kokkos::createDynRankView(_residuals[0].get_static_view(),"tmp",_vector.dimension(0), _num_quadrature_points, _num_dims);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
Integrator_GradBasisCrossVector<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 

  // Zero the residuals
  for (int i(0); i < static_cast<int>(_num_dims); ++i)
    Kokkos::deep_copy(_residuals[i].get_static_view(), ScalarT(0.0));

  // The curl operation will only do something if _num_dims == 3
  if(_num_dims != 3){
    return;
  }

  // do a scaled copy to initialize _tmp
  for (int i=0; i < _vector.extent_int(0); ++i)
    for (int j=0; j < _vector.extent_int(1); ++j)
      for (int k=0; k < _vector.extent_int(2); ++k)
        _tmp(i,j,k) = _multiplier * _vector(i,j,k);

  // Apply the field multipliers
  for(auto & field_data : _field_multipliers){
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      for (std::size_t qp = 0; qp < _num_quadrature_points; ++qp) {
        const ScalarT & tmpVar = field_data(cell,qp);
        for (std::size_t dim = 0; dim < _num_dims; ++dim){
          _tmp(cell,qp,dim) *= tmpVar;
        }
      }
    }
  }

  // const Kokkos::DynRankView<double,PHX::Device> & weighted_grad_basis = this->wda(workset).bases[basis_index]->weighted_grad_basis;
  const BasisValues2<double> & bv = *this->wda(workset).bases[_basis_index];

  ScalarT r_0,r_1,r_2;

  // this part gets annoying since dimensionality of the weighted_grad_basis is defined by the mesh
  // There might be a fix, but for now we will just make it ugly.
  if(_num_grad_dims == 1){

    // Curl only exists if the vector is three dimensional
    if(_num_dims == 3){
      for (index_t cell = 0; cell < workset.num_cells; ++cell){
        for (std::size_t basis = 0; basis < _num_basis_nodes; ++basis){
          r_1=r_2=0;
          for (std::size_t qp = 0; qp < _num_quadrature_points; ++qp){
            r_1 += _tmp(cell,qp,2)*bv.weighted_grad_basis(cell,basis,qp,0);
            r_2 += - _tmp(cell,qp,1)*bv.weighted_grad_basis(cell,basis,qp,0);
          }
          _residuals[1](cell,basis) += r_1;
          _residuals[2](cell,basis) += r_2;
        }
      }
    }

  } else if (_num_grad_dims == 2){

    // Curl only exists if the vector is three dimensional
    if(_num_dims == 3){
      for (index_t cell = 0; cell < workset.num_cells; ++cell){
        for (std::size_t basis = 0; basis < _num_basis_nodes; ++basis){
          r_0=r_1=r_2=0;
          for (std::size_t qp = 0; qp < _num_quadrature_points; ++qp){
            r_0 += - _tmp(cell,qp,2)*bv.weighted_grad_basis(cell,basis,qp,1);
            r_1 += _tmp(cell,qp,2)*bv.weighted_grad_basis(cell,basis,qp,0);
            r_2 += _tmp(cell,qp,0)*bv.weighted_grad_basis(cell,basis,qp,1) - _tmp(cell,qp,1)*bv.weighted_grad_basis(cell,basis,qp,0);
          }
          _residuals[0](cell,basis) += r_0;
          _residuals[1](cell,basis) += r_1;
          _residuals[2](cell,basis) += r_2;
        }
      }
    }

  } else if (_num_grad_dims == 3){

    // Will have thrown an error if _num_grad_dims != _num_dims != 3
    for (index_t cell = 0; cell < workset.num_cells; ++cell){
      for (std::size_t basis = 0; basis < _num_basis_nodes; ++basis){
        r_0=r_1=r_2=0;
        for (std::size_t qp = 0; qp < _num_quadrature_points; ++qp){
          r_0 += _tmp(cell,qp,1)*bv.weighted_grad_basis(cell,basis,qp,2) - _tmp(cell,qp,2)*bv.weighted_grad_basis(cell,basis,qp,1);
          r_1 += _tmp(cell,qp,2)*bv.weighted_grad_basis(cell,basis,qp,0) - _tmp(cell,qp,0)*bv.weighted_grad_basis(cell,basis,qp,2);
          r_2 += _tmp(cell,qp,0)*bv.weighted_grad_basis(cell,basis,qp,1) - _tmp(cell,qp,1)*bv.weighted_grad_basis(cell,basis,qp,0);
        }
        _residuals[0](cell,basis) += r_0;
        _residuals[1](cell,basis) += r_1;
        _residuals[2](cell,basis) += r_2;
      }
    }
  }

}

//**********************************************************************

}

#endif
