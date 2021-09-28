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

#ifndef PANZER_NEUMANN_RESIDUAL_IMPL_HPP
#define PANZER_NEUMANN_RESIDUAL_IMPL_HPP

#include <cstddef>
#include <string>
#include <vector>
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Teuchos_RCP.hpp"

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
NeumannResidual<EvalT, Traits>::
NeumannResidual(
  const Teuchos::ParameterList& p)
{
  std::string residual_name = p.get<std::string>("Residual Name");
  std::string flux_name = p.get<std::string>("Flux Name");
  std::string normal_name = p.get<std::string>("Normal Name");
  std::string normal_dot_flux_name = normal_name + " dot " + flux_name;
  
  const Teuchos::RCP<const panzer::PureBasis> basis =
    p.get< Teuchos::RCP<const panzer::PureBasis> >("Basis");

  const Teuchos::RCP<const panzer::IntegrationRule> ir = 
    p.get< Teuchos::RCP<const panzer::IntegrationRule> >("IR");


  residual = PHX::MDField<ScalarT>(residual_name, basis->functional);
  normal_dot_flux = PHX::MDField<ScalarT>(normal_dot_flux_name, ir->dl_scalar);
  flux = PHX::MDField<const ScalarT>(flux_name, ir->dl_vector);
  normal = PHX::MDField<const ScalarT>(normal_name, ir->dl_vector);

  this->addEvaluatedField(residual);
  this->addEvaluatedField(normal_dot_flux);
  this->addDependentField(normal);
  this->addDependentField(flux);
 
  basis_name = panzer::basisIRLayout(basis,*ir)->name();

  std::string n = "Neumann Residual Evaluator";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
NeumannResidual<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  num_ip = flux.extent(1);
  num_dim = flux.extent(2);

  TEUCHOS_ASSERT(flux.extent(1) == normal.extent(1));
  TEUCHOS_ASSERT(flux.extent(2) == normal.extent(2));

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0], this->wda);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
NeumannResidual<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 
  residual.deep_copy(ScalarT(0.0));
  auto normal_dot_flux_v = normal_dot_flux.get_static_view();
  auto normal_v = normal.get_static_view();
  auto flux_v = flux.get_static_view();

  std::size_t l_num_ip = num_ip, l_num_dim = num_dim;
  Kokkos::parallel_for("NeumannResidual", workset.num_cells, KOKKOS_LAMBDA (const index_t cell) {
      for (std::size_t ip = 0; ip < l_num_ip; ++ip) {
	normal_dot_flux_v(cell,ip) = ScalarT(0.0);
	for (std::size_t dim = 0; dim < l_num_dim; ++dim) {
	  normal_dot_flux_v(cell,ip) += normal_v(cell,ip,dim) * flux_v(cell,ip,dim); 
	}
      }
    });

  // const Kokkos::DynRankView<double,PHX::Device> & weighted_basis = this->wda(workset).bases[basis_index]->weighted_basis;
  const Teuchos::RCP<const BasisValues2<double> > bv = this->wda(workset).bases[basis_index];
  auto weighted_basis_scalar = this->wda(workset).bases[basis_index]->weighted_basis_scalar.get_static_view();
  auto residual_v = residual.get_static_view();
  Kokkos::parallel_for("NeumannResidual", workset.num_cells, KOKKOS_LAMBDA (const index_t cell) {
      for (std::size_t basis = 0; basis < residual_v.extent(1); ++basis) {
	for (std::size_t qp = 0; qp < l_num_ip; ++qp) {
	  residual_v(cell,basis) += normal_dot_flux_v(cell,qp)*weighted_basis_scalar(cell,basis,qp);
	}
      }
    });

  if(workset.num_cells>0)
    Intrepid2::FunctionSpaceTools<PHX::exec_space>::
      integrate<ScalarT>(residual.get_view(),
                         normal_dot_flux.get_view(), 
			 (this->wda(workset).bases[basis_index])->weighted_basis_scalar.get_view());

  Kokkos::fence();

}

//**********************************************************************

}

#endif
