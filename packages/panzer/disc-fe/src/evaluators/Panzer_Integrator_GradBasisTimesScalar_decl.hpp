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

#ifndef PANZER_EVALUATOR_GRADBASISTIMESSCALAR_DECL_HPP
#define PANZER_EVALUATOR_GRADBASISTIMESSCALAR_DECL_HPP

#include <string>
#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Kokkos_DynRankView.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
template<typename EvalT, typename Traits>
class Integrator_GradBasisTimesScalar
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Integrator_GradBasisTimesScalar(
      const Teuchos::ParameterList& p);

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;
  
  std::vector<PHX::MDField<ScalarT,Cell,BASIS> > _residuals;
    
  PHX::MDField<const ScalarT,Cell,IP> _scalar;
  std::vector<PHX::MDField<const ScalarT,Cell,IP> > _field_multipliers;

  // Number of nodes in the basis
  std::size_t _num_basis_nodes;

  // Number of nodes in the integration rule
  std::size_t _num_quadrature_points;

  // Number of dimensions associated with the gradient
  std::size_t _num_dims;

  // Scalar multiplier against vector
  ScalarT _multiplier;

  // Name of the basis
  std::string _basis_name;

  // Index of basis in workset bases
  std::size_t _basis_index;

  // Temporary variable to store tmp = multipliers * vector
  Kokkos::DynRankView<ScalarT,PHX::Device> _tmp;

}; // end of class Integrator_GradBasisTimesScalar


}

#endif
