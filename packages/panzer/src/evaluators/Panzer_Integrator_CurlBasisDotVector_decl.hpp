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

#ifndef PANZER_EVALUATOR_CURLBASISDOTVECTOR_DECL_HPP
#define PANZER_EVALUATOR_CURLBASISDOTVECTOR_DECL_HPP

#include <string>
#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Intrepid_FieldContainer.hpp"

namespace panzer {
    
/** In 3D this computes
  * 
  *  \f$\int \nabla\times \phi \cdot v \f$
  *
  * however the name can be misleading. The curl of a vector
  * in 2D is simply a scalar, here the evaluators handles
  * both cases.
  */
PHX_EVALUATOR_CLASS(Integrator_CurlBasisDotVector)
  
  PHX::MDField<ScalarT,Cell,BASIS> residual;
  PHX::MDField<ScalarT,Cell,BASIS> dof_orientation;
    
  PHX::MDField<ScalarT> flux; // note that this is used for both vector and scalar fields
                              // but this follows the distinction between 2D and 3D curls
    
  std::vector<PHX::MDField<ScalarT,Cell,IP> > field_multipliers;

  std::size_t num_nodes;
  std::size_t num_qp;
  std::size_t num_dim;

  double multiplier;

  std::string basis_name;
  std::size_t basis_index;

  bool useScalarField;

  Intrepid::FieldContainer<ScalarT> tmp;

private:
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;
PHX_EVALUATOR_CLASS_END

}

#endif
