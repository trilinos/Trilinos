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

#ifndef PANZER_EVALUATOR_DIRICHLET_RESIDUAL_EDGEBASIS_HPP
#define PANZER_EVALUATOR_DIRICHLET_RESIDUAL_EDGEBASIS_HPP

#include "Teuchos_RCP.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_config.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_PointRule.hpp"
#include "Panzer_PointValues.hpp"

#include "Intrepid_FieldContainer.hpp"

namespace panzer {
    
/** Evaluates a Dirichlet BC residual corresponding to a field value
  * at a set of points defined by a point rule. Note that this assumes
  * a vector basis is used.
  */
PHX_EVALUATOR_CLASS(DirichletResidual_EdgeBasis)
  
  PHX::MDField<ScalarT,Cell,BASIS> residual;
  PHX::MDField<ScalarT,Cell,Point,Dim> dof;
  PHX::MDField<ScalarT,Cell,Point,Dim> value;
  PHX::MDField<ScalarT,Cell,BASIS> dof_orientation; // will scale residual
                                                    // by orientation to ensure
                                                    // parallel consistency

  Teuchos::RCP<const panzer::PureBasis> basis; 
  Teuchos::RCP<const panzer::PointRule> pointRule; 
  Intrepid::FieldContainer<ScalarT> edgeTan; // edge tangents
  Intrepid::FieldContainer<ScalarT> refEdgeTan; // reference edge tangents

  PointValues<ScalarT,PHX::MDField<ScalarT> > pointValues;

PHX_EVALUATOR_CLASS_END

}

#endif
