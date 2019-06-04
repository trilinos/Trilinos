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

#ifndef PANZER_DOF_BASIS_TO_BASIS_IMPL_HPP
#define PANZER_DOF_BASIS_TO_BASIS_IMPL_HPP

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Intrepid2_Basis.hpp"
#include "Teuchos_Assert.hpp"

namespace panzer {

//**********************************************************************
template <typename EvalT, typename TRAITST>
DOF_BasisToBasis<EvalT,TRAITST>::
DOF_BasisToBasis(const std::string & fieldName,
		 const PureBasis & sourceBasis,
		 const PureBasis & targetBasis)
{
  TEUCHOS_ASSERT(sourceBasis.numCells() == targetBasis.numCells());
  
  // **************
  // Declare fields
  // **************
  dof_source_coeff = PHX::MDField<const ScalarT>(fieldName,sourceBasis.functional);
  dof_target_coeff = PHX::MDField<ScalarT>(fieldName,targetBasis.functional);

  this->addDependentField(dof_source_coeff);
  this->addEvaluatedField(dof_target_coeff);

  // **************
  // Get coordinate points for reference cell on target basis 
  // **************
  Kokkos::DynRankView<double,PHX::Device>intrpCoords =
    Kokkos::DynRankView<double,PHX::Device>("intrpCoords",targetBasis.cardinality(),targetBasis.dimension());
  
  targetBasis.getIntrepid2Basis<PHX::exec_space,double,double>()->getDofCoords(intrpCoords);

  // **************
  // Evaluate source basis values at target basis coordinates
  // **************
  Kokkos::DynRankView<double,PHX::Device> basisRef = 
    Kokkos::DynRankView<double,PHX::Device>("basisRef",sourceBasis.cardinality(),targetBasis.cardinality());

  sourceBasis.getIntrepid2Basis()->getValues(basisRef, intrpCoords, Intrepid2::OPERATOR_VALUE);
  
  // **************
  // Copy the reference basis values for all cells in workset
  // **************
  basis = Kokkos::DynRankView<double,PHX::Device>("basis",sourceBasis.numCells(),sourceBasis.cardinality(),targetBasis.cardinality());
  Intrepid2::FunctionSpaceTools<PHX::exec_space>::HGRADtransformVALUE(basis,basisRef);
    
  std::string n = "DOF_BasisToBasis: " + dof_target_coeff.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template <typename EvalT, typename TRAITST>
void DOF_BasisToBasis<EvalT,TRAITST>::evaluateFields(typename TRAITST::EvalData workset)
{ 
  // Zero out arrays (intrepid does a sum!)
  dof_target_coeff.deep_copy(ScalarT(0.0));

  if(workset.num_cells>0) {

    // evaluate function at specified points
    Intrepid2::FunctionSpaceTools<PHX::exec_space>::evaluate(dof_target_coeff.get_view(),
                                                             dof_source_coeff.get_view(),
                                                             basis);
  }
}

//**********************************************************************

}

#endif
