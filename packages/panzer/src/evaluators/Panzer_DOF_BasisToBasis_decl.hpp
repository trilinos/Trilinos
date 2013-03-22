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

#ifndef PANZER_EVALUATOR_DOF_BASIS_TO_BASIS_DECL_HPP
#define PANZER_EVALUATOR_DOF_BASIS_TO_BASIS_DECL_HPP

#include <string>

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_config.hpp"
#include "Panzer_PureBasis.hpp"

#include "Intrepid_FieldContainer.hpp"

namespace panzer {
    
//! Interpolates DOF coefficients on one basis to points on another basis.  This is used with nodal bases to map DOF coefficient values from one nodal basis to dof coefficients on another basis. 
template <typename EvalT, typename TraitsT>
class DOF_BasisToBasis 
  : public PHX::EvaluatorWithBaseImpl<TraitsT>,
    public PHX::EvaluatorDerived<EvalT, TraitsT> {
public:

  /** \brief Ctor
    *
    * \param[in] fieldName Name of the field in the field manager (used for both source and target fields
    * \param[in] sourceBasis Basis that the source DOF coefficients are defined on
    * \param[in] targetBasis Basis that provides the target coordinate points for the field to be interpolated to
    */
  DOF_BasisToBasis(const std::string & fieldName,
		   const PureBasis & sourceBasis,
		   const PureBasis & targetBasis);

  void postRegistrationSetup(typename TraitsT::SetupData d,
			     PHX::FieldManager<TraitsT>& vm);

  void evaluateFields(typename TraitsT::EvalData workset);

private:
  typedef typename EvalT::ScalarT ScalarT;

  //! Dependent field: DOF coefficient values at source basis
  PHX::MDField<ScalarT,Cell,BASIS> dof_source_coeff;
  
  //! Vealuated field: DOF coefficient values at target basis
  PHX::MDField<ScalarT,Cell,BASIS> dof_target_coeff;

  //! Reference cell basis values at target points, replicated for each cell in workset
  Intrepid::FieldContainer<double> basis;
};

}

#endif
