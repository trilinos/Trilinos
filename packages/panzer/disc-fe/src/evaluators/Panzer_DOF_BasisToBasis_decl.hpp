// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_DOF_BASIS_TO_BASIS_DECL_HPP
#define PANZER_EVALUATOR_DOF_BASIS_TO_BASIS_DECL_HPP

#include <string>

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "PanzerDiscFE_config.hpp"
#include "Panzer_PureBasis.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {
    
//! Interpolates DOF coefficients on one basis to points on another basis.  This is used with nodal bases to map DOF coefficient values from one nodal basis to dof coefficients on another basis. 
template <typename EvalT, typename TRAITST>
class DOF_BasisToBasis 
  : public panzer::EvaluatorWithBaseImpl<TRAITST>,
    public PHX::EvaluatorDerived<EvalT, TRAITST> {
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

  void evaluateFields(typename TRAITST::EvalData workset);

private:
  typedef typename EvalT::ScalarT ScalarT;

  //! Dependent field: DOF coefficient values at source basis
  PHX::MDField<const ScalarT> dof_source_coeff;
  
  //! Evaluated field: DOF coefficient values at target basis
  PHX::MDField<ScalarT> dof_target_coeff;

  //! Reference cell basis values at target points, replicated for each cell in workset
  Kokkos::DynRankView<double,PHX::Device> basis;
};

}

#endif
