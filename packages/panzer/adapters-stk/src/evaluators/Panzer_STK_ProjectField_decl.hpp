// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_STK_PROJECT_FIELD_DECL_HPP
#define PANZER_STK_PROJECT_FIELD_DECL_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Traits.hpp"
#include "Panzer_STK_Interface.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer_stk {

/** \brief Given a field, perform a local L2 projection onto the desired basis.
 * 
 * \note Requires that orientations be given in the \c postRegistrationSetup phase.
*/
template<typename EvalT, typename Traits> 
class ProjectField
  : public panzer::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits>
    {
   
public:

  /**
   * Constructor for the ProjectField evaluator. 
   * 
   * \param[in] inName Name of the source MDField
   * \param[in] src Basis of the source field
   * \param[in] dst Target basis
   * \param[in] outname (Optional) Name for the projected MDField
   */
  
  ProjectField(const std::string & inName, Teuchos::RCP<panzer::PureBasis> src,
               Teuchos::RCP<panzer::PureBasis> dst, std::string outname = "");
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& fm);
  
  void evaluateFields(typename Traits::EvalData d);

private:

  using ScalarT = typename EvalT::ScalarT;

  const std::string field_name_;
  Teuchos::RCP<const panzer::PureBasis> srcBasis_;
  Teuchos::RCP<const panzer::PureBasis> dstBasis_;
  Teuchos::RCP<const std::vector<Intrepid2::Orientation> > orientations_;

  PHX::MDField<double,panzer::Cell,panzer::BASIS> result_;
  PHX::MDField<double,panzer::Cell,panzer::BASIS> source_;

  Kokkos::DynRankView<Intrepid2::Orientation,PHX::Device> local_orts_;

};

}

// **************************************************************
#endif
