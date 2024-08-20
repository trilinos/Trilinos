// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_GLOBAL_STATISTICS_DECL_HPP
#define PANZER_GLOBAL_STATISTICS_DECL_HPP

#include <iostream>
#include <string>
#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Teuchos_Comm.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {

struct GlobalData;
    
template<typename EvalT, typename Traits>
class GlobalStatistics
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    GlobalStatistics(
      const Teuchos::ParameterList& p);

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d);

    void
    preEvaluate(
      typename Traits::PreEvalData d);

    void
    postEvaluate(
      typename Traits::PostEvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;
  
  PHX::MDField<ScalarT,Cell> volumes;
    
  PHX::MDField<ScalarT,Cell> tmp;

  PHX::MDField<ScalarT,Cell,IP> ones;

  std::vector<PHX::MDField<const ScalarT,Cell,IP> > field_values;

  ScalarT total_volume;
  std::vector<ScalarT> averages;
  std::vector<ScalarT> maxs;
  std::vector<ScalarT> mins;
  ScalarT global_total_volume;
  std::vector<ScalarT> global_averages;
  std::vector<ScalarT> global_maxs;
  std::vector<ScalarT> global_mins;

  int ir_order;
  std::size_t ir_index;

  Teuchos::RCP<const Teuchos::Comm<int> > comm;

  Teuchos::RCP<panzer::GlobalData> global_data;

  void postprocess(std::ostream& os);

public:
  const PHX::FieldTag& getRequiredFieldTag();

}; // end of class GlobalStatistics


}

#endif
