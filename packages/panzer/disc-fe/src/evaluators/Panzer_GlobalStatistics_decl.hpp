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
