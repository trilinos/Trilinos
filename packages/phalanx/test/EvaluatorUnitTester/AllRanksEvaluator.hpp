// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#ifndef PHX_ALL_RANKS_EVALUATOR_HPP
#define PHX_ALL_RANKS_EVALUATOR_HPP

#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

/** \brief Evaluator to unit test the EvaluatorUnitTester for field checks for all ranks

    This class is used to test the EvaluatorUnitTester object.
*/
template<typename EvalT, typename Traits>
class AllRanksEvaluator : public PHX::EvaluatorWithBaseImpl<Traits>,
                          public PHX::EvaluatorDerived<EvalT, Traits>  {

  using ScalarT = typename EvalT::ScalarT;
  PHX::MDField<const ScalarT,R> f1;
  PHX::MDField<const ScalarT,R,R> f2;
  PHX::MDField<const ScalarT,R,R,R> f3;
  PHX::MDField<const ScalarT,R,R,R,R> f4;
  PHX::MDField<const ScalarT,R,R,R,R,R> f5;
  PHX::MDField<const ScalarT,R,R,R,R,R,R> f6;
  PHX::MDField<ScalarT,R> x1;
  PHX::MDField<ScalarT,R,R> x2;
  PHX::MDField<ScalarT,R,R,R> x3;
  PHX::MDField<ScalarT,R,R,R,R> x4;
  PHX::MDField<ScalarT,R,R,R,R,R> x5;
  PHX::MDField<ScalarT,R,R,R,R,R,R> x6;
  
public:
  AllRanksEvaluator(const Teuchos::RCP<PHX::DataLayout>& dl1,
                    const Teuchos::RCP<PHX::DataLayout>& dl2,
                    const Teuchos::RCP<PHX::DataLayout>& dl3,
                    const Teuchos::RCP<PHX::DataLayout>& dl4,
                    const Teuchos::RCP<PHX::DataLayout>& dl5,
                    const Teuchos::RCP<PHX::DataLayout>& dl6);
  void evaluateFields(typename Traits::EvalData workset) override;
  KOKKOS_INLINE_FUNCTION
  void operator () (const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const;
};

#include "AllRanksEvaluator_Def.hpp"

#endif
