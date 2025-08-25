// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef POINT_EVALUATOR
#define POINT_EVALUATOR

#include <string>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Kokkos_Random.hpp"

template <typename ScalarT>
class PointEvaluation {
public:
   virtual void evaluateContainer(const Kokkos::DynRankView<double,PHX::Device> & points,
                                  PHX::MDField<ScalarT> & field) const = 0;
};

template<typename EvalT, typename Traits>
class RandomFieldEvaluator
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:
  void evaluateFields(typename Traits::EvalData d);

private:
  using ScalarT = typename EvalT::ScalarT;
  PHX::MDField<ScalarT> field_;
  Kokkos::Random_XorShift64_Pool<PHX::Device> random_pool_;

public:
  RandomFieldEvaluator(const std::string & name,
                       const Teuchos::RCP<PHX::DataLayout> & dl)
    : field_(name,dl),random_pool_(12345)
  { this->addEvaluatedField(field_); }
};

//**********************************************************************
template<typename EvalT, typename Traits>
void
RandomFieldEvaluator<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  auto field = field_;
  auto random_pool = random_pool_;

  Kokkos::parallel_for("panzer::RandomFieldEvalautor",workset.num_cells,KOKKOS_LAMBDA(const int i){
    auto generator = random_pool.get_state();
    for (size_t j=0; j < field.extent(1); ++j) {
      field(i,j) = generator.drand(0.0,1.0);
    }
    random_pool.free_state(generator);
  });
}

//**********************************************************************

#endif
