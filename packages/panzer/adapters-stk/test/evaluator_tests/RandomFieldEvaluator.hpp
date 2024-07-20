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

template <typename ScalarT>
class PointEvaluation {
public:
   virtual void evaluateContainer(const Kokkos::DynRankView<double,PHX::Device> & points,
                                  PHX::MDField<ScalarT> & field) const = 0;
};

template<typename EvalT, typename Traits>
class RandomFieldEvaluator
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    RandomFieldEvaluator(
      const Teuchos::ParameterList& p);

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

  PHX::MDField<ScalarT> field;

public:
  RandomFieldEvaluator(const std::string & name,
                       const Teuchos::RCP<PHX::DataLayout> & dl)
     : field(name,dl) { this->addEvaluatedField(field); }
}; // end of class RandomFieldEvaluator


//**********************************************************************
template<typename EvalT, typename Traits>
RandomFieldEvaluator<EvalT, Traits>::
RandomFieldEvaluator(
  const Teuchos::ParameterList& p)
{
   TEUCHOS_ASSERT(false); // don't do this
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
RandomFieldEvaluator<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData  /* sd */,
  PHX::FieldManager<Traits>&  fm)
{
  this->utils.setFieldData(field,fm);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
RandomFieldEvaluator<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  /* workset */)
{
  auto field_h = Kokkos::create_mirror_view(field.get_view());
   for(int i=0;i<static_cast<int>(field.size());i++)
      field_h[i] = double(std::rand())/double(RAND_MAX);
   Kokkos::deep_copy(field.get_view(),field_h);
}

//**********************************************************************

#endif
