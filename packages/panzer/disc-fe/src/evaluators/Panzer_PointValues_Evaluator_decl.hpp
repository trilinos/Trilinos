// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_POINT_VALUES_EVALUATOR_DECL_HPP
#define PANZER_POINT_VALUES_EVALUATOR_DECL_HPP

#include <string>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Panzer_PointValues2.hpp"
#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
//! Interpolates basis DOF values to IP DOF values
template<typename EvalT, typename Traits>
class PointValues_Evaluator
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    PointValues_Evaluator(
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

  // is anything other than ScalarT really needed here?
  PointValues2<double> pointValues;
 
  PHX::MDField<double,NODE,Dim> refPointArray;

  bool useBasisValuesRefArray; // if true then basis is non-null
  Teuchos::RCP<const panzer::PureBasis> basis;
  std::size_t basis_index;

  //! Initialization method to unify the constructors.
  template <typename ArrayT>
  void initialize(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                  const Teuchos::Ptr<const ArrayT> & userArray,
                  // const Teuchos::Ptr<const Kokkos::DynRankView<double,PHX::Device> > & userArray,
                  const Teuchos::RCP<const panzer::PureBasis> & pureBasis);

public:
  PointValues_Evaluator(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                        const Kokkos::DynRankView<double,PHX::Device> & userArray);

  PointValues_Evaluator(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                        const PHX::MDField<double, panzer::IP, panzer::Dim> & userArray);

  //! This builds a point rule from the basis function reference points in the workset
  PointValues_Evaluator(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                        const Teuchos::RCP<const panzer::PureBasis> & pureBasis);

}; // end of class PointValues_Evaluator


}

#endif
