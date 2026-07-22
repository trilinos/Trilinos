// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_INTEGRATE_SOURCE_TERM_HPP
#define PHX_INTEGRATE_SOURCE_TERM_HPP

#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"
#include "Dimension.hpp"

template<typename EvalT, typename Traits>
class IntegrateSourceTerm : public PHX::EvaluatorWithBaseImpl<Traits>,
                            public PHX::EvaluatorDerived<EvalT, Traits>  {

  using ScalarT = typename EvalT::ScalarT;
  PHX::MDField<const ScalarT,CELL,QP> source;
  PHX::MDField<ScalarT,CELL,BASIS> residual;
#ifdef PHX_ENABLE_KOKKOS_AMT
  // Make residual atomic so that AMT mode can sum diffusion and source terms at same time
  Kokkos::View<ScalarT**,typename PHX::DevLayout<ScalarT>::type,PHX::Device,Kokkos::MemoryTraits<Kokkos::Atomic>> residual_atomic;
#endif
  Kokkos::View<const double**,PHX::Device> basis_view;
  Kokkos::View<const double*,PHX::Device> weights;
  Kokkos::View<const double**,PHX::Device> cell_measure;
  
public:
  IntegrateSourceTerm(const std::string& source_name,
                      const Teuchos::RCP<PHX::DataLayout>& source_layout,
                      const std::string& residual_name,
                      const Teuchos::RCP<PHX::DataLayout>& residual_layout);
  void evaluateFields(typename Traits::EvalData workset) override;
  KOKKOS_INLINE_FUNCTION
  void operator () (const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const;
};

#endif
