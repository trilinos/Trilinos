// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MINIEM_RANDOM_FORCING_DECL_HPP
#define MINIEM_RANDOM_FORCING_DECL_HPP

#include "PanzerAdaptersSTK_config.hpp"

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_FieldLibrary.hpp"

#include <string>

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace mini_em {

  using panzer::Cell;
  using panzer::Point;
  using panzer::Dim;

/** Random current source
  */
template<typename EvalT, typename Traits>
class RandomForcing : public panzer::EvaluatorWithBaseImpl<Traits>,
                      public PHX::EvaluatorDerived<EvalT, Traits>  {

public:
    RandomForcing(const std::string & name,
                  const panzer::IntegrationRule & ir,
                  const panzer::FieldLayoutLibrary & fl,
                  const unsigned int & seed,
                  const double & rangeMin,
                  const double & rangeMax,
                  const std::string& basisName="E_edge");

    void evaluateFields(typename Traits::EvalData d);


private:
  typedef typename EvalT::ScalarT ScalarT;

  // Simulation source
  PHX::MDField<ScalarT,Cell,Point,Dim> forcingVector;
  PHX::MDField<ScalarT,Cell,Point> forcingScalar;
  int ir_degree, ir_index, ir_dim;
  double rangeMin_, rangeMax_;
  bool vectorBasis_;

  using device_type = PHX::Device;
  using pool_type = Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> ;
  pool_type rand_pool_;
};

}

#include "MiniEM_RandomForcing_impl.hpp"

#endif
