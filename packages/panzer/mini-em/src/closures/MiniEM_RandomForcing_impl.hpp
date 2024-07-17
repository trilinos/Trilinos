// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MINIEM_RANDOM_FORCING_IMPL_HPP
#define MINIEM_RANDOM_FORCING_IMPL_HPP

#include <cmath>

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_GatherBasisCoordinates.hpp"

#include "Panzer_Traits.hpp"
#include "Kokkos_Random.hpp"
#include "MiniEM_Sacado_Kokkos_Random.hpp"
#include "Kokkos_ArithTraits.hpp"


namespace mini_em {

//**********************************************************************
template <typename EvalT,typename Traits>
RandomForcing<EvalT,Traits>::RandomForcing(const std::string & name,
                                           const panzer::IntegrationRule & ir,
                                           const panzer::FieldLayoutLibrary & fl,
                                           const unsigned int & seed,
                                           const double & rangeMin,
                                           const double & rangeMax,
                                           const std::string& basisName)
{
  using Teuchos::RCP;

  Teuchos::RCP<const panzer::PureBasis> basis = fl.lookupBasis(basisName);

  vectorBasis_ = basis->isVectorBasis();
  if (vectorBasis_) {
    Teuchos::RCP<PHX::DataLayout> data_layout = ir.dl_vector;
    ir_degree = ir.cubature_degree;
    ir_dim = ir.spatial_dimension;

    forcingVector = PHX::MDField<ScalarT,Cell,Point,Dim>(name, data_layout);
    this->addEvaluatedField(forcingVector);
  } else {
    Teuchos::RCP<PHX::DataLayout> data_layout = ir.dl_scalar;
    ir_degree = ir.cubature_degree;
    ir_dim = ir.spatial_dimension;

    forcingScalar = PHX::MDField<ScalarT,Cell,Point>(name, data_layout);
    this->addEvaluatedField(forcingScalar);
  }

  std::string n = "Random Forcing";
  this->setName(n);

  rangeMin_ = rangeMin;
  rangeMax_ = rangeMax;

  rand_pool_ = pool_type(seed);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void RandomForcing<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  using panzer::index_t;

  // double time = workset.time;

  using IST = typename Kokkos::ArithTraits<ScalarT>::val_type;

  const IST max = static_cast<IST>(rangeMin_);
  const IST min = static_cast<IST>(rangeMax_);

  if (vectorBasis_) {
    auto tmp_forcing = forcingVector.get_static_view();

    Kokkos::fill_random(tmp_forcing, rand_pool_, min, max);
  } else {
    auto tmp_forcing = forcingScalar.get_static_view();

    Kokkos::fill_random(tmp_forcing, rand_pool_, min, max);
  }
}

//**********************************************************************
}

#endif
