// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_PARAMETER_LIBRARY_HPP
#define PANZER_PARAMETER_LIBRARY_HPP

#include "Sacado_ScalarParameterLibrary.hpp"
#include "Sacado_ScalarParameterVector.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_EvaluationTraits.hpp"

namespace panzer {

  /// \brief Panzer's registry of named scalar parameters (see panzer::ScalarParameterEntry), keyed and evaluated using panzer::EvaluationTraits.
  typedef Sacado::ScalarParameterLibrary<panzer::EvaluationTraits> ParamLib;
  /// \brief A vector of named scalar parameter entries drawn from a ParamLib, e.g. for use as LOCA continuation/bifurcation parameters.
  typedef Sacado::ScalarParameterVector<panzer::EvaluationTraits> ParamVec;

}

#endif
