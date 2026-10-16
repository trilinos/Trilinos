// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Example_AnisotropicTensor_hpp__
#define __Example_AnisotropicTensor_hpp__

#include "PanzerAdaptersSTK_config.hpp"

#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_FieldLibrary.hpp"

#include <string>

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace Example {

  using panzer::Cell;
  using panzer::Point;
  using panzer::Dim;

/** Entries of the constant, symmetric positive definite tensor K. The
  * off-diagonal entries are nonzero so every entry of K contributes to
  * the source, \f$\nabla\cdot(K\nabla\phi)\f$.
  */
struct AnisotropicTensorEntries {
  double xx = 2.0, yy = 1.5, zz = 1.0;
  double xy = 0.5, xz = 0.25, yz = 0.3;
};

/** Evaluates the constant tensor K at the integration points.
  */
template<typename EvalT, typename Traits>
class AnisotropicTensor : public panzer::EvaluatorWithBaseImpl<Traits>,
                          public PHX::EvaluatorDerived<EvalT, Traits>  {

public:
    AnisotropicTensor(const std::string & name,
                      const panzer::IntegrationRule & ir);

    void evaluateFields(typename Traits::EvalData d);


private:
  typedef typename EvalT::ScalarT ScalarT;

  PHX::MDField<ScalarT,Cell,Point,Dim,Dim> tensor;
};

}

#include "Example_AnisotropicTensor_impl.hpp"

#endif
