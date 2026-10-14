// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Example_SineSource_hpp__
#define __Example_SineSource_hpp__

#include "PanzerAdaptersSTK_config.hpp"

#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_FieldLibrary.hpp"

#include <string>

#include "Panzer_Evaluator_WithBaseImpl.hpp"
#include "Example_AnisotropicTensor.hpp"

namespace Example {
    
  using panzer::Cell;
  using panzer::Point;
  using panzer::Dim;

/** A source for the curl Laplacian that results in the solution. If
  * anisotropic, the source is \f$\nabla\cdot(K\nabla\phi)\f$ for the
  * AnisotropicTensor K. If weighted, the source is scaled by
  * (1+x)(1+y)(1+z).
  */
template<typename EvalT, typename Traits>
class SineSource : public panzer::EvaluatorWithBaseImpl<Traits>,
                        public PHX::EvaluatorDerived<EvalT, Traits>  {

public:
    SineSource(const std::string & name,
                       const panzer::IntegrationRule & ir,
                       const bool weighted = false,
                       const bool anisotropic = false);
                                                                        
    void postRegistrationSetup(typename Traits::SetupData d,           
                               PHX::FieldManager<Traits>& fm);        
                                                                     
    void evaluateFields(typename Traits::EvalData d);               


private:
  typedef typename EvalT::ScalarT ScalarT;

  // Simulation source
  PHX::MDField<ScalarT,Cell,Point> source;
  int ir_degree, ir_index;
  bool weighted;
  AnisotropicTensorEntries K;
};

}

#include "Example_SineSource_impl.hpp"

#endif
