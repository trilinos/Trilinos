// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __PANZER_STK_GatherRefCoords_HPP__
#define __PANZER_STK_GatherRefCoords_HPP__

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_STK_Interface.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer_stk {

/** This gathers the physical (mesh) coordinates of a basis's nodes
  * directly from the STK mesh, into a field with one entry per
  * cell/node/spatial dimension.
  */
template<typename EvalT, typename Traits>
class GatherRefCoords
  : public panzer::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits> {

public:
  /// \brief Construct from the mesh, the basis whose node coordinates to gather, and the name to give the output field.
  GatherRefCoords(const Teuchos::RCP<const STK_Interface> & mesh,
                  const panzer::BasisIRLayout & basis,
                  const std::string & fieldName);

  /// \brief Gathers the node coordinate field for this workset per the class description.
  void evaluateFields(typename Traits::EvalData d);

private:
  typedef typename EvalT::ScalarT ScalarT;

  PHX::MDField<ScalarT,panzer::Cell,panzer::NODE,panzer::Dim> coordField_;
 
  Teuchos::RCP<const STK_Interface> mesh_;

  GatherRefCoords();

};

}

// **************************************************************
#endif
