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

/** This class is a gather operation from the mesh. It
  * takes a set of field names and basis objects and
  * then reads them in from the mesh object.
  *
  * The constructor takes a STK_Interface RCP and parameter list
  * that is required to contain the following two fields
  * "Field Names" of type <code>Teuchos::RCP<std::vector<std::string> ></code>
  * and "Basis" of type <code>Teuchos::RCP<panzer::Basis></code>.
  */
template<typename EvalT, typename Traits> 
class GatherRefCoords
  : public panzer::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits> {

public:
  GatherRefCoords(const Teuchos::RCP<const STK_Interface> & mesh,
                  const panzer::BasisIRLayout & basis,
                  const std::string & fieldName);
  
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
