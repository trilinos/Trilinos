// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_CellExtreme_hpp__
#define __Panzer_CellExtreme_hpp__

#include <string>
#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Kokkos_DynRankView.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
/** This finds the maximum or minimum value of a scalar quantity, over
  * the integration points of each cell.
  * It is useful for computing extreme-value responses.

  \verbatim
    <ParameterList>
      <Parameter name="Extreme Name" type="string" value="<Name to give to the extreme field>"/>
      <Parameter name="Field Name" type="string" value="<Name of field to find extreme of>"/>
      <Parameter name="IR" type="RCP<IntegrationRule>" value="<user specified IntegrationRule>"/>
      <Parameter name="Use Max" type="bool" value="<Compute maximum (default = true) or minimum (false)>"/>
      <Parameter name="Multiplier" type="double" value="<Scaling factor, default=1>"/>
      <Parameter name="Field Multipliers" type="RCP<const vector<string> >" value="<Other scalar multiplier fields>"/>
    </ParameterList>
  \endverbatim
  */
template<typename EvalT, typename Traits>
class CellExtreme
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    /// \brief Construct from a ParameterList as documented in the class description.
    CellExtreme(
      const Teuchos::ParameterList& p);

    /// \brief Looks up the quadrature index/order needed by evaluateFields(). Called once before the first evaluation.
    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    /// \brief Computes the "Extreme Name" field for this workset per the class description.
    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;
  
  PHX::MDField<ScalarT> extreme;  // result
    
  PHX::MDField<const ScalarT,Cell,IP> scalar; // function to be integrated

  std::vector<PHX::MDField<const ScalarT,Cell,IP> > field_multipliers;
  double multiplier;

  std::size_t num_qp;
  std::size_t quad_index;
  int quad_order;
 
  bool use_max; // true ... if false then this is a "min"

public:
  // for testing purposes
  /// \brief Returns the FieldTag of the extreme-value output field. Used only for unit testing.
  const PHX::FieldTag & getFieldTag() const
  { return extreme.fieldTag(); }

private:
  /// \brief Returns the ParameterList of valid parameters/defaults for this evaluator's constructor.
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class CellExtreme


}

#endif
