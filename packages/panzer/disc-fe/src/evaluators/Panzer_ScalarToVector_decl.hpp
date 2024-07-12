// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_SCALAR_TO_VECTOR_DECL_HPP
#define PANZER_EVALUATOR_SCALAR_TO_VECTOR_DECL_HPP

#include <vector>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
//! Interpolates basis DOF values to IP DOF values
template<typename EvalT, typename Traits>
class ScalarToVector
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    ScalarToVector(
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
  std::vector< PHX::MDField<const ScalarT,Cell,Point> > scalar_fields;
  PHX::MDField<ScalarT,Cell,Point,Dim> vector_field;

protected:
  typedef PHX::View<const ScalarT**> KokkosScalarFields_t;
  PHX::View<KokkosScalarFields_t*> internal_scalar_fields;
public:

  /**
   * \brief Tag only constructor for this class.
   */
  ScalarToVector(const std::vector<PHX::Tag<ScalarT>> & input,
                 const PHX::FieldTag & output);

}; // end of class ScalarToVector


}

#endif
