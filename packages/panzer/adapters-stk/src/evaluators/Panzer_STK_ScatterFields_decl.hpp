// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __PANZER_STK_ScatterFields_decl_HPP__
#define __PANZER_STK_ScatterFields_decl_HPP__

#include "Phalanx_config.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_STK_Interface.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer_stk {

/** This class is a scatter operation to the mesh. It
  * takes a set of field names and basis objects and
  * then writes them to the mesh object. Note that <code>scaling</code> vector
  * must be the same length as the <code>names</code> vector. The scaling
  * is applied to each field.
  */
template <typename EvalT,typename TraitsT>
class ScatterFields : public panzer::EvaluatorWithBaseImpl<TraitsT>,
                      public PHX::EvaluatorDerived<EvalT, TraitsT>  { 
  typedef typename EvalT::ScalarT ScalarT;
  typedef panzer_stk::STK_Interface::SolutionFieldType VariableField;

  std::vector< PHX::MDField<const ScalarT,panzer::Cell,panzer::NODE> > scatterFields_;
  Teuchos::RCP<STK_Interface> mesh_;

  std::vector<double> scaling_;

  bool cellFields_;

  /** \brief Shared setup logic for both constructors: sets up dependent fields and per-field scaling.
    *
    * \param[in] scatterName name to give this evaluator.
    * \param[in] mesh the STK mesh to write field values to.
    * \param[in] basis the basis whose node values are gathered from the fields and written to the mesh.
    * \param[in] names names of the fields to scatter.
    * \param[in] scaling per-field scale factor applied before writing; must be the same length as names.
    */
  void initialize(const std::string & scatterName,
                  const Teuchos::RCP<STK_Interface> mesh,
                  const Teuchos::RCP<const panzer::PureBasis> & basis,
                  const std::vector<std::string> & names,
                  const std::vector<double> & scaling);

public:

  /** \brief Construct from the mesh, basis, field names to scatter, and a per-field scaling factor (must be the same length as names).
    *
    * \param[in] scatterName name to give this evaluator.
    * \param[in] mesh the STK mesh to write field values to.
    * \param[in] basis the basis whose node values are gathered from the fields and written to the mesh.
    * \param[in] names names of the fields to scatter.
    */
  ScatterFields(const std::string & scatterName,
                const Teuchos::RCP<STK_Interface> mesh,
                const Teuchos::RCP<const panzer::PureBasis> & basis,
                const std::vector<std::string> & names);

  /** \brief Construct from the mesh, basis, field names to scatter, and a per-field scaling factor (must be the same length as names).
    *
    * \param[in] scatterName name to give this evaluator.
    * \param[in] mesh the STK mesh to write field values to.
    * \param[in] basis the basis whose node values are gathered from the fields and written to the mesh.
    * \param[in] names names of the fields to scatter.
    * \param[in] scaling per-field scale factor applied before writing; must be the same length as names.
    */
  ScatterFields(const std::string & scatterName,
                const Teuchos::RCP<STK_Interface> mesh,
                const Teuchos::RCP<const panzer::PureBasis> & basis,
                const std::vector<std::string> & names,
                const std::vector<double> & scaling);

  /// \brief Looks up the STK field pointers needed by evaluateFields(). Called once before the first evaluation.
  void postRegistrationSetup(typename TraitsT::SetupData d,
                             PHX::FieldManager<TraitsT>& fm);

  /// \brief Writes this workset's (scaled) field values to the STK mesh per the class description.
  void evaluateFields(typename TraitsT::EvalData d);
}; 

}

// **************************************************************
#endif
