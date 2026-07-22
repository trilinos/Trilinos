// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_DOF_PointField_IMPL_HPP
#define PANZER_EVALUATOR_DOF_PointField_IMPL_HPP

#include <string>

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "PanzerDiscFE_config.hpp"
#include "Panzer_PureBasis.hpp"
#include "Intrepid2_Basis.hpp"
#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {
    
//! Interpolates basis DOF using reference coordinates defined by a field
template <typename EvalT, typename TRAITST>
class DOF_PointField 
  : public panzer::EvaluatorWithBaseImpl<TRAITST>,
    public PHX::EvaluatorDerived<EvalT, TRAITST> {
public:

  /** \basic Constructor that allows user to specify a postfix for the
    * field.
    *
    * This constructor builds an evaluator from coordinates defined on
    * the reference element. The name of the evaluated field is flexible,
    * the name being <code>fieldName+postfixFieldName</code>.
    *
    * \param[in] postfixFieldName Postfix string to modify field name
    * \param[in] fieldName Name of DOF field (dimensioned number cells
    *                      by number of basis functions)
    * \param[in] fieldBasis Datalayout describing DOF field
    * \param[in] coordinateName Name of reference coordinates (sized
    *                           number of points by dimension)
    * \param[in] coordLayout Layout for coordinates
    */
  DOF_PointField(const std::string & postfixFieldName,
                 const std::string & fieldName,
                 const PureBasis & fieldBasis,
                 const std::string & coordinateName,
                 const Teuchos::RCP<PHX::DataLayout> & coordLayout,
                 const Teuchos::RCP<PHX::DataLayout> & quadLayout)
  { initialize(fieldName,fieldBasis,coordinateName,coordLayout,quadLayout,postfixFieldName); }

  /** \basic Constructor that appends (or not) the coordinate name to the
    * field.
    *
    * This constructor builds an evaluator from coordinates defined on
    * the reference element. The name of the evaluated field is either
    * <code>fieldName+coordinateName</code> or simply <code>fieldName</code>.
    *
    * \param[in] fieldName Name of DOF field (dimensioned number cells
    *                      by number of basis functions)
    * \param[in] fieldBasis Datalayout describing DOF field
    * \param[in] coordinateName Name of reference coordinates (sized
    *                           number of points by dimension)
    * \param[in] coordLayout Layout for coordinates
    * \param[in] useCoordPostfix Postfix field name with coordinate name.
    */
  DOF_PointField(const std::string & fieldName,
                 const PureBasis & fieldBasis,
                 const std::string & coordinateName,
                 const Teuchos::RCP<PHX::DataLayout> & coordLayout,
                 const Teuchos::RCP<PHX::DataLayout> & quadLayout,
                 bool useCoordPostfix)
  { std::string postfixFieldName = (useCoordPostfix ? coordinateName : ""); 
    initialize(fieldName,fieldBasis,coordinateName,coordLayout,quadLayout,postfixFieldName); }
  
  void evaluateFields(typename TRAITST::EvalData workset);

private:
  typedef typename EvalT::ScalarT ScalarT;

  //! Convenience initialization routine, see constructor above.
  void initialize(const std::string & fieldName,
                  const PureBasis & fieldBasis,
                  const std::string & coordinateName,
                  const Teuchos::RCP<PHX::DataLayout> & coordLayout,
                  const Teuchos::RCP<PHX::DataLayout> & quadLayout,
                  const std::string & postfixFieldName);

  PHX::MDField<const ScalarT,Point,Dim> coordinates; // reference coordinates
  PHX::MDField<const ScalarT> dof_coeff;  // coefficient values   
  PHX::MDField<ScalarT> dof_field;  // evaluate field

  Teuchos::RCP<Intrepid2::Basis<PHX::exec_space,double,double>> intrepidBasis;
  Kokkos::DynRankView<double,PHX::Device> intrpCoords, basisRef, basis;
};

}

#endif
