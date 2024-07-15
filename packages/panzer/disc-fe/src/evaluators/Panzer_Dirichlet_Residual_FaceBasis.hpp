// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_DIRICHLET_RESIDUAL_FACEBASIS_HPP
#define PANZER_EVALUATOR_DIRICHLET_RESIDUAL_FACEBASIS_HPP

#include "Teuchos_RCP.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_PointRule.hpp"
#include "Panzer_PointValues2.hpp"

#include "Kokkos_DynRankView.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
/** Evaluates a Dirichlet BC residual corresponding to a field value
  * at a set of points defined by a point rule. Note that this assumes
  * a vector basis is used.
  */
template<typename EvalT, typename Traits>
class DirichletResidual_FaceBasis
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    DirichletResidual_FaceBasis(
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
  
  PHX::MDField<ScalarT,Cell,BASIS> residual;
  PHX::MDField<const ScalarT,Cell,Point,Dim> dof;
  PHX::MDField<const ScalarT,Cell,Point,Dim> value;

  Teuchos::RCP<const panzer::PureBasis> basis; 
  Teuchos::RCP<const panzer::PointRule> pointRule; 
  Kokkos::DynRankView<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device> faceNormal; // face normals
  Kokkos::DynRankView<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device> refFaceNormal; // reference face normals

  PointValues2<double> pointValues;
  PHX::MDField<const double,Cell,IP,Dim,Dim> constJac_;

  Teuchos::RCP<const std::vector<Intrepid2::Orientation> > orientations;

}; // end of class DirichletResidual_FaceBasis


}

#endif
