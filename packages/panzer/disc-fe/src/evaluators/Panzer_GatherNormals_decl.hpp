// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_GATHER_NORMALS_DECL_HPP
#define PANZER_EVALUATOR_GATHER_NORMALS_DECL_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"
#include "Panzer_PointValues2.hpp"


#include "Intrepid2_CellTools.hpp"

namespace panzer {

/** \brief Gathers tangent vectors per field from the global indexer and
    stores them in the field manager.
*/
template<typename EvalT, typename Traits> 
class GatherNormals
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits>,
    public CloneableEvaluator  {
   
public:
  
  GatherNormals(const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  {return Teuchos::rcp(new GatherNormals<EvalT,Traits>(pl));}
  
private:

  typedef typename EvalT::ScalarT ScalarT;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  std::string dof_name_;

  PHX::MDField<ScalarT,Cell,NODE,Dim> gatherFieldNormals_;

  GatherNormals();

  Teuchos::RCP<const PureBasis> basis_;
  Teuchos::RCP<const PointRule> pointRule_;
  Intrepid2::RefSubcellParametrization<PHX::Device>::ConstViewType sideParam_;
  PointValues2<double> pointValues_;
  PHX::MDField<const double,Cell,IP,Dim,Dim> constJac_;

  Kokkos::View<Intrepid2::Orientation*> orientations_;
  Kokkos::View<unsigned int*> keys_;
  
  // Temporaries
  Kokkos::DynRankView<ScalarT,PHX::Device> refEdges_;
  Kokkos::DynRankView<ScalarT,PHX::Device> phyEdges_;
};

}

// **************************************************************
#endif
