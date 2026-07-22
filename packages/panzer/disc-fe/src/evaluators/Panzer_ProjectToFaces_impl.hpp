// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_PROJECT_TO_FACES_IMPL_HPP
#define PANZER_PROJECT_TO_FACES_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Intrepid2_Cubature.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_Kernels.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_OrientationTools.hpp"

#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_CommonArrayFactories.hpp"
#include "Panzer_HierarchicParallelism.hpp"
#include "Kokkos_ViewFactory.hpp"

#include "Teuchos_FancyOStream.hpp"

#include <cstring>

template<typename EvalT,typename Traits>
panzer::ProjectToFaces<EvalT, Traits>::
ProjectToFaces(const Teuchos::ParameterList& p)
{ 
  dof_name_ = (p.get< std::string >("DOF Name"));

  if(p.isType< Teuchos::RCP<PureBasis> >("Basis"))
    basis_ = p.get< Teuchos::RCP<PureBasis> >("Basis");
  else
    basis_ = p.get< Teuchos::RCP<const PureBasis> >("Basis");

  Teuchos::RCP<PHX::DataLayout> basis_layout  = basis_->functional;
  Teuchos::RCP<PHX::DataLayout> vector_layout  = basis_->functional_grad;

  // some sanity checks
  TEUCHOS_ASSERT(basis_->isVectorBasis());

  result_ = PHX::MDField<ScalarT,Cell,BASIS>(dof_name_,basis_layout);
  this->addEvaluatedField(result_);

  normals_ = PHX::MDField<const ScalarT,Cell,BASIS,Dim>(dof_name_+"_Normals",vector_layout);
  this->addDependentField(normals_);

  vector_values_ = PHX::MDField<const ScalarT,Cell,BASIS,Dim>(dof_name_+"_Vector",vector_layout);
  this->addDependentField(vector_values_);
  
  this->setName("Project To Faces");
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::ProjectToFaces<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData d, 
		      PHX::FieldManager<Traits>& /* fm */)
{
  num_faces_  = result_.extent(1);
  num_dim_  = vector_values_.extent(2);

  TEUCHOS_ASSERT(result_.extent(1) == normals_.extent(1));
  TEUCHOS_ASSERT(vector_values_.extent(2) == normals_.extent(2));
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::ProjectToFaces<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 

  // Restricting HDiv field, multiplied by the normal to the faces, into HVol on the faces.
  // This code assumes affine mapping and the projection into 1 quadrature point for each face,
  // which is identified with the face. This makes sense only for low order bases, for which
  // HVol is constant

  //TODO: make this work w/ high order basis
  const int intDegree = basis_->order();
  TEUCHOS_ASSERT(intDegree == 1);

  auto result = result_.get_static_view();
  auto vector_values = vector_values_.get_static_view();
  auto normals = normals_.get_static_view();
  auto num_faces = num_faces_;
  auto num_dim = num_dim_;
    
  auto policy = panzer::HP::inst().teamPolicy<ScalarT,PHX::exec_space>(workset.num_cells);
  Kokkos::parallel_for("panzer::ProjectToFaces",policy,KOKKOS_LAMBDA(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) {
    const auto cell = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_faces),[&] (const int p) {
      result(cell,p) = ScalarT(0.0);
      for (int dim = 0; dim < num_dim; ++dim)
        result(cell,p) += vector_values(cell,p,dim) * normals(cell,p,dim);
    });
  });
}


#endif
