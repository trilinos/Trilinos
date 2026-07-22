// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_GATHER_BASIS_COORDINATES_IMPL_HPP
#define PANZER_GATHER_BASIS_COORDINATES_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_Workset_Utilities.hpp"

#include "Teuchos_FancyOStream.hpp"

template<typename EvalT,typename TRAITS>
std::string 
panzer::GatherBasisCoordinates<EvalT, TRAITS>::
fieldName(const std::string & basisName)
{
   std::stringstream ss; 
   ss << "Basis_" << basisName << " BasisCoordinates";
   return ss.str();
}

template<typename EvalT,typename TRAITS>
panzer::GatherBasisCoordinates<EvalT, TRAITS>::
GatherBasisCoordinates(const panzer::PureBasis & basis)
{ 
  basisName_ = basis.name();

  basisCoordinates_ = PHX::MDField<ScalarT,Cell,BASIS,Dim>(fieldName(basisName_),basis.coordinates);

  this->addEvaluatedField(basisCoordinates_);
  this->addUnsharedField(Teuchos::rcp_const_cast<PHX::FieldTag>(basisCoordinates_.fieldTagPtr()));

  this->setName("GatherBasisCoordinates: "+fieldName(basisName_));
}

// **********************************************************************
template<typename EvalT,typename TRAITS>
void panzer::GatherBasisCoordinates<EvalT, TRAITS>::
postRegistrationSetup(typename TRAITS::SetupData sd, 
		      PHX::FieldManager<TRAITS>& /* fm */)
{
  basisIndex_ = panzer::getPureBasisIndex(basisName_, (*sd.worksets_)[0], this->wda);
  // make sure to zero out derivative array as we only set the value for AD types in the loop below
  Kokkos::deep_copy(basisCoordinates_.get_static_view(),0.0);
}

// **********************************************************************
template<typename EvalT,typename TRAITS> 
void panzer::GatherBasisCoordinates<EvalT, TRAITS>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
  // const Kokkos::DynRankView<double,PHX::Device> & basisCoords = this->wda(workset).bases[basisIndex_]->basis_coordinates;  
  const Teuchos::RCP<const BasisValues2<double> > bv = this->wda(workset).bases[basisIndex_];

  // just copy the array
  auto d_basisCoordinates = basisCoordinates_.get_static_view();
  auto s_basis_coordinates = bv->basis_coordinates.get_static_view();

  Kokkos::MDRangePolicy<PHX::Device,Kokkos::Rank<3>> policy({0,0,0},{int(workset.num_cells),s_basis_coordinates.extent_int(1),s_basis_coordinates.extent_int(2)});
  Kokkos::parallel_for("GatherBasisCoords",policy, KOKKOS_LAMBDA(const int i, const int j, const int k) {
    auto d_basisCoordinates_tmp = d_basisCoordinates;
    auto s_basis_coordinates_tmp = s_basis_coordinates;
    if constexpr(Sacado::IsADType<typename EvalT::ScalarT>::value) {
      d_basisCoordinates_tmp(i,j,k).val() = s_basis_coordinates_tmp(i,j,k);
    }
    else {
      d_basisCoordinates_tmp(i,j,k) = s_basis_coordinates_tmp(i,j,k);
    }
  });
  Kokkos::fence();
}

#endif
