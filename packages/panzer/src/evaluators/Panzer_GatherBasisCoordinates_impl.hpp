#ifndef PANZER_GATHER_BASIS_COORDINATES_IMPL_HPP
#define PANZER_GATHER_BASIS_COORDINATES_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_Workset_Utilities.hpp"

#include "Teuchos_FancyOStream.hpp"

template<typename EvalT,typename Traits>
std::string 
panzer::GatherBasisCoordinates<EvalT, Traits>::
fieldName(const std::string & basisName)
{
   std::stringstream ss; 
   ss << "Basis_" << basisName << " BasisCoordinates";
   return ss.str();
}

template<typename EvalT,typename Traits>
panzer::GatherBasisCoordinates<EvalT, Traits>::
GatherBasisCoordinates(const panzer::PureBasis & basis)
{ 
  basisName_ = basis.name();

  basisCoordinates_ = PHX::MDField<ScalarT,Cell,BASIS,Dim>(fieldName(basisName_),basis.coordinates);

  this->addEvaluatedField(basisCoordinates_);

  this->setName("Gather "+fieldName(basisName_));
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::GatherBasisCoordinates<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData sd, 
		      PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(basisCoordinates_,fm);

  basisIndex_ = panzer::getBasisIndex(basisName_, (*sd.worksets_)[0]);
}

// **********************************************************************
template<typename EvalT,typename Traits> 
void panzer::GatherBasisCoordinates<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
  const Intrepid::FieldContainer<double> & basisCoords = workset.bases[basisIndex_]->basis_coordinates;  

  // just copy the array
  for(int i=0;i<basisCoords.size();i++)
     basisCoordinates_[i] = basisCoords[i];
}

#endif
