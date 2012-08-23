// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

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
