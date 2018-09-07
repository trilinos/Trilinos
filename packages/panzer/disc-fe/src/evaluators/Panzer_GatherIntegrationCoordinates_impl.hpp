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

#ifndef PANZER_GATHER_INTEGRATION_COORDINATES_IMPL_HPP
#define PANZER_GATHER_INTEGRATION_COORDINATES_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset_Utilities.hpp"

#include "Teuchos_FancyOStream.hpp"

template<typename EvalT,typename TRAITS>
std::string 
panzer::GatherIntegrationCoordinates<EvalT, TRAITS>::
fieldName(int degree)
{
   std::stringstream ss; 
   ss << "IR_" << degree << " IntegrationCoordinates";
   return ss.str();
}

template<typename EvalT,typename TRAITS>
panzer::GatherIntegrationCoordinates<EvalT, TRAITS>::
GatherIntegrationCoordinates(const panzer::IntegrationRule & quad)
{ 
  quadDegree_ = quad.cubature_degree;

  quadCoordinates_ = PHX::MDField<ScalarT,Cell,Point,Dim>(fieldName(quadDegree_),quad.dl_vector);

  this->addEvaluatedField(quadCoordinates_);

  this->setName("Gather "+fieldName(quadDegree_));
}

// **********************************************************************
template<typename EvalT,typename TRAITS>
void panzer::GatherIntegrationCoordinates<EvalT, TRAITS>::
postRegistrationSetup(typename TRAITS::SetupData sd, 
		      PHX::FieldManager<TRAITS>& /* fm */)
{
  quadIndex_ = panzer::getIntegrationRuleIndex(quadDegree_, (*sd.worksets_)[0], this->wda);
}

// **********************************************************************
template<typename EvalT,typename TRAITS> 
void panzer::GatherIntegrationCoordinates<EvalT, TRAITS>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
  // const Kokkos::DynRankView<double,PHX::Device> & quadCoords = this->wda(workset).int_rules[quadIndex_]->ip_coordinates;  
  const IntegrationValues2<double> & iv = *this->wda(workset).int_rules[quadIndex_];

  // just copy the array
  for(int i=0;i<iv.ip_coordinates.extent_int(0);i++)
    for(int j=0;j<iv.ip_coordinates.extent_int(1);j++)
      for(int k=0;k<iv.ip_coordinates.extent_int(2);k++)
	quadCoordinates_(i,j,k) = iv.ip_coordinates(i,j,k);
}

// **********************************************************************
#endif
