// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid2_ProjectionToolsDefHVOL.hpp
    \brief  Header file for the Intrepid2::Experimental::ProjectionTools
            containing definitions for HVOL projections.
    \author Created by Mauro Perego
 */

#ifndef __INTREPID2_PROJECTIONTOOLSDEFHVOL_HPP__
#define __INTREPID2_PROJECTIONTOOLSDEFHVOL_HPP__

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"


namespace Intrepid2 {
namespace Experimental {

template<typename SpT>
template<typename BasisType,
typename ortValueType,       class ...ortProperties>
void
ProjectionTools<SpT>::getHVolEvaluationPoints(typename BasisType::ScalarViewType ePoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  /*orts*/,
    const BasisType* cellBasis,
    ProjectionStruct<SpT, typename BasisType::scalarType> * projStruct,
    const EvalPointsType ePointType) {
  ordinal_type dim = cellBasis->getBaseCellTopology().getDimension();
  auto refEPoints = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getEvalPoints(dim,0,ePointType));
  auto ePointsRange  = projStruct->getPointsRange(ePointType);
  RealSpaceTools<SpT>::clone(Kokkos::subview(ePoints, Kokkos::ALL(), ePointsRange(dim, 0), Kokkos::ALL()), refEPoints);
}


template<typename SpT>
template<typename basisCoeffsValueType, class ...basisCoeffsProperties,
typename funValsValueType, class ...funValsProperties,
typename BasisType,
typename ortValueType,class ...ortProperties>
void
ProjectionTools<SpT>::getHVolBasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtTargetEPoints,
    const typename BasisType::ScalarViewType targetEPoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts,
    const BasisType* cellBasis,
    ProjectionStruct<SpT, typename BasisType::scalarType> * projStruct){

  typedef typename BasisType::scalarType scalarType;
  typedef Kokkos::DynRankView<scalarType,SpT> ScalarViewType;
  ordinal_type dim = cellBasis->getBaseCellTopology().getDimension();

  ordinal_type basisCardinality = cellBasis->getCardinality();

  ordinal_type numCells = targetAtTargetEPoints.extent(0);

  auto refTargetEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetEvalWeights(dim,0));
  auto targetEPointsRange  = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getTargetPointsRange());

  auto refBasisEWeights = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisEvalWeights(dim,0));
  auto basisEPointsRange  = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(),projStruct->getBasisPointsRange());

  ordinal_type numTargetEPoints = range_size(targetEPointsRange(dim,0));
  ordinal_type numBasisEPoints = range_size(basisEPointsRange(dim,0));

  ScalarViewType basisAtBasisEPoints("basisAtBasisEPoints", 1, basisCardinality, numBasisEPoints);
  ScalarViewType basisAtTargetEPoints("basisAtTargetEPoints", basisCardinality, numTargetEPoints);

  ScalarViewType basisEPoints("basisEPoints",numCells,projStruct->getNumBasisEvalPoints(), dim);
  getHVolEvaluationPoints(basisEPoints, orts, cellBasis, projStruct, EvalPointsType::BASIS);

  cellBasis->getValues(Kokkos::subview(basisAtBasisEPoints, 0, Kokkos::ALL(), Kokkos::ALL()), Kokkos::subview(basisEPoints,0, Kokkos::ALL(), Kokkos::ALL()));
  if(targetEPoints.rank()==3)
    cellBasis->getValues(basisAtTargetEPoints, Kokkos::subview(targetEPoints, 0, Kokkos::ALL(), Kokkos::ALL()));
  else
    cellBasis->getValues(basisAtTargetEPoints, targetEPoints);

  ScalarViewType weightedBasisAtTargetEPoints("weightedBasisAtTargetEPoints_",numCells, basisCardinality, numTargetEPoints);
  ScalarViewType weightedBasisAtBasisEPoints("weightedBasisAtBasisEPoints", 1, basisCardinality, numBasisEPoints);

  auto tagToOrdinal = Kokkos::create_mirror_view_and_copy(typename SpT::memory_space(), cellBasis->getAllDofOrdinal());
  auto cellDofs = Kokkos::subview(tagToOrdinal, dim, 0, Kokkos::ALL());

  ScalarViewType
  massMat0("massMat0", 1, basisCardinality, basisCardinality),
  massMat("massMat", numCells, basisCardinality, basisCardinality),
  rhsMat("rhsMat", numCells, basisCardinality  );

  ordinal_type offsetBasis = basisEPointsRange(dim,0).first;
  ordinal_type offsetTarget = targetEPointsRange(dim,0).first;
  for(ordinal_type j=0; j <basisCardinality; ++j) {
    ordinal_type idof = cellDofs(j);
      for(ordinal_type iq=0; iq <ordinal_type(refBasisEWeights.extent(0)); ++iq)
        weightedBasisAtBasisEPoints(0,j,iq) = basisAtBasisEPoints(0,idof,offsetBasis+iq) * refBasisEWeights(iq);
      for(ordinal_type iq=0; iq <ordinal_type(refTargetEWeights.extent(0)); ++iq)
        weightedBasisAtTargetEPoints(0,j,iq) = basisAtTargetEPoints(idof,offsetTarget+iq)* refTargetEWeights(iq);
  }

  FunctionSpaceTools<SpT >::integrate(massMat0, basisAtBasisEPoints, weightedBasisAtBasisEPoints);
  RealSpaceTools<SpT>::clone(massMat, Kokkos::subview(massMat0,0,Kokkos::ALL(), Kokkos::ALL()));
  RealSpaceTools<SpT>::clone(weightedBasisAtTargetEPoints,  Kokkos::subview(weightedBasisAtTargetEPoints,0,Kokkos::ALL(), Kokkos::ALL()));
  FunctionSpaceTools<SpT >::integrate(rhsMat, targetAtTargetEPoints, weightedBasisAtTargetEPoints);

  typedef Kokkos::DynRankView<scalarType, Kokkos::LayoutRight, SpT> WorkArrayViewType;
  ScalarViewType t_("t",numCells, basisCardinality);
  WorkArrayViewType w_("w",numCells,basisCardinality);

  ElemSystem cellSystem("cellSystem", true);
  cellSystem.solve(basisCoeffs, massMat, rhsMat, t_, w_, cellDofs, basisCardinality);
}

}
}

#endif

