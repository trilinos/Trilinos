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
ProjectionTools<SpT>::getHVolEvaluationPoints(typename BasisType::scalarViewType evaluationPoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  /*orts*/,
    const BasisType* cellBasis,
    ProjectionStruct<SpT, typename BasisType::scalarType> * projStruct,
    const EvalPointsType evalPointType) {
  typedef typename BasisType::scalarType scalarType;
  typedef Kokkos::DynRankView<scalarType,SpT> scalarViewType;
  ordinal_type dim = cellBasis->getBaseCellTopology().getDimension();

  scalarViewType cubPoints;
  if(evalPointType == TARGET) {
    cubPoints = projStruct->getTargetEvalPoints(dim, 0);
  } else {
    cubPoints = projStruct->getBasisEvalPoints(dim, 0);
  }
  RealSpaceTools<SpT>::clone(evaluationPoints,cubPoints);
}


template<typename SpT>
template<typename basisCoeffsValueType, class ...basisCoeffsProperties,
typename funValsValueType, class ...funValsProperties,
typename BasisType,
typename ortValueType,class ...ortProperties>
void
ProjectionTools<SpT>::getHVolBasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtEvalPoints,
    const typename BasisType::scalarViewType evaluationPoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts,
    const BasisType* cellBasis,
    ProjectionStruct<SpT, typename BasisType::scalarType> * projStruct){

  typedef typename Kokkos::Impl::is_space<SpT>::host_mirror_space::execution_space host_space_type;
  typedef typename BasisType::scalarType scalarType;
  typedef Kokkos::DynRankView<scalarType,SpT> scalarViewType;
  ordinal_type dim = cellBasis->getBaseCellTopology().getDimension();

  ordinal_type basisCardinality = cellBasis->getCardinality();

  ordinal_type numCubPoints = projStruct->getNumBasisEvalPoints(dim, 0);
  ordinal_type numTargetCubPoints = projStruct->getNumTargetEvalPoints(dim, 0);
  scalarViewType cubPoints = projStruct->getBasisEvalPoints(dim, 0);
  scalarViewType cubWeights = projStruct->getBasisEvalWeights(dim, 0);
  scalarViewType cubTargetWeights = projStruct->getTargetEvalWeights(dim, 0);

  ordinal_type numCells = targetAtEvalPoints.extent(0);

  scalarViewType basisAtCubPoints("basisAtcubPoints", basisCardinality, numCubPoints);
  scalarViewType basisAtcubTargetPoints("basisAtcubTargetPoints", basisCardinality, numTargetCubPoints);

  cellBasis->getValues(basisAtCubPoints, cubPoints);
  if(evaluationPoints.rank()==3)
    cellBasis->getValues(basisAtcubTargetPoints, Kokkos::subview(evaluationPoints,0,Kokkos::ALL(),Kokkos::ALL()));
  else
    cellBasis->getValues(basisAtcubTargetPoints, evaluationPoints);


  scalarViewType weightedBasisAtcubTargetPoints_("weightedBasisAtcubTargetPoints_",numCells, basisCardinality, numTargetCubPoints);
  scalarViewType cubWeights_(cubWeights.data(),1,numCubPoints);
  scalarViewType evaluationWeights_(cubTargetWeights.data(),1,numTargetCubPoints);
  scalarViewType basisAtcubTargetPoints_(basisAtcubTargetPoints.data(),1, basisCardinality, numTargetCubPoints);
  scalarViewType basisAtCubPoints_(basisAtCubPoints.data(),1, basisCardinality, numCubPoints);
  scalarViewType weightedBasisAtCubPoints("weightedBasisAtCubPoints",1,basisCardinality, numCubPoints);
  scalarViewType weightedBasisAtcubTargetPoints("weightedBasisAtcubTargetPoints",1, basisCardinality, numTargetCubPoints);
  ArrayTools<SpT>::scalarMultiplyDataField( weightedBasisAtCubPoints, cubWeights_, basisAtCubPoints_, false);
  ArrayTools<SpT>::scalarMultiplyDataField( weightedBasisAtcubTargetPoints, evaluationWeights_, basisAtcubTargetPoints, false);
  RealSpaceTools<SpT>::clone(weightedBasisAtcubTargetPoints_,Kokkos::subview(weightedBasisAtcubTargetPoints,0,Kokkos::ALL(), Kokkos::ALL()));

  Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type>
  massMat("massMat", basisCardinality, basisCardinality),
  rhsMat("rhsMat", basisCardinality, numCells );

  Kokkos::DynRankView<funValsValueType,Kokkos::LayoutLeft,host_space_type> massMat_(massMat.data(),1,basisCardinality,basisCardinality);
  Kokkos::DynRankView<funValsValueType,Kokkos::LayoutLeft,host_space_type> rhsMatTrans("rhsMatTrans",numCells,basisCardinality);

  FunctionSpaceTools<SpT >::integrate(massMat_, basisAtCubPoints_, weightedBasisAtCubPoints);
  FunctionSpaceTools<SpT >::integrate(rhsMatTrans, targetAtEvalPoints, weightedBasisAtcubTargetPoints_);

  for(ordinal_type i=0; i<basisCardinality; ++i)
    for(ordinal_type j=0; j<numCells; ++j)
      rhsMat(i,j) = rhsMatTrans(j,i);

  Teuchos::LAPACK<ordinal_type,funValsValueType> lapack;
  ordinal_type info = 0;

  lapack.POSV('U', basisCardinality, numCells,
      massMat.data(),
      massMat.stride_1(),
      rhsMat.data(),
      rhsMat.stride_1(),
      &info);

  for(ordinal_type i=0; i<basisCardinality; ++i)
    for(ordinal_type j=0; j<numCells; ++j) {
      basisCoeffs(j,i) = rhsMat(i,j);
    }

  if (info) {
    std::stringstream ss;
    ss << ">>> ERROR (Intrepid::ProjectionTools::getBasisCoeffs): "
        << "LAPACK return with error code: "
        << info;
    INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str().c_str() );
  }
}



}
}

#endif

