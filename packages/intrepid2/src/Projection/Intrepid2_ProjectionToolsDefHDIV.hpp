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

/** \file   Intrepid2_ProjectionToolsDefHDIV.hpp
    \brief  Header file for the Intrepid2::Experimental::ProjectionTools
            containing definitions for HDIV projections.
    \author Created by Mauro Perego
 */

#ifndef __INTREPID2_PROJECTIONTOOLSDEFHDIV_HPP__
#define __INTREPID2_PROJECTIONTOOLSDEFHDIV_HPP__

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"


namespace Intrepid2 {
namespace Experimental {


template<typename SpT>
template<typename BasisType,
typename ortValueType,       class ...ortProperties>
void
ProjectionTools<SpT>::getHDivEvaluationPoints(typename BasisType::scalarViewType evaluationPoints,
    typename BasisType::scalarViewType extDerivEvaluationPoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts,
    const BasisType* cellBasis,
    ProjectionStruct<SpT, typename BasisType::scalarType> * projStruct,
    const EvalPointsType evalPointType){
  typedef typename BasisType::scalarType scalarType;
  typedef Kokkos::DynRankView<scalarType,SpT> scalarViewType;
  typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
  ordinal_type dim = cellBasis->getBaseCellTopology().getDimension();
  ordinal_type sideDim = dim-1;

  ordinal_type numSides = cellBasis->getBaseCellTopology().getSideCount();

  ordinal_type numCells = orts.extent(0);
  Kokkos::DynRankView<ordinal_type> sOrt("sOrt", numSides);

  for(ordinal_type is=0; is<numSides; ++is)  {
    range_type sidePointsRange;
    scalarViewType sideCubPoints;
    if(evalPointType == TARGET) {
      sidePointsRange = projStruct->getTargetPointsRange(sideDim, is);
      sideCubPoints = projStruct->getTargetEvalPoints(sideDim, is);
    }
    else {
      sidePointsRange = projStruct->getBasisPointsRange(sideDim, is);
      sideCubPoints = projStruct->getBasisEvalPoints(sideDim, is);
    }

    scalarViewType orientedTargetCubPoints("orientedTargetCubPoints", sideCubPoints.extent(0),sideDim);

    const auto topoKey = projStruct->getTopologyKey(sideDim,is);

    for(ordinal_type ic=0; ic<numCells; ++ic)  {
      if(dim == 3)
        orts(ic).getFaceOrientation(sOrt.data(), numSides);
      else
        orts(ic).getEdgeOrientation(sOrt.data(), numSides);
      ordinal_type ort = sOrt(is);
      Impl::OrientationTools::mapToModifiedReference(orientedTargetCubPoints,sideCubPoints,topoKey,ort);
      CellTools<SpT>::mapToReferenceSubcell(Kokkos::subview(evaluationPoints,ic,sidePointsRange,Kokkos::ALL()), orientedTargetCubPoints, sideDim, is, cellBasis->getBaseCellTopology());
    }
  }

  if(cellBasis->getDofCount(dim,0) <= 0)
    return;

  range_type cellDivPointsRange;
  scalarViewType divCubPoints;
  if(evalPointType == TARGET) {
    divCubPoints = projStruct->getTargetDerivEvalPoints(dim, 0);
    cellDivPointsRange = projStruct->getTargetDerivPointsRange(dim, 0);
  } else {
    divCubPoints = projStruct->getBasisDerivEvalPoints(dim, 0);
    cellDivPointsRange = projStruct->getBasisDerivPointsRange(dim, 0);
  }
  RealSpaceTools<SpT>::clone(Kokkos::subview(extDerivEvaluationPoints, Kokkos::ALL(), cellDivPointsRange, Kokkos::ALL()), divCubPoints);

  if(projStruct->getTargetEvalPoints(dim, 0).data() != NULL)
  {
    range_type cellPointsRange;
    scalarViewType cubPoints;
    if(evalPointType == TARGET) {
      cubPoints = projStruct->getTargetEvalPoints(dim, 0);
      cellPointsRange = projStruct->getTargetPointsRange(dim, 0);
    } else {
      cubPoints = projStruct->getBasisEvalPoints(dim, 0);
      cellPointsRange = projStruct->getBasisPointsRange(dim, 0);
    }
    RealSpaceTools<SpT>::clone(Kokkos::subview(evaluationPoints, Kokkos::ALL(), cellPointsRange, Kokkos::ALL()), cubPoints);
  }
}


template<typename SpT>
template<typename basisCoeffsValueType, class ...basisCoeffsProperties,
typename funValsValueType, class ...funValsProperties,
typename BasisType,
typename ortValueType,class ...ortProperties>
void
ProjectionTools<SpT>::getHDivBasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtEvalPoints,
    const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetDivAtDivEvalPoints,
    const typename BasisType::scalarViewType evaluationPoints,
    const typename BasisType::scalarViewType extDerivEvaluationPoints,
    const Kokkos::DynRankView<ortValueType,   ortProperties...>  orts,
    const BasisType* cellBasis,
    ProjectionStruct<SpT, typename BasisType::scalarType> * projStruct){
  typedef typename Kokkos::Impl::is_space<SpT>::host_mirror_space::execution_space host_space_type;
  typedef typename BasisType::scalarType scalarType;
  typedef Kokkos::DynRankView<scalarType,SpT> scalarViewType;
  typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
  const auto cellTopo = cellBasis->getBaseCellTopology();
  ordinal_type dim = cellTopo.getDimension();
  ordinal_type numTotalEvaluationPoints(targetAtEvalPoints.extent(1)),
      numTotalDivEvaluationPoints(targetDivAtDivEvalPoints.extent(1));
  ordinal_type basisCardinality = cellBasis->getCardinality();
  ordinal_type numCells = targetAtEvalPoints.extent(0);
  const ordinal_type sideDim = dim-1;

  const std::string& name = cellBasis->getName();

  ordinal_type numSides = cellBasis->getBaseCellTopology().getSideCount();
  scalarViewType refSideNormal("refSideNormal", dim);

  ordinal_type numSideDofs(0);
  for(ordinal_type is=0; is<numSides; ++is)
    numSideDofs += cellBasis->getDofCount(sideDim,is);

  Kokkos::View<ordinal_type*> computedDofs("computedDofs",numSideDofs);

  ordinal_type computedDofsCount = 0;

  ordinal_type numTotalCubPoints = projStruct->getNumBasisEvalPoints(), numTotalDivCubPoints = projStruct->getNumBasisDerivEvalPoints();
  scalarViewType cubPoints_("cubPoints",numCells,numTotalCubPoints, dim);
  scalarViewType divCubPoints("divCubPoints",numCells,numTotalDivCubPoints, dim);
  getHDivEvaluationPoints(cubPoints_, divCubPoints, orts, cellBasis, projStruct, BASIS);

  scalarViewType basisAtCubPoints("basisAtCubPoints",numCells,basisCardinality, numTotalCubPoints, dim);
  scalarViewType basisAtTargetCubPoints("basisAtTargetCubPoints",numCells,basisCardinality, numTotalEvaluationPoints, dim);
  {
    scalarViewType nonOrientedBasisAtCubPoints("nonOrientedBasisAtCubPoints",numCells,basisCardinality, numTotalCubPoints, dim);
    scalarViewType nonOrientedBasisAtTargetCubPoints("nonOrientedBasisAtTargetCubPoints",numCells,basisCardinality, numTotalEvaluationPoints, dim);
    for(ordinal_type ic=0; ic<numCells; ++ic) {
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtTargetCubPoints,ic,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(evaluationPoints, ic, Kokkos::ALL(), Kokkos::ALL()));
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisAtCubPoints,ic,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(cubPoints_, ic, Kokkos::ALL(), Kokkos::ALL()));
    }

    OrientationTools<SpT>::modifyBasisByOrientation(basisAtCubPoints, nonOrientedBasisAtCubPoints, orts, cellBasis);
    OrientationTools<SpT>::modifyBasisByOrientation(basisAtTargetCubPoints, nonOrientedBasisAtTargetCubPoints, orts, cellBasis);
  }

  scalarViewType basisDivAtDivCubPoints;
  scalarViewType basisDivAtTargetDivCubPoints;
  if(numTotalDivEvaluationPoints>0) {
    scalarViewType nonOrientedBasisDivAtTargetDivCubPoints, nonOrientedBasisDivAtDivCubPoints;
    basisDivAtDivCubPoints = scalarViewType ("basisDivAtDivCubPoints",numCells,basisCardinality, numTotalDivCubPoints);
    nonOrientedBasisDivAtDivCubPoints = scalarViewType ("nonOrientedBasisDivAtDivCubPoints",numCells,basisCardinality, numTotalDivCubPoints);
    basisDivAtTargetDivCubPoints = scalarViewType("basisDivAtTargetDivCubPoints",numCells,basisCardinality, numTotalDivEvaluationPoints);
    nonOrientedBasisDivAtTargetDivCubPoints = scalarViewType("nonOrientedBasisDivAtTargetDivCubPoints",numCells,basisCardinality, numTotalDivEvaluationPoints);

    for(ordinal_type ic=0; ic<numCells; ++ic) {
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisDivAtDivCubPoints,ic,Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(divCubPoints, ic, Kokkos::ALL(), Kokkos::ALL()),OPERATOR_DIV);
      cellBasis->getValues(Kokkos::subview(nonOrientedBasisDivAtTargetDivCubPoints,ic,Kokkos::ALL(),Kokkos::ALL()), Kokkos::subview(extDerivEvaluationPoints, ic, Kokkos::ALL(), Kokkos::ALL()),OPERATOR_DIV);
    }
    OrientationTools<SpT>::modifyBasisByOrientation(basisDivAtDivCubPoints, nonOrientedBasisDivAtDivCubPoints, orts, cellBasis);
    OrientationTools<SpT>::modifyBasisByOrientation(basisDivAtTargetDivCubPoints, nonOrientedBasisDivAtTargetDivCubPoints, orts, cellBasis);
  }


  for(ordinal_type is=0; is<numSides; ++is)  {

    ordinal_type sideCardinality = cellBasis->getDofCount(sideDim,is);
    ordinal_type numCubPoints = projStruct->getNumBasisEvalPoints(sideDim,is);
    ordinal_type numTargetCubPoints = projStruct->getNumTargetEvalPoints(sideDim,is);

    for(ordinal_type i=0; i<sideCardinality; ++i)
      computedDofs(computedDofsCount++) = cellBasis->getDofOrdinal(sideDim, is, i);

    CellTools<SpT>::getReferenceSideNormal(refSideNormal, is, cellBasis->getBaseCellTopology());

    scalarViewType normalBasisAtElemcubPoints("normalBasisAtElemcubPoints",numCells,sideCardinality, numCubPoints);
    scalarViewType normalBasisAtTargetcubPoints("normalBasisAtTargetcubPoints",numCells,sideCardinality, numTargetCubPoints);
    scalarViewType weightedNormalBasisAtElemcubPoints("weightedNormalBasisAtElemcubPoints",numCells,sideCardinality, numCubPoints);
    scalarViewType weightedNormalBasisAtTargetcubPoints("weightedNormalBasisAtTargetcubPoints",numCells,sideCardinality, numTargetCubPoints);
    scalarViewType normalTargetAtTargetcubPoints("normalTargetAtTargetcubPoints",numCells, numTargetCubPoints);

    scalarViewType targetEvalWeights = projStruct->getTargetEvalWeights(sideDim, is);
    scalarViewType basisEvalWeights = projStruct->getBasisEvalWeights(sideDim, is);

    //Note: we are not considering the jacobian of the orientation map since it is simply a scalar term for the integrals and it does not affect the projection
    ordinal_type offsetBasis = projStruct->getBasisPointsRange(sideDim, is).first;
    ordinal_type offsetTarget = projStruct->getTargetPointsRange(sideDim, is).first;
    for(ordinal_type ic=0; ic<numCells; ++ic)  {
      for(ordinal_type j=0; j <sideCardinality; ++j) {
        ordinal_type side_dof = cellBasis->getDofOrdinal(sideDim, is, j);
        for(ordinal_type iq=0; iq <numCubPoints; ++iq) {
          for(ordinal_type d=0; d <dim; ++d)
            normalBasisAtElemcubPoints(ic,j,iq) += refSideNormal(d)*basisAtCubPoints(ic,side_dof,offsetBasis+iq,d);
          weightedNormalBasisAtElemcubPoints(ic,j,iq) = normalBasisAtElemcubPoints(ic,j,iq)*basisEvalWeights(iq);
        }
        for(ordinal_type iq=0; iq <numTargetCubPoints; ++iq) {
          for(ordinal_type d=0; d <dim; ++d)
            normalBasisAtTargetcubPoints(ic,j,iq) += refSideNormal(d)*basisAtTargetCubPoints(ic,side_dof,offsetTarget+iq,d);
          weightedNormalBasisAtTargetcubPoints(ic,j,iq) = normalBasisAtTargetcubPoints(ic,j,iq)*targetEvalWeights(iq);
        }
      }
      for(ordinal_type iq=0; iq <numTargetCubPoints; ++iq)
        for(ordinal_type d=0; d <dim; ++d)
          normalTargetAtTargetcubPoints(ic,iq) += refSideNormal(d)*targetAtEvalPoints(ic,offsetTarget+iq,d);
    }


    scalarViewType sideMassMat_("sideMassMat_", numCells, sideCardinality+1, sideCardinality+1),
        sideRhsMat_("rhsMat_", numCells, sideCardinality+1);

    scalarViewType targetEvalWeights_("targetEvalWeights", numCells, 1, targetEvalWeights.extent(0));
    RealSpaceTools<SpT>::clone(targetEvalWeights_, targetEvalWeights);

    range_type range_H(0, sideCardinality);
    range_type range_B(sideCardinality, sideCardinality+1);
    scalarViewType ones("ones",numCells,1,numCubPoints);
    Kokkos::deep_copy(ones,1);

    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(sideMassMat_, Kokkos::ALL(), range_H, range_H), normalBasisAtElemcubPoints, weightedNormalBasisAtElemcubPoints);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(sideMassMat_, Kokkos::ALL(), range_H, range_B), weightedNormalBasisAtElemcubPoints, ones);

    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(sideRhsMat_, Kokkos::ALL(), range_H), normalTargetAtTargetcubPoints, weightedNormalBasisAtTargetcubPoints);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(sideRhsMat_, Kokkos::ALL(), range_B), normalTargetAtTargetcubPoints, targetEvalWeights_);

    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> sideMassMat("sideMassMat", sideCardinality+1,sideCardinality+1);
    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> sideRhsMat("sideRhsMat",sideCardinality+1, 1);

    Teuchos::LAPACK<ordinal_type,funValsValueType> lapack;
    ordinal_type info = 0;
    Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> pivVec("pivVec", sideCardinality+1, 1);

    for(ordinal_type ic=0; ic<numCells; ++ic)  {
      Kokkos::deep_copy(sideMassMat,funValsValueType(0));  //LAPACK might overwrite the matrix
      for(ordinal_type i=0; i<sideCardinality; ++i) {
        sideRhsMat(i,0) = sideRhsMat_(ic,i);
        for(ordinal_type j=0; j<sideCardinality+1; ++j){
          sideMassMat(i,j) = sideMassMat_(ic,i,j);
        }
        sideMassMat(sideCardinality,i) =  sideMassMat_(ic,i,sideCardinality);
      }
      sideRhsMat(sideCardinality,0) = sideRhsMat_(ic,sideCardinality);


      lapack.GESV(sideCardinality+1, 1,
          sideMassMat.data(),
          sideMassMat.stride_1(),
          (ordinal_type*)pivVec.data(),
          sideRhsMat.data(),
          sideRhsMat.stride_1(),
          &info);

      for(ordinal_type i=0; i<sideCardinality; ++i){
        ordinal_type facet_dof = cellBasis->getDofOrdinal(dim-1, is, i);
        basisCoeffs(ic,facet_dof) = sideRhsMat(i,0);
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


  //elem
  ordinal_type numElemDofs = cellBasis->getDofCount(dim,0);
  if(numElemDofs==0)
    return;

  Basis<host_space_type,scalarType,scalarType> *hcurlBasis = NULL;
  if(name.find("HEX")!=std::string::npos)
    hcurlBasis = new Basis_HCURL_HEX_In_FEM<host_space_type,scalarType,scalarType>(cellBasis->getDegree());
  else if(name.find("TET")!=std::string::npos)
    hcurlBasis = new Basis_HCURL_TET_In_FEM<host_space_type,scalarType,scalarType>(cellBasis->getDegree());
  else if(name.find("QUAD")!=std::string::npos)
    hcurlBasis = new Basis_HGRAD_QUAD_Cn_FEM<host_space_type,scalarType,scalarType>(cellBasis->getDegree());
  else if(name.find("TRI")!=std::string::npos)
    hcurlBasis = new Basis_HGRAD_TRI_Cn_FEM<host_space_type,scalarType,scalarType>(cellBasis->getDegree());
  else  {
    std::stringstream ss;
    ss << ">>> ERROR (Intrepid2::ProjectionTools::getHDivEvaluationPoints): "
        << "Method not implemented for basis " << name;
    INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str().c_str() );
  }
  if(hcurlBasis == NULL) return;


  ordinal_type numTargetDivCubPoints = projStruct->getNumTargetDerivEvalPoints(dim,0);
  ordinal_type numDivCubPoints = projStruct->getNumBasisDerivEvalPoints(dim,0);

  scalarViewType weightedBasisDivAtcubPoints("weightedBasisDivAtcubPoints",numCells,numElemDofs, numDivCubPoints);
  scalarViewType weightedBasisDivAtcubTargetPoints("weightedBasisDivAtcubTargetPoints",numCells, numElemDofs, numTargetDivCubPoints);

  scalarViewType internalBasisDivAtcubPoints("basisDivAtcubPoints",numCells,numElemDofs, numDivCubPoints);

  scalarViewType targetDivEvalWeights = projStruct->getTargetDerivEvalWeights(dim, 0);
  scalarViewType divEvalWeights = projStruct->getBasisDerivEvalWeights(dim, 0);
  ordinal_type offsetBasisDiv = projStruct->getBasisDerivPointsRange(dim, 0).first;
  ordinal_type offsetTargetDiv = projStruct->getTargetDerivPointsRange(dim, 0).first;

  for(ordinal_type ic=0; ic<numCells; ++ic) {
    for(ordinal_type i=0; i<numElemDofs; ++i) {
      ordinal_type idof = cellBasis->getDofOrdinal(dim, 0, i);
      for(ordinal_type iq=0; iq<numDivCubPoints; ++iq) {
        internalBasisDivAtcubPoints(ic,i,iq) = basisDivAtDivCubPoints(ic,idof,offsetBasisDiv+iq);
        weightedBasisDivAtcubPoints(ic,i,iq) = internalBasisDivAtcubPoints(ic,i,iq) * divEvalWeights(iq);
      }
    }
    for(ordinal_type i=0; i<numElemDofs; ++i) {
      ordinal_type idof = cellBasis->getDofOrdinal(dim, 0, i);
      for(ordinal_type iq=0; iq<numTargetDivCubPoints; ++iq)
        weightedBasisDivAtcubTargetPoints(ic,i,iq) = targetDivEvalWeights(iq)*basisDivAtTargetDivCubPoints(ic,idof,offsetTargetDiv+iq);
    }
  }

  ordinal_type hcurlBasisCardinality = hcurlBasis->getCardinality();
  ordinal_type numCurlInteriorDOFs = hcurlBasis->getDofCount(dim,0);

  range_type range_H(0, numElemDofs);
  range_type range_B(numElemDofs, numElemDofs+numCurlInteriorDOFs);

  Kokkos::DynRankView<funValsValueType,Kokkos::LayoutLeft,host_space_type> massMat_("massMat_",numCells,numElemDofs+numCurlInteriorDOFs,numElemDofs+numCurlInteriorDOFs);
  Kokkos::DynRankView<funValsValueType,Kokkos::LayoutLeft,host_space_type> rhsMatTrans("rhsMatTrans",numCells,numElemDofs+numCurlInteriorDOFs);

  scalarViewType targetSideDivAtcubPoints("targetSideDivAtcubPoints",numCells, numDivCubPoints);
  for(ordinal_type i=0; i<numSideDofs; ++i) {
    ordinal_type idof = computedDofs(i);
    for(ordinal_type ic=0; ic<numCells; ++ic)
      for(ordinal_type iq=0; iq<numDivCubPoints; ++iq){
        targetSideDivAtcubPoints(ic,iq) -= basisCoeffs(ic,idof)*basisDivAtDivCubPoints(ic,idof,offsetBasisDiv+iq);
      }
  }

  FunctionSpaceTools<SpT >::integrate(Kokkos::subview(massMat_, Kokkos::ALL(), range_H,range_H), internalBasisDivAtcubPoints, weightedBasisDivAtcubPoints);
  FunctionSpaceTools<SpT >::integrate(Kokkos::subview(rhsMatTrans, Kokkos::ALL(), range_H), targetDivAtDivEvalPoints, weightedBasisDivAtcubTargetPoints);
  FunctionSpaceTools<SpT >::integrate(Kokkos::subview(rhsMatTrans, Kokkos::ALL(), range_H), targetSideDivAtcubPoints, weightedBasisDivAtcubPoints,true);

  if(numCurlInteriorDOFs>0){
    scalarViewType cubPoints = projStruct->getBasisEvalPoints(dim,0);
    ordinal_type numCubPoints = projStruct->getNumBasisEvalPoints(dim,0);
    ordinal_type numTargetCubPoints = projStruct->getNumTargetEvalPoints(dim,0);

    scalarViewType targetSideApproxAtcubPoints("targetSideAtcubPoints",numCells, numCubPoints, dim);
    scalarViewType internalBasisAtcubPoints("basisAtcubPoints",numCells,numElemDofs, numCubPoints, dim);
    scalarViewType hcurlBasisCurlAtcubPoints("hcurlBasisCurlAtcubPoints",hcurlBasisCardinality, numCubPoints,dim);
    scalarViewType internalHcurlBasisCurlAtcubPoints("internalHcurlBasisCurlAtcubPoints",numCells,numCurlInteriorDOFs, numCubPoints,dim);
    scalarViewType hcurlBasisCurlAtcubTargetPoints("hcurlBasisCurlAtcubTargetPoints", hcurlBasisCardinality,numTargetCubPoints, dim);
    scalarViewType internalHcurlBasisCurlAtcubTargetPoints("internalHcurlBasisCurlAtcubTargetPoints",numCells, numCurlInteriorDOFs, numTargetCubPoints, dim);
    scalarViewType weightedHcurlBasisCurlAtcubPoints("weightedHcurlBasisHcurlAtcubPoints", numCells, numCurlInteriorDOFs, numCubPoints,dim);
    scalarViewType weightedHcurlBasisCurlAtcubTargetPoints("weightedHcurlBasisHcurlAtcubTargetPoints",numCells, numCurlInteriorDOFs, numTargetCubPoints,dim);

    hcurlBasis->getValues(hcurlBasisCurlAtcubPoints, cubPoints, OPERATOR_CURL);

    ordinal_type offsetBasis = projStruct->getBasisPointsRange(dim, 0).first;
    range_type targetPointsRange = projStruct->getTargetPointsRange(dim, 0);

    scalarViewType targetEvalWeights = projStruct->getTargetEvalWeights(dim, 0);
    scalarViewType basisEvalWeights = projStruct->getBasisEvalWeights(dim, 0);


    for(ordinal_type ic=0; ic<numCells; ++ic) {

      for(ordinal_type i=0; i<numSideDofs; ++i) {
        ordinal_type idof = computedDofs(i);
        for(ordinal_type iq=0; iq<numCubPoints; ++iq){
          for(ordinal_type d=0; d <dim; ++d)
            targetSideApproxAtcubPoints(ic,iq,d) -= basisCoeffs(ic,idof)*basisAtCubPoints(ic,idof,offsetBasis+iq,d);
        }
      }

      for(ordinal_type i=0; i<numElemDofs; ++i) {
        ordinal_type idof = cellBasis->getDofOrdinal(dim, 0, i);
        for(ordinal_type iq=0; iq<numCubPoints; ++iq) {
          for(ordinal_type d=0; d<dim; ++d)
            internalBasisAtcubPoints(ic,i,iq,d) = basisAtCubPoints(ic,idof,offsetBasis+iq,d);
        }
      }

      for(ordinal_type i=0; i<numCurlInteriorDOFs; ++i) {
        ordinal_type idof = hcurlBasis->getDofOrdinal(dim, 0, i);
        for(ordinal_type d=0; d<dim; ++d)
          for(ordinal_type iq=0; iq<numCubPoints; ++iq) {
            internalHcurlBasisCurlAtcubPoints(ic,i,iq,d) = hcurlBasisCurlAtcubPoints(idof,iq,d);
            weightedHcurlBasisCurlAtcubPoints(ic,i,iq,d) = internalHcurlBasisCurlAtcubPoints(ic,i,iq,d)*basisEvalWeights(iq);
          }
      }

      hcurlBasis->getValues(hcurlBasisCurlAtcubTargetPoints, Kokkos::subview(evaluationPoints,ic,targetPointsRange,Kokkos::ALL()), OPERATOR_CURL);
      for(ordinal_type i=0; i<numCurlInteriorDOFs; ++i) {
        ordinal_type idof = hcurlBasis->getDofOrdinal(dim, 0, i);
        for(ordinal_type d=0; d<dim; ++d)
          for(ordinal_type iq=0; iq<numTargetCubPoints; ++iq) {
            internalHcurlBasisCurlAtcubTargetPoints(ic,i,iq,d) = hcurlBasisCurlAtcubTargetPoints(idof,iq,d);
            weightedHcurlBasisCurlAtcubTargetPoints(ic,i,iq,d) = internalHcurlBasisCurlAtcubTargetPoints(ic,i,iq,d)*targetEvalWeights(iq);
          }
      }
    }

    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(massMat_, Kokkos::ALL(), range_H,range_B), internalBasisAtcubPoints, weightedHcurlBasisCurlAtcubPoints);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(rhsMatTrans, Kokkos::ALL(), range_B), Kokkos::subview(targetAtEvalPoints, Kokkos::ALL(), targetPointsRange, Kokkos::ALL()), weightedHcurlBasisCurlAtcubTargetPoints);
    FunctionSpaceTools<SpT >::integrate(Kokkos::subview(rhsMatTrans, Kokkos::ALL(), range_B), targetSideApproxAtcubPoints, weightedHcurlBasisCurlAtcubPoints,true);


  }
  delete hcurlBasis;

  Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type>
  massMat("massMat", numElemDofs+numCurlInteriorDOFs, numElemDofs+numCurlInteriorDOFs),
  rhsMat("rhsMat", numElemDofs+numCurlInteriorDOFs, 1 );

  Teuchos::LAPACK<ordinal_type,funValsValueType> lapack;
  ordinal_type info = 0;
  Kokkos::View<funValsValueType**,Kokkos::LayoutLeft,host_space_type> pivVec("pivVec", numElemDofs+numCurlInteriorDOFs, 1);

  for(ordinal_type ic=0; ic<numCells; ++ic) {
    Kokkos::deep_copy(massMat,funValsValueType(0));  //LAPACK might overwrite the matrix

    for(ordinal_type i=0; i<numElemDofs+numCurlInteriorDOFs; ++i) {
      rhsMat(i,0) = rhsMatTrans(ic,i);
      for(ordinal_type j=0; j<numElemDofs; ++j)
        massMat(j,i) = massMat_(ic,j,i);
    }

    for(ordinal_type i=numElemDofs; i<numElemDofs+numCurlInteriorDOFs; ++i)
      for(ordinal_type j=0; j<numElemDofs; ++j)
        massMat(i,j) = massMat(j,i);

    lapack.GESV(numElemDofs+numCurlInteriorDOFs, 1,
        massMat.data(),
        massMat.stride_1(),
        (ordinal_type*)pivVec.data(),
        rhsMat.data(),
        rhsMat.stride_1(),
        &info);

    if (info) {
      std::stringstream ss;
      ss << ">>> ERROR (Intrepid::ProjectionTools::getBasisCoeffs): "
          << "LAPACK return with error code: "
          << info;
      INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str().c_str() );
    }

    for(ordinal_type i=0; i<numElemDofs; ++i) {
      ordinal_type idof = cellBasis->getDofOrdinal(dim, 0, i);
      basisCoeffs(ic,idof) = rhsMat(i,0);
    }
  }
}



}
}

#endif

