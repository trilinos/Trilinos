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

/** \file   Intrepid2_ProjectionStruct.hpp
    \brief  Header file for the Intrepid2::Experimental::ProjectionStruct.
    \author Created by Mauro Perego
 */
#ifndef __INTREPID2_PROJECTIONSTRUCT_HPP__
#define __INTREPID2_PROJECTIONSTRUCT_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Types.hpp"

#include <array>

namespace Intrepid2 {

namespace Experimental {

/** \class  Intrepid2::Experimental::ProjectionStruct
    \brief  An helper class to compute the evaluation points and weights needed for performing projections

    In order to perform projections, the basis functions and the target function need to be evaluated
    at several sets of evaluation points (cubature points) defined on subcell entities (edges, faces,
    volumes).
    Depending on the projection, the evaluation of derivatives of the basis functions and of the target
    function may be needed as well.
    This class provides a struct to store the evaluation points/weights on the reference cell.

    Use: create the proper ProjectionStruct rule by calling one of the functions:
    <var><b>createL2ProjectionStruct</b></var>,
    <var><b>createHGradProjectionStruct</b></var>,
    <var><b>createHCurlProjectionStruct</b></var>,
    <var><b>createHDivProjectionStruct</b></var>,
    <var><b>createHVolProjectionStruct</b></var>,
    depending on the type of projection wanted.

    The created class is then used with the Projection Tools. See ProjectionTools class for more info.
 */

ordinal_type
KOKKOS_INLINE_FUNCTION
range_size(const Kokkos::pair<ordinal_type, ordinal_type>& range) {
  return range.second - range.first;
}

template<typename SpT, typename ValueType>
class ProjectionStruct {
public:

  enum EvalPointsType {BASIS, TARGET};

  typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
  typedef typename Kokkos::Impl::is_space<SpT>::host_mirror_space::execution_space host_space_type;
  typedef Kokkos::DynRankView<ValueType,host_space_type> view_type;
  typedef Kokkos::View<range_type**,host_space_type> range_tag;
  static constexpr int numberSubCellDims = 4; //{0 for vertex, 1 for edges, 2 for faces, 3 for volumes}
  //max of numVertices, numEdges, numFaces for a reference cell.
  //12 is the number of edges in a Hexahderon.
  //We'll need to change this if we consider generic polyhedra
  static constexpr int maxSubCellsCount = 12;
  typedef std::array<std::array<view_type, maxSubCellsCount>, numberSubCellDims> view_tag;
  typedef Kokkos::View<unsigned**,host_space_type> key_tag;

  /** \brief  Returns number of basis evaluation points
   */
  ordinal_type getNumBasisEvalPoints() {
    return numBasisEvalPoints;
  }

  /** \brief  Returns number of evaluation points for basis derivatives
   */
  ordinal_type getNumBasisDerivEvalPoints() {
    return numBasisDerivEvalPoints;
  }

  /** \brief  Returns number of points where to evaluate the target function
   */
  ordinal_type getNumTargetEvalPoints() {
    return numTargetEvalPoints;
  }

  /** \brief  Returns number of points where to evaluate the derivatives of the target function
   */
  ordinal_type getNumTargetDerivEvalPoints() {
    return numTargetDerivEvalPoints;
  }

  /** \brief  Returns the maximum number of derivative evaluation points across all the subcells
      \param  evalPointType [in] - enum selecting whether the points should be computed for the basis
                               functions or for the target function

      \return the maximum number of the derivative evaluation points across all the subcells
   */
  ordinal_type getMaxNumDerivPoints(const EvalPointsType type) const {
    if(type == BASIS)
      return maxNumBasisDerivEvalPoints;
    else
      return maxNumTargetDerivEvalPoints;
  }

  /** \brief  Returns the maximum number of evaluation points across all the subcells
      \param  evalPointType [in] - enum selecting whether the points should be computed for the basis
                               functions or for the target function

      \return the maximum number of the evaluation points across all the subcells
   */
  ordinal_type getMaxNumEvalPoints(const EvalPointsType type) const {
    if(type == BASIS)
      return maxNumBasisEvalPoints;
    else
      return maxNumTargetEvalPoints;
  }

  /** \brief  Returns the basis evaluation points on a subcell

      \code
      P  - num. evaluation points
      D  - spatial dimension
      \endcode

      \param  subCellDim  [in]  - dimension of the subcell
      \param  subCellId   [in]  - ordinal of the subcell defined by cell topology

      \return a rank-2 view (P,D) containing the basis evaluation points on the selected subcell
   */
  view_type getBasisEvalPoints(const ordinal_type subCellDim, const ordinal_type subCellId) {
    return basisCubPoints[subCellDim][subCellId];
  }


  /** \brief  Returns the evaluation points for basis derivatives on a subcell

      \code
      P  - num. evaluation points
      D  - spatial dimension
      \endcode

      \param  subCellDim  [in]   - dimension of the subcell
      \param  subCellId   [in]   - ordinal of the subcell defined by cell topology

      \return a rank-2 view (P,D) containing the basis derivatives evaluation points on the selected subcell
   */
  view_type getBasisDerivEvalPoints(const ordinal_type subCellDim, const ordinal_type subCellId) {
    return basisDerivCubPoints[subCellDim][subCellId];
  }


  /** \brief  Returns the points where to evaluate the target function on a subcell

      \code
      P  - num. evaluation points
      D  - spatial dimension
      \endcode

      \param  subCellDim  [in]  - dimension of the subcell
      \param  subCellId   [in]  - ordinal of the subcell defined by cell topology

      \return a rank-2 view (P,D) containing the target evaluation points on the selected subcell
   */
  view_type getTargetEvalPoints(const ordinal_type subCellDim, const ordinal_type subCellId) {
    return targetCubPoints[subCellDim][subCellId];
  }


  /** \brief  Returns the points where to evaluate the derivatives of the target function on a subcell

      \code
      P  - num. evaluation points
      D  - spatial dimension
      \endcode

      \param  subCellDim  [in]  - dimension of the subcell
      \param  subCellId   [in]  - ordinal of the subcell defined by cell topology

      \return a rank-2 view (P,D) containing the target derivatives evaluation points on the selected subcell
   */
  view_type getTargetDerivEvalPoints(const ordinal_type subCellDim, const ordinal_type subCellId) {
    return targetDerivCubPoints[subCellDim][subCellId];
  }


  /** \brief  Returns the basis/target evaluation points on a subcell

      \code
      P  - num. evaluation points
      D  - spatial dimension
      \endcode

      \param  subCellDim    [in]  - dimension of the subcell
      \param  subCellId     [in]  - ordinal of the subcell defined by cell topology
      \param  evalPointType [in]  - enum selecting whether the points should be computed for the basis
                                    functions or for the target function

      \return a rank-2 view (P,D) containing the basis/target function evaluation points on the selected subcell
   */
  view_type getEvalPoints(const ordinal_type subCellDim, const ordinal_type subCellId, EvalPointsType type) const{
    if(type == BASIS)
      return basisCubPoints[subCellDim][subCellId];
    else
      return targetCubPoints[subCellDim][subCellId];
  }

  /** \brief  Returns the evaluation points for basis/target derivatives on a subcell

      \code
      P  - num. evaluation points
      D  - spatial dimension
      \endcode

      \param  subCellDim    [in]  - dimension of the subcell
      \param  subCellId     [in]  - ordinal of the subcell defined by cell topology
      \param  evalPointType [in] - enum selecting whether the points should be computed for the basis
                                   functions or for the target function

      \return a rank-2 view (P,D) containing the basis/target  derivatives evaluation points on the selected subcell
   */
  view_type getDerivEvalPoints(const ordinal_type subCellDim, const ordinal_type subCellId, EvalPointsType type) const{
    if(type == BASIS)
      return basisDerivCubPoints[subCellDim][subCellId];
    else
      return targetDerivCubPoints[subCellDim][subCellId];
  }



  /** \brief  Returns the basis evaluation weights on a subcell

      \code
      P  - num. evaluation points
      \endcode

      \param  subCellDim  [in]  - dimension of the subcell
      \param  subCellId   [in]  - ordinal of the subcell defined by cell topology

      \return a rank-1 view (P) containing the basis evaluation weights on the selected subcell
   */
  view_type getBasisEvalWeights(const ordinal_type subCellDim, const ordinal_type subCellId) {
    return basisCubWeights[subCellDim][subCellId];
  }


  /** \brief  Returns the basis derivatives evaluation weights on a subcell

      \code
      P  - num. evaluation points
      \endcode

      \param  subCellDim  [in]  - dimension of the subcell
      \param  subCellId   [in]  - ordinal of the subcell defined by cell topology

      \return a rank-1 view (P) containing the basis derivatives evaluation weights on the selected subcell
   */
  view_type getBasisDerivEvalWeights(const ordinal_type subCellDim, const ordinal_type subCellId) {
    return basisDerivCubWeights[subCellDim][subCellId];
  }


  /** \brief  Returns the function evaluation weights on a subcell

      \code
      P  - num. evaluation points
      \endcode

      \param  subCellDim  [in]  - dimension of the subcell
      \param  subCellId   [in]  - ordinal of the subcell defined by cell topology

      \return a rank-1 view (P) containing the target evaluation weights on the selected subcell
   */
  view_type getTargetEvalWeights(const ordinal_type subCellDim, const ordinal_type subCellId) {
    return targetCubWeights[subCellDim][subCellId];
  }


  /** \brief  Returns the function derivatives evaluation weights on a subcell

      \code
      P  - num. evaluation points
      \endcode

      \param  subCellDim  [in]  - dimension of the subcell
      \param  subCellId   [in]  - ordinal of the subcell defined by cell topology

      \return a rank-1 view (P) containing the target derivatives evaluation weights on the selected subcell
   */
  view_type getTargetDerivEvalWeights(const ordinal_type subCellDim, const ordinal_type subCellId) {
    return targetDerivCubWeights[subCellDim][subCellId];
  }


  /** \brief  Returns the range tag of the basis evaluation points subcells

      \return the range tag of the basis evaluation points on subcells
   */
  const range_tag getBasisPointsRange() const {
    return basisPointsRange;
  }


  /** \brief  Returns the range tag of the basis/target evaluation points in subcells
      \param  evalPointType [in] - enum selecting whether the points should be computed for the basis
                               functions or for the target function
      \return the range tag of the basis/target evaluation points on subcells
   */
  const range_tag getPointsRange(const EvalPointsType type) const {
    if(type == BASIS)
      return basisPointsRange;
    else
      return targetPointsRange;
  }


  /** \brief  Returns the range tag of the derivative evaluation points on subcell

      \return the range tag of the basis derivative evaluation points corresponding on subcell
   */
  const range_tag getBasisDerivPointsRange() const {
    return basisDerivPointsRange;
  }

  /** \brief  Returns the range tag of the basis/target derivative evaluation points on subcells
      \param  evalPointType [in] - enum selecting whether the points should be computed for the basis
                               functions or for the target function

      \return the range tag of the basis/target derivative evaluation points on subcells
   */
  const range_tag getDerivPointsRange(const EvalPointsType type) const {
    if(type == BASIS)
      return basisDerivPointsRange;
    else
      return targetDerivPointsRange;
  }


  /** \brief  Returns the range of the target function evaluation points on subcells

      \return the range of the target function evaluation points corresponding on subcells
   */
  const range_tag getTargetPointsRange() const {
    return targetPointsRange;
  }


  /** \brief  Returns the range tag of the target function derivative evaluation points on subcells

      \return the range of the target function derivative evaluation points corresponding on subcells
   */
  const range_tag getTargetDerivPointsRange() const {
    return targetDerivPointsRange;
  }

  /** \brief  Returns the key tag for subcells

      \return the key tag of the selected subcells
   */
  const key_tag getTopologyKey() const {
    return subCellTopologyKey;
  }




  /** \brief  Initialize the ProjectionStruct for L2 projections
      \param  cellBasis       [in]  - basis functions for the projection
      \param  targetCubDegree [in]  - degree of the cubature needed to integrate the target function
   */
  template<typename BasisPtrType>
  void createL2ProjectionStruct(const BasisPtrType cellBasis,
      const ordinal_type targetCubDegree);


  /** \brief  Initialize the ProjectionStruct for HGRAD projections
      \param  cellBasis          [in]  - HGRAD basis functions for the projection
      \param  targetCubDegree    [in]  - degree of the cubature needed to integrate the target function
      \param  targetGradCubDegre [in]  - degree of the cubature needed to integrate the derivative of target function
   */
  template<typename BasisPtrType>
  void createHGradProjectionStruct(const BasisPtrType cellBasis,
      const ordinal_type targetCubDegree,
      const ordinal_type targetGradCubDegre);


  /** \brief  Initialize the ProjectionStruct for HCURL projections
      \param  cellBasis          [in]  - HCURL basis functions for the projection
      \param  targetCubDegree    [in]  - degree of the cubature needed to integrate the target function
      \param  targetGradCubDegre [in]  - degree of the cubature needed to integrate the derivative of target function
   */
  template<typename BasisPtrType>
  void createHCurlProjectionStruct(const BasisPtrType cellBasis,
      const ordinal_type targetCubDegree,
      const ordinal_type targetCurlCubDegre);


  /** \brief  Initialize the ProjectionStruct for HDIV projections
      \param  cellBasis          [in]  - HDIV basis functions for the projection
      \param  targetCubDegree    [in]  - degree of the cubature needed to integrate the target function
      \param  targetGradCubDegre [in]  - degree of the cubature needed to integrate the derivative of target function
   */
  template<typename BasisPtrType>
  void createHDivProjectionStruct(const BasisPtrType cellBasis,
      const ordinal_type targetCubDegree,
      const ordinal_type targetDivCubDegre);

  /** \brief  Initialize the ProjectionStruct for HVOL (local-L2) projection
      \param  cellBasis        [in]  - HVOL basis functions for the projection
      \param  targetCubDegree  [in]  - degree of the cubature needed to integrate the target function
   */
  template<typename BasisPtrType>
  void createHVolProjectionStruct(const BasisPtrType cellBasis,
      const ordinal_type targetCubDegree);

  key_tag subCellTopologyKey;
  range_tag basisPointsRange;
  range_tag basisDerivPointsRange;
  range_tag targetPointsRange;
  range_tag targetDerivPointsRange;
  view_tag basisCubPoints;
  view_tag basisCubWeights;
  view_tag basisDerivCubPoints;
  view_tag basisDerivCubWeights;
  view_tag targetCubPoints;
  view_tag targetCubWeights;
  view_tag targetDerivCubPoints;
  view_tag targetDerivCubWeights;
  ordinal_type numBasisEvalPoints;
  ordinal_type numBasisDerivEvalPoints;
  ordinal_type numTargetEvalPoints;
  ordinal_type numTargetDerivEvalPoints;
  ordinal_type maxNumBasisEvalPoints;
  ordinal_type maxNumTargetEvalPoints;
  ordinal_type maxNumBasisDerivEvalPoints;
  ordinal_type maxNumTargetDerivEvalPoints;
};

}
}
#include "Intrepid2_ProjectionStructDef.hpp"
#endif





