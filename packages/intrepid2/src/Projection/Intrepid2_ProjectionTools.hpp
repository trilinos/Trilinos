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

/** \file   Intrepid2_ProjectionTools.hpp
    \brief  Header file for the Intrepid2::Experimental::ProjectionTools.
    \author Created by Mauro Perego
 */
#ifndef __INTREPID2_PROJECTIONTOOLS_HPP__
#define __INTREPID2_PROJECTIONTOOLS_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Shards_CellTopology.hpp"
#include "Shards_BasicTopologies.hpp"

#include "Intrepid2_PointTools.hpp"

#include "Intrepid2_Basis.hpp"

// -- HGRAD family
#include "Intrepid2_HGRAD_LINE_Cn_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_Cn_FEM.hpp"

#include "Intrepid2_HGRAD_TRI_Cn_FEM.hpp"
#include "Intrepid2_HGRAD_TET_Cn_FEM.hpp"

// -- HCURL family
#include "Intrepid2_HCURL_QUAD_In_FEM.hpp"
#include "Intrepid2_HCURL_HEX_In_FEM.hpp"

#include "Intrepid2_HCURL_TRI_In_FEM.hpp"
#include "Intrepid2_HCURL_TET_In_FEM.hpp"
#include "Intrepid2_HVOL_LINE_Cn_FEM.hpp"

// -- HDIV family
#include "Intrepid2_HDIV_QUAD_In_FEM.hpp"
#include "Intrepid2_HDIV_HEX_In_FEM.hpp"

#include "Intrepid2_HDIV_TRI_In_FEM.hpp"
#include "Intrepid2_HDIV_TET_In_FEM.hpp"
#include "Intrepid2_HVOL_TRI_Cn_FEM.hpp"

// -- Lower order family
#include "Intrepid2_HCURL_QUAD_I1_FEM.hpp"
#include "Intrepid2_HCURL_TRI_I1_FEM.hpp"

#include "Intrepid2_HDIV_QUAD_I1_FEM.hpp"
#include "Intrepid2_HDIV_TRI_I1_FEM.hpp"

#include "Intrepid2_HCURL_HEX_I1_FEM.hpp"
#include "Intrepid2_HCURL_TET_I1_FEM.hpp"
#include "Intrepid2_HCURL_WEDGE_I1_FEM.hpp"

#include "Intrepid2_HDIV_HEX_I1_FEM.hpp"
#include "Intrepid2_HDIV_TET_I1_FEM.hpp"
#include "Intrepid2_HDIV_WEDGE_I1_FEM.hpp"

#include "Teuchos_LAPACK.hpp"
#include "Intrepid2_OrientationTools.hpp"

#include "Intrepid2_ProjectionStruct.hpp"

#include "KokkosBatched_QR_Serial_Internal.hpp"
#include "KokkosBatched_ApplyQ_Serial_Internal.hpp"
#include "KokkosBatched_Trsv_Serial_Internal.hpp"
#include "KokkosBatched_Util.hpp"

namespace Intrepid2 {

namespace Experimental {



/** \class  Intrepid2::Experimental::ProjectionTools
    \brief  A class providing static members to perform projection-based interpolations:

    This class provides tools to perform projection-based interpolations of a target function
    \f$f\f$ that satisfy the commutativity property of the De Rham complex and that have optimal
    accuracy w.r.t. the corresponding norms (\f$H^1\f$, \f$H^{\text{div}}\f$, \f$H^{\text{curl}}\f$, \f$L^2\f$).

    The projection-based interpolation \f$\Pi_h^X\f$ of the target function \f$f\f$ reads
    \f[
    \Pi_h^X(f) = \sum_i \alpha^f_i \phi_i, \quad X \in \{\text{grad},\, \text{div},\, \text{curl},\, \text{vol}\},
    \f]
    where \f$\{\phi_i\}\f$ is the basis of the finite element, \f$\alpha_i^f\f$ are the
    <var><b>basisCoeffs</b></var>.



    It also provides tools to perform a local L2 projection into HGrad, HCurl, HDiv and L2 fields.
    This projection does not satisfy the properties of the projection-based interpolations, but it
    is simpler and does not require to evaluate the derivatives of the target functions.

    Use:
    1. create a ProjectionStruct object
    2. allocate views for storing the points where to evaluate the target function and its derivatives
    3. evaluate the points/weights using one of the methods
       <var><b>getHGradEvaluationPoints</b></var>,
       <var><b>getHCURLEvaluationPoints</b></var>,
       <var><b>getHDivEvaluationPoints</b></var>,
       <var><b>getHVolEvaluationPoints</b></var>, or
       <var><b>getL2EvaluationPoints</b></var>
    4. Map to the physical elements the evaluation points,
       evaluate the target function and its derivatives at these points and
       map them back (inverse of pullback operator) to the reference points.
       Note: if the target function is defined at reference element (e.g. in case the target function is a
       combination of basis functions) there is no need to map points and functions between the reference
       and physical spaces, but one needs to simply evaluate the target function and its derivatives at the
       evaluation points
    5. Evaluate the basis coefficients using one of the methods
       <var><b>getHGradBasisCoeffs</b></var>,
       <var><b>getHCURLBasisCoeffs</b></var>,
       <var><b>getHDivBasisCoeffs</b></var>,
       <var><b>getHVolBasisCoeffs</b></var>, or
       <var><b>getL2BasisCoeffs</b></var>

    \remark The projections are performed at the oriented reference element. Therefore, the target function \f$f\f$,
            which is contravariant, needs to be mapped back to the reference element (using inverse of pullback operator)
            before calling <var><b>getXXXBasisCoeffs</b></var>.

    \remark The projections are implemented as in "Demkowics et al., Computing with Hp-Adaptive Finite Elements,
            Vol. 2, Chapman & Hall/CRC 2007". However, the projections along edges for HGrad and HCurl elements have been
            performed on the \f$H^1\f$ seminorm and the \f$L^2\f$ norm respectively, instead of on the \f$L^2\f$  and
            \f$H^{-1}\f$ and norms. This requires more regularity of the target function.

    \todo  There is room for significant improvement.
           One could separate the computation of the basis function values and derivatives from the functions getXXXBasisCoeffs,
           so that they can be stored and reused for projecting other target functions.
           Similarly one could store all the QR factorizations and reuse them for other target functions.
           For internal evaluation points (that are not affected by orientation) one could compute the QR factorization on the reference cell
           and then use on all the cells.

           Note: Other algorithmic improvements could be enabled by accessing the implementation of the orientation tools,
           however, we preferred the projections to work with any orientation, and assuming only that internal basis functions are not affected by
           the orientation.
 */

template<typename ExecSpaceType>
class ProjectionTools {
public:

  using EvalPointsType = typename ProjectionStruct<ExecSpaceType, double>::EvalPointsType;


  /** \brief  Computes evaluation points for L2 projection

      \code
      C  - num. cells
      P  - num. evaluation points
      D  - spatial dimension
      \endcode

      \param  evaluationPoints [out] - rank-3 view (C,P,D) containing the coordinates of the evaluation
                                       points for the projection at each cell
      \param  cellOrientations [in]  - rank-1 view (C) containing the Orientation objects at each cell
      \param  cellBasis        [in]  - pointer to the basis for the projection
      \param  projStruct       [in]  - pointer to ProjectionStruct object
      \param  evalPointType    [in]  - enum selecting whether the points should be computed for the basis
                                       functions or for the target function
   */
  template<typename BasisType,
  typename ortValueType,       class ...ortProperties>
  static void
  getL2EvaluationPoints(typename BasisType::ScalarViewType evaluationPoints,
      const Kokkos::DynRankView<ortValueType,   ortProperties...>  cellOrientations,
      const BasisType* cellBasis,
      ProjectionStruct<ExecSpaceType, typename BasisType::scalarType> * projStruct,
      const EvalPointsType evalPointType = EvalPointsType::TARGET
  );

  /** \brief  Computes the basis coefficients of the L2 projection of the target function

      \code
      C  - num. cells
      F  - num. fields
      P  - num. evaluation points
      D  - spatial dimension
      \endcode

      \param  basisCoeffs         [out] - rank-2 view (C,F) containing the basis coefficients
      \param  targetAtEvalPoints  [in]  - variable rank view containing the values of the target function
                                          evaluated at the evaluation points
      \param  cellOrientations    [in]  - 1-rank view (C) containing the Orientation objects at each cell
      \param  cellBasis           [in]  - pointer to the basis for the projection
      \param  projStruct          [in]  - pointer to ProjectionStruct object

      \remark targetAtEvalPoints has rank 2 (C,P) for HGRAD and HVOL elements, and rank 3 (C,P,D)
              for HCURL and HDIV elements
   */
  template<typename basisCoeffsValueType, class ...basisCoeffsProperties,
  typename funValsValueType, class ...funValsProperties,
  typename BasisType,
  typename ortValueType,       class ...ortProperties>
  static void
  getL2BasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
      const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtEvalPoints,
      const typename BasisType::ScalarViewType evaluationPoints,
      const Kokkos::DynRankView<ortValueType,   ortProperties...>  cellOrientations,
      const BasisType* cellBasis,
      ProjectionStruct<ExecSpaceType, typename BasisType::scalarType> * projStruct);


  /** \brief  Computes evaluation points for HGrad projection

      \code
      C  - num. cells
      P1 - num. evaluation points
      P2 - num. evaluation points for derivatives
      D  - spatial dimension
      \endcode

      \param  evaluationPoints [out] - rank-3 view (C,P1,D) containing the coordinates of the evaluation
                                       points, at each cell
      \param  gradEvalPoints   [in]  - rank-3 view (C,P2,D) containing the coordinates of the points
                                       where to evaluate the function gradients, at each cell
      \param  cellOrientations [in]  - rank-1 view (C) containing the Orientation objects at each cell
      \param  cellBasis        [in]  - pointer to the HGRAD basis for the projection
      \param  projStruct       [in]  - pointer to ProjectionStruct object
      \param  evalPointType    [in]  - enum selecting whether the points should be computed for the basis
                                       functions or for the target function
   */
  template<typename BasisType,
  typename ortValueType,       class ...ortProperties>
  static void
  getHGradEvaluationPoints(typename BasisType::ScalarViewType evaluationPoints,
      typename BasisType::ScalarViewType gradEvalPoints,
      const Kokkos::DynRankView<ortValueType,   ortProperties...>  cellOrientations,
      const BasisType* cellBasis,
      ProjectionStruct<ExecSpaceType, typename BasisType::scalarType> * projStruct,
      const EvalPointsType evalPointType = EvalPointsType::TARGET
  );

  /** \brief  Computes the basis coefficients of the HGrad projection of the target function

      \code
      C  - num. cells
      F  - num. fields
      P1 - num. evaluation points
      P2 - num. evaluation points for derivatives
      D  - spatial dimension
      \endcode

      \param  basisCoeffs                [out] - rank-2 view (C,F) containing the basis coefficients
      \param  targetAtEvalPoints         [in]  - rank-2 view (C,P1) containing the values of the target function
                                                 evaluated at the evaluation points
      \param  targetGradAtGradEvalPoints [in]  - rank-3 view (C,P2,D) view containing the values of the gradient
                                                 of the target function evaluated at the evaluation points
      \param  evaluationPoints           [in]  - rank-3 view (C,P1,D) containing the coordinates of the evaluation
                                                 points, at each cell
      \param  gradEvalPoints             [in]  - rank-3 view (C,P2,D) containing the coordinates of the points
                                                 where to evaluate the function gradients, at each cell
      \param  cellOrientations           [in]  - 1-rank view (C) containing the Orientation objects at each cell
      \param  cellBasis                  [in]  - pointer to the HGRAD basis for the projection
      \param  projStruct                 [in]  - pointer to ProjectionStruct object
   */
  template<typename basisCoeffsValueType, class ...basisCoeffsProperties,
  typename funValsValueType, class ...funValsProperties,
  typename BasisType,
  typename ortValueType,       class ...ortProperties>
  static void
  getHGradBasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
      const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtEvalPoints,
      const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetGradAtGradEvalPoints,
      const typename BasisType::ScalarViewType evaluationPoints,
      const typename BasisType::ScalarViewType gradEvalPoints,
      const Kokkos::DynRankView<ortValueType,   ortProperties...>  cellOrientations,
      const BasisType* cellBasis,
      ProjectionStruct<ExecSpaceType, typename BasisType::scalarType> * projStruct);


  /** \brief  Computes evaluation points for HCurl projection

      \code
      C  - num. cells
      P1 - num. evaluation points
      P2 - num. evaluation points for derivatives
      D  - spatial dimension
      \endcode

      \param  evaluationPoints [out] - rank-3 view (C,P1,D) containing the coordinates of the evaluation
                                       points for the projection at each cell
      \param  curlEvalPoints   [in]  - rank-3 view (C,P2,D) containing the coordinates of the points
                                       where to evaluate the function curls, at each cell
      \param  cellOrientations [in]  - rank-1 view (C) containing the Orientation objects at each cell
      \param  cellBasis        [in]  - pointer to the HCURL basis for the projection
      \param  projStruct       [in]  - pointer to ProjectionStruct object
      \param  evalPointType    [in]  - enum selecting whether the points should be computed for the basis
                                       functions or for the target function
   */
  template<typename BasisType,
  typename ortValueType,       class ...ortProperties>
  static void
  getHCurlEvaluationPoints(typename BasisType::ScalarViewType evaluationPoints,
      typename BasisType::ScalarViewType curlEvalPoints,
      const Kokkos::DynRankView<ortValueType,   ortProperties...>  cellOrientations,
      const BasisType* cellBasis,
      ProjectionStruct<ExecSpaceType, typename BasisType::scalarType> * projStruct,
      const EvalPointsType evalPointType = EvalPointsType::TARGET
  );

  /** \brief  Computes the basis coefficients of the HCurl projection of the target function

      \code
      C  - num. cells
      F  - num. fields
      P1 - num. evaluation points
      P2 - num. evaluation points for derivatives
      D  - spatial dimension
      \endcode

      \param  basisCoeffs                [out] - rank-2 view (C,F) containing the basis coefficients
      \param  targetAtEvalPoints         [in]  - rank-3 view (C,P1,D) containing the values of the target function
                                                 evaluated at the evaluation points
      \param  targetcurlAtCurlEvalPoints [in]  - variable rank view containing the values of the curl of the target
                                                 function evaluated at the evaluation points
      \param  evaluationPoints           [in]  - rank-3 view (C,P1,D) containing the coordinates of the evaluation
                                                 points for the projection at each cell
      \param  curlEvalPoints             [in]  - rank-3 view (C,P2,D) containing the coordinates of the points
                                                 where to evaluate the function curls, at each cell
      \param  cellOrientations           [in]  - 1-rank view (C) containing the Orientation objects at each cell
      \param  cellBasis                  [in]  - pointer to the HCURL basis for the projection
      \param  projStruct                 [in]  - pointer to ProjectionStruct object

      \remark targetAtCurlEvalPoints has rank 2 (C,P2) in 2D, and rank 3 (C,P2,D) in 3D
   */
  template<typename basisCoeffsValueType, class ...basisCoeffsProperties,
  typename funValsValueType, class ...funValsProperties,
  typename BasisType,
  typename ortValueType,       class ...ortProperties>
  static void
  getHCurlBasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
      const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtEvalPoints,
      const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetCurlAtCurlEvalPoints,
      const typename BasisType::ScalarViewType evaluationPoints,
      const typename BasisType::ScalarViewType curlEvalPoints,
      const Kokkos::DynRankView<ortValueType,   ortProperties...>  cellOrientations,
      const BasisType* cellBasis,
      ProjectionStruct<ExecSpaceType, typename BasisType::scalarType> * projStruct);


  /** \brief  Computes evaluation points for HDiv projection

      \code
      C  - num. cells
      P1 - num. evaluation points
      P2 - num. evaluation points for derivatives
      D  - spatial dimension
      \endcode

      \param  evaluationPoints [out] - rank-3 view (C,P1,D) containing the coordinates of the evaluation
                                       points for the projection at each cell
      \param  divEvalPoints    [in]  - rank-3 view (C,P2,D) containing the coordinates of the points
                                       where to evaluate the function divergence, at each cell
      \param  cellOrientations [in]  - rank-1 view (C) containing the Orientation objects at each cell
      \param  cellBasis        [in]  - pointer to the HDIV basis for the projection
      \param  projStruct       [in]  - pointer to ProjectionStruct object
      \param  evalPointType    [in]  - enum selecting whether the points should be computed for the basis
                                       functions or for the target function
   */
  template<typename BasisType,
  typename ortValueType,       class ...ortProperties>
  static void
  getHDivEvaluationPoints(typename BasisType::ScalarViewType evaluationPoints,
      typename BasisType::ScalarViewType divEvalPoints,
      const Kokkos::DynRankView<ortValueType,   ortProperties...>  cellOrientations,
      const BasisType* cellBasis,
      ProjectionStruct<ExecSpaceType, typename BasisType::scalarType> * projStruct,
      const EvalPointsType evalPointType = EvalPointsType::TARGET
  );

  /** \brief  Computes the basis coefficients of the HDiv projection of the target function

      \code
      C  - num. cells
      F  - num. fields
      P1 - num. evaluation points
      P2 - num. evaluation points for derivatives
      D  - spatial dimension
      \endcode

      \param  basisCoeffs              [out] - rank-2 view (C,F) containing the basis coefficients
      \param  targetAtEvalPoints       [in]  - rank-3 view (C,P1,D) containing the values of the target function
                                               evaluated at the evaluation points
      \param  targetDivAtDivEvalPoints [in]  - rank-2 view (C,P2) view containing the values of the divergence
                                               of the target function evaluated at the evaluation points
      \param  evaluationPoints         [in]  - rank-3 view (C,P1,D) containing the coordinates of the evaluation
                                               points, at each cell
      \param  divEvalPoints            [in]  - rank-3 view (C,P2,D) containing the coordinates of the points
                                               where to evaluate the function divergence, at each cell
      \param  cellOrientations         [in]  - 1-rank view (C) containing the Orientation objects at each cell
      \param  cellBasis                [in]  - pointer to the HDIV basis for the projection
      \param  projStruct               [in]  - pointer to ProjectionStruct object
   */
  template<typename basisCoeffsValueType, class ...basisCoeffsProperties,
  typename funValsValueType, class ...funValsProperties,
  typename BasisType,
  typename ortValueType,       class ...ortProperties>
  static void
  getHDivBasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
      const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtEvalPoints,
      const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetDivAtDivEvalPoints,
      const typename BasisType::ScalarViewType evaluationPoints,
      const typename BasisType::ScalarViewType divEvalPoints,
      const Kokkos::DynRankView<ortValueType,   ortProperties...>  cellOrientations,
      const BasisType* cellBasis,
      ProjectionStruct<ExecSpaceType, typename BasisType::scalarType> * projStruct);

  /** \brief  Computes evaluation points for HVol projection

      \code
      C  - num. cells
      P  - num. evaluation points
      D  - spatial dimension
      \endcode

      \param  evaluationPoints [out] - rank-3 view (C,P,D) containing the coordinates of the evaluation
                                       points, at each cell
      \param  cellOrientations [in]  - rank-1 view (C) containing the Orientation objects at each cell
      \param  cellBasis        [in]  - pointer to the HVOL basis for the projection
      \param  projStruct       [in]  - pointer to ProjectionStruct object
      \param  evalPointType    [in]  - enum selecting whether the points should be computed for the basis
                                       functions or for the target function
   */
  template<typename BasisType,
  typename ortValueType,       class ...ortProperties>
  static void
  getHVolEvaluationPoints(typename BasisType::ScalarViewType evaluationPoints,
      const Kokkos::DynRankView<ortValueType, ortProperties...>  cellOrientations,
      const BasisType* cellBasis,
      ProjectionStruct<ExecSpaceType, typename BasisType::scalarType> * projStruct,
      const EvalPointsType evalPointType = EvalPointsType::TARGET
  );

  /** \brief  Computes the basis coefficients of the HVol projection of the target function

      \code
      C  - num. cells
      F  - num. fields
      P  - num. evaluation points
      D  - spatial dimension
      \endcode

      \param  basisCoeffs           [out] - rank-2 view (C,F) containing the basis coefficients
      \param  targetAtEvalPoints    [in]  - rank-2 view (C,P) containing the values of the target function
                                            evaluated at the evaluation points
      \param  evaluationPoints      [in]  - rank-3 view (C,P,D) containing the coordinates of the evaluation
                                            points, at each cell
      \param  cellOrientations      [in]  - 1-rank view (C) containing the Orientation objects at each cell
      \param  cellBasis             [in]  - pointer to the HGRAD basis for the projection
      \param  projStruct            [in]  - pointer to ProjectionStruct object
   */
  template<typename basisCoeffsValueType, class ...basisCoeffsProperties,
  typename funValsValueType, class ...funValsProperties,
  typename BasisType,
  typename ortValueType,       class ...ortProperties>
  static void
  getHVolBasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
      const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtEvalPoints,
      const typename BasisType::ScalarViewType evaluationPoints,
      const Kokkos::DynRankView<ortValueType,   ortProperties...>  cellOrientations,
      const BasisType* cellBasis,
      ProjectionStruct<ExecSpaceType, typename BasisType::scalarType> * projStruct);



  /** \brief  Functor to solve a square system A x = b on each cell using QR method implemented in KokkosKernels
              A is expected to be saddle a point (KKT) matrix of the form [C B; B^T 0],
              where C has size nxn and B nxm, with n>0, m>=0.
              B^T is copied from B, so one does not have to define the B^T portion of A.
              b will contain the solution x.
              The first n-entries of x are copied into the provided basis coefficients using the provided indexing.
   */
  template<typename ViewType1, typename ViewType2, typename ViewType3, typename ViewType4>
  struct SolveSystem {
    ViewType1 basisCoeffs_; // rank-2 view (C,F) containing the basis coefficients on each cell
    ViewType2 elemMat_;     // rank-3 view (C,P,P) containing the element matrix on each cell
    ViewType2 elemRhs_;     // rank-2 view (C,P) containing the element rhs on each cell
    ViewType2 tau_;         // rank-2 view (C,P) used to store the QR factorization
    ViewType3 w_;           // rank-2 view (C,P) used has a workspace (needs to be of Layout Right)

    const ViewType4 elemDof_; // rank-1 view having dimension n, containing the basis numbering
    ordinal_type n_, m_;      // basis cardinality and dimension of the constraint of the KKT system

    /** \brief  Functor constructor

        \code
          C  - num. cells
          P  - num. evaluation points
        \endcode


        \param  basisCoeffs           [out]     - rank-2 view (C,F) containing the basis coefficients
        \param  elemMat               [in/out]  - rank-3 view (C,P,P) containing the element matrix of size
                                                  numCells x (n+m)x(n+m) on each cell
                                                  it will be overwritten.
        \param  elemRhs               [in/out]  - rank-2 view (C,P) containing the element rhs on each cell
                                                  of size numCells x (n+m)
                                                  it will contain the solution of the system on output
        \param  tau                   [out]     - rank-2 view (C,P) used to store the QR factorization
                                                  size: numCells x (n+m)
        \param  w                     [out]     - rank-2 view (C,P) used has a workspace
                                                  Layout Right, size: numCells x (n+m)
        \param  elemDof               [in]      - rank-1 view having dimension n, containing the basis numbering
        \param  n                     [in]      - ordinal_type, basis cardinality
        \param  m                     [in]      - ordinal_type, dimension of the constraint of the KKT system
     */
    SolveSystem (ViewType1 basisCoeffs, ViewType2 elemMat, ViewType2 elemRhs, ViewType2 tau,
        ViewType3 w,const  ViewType4 elemDof, ordinal_type n, ordinal_type m=0) :
          basisCoeffs_(basisCoeffs), elemMat_(elemMat), elemRhs_(elemRhs),
          tau_(tau), w_(w), elemDof_(elemDof), n_(n), m_(m){};


    void
    KOKKOS_INLINE_FUNCTION
    operator()(const ordinal_type ic) const {
      auto A = Kokkos::subview(elemMat_, ic, Kokkos::ALL(), Kokkos::ALL());
      auto b = Kokkos::subview(elemRhs_, ic, Kokkos::ALL());
      auto tau = Kokkos::subview(tau_, ic, Kokkos::ALL());
      auto w = Kokkos::subview(w_, ic, Kokkos::ALL());

      for(ordinal_type i=n_; i<n_+m_; ++i)
        for(ordinal_type j=0; j<n_; ++j)
          A(i,j) = A(j,i);

      //computing QR of A. QR is saved in A and tau
      KokkosBatched::SerialQR_Internal::invoke(A.extent(0), A.extent(1),
          A.data(), A.stride_0(), A.stride_1(), tau.data(), tau.stride_0(), w.data());

      //b'*Q -> b
      KokkosBatched::SerialApplyQ_RightNoTransForwardInternal::invoke(
          1, A.extent(0), A.extent(1),
          A.data(),  A.stride_0(), A.stride_1(),
          tau.data(), tau.stride_0(),
          b.data(),  1, b.stride_0(),
          w.data());

      // R^{-1} b -> b
      KokkosBatched::SerialTrsvInternalUpper<KokkosBatched::Algo::Trsv::Unblocked>::invoke(false,
          A.extent(0),
          1.0,
          A.data(), A.stride_0(), A.stride_1(),
          b.data(),  b.stride_0());

      //scattering b into the basis coefficients
      for(ordinal_type i=0; i<n_; ++i){
        basisCoeffs_(ic,elemDof_(i)) = b(i);
      }
    }
  };
};
}
}


// include templated function definitions
#include "Intrepid2_ProjectionToolsDefL2.hpp"
#include "Intrepid2_ProjectionToolsDefHGRAD.hpp"
#include "Intrepid2_ProjectionToolsDefHCURL.hpp"
#include "Intrepid2_ProjectionToolsDefHDIV.hpp"
#include "Intrepid2_ProjectionToolsDefHVOL.hpp"

#endif





