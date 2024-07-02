// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_ProjectionTools.hpp
    \brief  Header file for the Intrepid2::ProjectionTools.
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

#include "Intrepid2_NodalBasisFamily.hpp"

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

#ifdef HAVE_INTREPID2_KOKKOSKERNELS
#include "KokkosBatched_QR_Serial_Internal.hpp"
#include "KokkosBatched_ApplyQ_Serial_Internal.hpp"
#include "KokkosBatched_Trsv_Serial_Internal.hpp"
#include "KokkosBatched_Util.hpp"
#endif

namespace Intrepid2 {

/** \class  Intrepid2::ProjectionTools
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
    2. get the evaluation points where to evaluate the target function and its derivative using
       the ProjectionStruct methods <var><b>getAllEvalPoints</b></var> 
       and <var><b>getAllDerivEvalPoints</b></var>
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

    \todo  There is room for improvement.
           One could separate the computation of the basis function values and derivatives from the functions getXXXBasisCoeffs,
           so that they can be stored and reused for projecting other target functions.
           Similarly one could store all the QR factorizations and reuse them for other target functions.
 */

template<typename DeviceType>
class ProjectionTools {
public:
  using ExecSpaceType = typename DeviceType::execution_space;
  using MemSpaceType = typename DeviceType::memory_space;
  using EvalPointsType = typename ProjectionStruct<DeviceType, double>::EvalPointsType;


  /** \brief  Computes the basis coefficients of the L2 projection of the target function

      \code
      C  - num. cells
      F  - num. fields
      P  - num. evaluation points
      D  - spatial dimension
      \endcode

      \param  basisCoeffs         [out] - rank-2 view (C,F) containing the basis coefficients
      \param  targetAtEvalPoints  [in]  - variable rank view containing the values of the target function
                                          evaluated at the evaluation points given by <var><b>projStruct->getAllEvalPoints()</var></b>
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
      const Kokkos::DynRankView<ortValueType,   ortProperties...>  cellOrientations,
      const BasisType* cellBasis,
      ProjectionStruct<DeviceType, typename BasisType::scalarType> * projStruct);


  /** \brief  Computes evaluation points for local L2 projection
     for broken HGRAD HCURL HDIV and HVOL spaces

     This function simply perform an L2 projection in each cell with no guarantee
     of preserving continuity across cells

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
  getL2DGBasisCoeffs(Kokkos::DynRankView<basisCoeffsValueType,basisCoeffsProperties...> basisCoeffs,
      const Kokkos::DynRankView<funValsValueType,funValsProperties...> targetAtEvalPoints,
      const Kokkos::DynRankView<ortValueType,   ortProperties...>  cellOrientations,
      const BasisType* cellBasis,
      ProjectionStruct<DeviceType, typename BasisType::scalarType> * projStruct);


  /** \brief  Computes evaluation points for local L2 projection
     for broken HGRAD HCURL HDIV and HVOL spaces

     This function simply perform an L2 projection in each cell with no guarantee
     of preserving continuity across cells. It also does not account for orientation.

      \code
      C  - num. cells
      F  - num. fields
      P  - num. evaluation points
      D  - spatial dimension
      \endcode

      \param  basisCoeffs         [out] - rank-2 view (C,F) containing the basis coefficients
      \param  targetAtEvalPoints  [in]  - variable rank view containing the values of the target function
                                          evaluated at the evaluation points
      \param  cellBasis           [in]  - pointer to the basis for the projection
      \param  projStruct          [in]  - pointer to ProjectionStruct object

      \remark targetAtEvalPoints has rank 2 (C,P) for HGRAD and HVOL elements, and rank 3 (C,P,D)
              for HCURL and HDIV elements
   */
  template<typename basisViewType, typename targetViewType, typename BasisType>
  static void
  getL2DGBasisCoeffs(basisViewType basisCoeffs,
      const targetViewType targetAtTargetEPoints,
      const BasisType* cellBasis,
      ProjectionStruct<DeviceType, typename BasisType::scalarType> * projStruct);

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
                                                 evaluated at the evaluation points given by
                                                 <var><b>projStruct->getAllEvalPoints()</var></b>
      \param  targetGradAtGradEvalPoints [in]  - rank-3 view (C,P2,D) view containing the values of the gradient
                                                 of the target function evaluated at the evaluation points given by
                                                 <var><b>projStruct->getAllDerivEvalPoints()</var></b>
      \param  cellOrientations           [in]  - 1-rank view (C) containing the Orientation objects at each cell
      \param  cellBasis                  [in]  - pointer to the HGRAD basis for the projection
      \param  projStruct                 [in]  - pointer to ProjectionStruct object
  */
  template<class BasisCoeffsViewType, class TargetValueViewType, class TargetGradViewType,
           class BasisType, class OrientationViewType>
  static void
  getHGradBasisCoeffs(BasisCoeffsViewType basisCoeffs,
                      const TargetValueViewType targetAtEvalPoints,
                      const TargetGradViewType targetGradAtGradEvalPoints,
                      const OrientationViewType cellOrientations,
                      const BasisType* cellBasis,
                      ProjectionStruct<DeviceType, typename BasisType::scalarType> * projStruct);


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
                                                 evaluated at the evaluation points given by
                                                 <var><b>projStruct->getAllEvalPoints()</var></b>
      \param  targetcurlAtCurlEvalPoints [in]  - variable rank view containing the values of the curl of the target
                                                 function evaluated at the evaluation points given by
                                                 <var><b>projStruct->getAllDerivEvalPoints()</var></b>
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
      const Kokkos::DynRankView<ortValueType,   ortProperties...>  cellOrientations,
      const BasisType* cellBasis,
      ProjectionStruct<DeviceType, typename BasisType::scalarType> * projStruct);
  
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
                                               evaluated at the evaluation points given by
                                                 <var><b>projStruct->getAllEvalPoints()</var></b>
      \param  targetDivAtDivEvalPoints [in]  - rank-2 view (C,P2) view containing the values of the divergence
                                               of the target function evaluated at the evaluation points given by
                                                 <var><b>projStruct->getAllDerivEvalPoints()</var></b>
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
      const Kokkos::DynRankView<ortValueType,   ortProperties...>  cellOrientations,
      const BasisType* cellBasis,
      ProjectionStruct<DeviceType, typename BasisType::scalarType> * projStruct);


  /** \brief  Computes the basis coefficients of the HVol projection of the target function

      \code
      C  - num. cells
      F  - num. fields
      P  - num. evaluation points
      D  - spatial dimension
      \endcode

      \param  basisCoeffs           [out] - rank-2 view (C,F) containing the basis coefficients
      \param  targetAtEvalPoints    [in]  - rank-2 view (C,P) containing the values of the target function
                                            evaluated at the evaluation points given by
                                            <var><b>projStruct->getAllEvalPoints()</var></b>
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
      [[maybe_unused]] const Kokkos::DynRankView<ortValueType,   ortProperties...>  cellOrientations,
      const BasisType* cellBasis,
      ProjectionStruct<DeviceType, typename BasisType::scalarType> * projStruct);


  
/** \brief  Computes the L2 projection of a finite element field into a compatible finite element space

      \code
      C  - num. cells
      F  - num. fields
      P  - num. evaluation points
      D  - spatial dimension
      \endcode

      \param  dstCoeffs          [out] - rank-2 view (C,F) containing the basis coefficients of the projected field over the space generated by dstBasis
      \param  dstBasis           [in]  - pointer to the basis we want to project to
      \param  srcCoeffs          [in]  - rank-2 view (C,F) containing the basis coefficients of the field w.r.t. the original basis
      \param  srcBasis           [in]  - original basis
      \param  cellOrientations   [in]  - rank-1 view (C) containing the Orientation objects at each cell
   */
  template<typename dstViewType,
  typename dstBasisType,
  typename srcViewType,
  typename srcBasisType,
  typename OrientationViewType>
  static void
  projectField(dstViewType dstCoeffs,
      const dstBasisType* dstBasis,
      const srcViewType srcCoeffs,
      const srcBasisType* srcBasis,
      const OrientationViewType cellOrientations){


    INTREPID2_TEST_FOR_EXCEPTION(dstBasis->getFunctionSpace() != srcBasis->getFunctionSpace(), std::runtime_error, 
      "The source and destination bases are not compatible. They need to belong to the same function space"); 
    INTREPID2_TEST_FOR_EXCEPTION(dstBasis->getBaseCellTopology().getKey() != srcBasis->getBaseCellTopology().getKey(), std::runtime_error, 
      "The source and destination bases are not compatible. They do not have the same basic cell topology");  
      
    ProjectionStruct<DeviceType,typename srcBasisType::scalarType> projStruct;
    projStruct.createL2ProjectionStruct(dstBasis, srcBasis->getDegree());

    
    ordinal_type numCells = cellOrientations.extent(0);
    ordinal_type dim = srcBasis->getBaseCellTopology().getDimension();
    ordinal_type srcBasisCardinality = srcBasis->getCardinality();
    ordinal_type fieldDimension = (srcBasis->getFunctionSpace() == Intrepid2::FUNCTION_SPACE_HCURL || srcBasis->getFunctionSpace() == Intrepid2::FUNCTION_SPACE_HDIV) ? dim : 1;

    auto evaluationPoints = projStruct.getAllEvalPoints();
    ordinal_type numPoints = evaluationPoints.extent(0);

    using outViewType = Kokkos::DynRankView<typename srcBasisType::OutputValueType, DeviceType>;
    outViewType srcAtEvalPoints, refBasisAtEvalPoints, basisAtEvalPoints;
    if(fieldDimension == dim) {
      srcAtEvalPoints = outViewType("srcAtEvalPoints", numCells, numPoints, dim);
      refBasisAtEvalPoints = outViewType("refBasisAtEvalPoints", srcBasisCardinality, numPoints, dim);
      basisAtEvalPoints = outViewType("basisAtEvalPoints", numCells, srcBasisCardinality, numPoints, dim);
    } else {
      srcAtEvalPoints = outViewType("srcAtEvalPoints", numCells, numPoints);
      refBasisAtEvalPoints = outViewType("refBasisAtEvalPoints", srcBasisCardinality, numPoints);
      basisAtEvalPoints = outViewType("basisAtEvalPoints", numCells, srcBasisCardinality, numPoints);
    }
    
    srcBasis->getValues(refBasisAtEvalPoints,evaluationPoints);

    // Modify basis values to account for orientations
    OrientationTools<DeviceType>::modifyBasisByOrientation(basisAtEvalPoints,
        refBasisAtEvalPoints,
        cellOrientations,
        srcBasis);

    Kokkos::parallel_for(Kokkos::RangePolicy<typename DeviceType::execution_space>(0,numCells),
        KOKKOS_LAMBDA (const int &ic) {
      for(int j=0; j<numPoints; ++j) {
        for(int k=0; k<srcBasisCardinality; ++k) {
          for(int d=0; d<fieldDimension; ++d)
            srcAtEvalPoints.access(ic,j,d) += srcCoeffs(ic,k)*basisAtEvalPoints.access(ic,k,j,d);
        }
      }
    });
    ExecSpaceType().fence();

    getL2BasisCoeffs(dstCoeffs,
        srcAtEvalPoints,
        cellOrientations,
        dstBasis,
        &projStruct);      
  }



  /** \brief  Class to solve a square system A x = b on each cell
              A is expected to be saddle a point (KKT) matrix of the form [C B; B^T 0],
              where C has size nxn and B nxm, with n>0, m>=0.
              B^T is copied from B, so one does not have to define the B^T portion of A.
              b will contain the solution x.
              The first n-entries of x are copied into the provided basis coefficients using the provided indexing.
              The system is solved either with a QR factorization implemented in KokkosKernels or
              with Lapack GELS function.
   */
  struct ElemSystem {


    std::string systemName_;
    bool matrixIndependentOfCell_;

    /** \brief  Functor constructor
        \param  systemName               [in]     - string containing the name of the system (passed to parallel for)
        \param  matrixIndependentOfCell  [in]     - bool: whether the local cell matrix of the system changes from cell to cell
                                                          if true, the matrix factorization is preformed only on the first cell
                                                          and reused on other cells.
    */

    ElemSystem (std::string systemName, bool matrixIndependentOfCell) :
      systemName_(systemName), matrixIndependentOfCell_(matrixIndependentOfCell){};



    /** \brief  Solve the system and returns the basis coefficients
                solve the system either using Kokkos Kernel QR or Lapack GELS
                depending on whether Kokkos Kernel is enabled.

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
    template<typename ViewType1, typename ViewType2, typename ViewType3, typename ViewType4>
    void solve(ViewType1 basisCoeffs, ViewType2 elemMat, ViewType2 elemRhs, ViewType2 tau,
        ViewType3 w,const  ViewType4 elemDof, ordinal_type n, ordinal_type m=0) {
#ifdef HAVE_INTREPID2_KOKKOSKERNELS
      solveDevice(basisCoeffs, elemMat, elemRhs, tau,
                  w, elemDof, n, m);
#else
      solveHost(basisCoeffs, elemMat, elemRhs, tau,
          w, elemDof, n, m);
#endif

    }

    /** \brief  Parallel implementation of solve, using Kokkos Kernels QR factoriation
     */
#ifdef HAVE_INTREPID2_KOKKOSKERNELS
    template<typename ViewType1, typename ViewType2, typename ViewType3, typename ViewType4>
    void solveDevice(ViewType1 basisCoeffs, ViewType2 elemMat, ViewType2 elemRhs, ViewType2 taul,
        ViewType3 work,const  ViewType4 elemDof, ordinal_type n, ordinal_type m) {
      using HostDeviceType = Kokkos::Device<Kokkos::DefaultHostExecutionSpace,Kokkos::HostSpace>;

      ordinal_type numCells = basisCoeffs.extent(0);

      if(matrixIndependentOfCell_) {
        auto A0 = Kokkos::subview(elemMat, 0, Kokkos::ALL(), Kokkos::ALL());
        auto tau0 = Kokkos::subview(taul, 0, Kokkos::ALL());

        Kokkos::DynRankView<typename ViewType2::value_type, HostDeviceType> A0_host("A0_host", A0.extent(0),A0.extent(1));
        auto A0_device = Kokkos::create_mirror_view(typename DeviceType::memory_space(), A0_host);
        Kokkos::deep_copy(A0_device, A0);
        Kokkos::deep_copy(A0_host, A0_device);

        for(ordinal_type i=n; i<n+m; ++i)
          for(ordinal_type j=0; j<n; ++j)
            A0_host(i,j) = A0_host(j,i);

        Kokkos::DynRankView<typename ViewType2::value_type, HostDeviceType> tau0_host("A0_host", tau0.extent(0));
        auto tau0_device = Kokkos::create_mirror_view(typename DeviceType::memory_space(), tau0_host);
        auto w0_host = Kokkos::create_mirror_view(Kokkos::subview(work, 0, Kokkos::ALL()));

        //computing QR of A0. QR is saved in A0 and tau0
        KokkosBatched::SerialQR_Internal::invoke(A0_host.extent(0), A0_host.extent(1),
            A0_host.data(), A0_host.stride_0(), A0_host.stride_1(),
            tau0_host.data(), tau0_host.stride_0(), w0_host.data());

        Kokkos::deep_copy(A0_device, A0_host);
        Kokkos::deep_copy(A0, A0_device);
        Kokkos::deep_copy(tau0_device, tau0_host);
        Kokkos::deep_copy(tau0, tau0_device);

        Kokkos::parallel_for (systemName_,
            Kokkos::RangePolicy<ExecSpaceType, int> (0, numCells),
            KOKKOS_LAMBDA (const size_t ic) {
          auto w = Kokkos::subview(work, ic, Kokkos::ALL());

          auto b = Kokkos::subview(elemRhs, ic, Kokkos::ALL());

          //b'*Q0 -> b
          KokkosBatched::SerialApplyQ_RightForwardInternal::invoke(
              1, A0.extent(0), A0.extent(1),
              A0.data(),  A0.stride_0(), A0.stride_1(),
              tau0.data(), tau0.stride_0(),
              b.data(),  1, b.stride_0(),
              w.data());

          // R0^{-1} b -> b
          KokkosBatched::SerialTrsvInternalUpper<KokkosBatched::Algo::Trsv::Unblocked>::invoke(false,
              A0.extent(0),
              1.0,
              A0.data(), A0.stride_0(), A0.stride_1(),
              b.data(),  b.stride_0());

          //scattering b into the basis coefficients
          for(ordinal_type i=0; i<n; ++i){
            basisCoeffs(ic,elemDof(i)) = b(i);
          }
        });

      } else {

        Kokkos::parallel_for (systemName_,
            Kokkos::RangePolicy<ExecSpaceType, int> (0, numCells),
            KOKKOS_LAMBDA (const size_t ic) {

          auto A = Kokkos::subview(elemMat, ic, Kokkos::ALL(), Kokkos::ALL());
          auto tau = Kokkos::subview(taul, ic, Kokkos::ALL());
          auto w = Kokkos::subview(work, ic, Kokkos::ALL());

          for(ordinal_type i=n; i<n+m; ++i)
            for(ordinal_type j=0; j<n; ++j)
              A(i,j) = A(j,i);

          //computing QR of A. QR is saved in A and tau
          KokkosBatched::SerialQR_Internal::invoke(A.extent(0), A.extent(1),
              A.data(), A.stride_0(), A.stride_1(), tau.data(), tau.stride_0(), w.data());

          auto b = Kokkos::subview(elemRhs, ic, Kokkos::ALL());

          //b'*Q -> b
          KokkosBatched::SerialApplyQ_RightForwardInternal::invoke(
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
          for(ordinal_type i=0; i<n; ++i){
            basisCoeffs(ic,elemDof(i)) = b(i);
          }
        });
      }
    }
#endif

    /** \brief  Serial implementation of solve, using Lapack GELS function
     */

    template<typename ViewType1, typename ViewType2, typename ViewType3, typename ViewType4>
    void solveHost(ViewType1 basisCoeffs, ViewType2 elemMat, ViewType2 elemRhs, ViewType2 ,
                   ViewType3, const  ViewType4 elemDof, ordinal_type n, ordinal_type m) {
      using value_type = typename ViewType2::value_type;
      using device_type = DeviceType;
      using host_exec_space = Kokkos::DefaultHostExecutionSpace;
      using host_memory_space = Kokkos::HostSpace;
      using host_device_type = Kokkos::Device<host_exec_space,host_memory_space>;
      using vector_host_type = Kokkos::View<value_type*,host_device_type>;
      using scratch_host_type = Kokkos::View<value_type*,host_exec_space::scratch_memory_space>;
      using matrix_host_type = Kokkos::View<value_type**,Kokkos::LayoutLeft,host_device_type>;
      using do_not_init_tag = Kokkos::ViewAllocateWithoutInitializing;
      using host_team_policy_type = Kokkos::TeamPolicy<host_exec_space>;
      using host_range_policy_type = Kokkos::RangePolicy<host_exec_space>;

      /// make sure all on-going kernels are done
      Kokkos::fence();

      /// const values
      const ordinal_type numCells = basisCoeffs.extent(0);
      const ordinal_type numRows = m+n, numCols = n;

      /// capture without this pointer
      Teuchos::LAPACK<ordinal_type,value_type> lapack;

      /// stride view copy
      Kokkos::View<ordinal_type*,host_device_type> elemDof_host(do_not_init_tag("elemDof_host"), elemDof.extent(0));
      {
        auto elemDof_device = Kokkos::create_mirror_view(typename device_type::memory_space(), elemDof_host);
        Kokkos::deep_copy(elemDof_device, elemDof); Kokkos::fence();
        Kokkos::deep_copy(elemDof_host, elemDof_device);
      }

      /// mirror view to host
      auto elemRhs_host = Kokkos::create_mirror_view_and_copy(host_memory_space(), elemRhs);
      auto elemMat_host = Kokkos::create_mirror_view_and_copy(host_memory_space(), elemMat);

      /// this in-out variable
      auto basisCoeffs_host = Kokkos::create_mirror_view_and_copy(host_memory_space(), basisCoeffs);

      if (matrixIndependentOfCell_) {
        /// invert the first matrix and apply for all
        matrix_host_type A(do_not_init_tag("A"), numRows, numRows);
        {
          for (ordinal_type j=0;j<numRows;++j)
            for (ordinal_type i=0;i<numRows;++i)
              A(i, j) = elemMat_host(0, i, j);

          for (ordinal_type j=0;j<numCols;++j)
            for (ordinal_type i=numCols;i<numRows;++i)
              A(i, j) = A(j, i);
        }
        
        ordinal_type lwork(-1);
        { /// workspace query
          ordinal_type info(0);
          value_type work[2];
          lapack.GELS('N', 
                      numRows, numRows, numCells,
                      nullptr, std::max(1,numRows),
                      nullptr, std::max(1,numRows),
                      &work[0], lwork,
                      &info);
          lwork = work[0];
        }
        
        matrix_host_type C(do_not_init_tag("C"), numRows, numCells);

        host_range_policy_type policy(0, numCells);
        {
          Kokkos::parallel_for
            ("ProjectionTools::solveHost::matrixIndependentOfCell::pack",
             policy, [=](const ordinal_type & ic) {
              for (ordinal_type i=0;i<numRows;++i)
                C(i,ic) = elemRhs_host(ic, i);
            });
        }
        {
          /// GELS does scaling and separating QR and apply Q is not stable
          vector_host_type work(do_not_init_tag("work"), lwork);
          ordinal_type info(0);
          lapack.GELS('N', 
                      numRows, numRows, numCells,
                      A.data(), std::max(1,numRows),
                      C.data(), std::max(1,numRows),
                      work.data(), lwork,
                      &info);
          INTREPID2_TEST_FOR_EXCEPTION
            (info != 0, std::runtime_error, "GELS return non-zero info code");          
        }
        {
          Kokkos::parallel_for
            ("ProjectionTools::solveHost::matrixIndependentOfCell::unpack",
             policy, [=](const ordinal_type & ic) {
              for (ordinal_type i=0;i<numCols;++i)
                basisCoeffs_host(ic,elemDof_host(i)) = C(i,ic);
            });
        }
      } else {
        const ordinal_type level(0);
        host_team_policy_type policy(numCells, 1, 1);

        /// workspace query
        ordinal_type lwork(-1);
        {
          ordinal_type info(0);
          value_type work[2];
          lapack.GELS('N', 
                      numRows, numRows, 1,
                      nullptr, std::max(1,numRows),
                      nullptr, std::max(1,numRows),
                      &work[0], lwork,
                      &info);
          lwork = work[0];
        }

        const ordinal_type per_team_extent = numRows*numRows + numRows + lwork;
        const ordinal_type per_team_scratch = scratch_host_type::shmem_size(per_team_extent);
        policy.set_scratch_size(level, Kokkos::PerTeam(per_team_scratch));

        /// solve for all
        Kokkos::parallel_for
          ("ProjectionTools::solveHost::matrixDependentOfCell",
           policy, [=](const typename host_team_policy_type::member_type& member) {
            const ordinal_type ic = member.league_rank();
            
            scratch_host_type scratch(member.team_scratch(level), per_team_extent);
            value_type * sptr = scratch.data();
            
            /// either A comes from host or device it is not column major; it needs repack
            matrix_host_type A(sptr, numRows, numRows); sptr += A.span();
            for (ordinal_type j=0;j<numRows;++j)
              for (ordinal_type i=0;i<numRows;++i)
                A(i, j) = elemMat_host(ic, i, j);

            for (ordinal_type j=0;j<numCols;++j)
              for (ordinal_type i=numCols;i<numRows;++i)
                A(i, j) = A(j, i);

            vector_host_type c(sptr, numRows); sptr += c.span();
            for (ordinal_type i=0;i<numRows;++i)
              c(i) = elemRhs_host(ic, i);

            vector_host_type work(sptr, lwork); sptr += work.span();
            ordinal_type info(0);
            lapack.GELS('N', 
                        numRows, numRows, 1,
                        A.data(), std::max(1,numRows),
                        c.data(), std::max(1,numRows),
                        work.data(), lwork,
                        &info);
            INTREPID2_TEST_FOR_EXCEPTION
              (info != 0, std::runtime_error, "GELS return non-zero info code");          

            /// scatter back to system
            for (ordinal_type i=0;i<numCols;++i) 
              basisCoeffs_host(ic,elemDof_host(i)) = c(i);
          });
      }      
      Kokkos::deep_copy(basisCoeffs, basisCoeffs_host);
    }
  };
  
};

} //Intrepid2


// include templated function definitions
#include "Intrepid2_ProjectionToolsDefL2.hpp"
#include "Intrepid2_ProjectionToolsDefHGRAD.hpp"
#include "Intrepid2_ProjectionToolsDefHCURL.hpp"
#include "Intrepid2_ProjectionToolsDefHDIV.hpp"
#include "Intrepid2_ProjectionToolsDefHVOL.hpp"

#endif





