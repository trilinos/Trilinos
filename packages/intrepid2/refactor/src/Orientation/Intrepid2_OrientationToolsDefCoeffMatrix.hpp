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


/** \file   Intrepid_OrientationToolsDef.hpp
    \brief  Definition file for the Intrepid2::OrientationTools class.
    \author Created by Kyungjoo Kim
*/
#ifndef __INTREPID2_ORIENTATIONTOOLS_DEF_COEFF_MATRIX_HPP__
#define __INTREPID2_ORIENTATIONTOOLS_DEF_COEFF_MATRIX_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

//#include "Teuchos_LAPACK.hpp"
namespace Intrepid2 {

  namespace Impl {

    template<typename VT>
    inline 
    Kokkos::View<ValueType**,Kokkos::HostSpace>
    OrientationTools::
    getEdgeCoeffMatrix_HGRAD(const Basis<Kokkos::HostSpace,VT,VT> lineBasis,
                             const Basis<Kokkos::HostSpace,VT,VT> cellBasis,
                             const ordinal_type edgeId,
                             const ordinal_type edgeOrt) {
      // populate points on a line and map to subcell
      const shards::CellTopology cellTopo = cellBasis.getBaseCellTopology();
      const shards::CellTopology lineTopo = lineBasis.getBaseCellTopology();

      const auto cellDim = cellTopo.getDimension();
      const auto lineDim = lineTopo.getDimension();

      const auto degree = cellBasis.getDegree();

      const auto numCellBasis = cellBasis.getCardinality();
      const auto numLineBasis = lineBasis.getCardinality();

      const auto ordEdge = cellBasis.getDofOrdinal(lineDim, edgeId, 0);
      const auto ndofEdge = cellBasis.getDofTag(ordEdge)(3);

#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_EXCEPTION( !(ndofEdge == PointTools::getLatticeSize(lineTopo, degree, 1)),
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getEdgeCoeffMatrix_HGRAD): " \
                                    "The number of DOFs does not match to the number of collocation points.");
#endif

      // reference points between (-1 , 1)
      Kokkos::DynRankView<VT,Kokkos::HostSpace> refPtsLine("refPtsLine", ndofEdge, lineDim);
      PointTools::getLattice(refPtsLine,
                             lineTopo,
                             degree,
                             1,
                             POINTTYPE_EQUISPACED);

      // modified points with orientation
      Kokkos::DynRankView<VT,Kokkos::HostSpace> ortPtsLine("ortPtsLine", ndofEdge, lineDim);
      Impl::OrientationTools::mapToModifiedReference(ortPtsLine,
                                                     refPtsLine,
                                                     lineTopo,
                                                     edgeOrt);
      
      // map to reference coordinates
      Kokkos::DynRankView<VT,Kokkos::HostSpace> refPtsCell("refPtsCell", ndofEdge, cellDim);
      CellTools<Kokkos::HostSpace>::mapToReferenceSubcell(refPtsCell,
                                                          refPtsLine,
                                                          lineDim,
                                                          edgeId,
                                                          cellTopo);

      // evaluate values on the reference cell
      Kokkos::DynRankView<VT,Kokkos::HostSpace> refValues("refValues", numCellBasis, ndofEdge);
      cellBasis.getValues(refValues, refPtsCell, OPERATOR_VALUE);

      // evaluate values on the modified cell
      Kokkos::DynRankView<VT,Kokkos::HostSpace> outValues("outValues", numLineBasis, ndofEdge);
      lineBasis.getValues(outValues, ortPtsLine, OPERATOR_VALUE);

      // construct collocation matrix
      Kokkos::View<VT**,Kokkos::HostSpace,Kokkos::LayoutLeft> 
        refMat("refMat", ndofEdge, ndofEdge),
        ortMat("ortMat", ndofEdge, ndofEdge),
        pivVec("pivVec", ndofEdge, 1);
      
      for (auto i=0;i<ndofEdge;++i) {
        const auto iref = cellBasis.getDofOrdinal(lineDim, edgeId, i);
        const auto iout = lineBasis.getDofOrdinal(lineDim, 0,      i);
        for (auto j=0;j<ndofEdge;++j) {
          refMat(j,i) = refValues(iref,j);
          ortMat(j,i) = outValues(iout,j);
        }
      }

      // solve the system
      {
        Teuchos::LAPACK<ordinal_type,VT> lapack;
        ordinal_type info = 0;
        lapack.GESV(ndofEdge, ndofEdge,
                    refMat.data(),
                    refMat.stride(0),
                    (ordinal_type*)pivVec.data(),
                    ortMat.data(),
                    ortMat.stride(0),
                    &info);
        
        if (info) {
          std::stringstream ss;
          ss << ">>> ERROR (Intrepid::OrientationTools::getEdgeCoeffMatrix_HGRAD): "
             << "LAPACK return with error code: "
             << info;
          INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str() );
        }
      }

      return ortMat;
    }


    template<typename VT>
    inline 
    Kokkos::View<ValueType**,Kokkos::HostSpace>
    OrientationTools::
    getTriangleCoeffMatrix_HGRAD(const Basis<Kokkos::HostSpace,VT,VT> faceBasis,
                                 const Basis<Kokkos::HostSpace,VT,VT> cellBasis,
                                 const ordinal_type faceId,
                                 const ordinal_type faceOrt) {
      // populate points on a line and map to subcell
      const shards::CellTopology cellTopo = cellBasis.getBaseCellTopology();
      const shards::CellTopology faceTopo = faceBasis.getBaseCellTopology();

      // if the face is left-handed system, the orientation should be re-enumerated
      const auto leftHanded = cellTopo.getNodeMap(2, faceId, 1) > cellTopo.getNodeMap(2, faceId, 2);
      const auto leftOrt[] = { 0, 2, 1, 3, 5, 4 };
      const auto ort = (leftHanded ? leftOrt[faceOrt] : faceOrt);

      const auto cellDim = cellTopo.getDimension();
      const auto faceDim = faceTopo.getDimension();

      const auto degree = cellBasis.getDegree();

      const auto numCellBasis = cellBasis.getCardinality();
      const auto numFaceBasis = faceBasis.getCardinality();

      const auto ordFace = cellBasis.getDofOrdinal(faceDim, faceId, 0);
      const auto ndofFace = cellBasis.getDofTag(ordFace)(3);

#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_ABORT( !(ndofFace == PointTools::getLatticeSize<Scalar>(faceTopo, degree, 1)),
                                  std::logic_error,
                                  ">>> ERROR (Intrepid::OrientationTools::getTriangleCoeffMatrix_HGRAD): " \
                                  "The number of DOFs does not match to the number of collocation points.");
#endif

      // reference points in triangle
      Kokkos::DynRankView<VT,Kokkos::HostSpace> refPtsFace("refPtsFace", ndofFace, faceDim);
      PointTools::getLattice(refPtsFace,
                             faceTopo,
                             degree,
                             1,
                             POINTTYPE_EQUISPACED);
      
      // modified points with orientation
      Kokkos::DynRankView<VT,Kokkos::HostSpace> ortPtsFace("ortPtsFace", ndofFace, faceDim);
      Impl::OrientationTools::mapToModifiedReference(ortPtsFace,
                                                     refPtsFace,
                                                     faceTopo,
                                                     ort);
      
      // map to reference coordinates
      Kokkos::DynRankView<VT,Kokkos::HostSpace> refPtsCell("refPtsCell", ndofFace, cellDim);
      CellTools<Kokkos::HostSpace>::mapToReferenceSubcell(refPtsCell,
                                                          refPtsFace,
                                                          faceDim,
                                                          faceId,
                                                          cellTopo);

      // evaluate values on the reference cell
      Kokkos::DynRankView<VT,Kokkos::HostSpace> refValues("refValues", numCellBasis, ndofFace);
      cellBasis.getValues(refValues, refPtsCell, OPERATOR_VALUE);

      // evaluate values on the modified cell
      Kokkos::DynRankView<VT,Kokkos::HostSpace> outValues("outValues", numFaceBasis, ndofFace);
      faceBasis.getValues(outValues, ortPtsFace, OPERATOR_VALUE);

      // construct collocation matrix
      Kokkos::View<VT**,Kokkos::HostSpace,Kokkos::LayoutLeft> 
        refMat("refMat", ndofFace, ndofFace),
        ortMat("ortMat", ndofFace, ndofFace),
        pivVec("pivVec", ndofFace, 1);

      for (auto i=0;i<ndofFace;++i) {
        const auto iref = cellBasis.getDofOrdinal(faceDim, faceId, i);
        const auto iout = faceBasis.getDofOrdinal(faceDim, 0,      i);
        for (auto j=0;j<ndofFace;++j) {
          refMat(j,i) = refValues(iref,j);
          ortMat(j,i) = outValues(iout,j);
        }
      }
      
      // solve the system
      {
        Teuchos::LAPACK<ordinal_type,VT> lapack;
        ordinal_type info = 0;
        lapack.GESV(ndofFace, ndofFace,
                    refMat.data(),
                    refMat.stride(0),
                    (ordinal_type*)pivVec.data(),
                    ortMat.data(),
                    ortMat.stride(0),
                    &info);
        
        if (info) {
          std::stringstream ss;
          ss << ">>> ERROR (Intrepid::OrientationTools::getTriangleCoeffMatrix_HGRAD): "
             << "LAPACK return with error code: "
             << info;
          INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str() );
        }
      }
      
      return ortMat;
    }


    template<typename VT>
    inline 
    Kokkos::View<ValueType**,Kokkos::HostSpace>
    OrientationTools::
    getQuadrilateralCoeffMatrix_HGRAD(const Basis<Kokkos::HostSpace,VT,VT> faceBasis,
                                      const Basis<Kokkos::HostSpace,VT,VT> cellBasis,
                                      const ordinal_type faceId,
                                      const ordinal_type faceOrt) {
      // populate points on a line and map to subcell
      const shards::CellTopology cellTopo = cellBasis.getBaseCellTopology();
      const shards::CellTopology faceTopo = faceBasis.getBaseCellTopology();

      // if the face is left-handed system, the orientation should be re-enumerated
      const auto leftHanded = cellTopo.getNodeMap(2, faceId, 1) > cellTopo.getNodeMap(2, faceId, 3);
      const auto leftOrt[] = { 0, 3, 2, 1, 4, 7, 6, 5 };
      const auto ort = (leftHanded ? leftOrt[faceOrt] : faceOrt);

      const auto cellDim = cellTopo.getDimension();
      const auto faceDim = faceTopo.getDimension();

      const auto degree = cellBasis.getDegree();

      const auto numCellBasis = cellBasis.getCardinality();
      const auto numFaceBasis = faceBasis.getCardinality();

      const auto ordFace = cellBasis.getDofOrdinal(faceDim, faceId, 0);
      const auto ndofFace = cellBasis.getDofTag(ordFace)(3);

#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_ABORT( !(ndofFace == PointTools::getLatticeSize(faceTopo, degree, 1)),
                                  std::logic_error,
                                  ">>> ERROR (Intrepid::OrientationTools::getQuadrilateralCoeffMatrix_HGRAD): " \
                                  "The number of DOFs does not match to the number of collocation points.");
#endif

      // reference points in quadrilateral
      Kokkos::DynRankView<VT,Kokkos::HostSpace> refPtsFace("refPtsFace", ndofFace, faceDim);
      PointTools::getLattice(refPtsFace,
                             faceTopo,
                             degree,
                             1,
                             POINTTYPE_EQUISPACED);
      
      // modified points with orientation
      Kokkos::DynRankView<VT,Kokkos::HostSpace> ortPtsFace("ortPtsFace", ndofFace, faceDim);
      Impl::OrientationTools::mapToModifiedReference(ortPtsFace,
                                                     refPtsFace,
                                                     faceTopo,
                                                     ort);
      
      // map to reference coordinates
      Kokkos::DynRankView<VT,Kokkos::HostSpace> refPtsCell("refPtsCell", ndofFace, cellDim);
      CellTools<Kokkos::HostSpace>::mapToReferenceSubcell(refPtsCell,
                                                          refPtsFace,
                                                          faceDim,
                                                          faceId,
                                                          cellTopo);

      // evaluate values on the reference cell
      Kokkos::DynRankView<VT,Kokkos::HostSpace> refValues("refValues", numCellBasis, ndofFace);
      cellBasis.getValues(refValues, refPtsCell, OPERATOR_VALUE);

      // evaluate values on the modified cell
      Kokkos::DynRankView<VT,Kokkos::HostSpace> outValues("outValues", numFaceBasis, ndofFace);
      faceBasis.getValues(outValues, ortPtsFace, OPERATOR_VALUE);

      // construct collocation matrix
      Kokkos::View<VT**,Kokkos::HostSpace,Kokkos::LayoutLeft> 
        refMat("refMat", ndofFace, ndofFace),
        ortMat("ortMat", ndofFace, ndofFace),
        pivVec("pivVec", ndofFace, 1);

      for (auto i=0;i<ndofFace;++i) {
        const auto iref = cellBasis.getDofOrdinal(faceDim, faceId, i);
        const auto iout = faceBasis.getDofOrdinal(faceDim, 0,      i);
        for (auto j=0;j<ndofFace;++j) {
          refMat(j,i) = refValues(iref,j);
          ortMat(j,i) = outValues(iout,j);
        }
      }
      
      // solve the system
      {
        Teuchos::LAPACK<ordinal_type,VT> lapack;
        ordinal_type info = 0;
        lapack.GESV(ndofFace, ndofFace,
                    refMat.data(),
                    refMat.stride(0),
                    (ordinal_type*)pivVec.data(),
                    ortMat.data(),
                    ortMat.stride(0),
                    &info);
        
        if (info) {
          std::stringstream ss;
          ss << ">>> ERROR (Intrepid::OrientationTools::getTriangleCoeffMatrix_HGRAD): "
             << "LAPACK return with error code: "
             << info;
          INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str() );
        }
      }
      
      return ortMat;
    }
  }
}
#endif
