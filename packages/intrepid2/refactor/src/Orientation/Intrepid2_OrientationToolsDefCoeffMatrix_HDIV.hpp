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
#ifndef __INTREPID2_ORIENTATIONTOOLS_DEF_COEFF_MATRIX_HDIV_HPP__
#define __INTREPID2_ORIENTATIONTOOLS_DEF_COEFF_MATRIX_HDIV_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {

  namespace Impl {

    template<typename outputViewType,
             typename subcellBasisType,
             typename cellBasisType>
    inline
    void
    OrientationTools::
    getCoeffMatrix_HDIV(outputViewType &output,
                        const subcellBasisType subcellBasis,
                        const cellBasisType cellBasis,
                        const ordinal_type subcellId,
                        const ordinal_type subcellOrt) {
      typedef typename outputViewType::execution_space space_type;
      typedef typename outputViewType::value_type value_type;
      
      // with shards, everything should be computed on host space
      typedef typename
        Kokkos::Impl::is_space<space_type>::host_mirror_space::execution_space host_space_type;

      typedef Kokkos::DynRankView<value_type,host_space_type> DynRankViewHostType;

      ///
      /// Function space
      ///
      
      {
        const std::string cellBasisName(cellBasis.getName());
        if (cellBasisName.find("HDIV") != std::string::npos) {
          const std::string subcellBasisName(subcellBasis.getName());
          // edge hcurl is hgrad with gauss legendre points
          // INTREPID2_TEST_FOR_EXCEPTION( subcellBasisName.find("HDIV") == std::string::npos,
          //                               std::logic_error,
          //                               ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): " 
          //                               "subcellBasis function space is not consistent to cellBasis.");
        }
        // degree does not have to match each other unlike hgrad functions
      }

      ///
      /// Topology
      ///

      // populate points on a subcell and map to subcell
      const shards::CellTopology cellTopo = cellBasis.getBaseCellTopology();
      const shards::CellTopology subcellTopo = subcellBasis.getBaseCellTopology();

#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_EXCEPTION( subcellTopo.getBaseKey() != shards::Line<>::key,
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): " \
                                    "subcellBasis must have line topology for now.");
      INTREPID2_TEST_FOR_EXCEPTION( cellTopo.getBaseKey() != shards::Quadrilateral<>::key,
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): " \
                                    "cellBasis must have quad topology for now.");
#endif

      // if node map has left handed system, orientation should be re-enumerated.
      ordinal_type ort = -1;
      switch (subcellTopo.getBaseKey()) {
      case shards::Line<>::key: {
        if (subcellOrt >= 0 && subcellOrt <  2)
          ort = subcellOrt;
        break;
      }
      // case shards::Triangle<>::key: {
      //   if (subcellOrt >= 0 && subcellOrt <  6) {
      //     const ordinal_type leftHanded = cellTopo.getNodeMap(2, subcellId, 1) > cellTopo.getNodeMap(2, subcellId, 2);
      //     const ordinal_type leftOrt[] = { 0, 2, 1, 3, 5, 4 };
      //     ort = (leftHanded ? leftOrt[subcellOrt] : subcellOrt);
      //   }
      //   break;
      // }
      // case shards::Quadrilateral<>::key: {
      //   if (subcellOrt >= 0 && subcellOrt <  8) {
      //     const ordinal_type leftHanded = cellTopo.getNodeMap(2, subcellId, 1) > cellTopo.getNodeMap(2, subcellId, 3);
      //     const ordinal_type leftOrt[] = { 0, 3, 2, 1, 4, 7, 6, 5 };
      //     ort = (leftHanded ? leftOrt[subcellOrt] : subcellOrt);
      //   }
      //   break;
      // }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( subcellTopo.getBaseKey() != shards::Line<>::key ||
                                      subcellTopo.getBaseKey() != shards::Quadrilateral<>::key ||
                                      subcellTopo.getBaseKey() != shards::Triangle<>::key,
                                      std::logic_error,
                                      ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): " \
                                      "subcellBasis must have line, quad, or triangle topology.");
      }
      }
      INTREPID2_TEST_FOR_EXCEPTION( ort == -1,
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): " \
                                    "Orientation is not properly setup.");
      
      ///
      /// Collocation points
      ///
      
      const ordinal_type cellDim = cellTopo.getDimension();
      const ordinal_type subcellDim = subcellTopo.getDimension();
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_EXCEPTION( subcellDim >= cellDim,
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): " \
                                    "cellDim must be greater than subcellDim.");
#endif

      const ordinal_type degree = cellBasis.getDegree();

      const ordinal_type numCellBasis = cellBasis.getCardinality();
      const ordinal_type numSubcellBasis = subcellBasis.getCardinality();

      const ordinal_type ordSubcell = cellBasis.getDofOrdinal(subcellDim, subcellId, 0);
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_EXCEPTION( ordSubcell == -1,
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): " \
                                    "Invalid subcellId returns -1 ordSubcell.");
#endif

      const ordinal_type ndofSubcell = cellBasis.getDofTag(ordSubcell)(3);
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_EXCEPTION( ndofSubcell < PointTools::getLatticeSize(subcellTopo, degree+1, 1),
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): " \
                                    "The number of DOFs should be less/equal to the number of collocation points.");
#endif

      // reference points on a subcell
      DynRankViewHostType refPtsSubcell("refPtsSubcell", ndofSubcell, subcellDim);
      PointTools::getLattice(refPtsSubcell,
                             subcellTopo,
                             degree+1,
                             1, // offset by 1 so the points are located inside
                             POINTTYPE_EQUISPACED);

      // modified points with orientation
      DynRankViewHostType ortPtsSubcell("ortPtsSubcell", ndofSubcell, subcellDim);
      Impl::OrientationTools::mapToModifiedReference(ortPtsSubcell,
                                                     refPtsSubcell,
                                                     subcellTopo,
                                                     subcellOrt);
      
      // map to reference coordinates
      DynRankViewHostType refPtsCell("refPtsCell", ndofSubcell, cellDim);
      CellTools<host_space_type>::mapToReferenceSubcell(refPtsCell,
                                                        refPtsSubcell,
                                                        subcellDim,
                                                        subcellId,
                                                        cellTopo);

      ///
      /// Basis evaluation on the collocation points
      ///

      // evaluate values on the reference cell
      DynRankViewHostType refValues("refValues", numCellBasis, ndofSubcell, cellDim);
      cellBasis.getValues(refValues, refPtsCell, OPERATOR_VALUE);

      // evaluate values on the modified cell
      DynRankViewHostType outValues("outValues", numSubcellBasis, ndofSubcell, subcellDim);
      if (subcellDim == 1) {
        subcellBasis.getValues(Kokkos::subview(outValues, Kokkos::ALL(), Kokkos::ALL(), 0),
                               ortPtsSubcell, OPERATOR_VALUE);
        // account orientation of 1D basis 
        if (subcellOrt == 1) {
          const ordinal_type
            iend = outValues.dimension(0),
            jend = outValues.dimension(1);
          for (ordinal_type i=0;i<iend;++i)
            for (ordinal_type j=0;j<jend;++j)
              outValues(i,j,0) *= (-1.0);
        }
      } else {
        subcellBasis.getValues(outValues, ortPtsSubcell, OPERATOR_VALUE);
      }

      ///
      /// Restrict vector valued basis functions on the subcell dimensions
      ///
      auto normalize = [&](DynRankViewHostType v) { 
        value_type norm = 0.0;
        const ordinal_type iend = v.dimension(0);
        for (ordinal_type i=0;i<iend;++i)
          norm += v(i)*v(i);
        norm = std::sqrt(norm);
        for (ordinal_type i=0;i<iend;++i)
          v(i) /= norm;
      };
      
      switch (subcellTopo.getBaseKey()) {
      case shards::Line<>::key: {
        DynRankViewHostType edgeTan("edgeTan", cellDim);
        CellTools<host_space_type>::getReferenceEdgeTangent(edgeTan, subcellId, cellTopo);

        normalize(edgeTan);

        // 90 degree rotation (for 2D)
        if (subcellDim == 2) {
          const auto tan_x = edgeTan(0);
          const auto tan_y = edgeTan(1);
          
          edgeTan(0) = -tan_y;
          edgeTan(1) =  tan_x;
        }
        
        DynRankViewHostType tmpValues("tmpValues", numCellBasis, ndofSubcell, subcellDim);
        for (ordinal_type i=0;i<numCellBasis;++i)
          for (ordinal_type j=0;j<ndofSubcell;++j)
            for (ordinal_type k=0;k<cellDim;++k)
              tmpValues(i,j,0) += refValues(i,j,k)*edgeTan(k);
        refValues = tmpValues;
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, "Should not come here" );        
      }
      }
      
      ///
      /// Construct collocation matrix and solve problems
      ///

      // construct collocation matrix; using lapack, it should be left layout
      const ordinal_type dofDim = outValues.dimension(2);
      Kokkos::View<value_type**,Kokkos::LayoutLeft,host_space_type>
        refMat("refMat", ndofSubcell*dofDim, ndofSubcell),
        ortMat("ortMat", ndofSubcell*dofDim, ndofSubcell),
        lwork("lwork", std::max(ndofSubcell*dofDim, 10), ndofSubcell);

      for (ordinal_type i=0;i<ndofSubcell;++i) {
        const ordinal_type iref = cellBasis.getDofOrdinal(subcellDim, subcellId, i);
        const ordinal_type iout = subcellBasis.getDofOrdinal(subcellDim, 0, i);

        for (ordinal_type j=0;j<ndofSubcell;++j) {
          const ordinal_type joff = j*dofDim;
          for (ordinal_type k=0;k<dofDim;++k) {
            refMat(joff+k,i) = refValues(iref,j,k);
            ortMat(joff+k,i) = outValues(iout,j,k);
          }
        }
      }

      // solve the system
      {
        Teuchos::LAPACK<ordinal_type,value_type> lapack;
        ordinal_type info = 0;
        lapack.GELS('N', // trans
                    ndofSubcell*dofDim, // m
                    ndofSubcell, // n
                    ndofSubcell, // nrhs
                    refMat.data(), refMat.stride_1(), // A
                    ortMat.data(), ortMat.stride_1(), // B
                    lwork.data(), lwork.span(), // work space
                    &info);

        if (info) {
          std::stringstream ss;
          ss << ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): "
             << "LAPACK return with error code: "
             << info;
          INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str().c_str() );
        }

        // clean up numerical noise
        {
          const double eps = threshold();
          const ordinal_type
            iend = ortMat.dimension(0),
            jend = ortMat.dimension(1);
          for (ordinal_type i=0;i<iend;++i)
            for (ordinal_type j=0;j<jend;++j)
              if (std::abs(ortMat(i,j)) < eps) ortMat(i,j) = 0;
        }
      }

      {
        // move the data to original device memory
        const Kokkos::pair<ordinal_type,ordinal_type> range(0, ndofSubcell);
        Kokkos::deep_copy(Kokkos::subview(output, range, range),
                          Kokkos::subview(ortMat, range, range));
      }
    }
  }

}
#endif
