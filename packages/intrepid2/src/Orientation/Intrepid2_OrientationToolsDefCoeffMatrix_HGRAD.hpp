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


/** \file   Intrepid2_OrientationToolsDefCoeffMatrix_HGRAD.hpp
    \brief  Creation of orientation matrix A of a face or edge for HGRAD elements
    \li     \f$\sum_k A_ik \psi_k(F_s (\eta_o (\xi_j))) = \phi_i (\xi_j) \f$ where
    \li     \f$\psi_k\f$ are the basis functions of the reference cell ,
    \li     \f$\phi_i\f$ are the basis function of the subcell
    \li     \f$F_s\f$ is the map from the reference subcell to the subcell s of the reference cell (s == subcellId)
    \li     \f$\eta_o\f$ is the orientation map o associated to the subcell s (o == subcellOrt)
    \li     \f$\xi_j\f$ are points in the reference subcell. We take them to be the dof points of the subcell basis so that \f$\phi_i (\xi_j) = \delta_ij\f$

    \author Created by Kyungjoo Kim
*/
#ifndef __INTREPID2_ORIENTATIONTOOLS_DEF_COEFF_MATRIX_HGRAD_HPP__
#define __INTREPID2_ORIENTATIONTOOLS_DEF_COEFF_MATRIX_HGRAD_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {

  namespace Impl {

    template<typename OutputViewType,
             typename subcellBasisType,
             typename cellBasisType>
    inline
    void
    OrientationTools::
    getCoeffMatrix_HGRAD(OutputViewType &output,
                         const subcellBasisType subcellBasis,
                         const cellBasisType cellBasis,
                         const ordinal_type subcellId,
                         const ordinal_type subcellOrt) {
      typedef typename OutputViewType::execution_space space_type;
      typedef typename OutputViewType::value_type value_type;

      // with shards, everything should be computed on host space
      typedef typename
        Kokkos::Impl::is_space<space_type>::host_mirror_space::execution_space host_space_type;

      typedef Kokkos::DynRankView<value_type,host_space_type> DynRankViewHostType;


      //
      // Topology
      //

      // populate points on a subcell and map to subcell
      const shards::CellTopology cellTopo = cellBasis.getBaseCellTopology();
      const shards::CellTopology subcellTopo = subcellBasis.getBaseCellTopology();

      const ordinal_type cellDim = cellTopo.getDimension();
      const ordinal_type subcellDim = subcellTopo.getDimension();

      INTREPID2_TEST_FOR_EXCEPTION( subcellDim >= cellDim,
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HGRAD): " \
                                    "cellDim must be greater than subcellDim.");

      const auto subcellBaseKey = subcellTopo.getBaseKey();

      INTREPID2_TEST_FOR_EXCEPTION( subcellBaseKey != shards::Line<>::key &&
                                    subcellBaseKey != shards::Quadrilateral<>::key &&
                                    subcellBaseKey != shards::Triangle<>::key,
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HGRAD): " \
                                    "subcellBasis must have line, quad, or triangle topology.");


      //
      // Function space
      //
      
      {
        const std::string cellBasisName(cellBasis.getName());
        if (cellBasisName.find("HGRAD") != std::string::npos) {
          const std::string subcellBasisName(subcellBasis.getName());
          INTREPID2_TEST_FOR_EXCEPTION( subcellBasisName.find("HGRAD") == std::string::npos,
                                        std::logic_error,
                                        ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HGRAD): " \
                                        "subcellBasis function space is not consistent to cellBasis.");
        }

        INTREPID2_TEST_FOR_EXCEPTION( subcellBasis.getDegree() != cellBasis.getDegree(),
                                      std::logic_error,
                                      ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HGRAD): " \
                                      "subcellBasis has a different polynomial degree from cellBasis' degree.");
      }


      //
      // Collocation points
      //

      const ordinal_type numCellBasis = cellBasis.getCardinality();
      const ordinal_type numSubcellBasis = subcellBasis.getCardinality();

      const ordinal_type ordSubcell = cellBasis.getDofOrdinal(subcellDim, subcellId, 0);
      INTREPID2_TEST_FOR_EXCEPTION( ordSubcell == -1,
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HGRAD): " \
                                    "Invalid subcellId returns -1 ordSubcell.");

      const ordinal_type ndofSubcell = cellBasis.getDofTag(ordSubcell)(3);

      // reference points on a subcell
      DynRankViewHostType refPtsSubcell("refPtsSubcell", ndofSubcell, subcellDim);
      DynRankViewHostType subcellDofCoords("subcellDofCoords", numSubcellBasis, subcellDim);
      subcellBasis.getDofCoords(subcellDofCoords);
      for(ordinal_type i=0; i<ndofSubcell; ++i)
        for(ordinal_type d=0; d <subcellDim; ++d)
          refPtsSubcell(i,d) = subcellDofCoords(subcellBasis.getDofOrdinal(subcellDim, 0, i),d);

      // modified points with orientation
      DynRankViewHostType ortPtsSubcell("ortPtsSubcell", ndofSubcell, subcellDim);
      Impl::OrientationTools::mapToModifiedReference(ortPtsSubcell,
                                                     refPtsSubcell,
                                                     subcellTopo,
                                                     subcellOrt);

      // map to reference coordinates
      DynRankViewHostType refPtsCell("refPtsCell", ndofSubcell, cellDim);
      CellTools<host_space_type>::mapToReferenceSubcell(refPtsCell,
                                                        ortPtsSubcell,
                                                        subcellDim,
                                                        subcellId,
                                                        cellTopo);

      //
      // Basis evaluation on the collocation points
      //

      DynRankViewHostType refValues("refValues", numCellBasis, ndofSubcell);
      cellBasis.getValues(refValues, refPtsCell, OPERATOR_VALUE);
      

      //
      // Construct collocation matrix and solve problems
      //

      // construct collocation matrix; using lapack, it should be left layout
      Kokkos::View<value_type**,Kokkos::LayoutLeft,host_space_type>
        refMat("refMat", ndofSubcell, ndofSubcell),
        ortMat("ortMat", ndofSubcell, ndofSubcell),
        pivVec("pivVec", ndofSubcell, 1);

      for (ordinal_type i=0;i<ndofSubcell;++i) {
        const ordinal_type iref = cellBasis.getDofOrdinal(subcellDim, subcellId, i);
        for (ordinal_type j=0;j<ndofSubcell;++j) {
          refMat(i,j) = refValues(iref,j);
          ortMat(i,j) = (i==j);      //thanks to the kronecher property of the basis functions
        }
      }

      // solve the system
      {
        Teuchos::LAPACK<ordinal_type,value_type> lapack;
        ordinal_type info = 0;

        lapack.GESV(ndofSubcell, ndofSubcell,
                    refMat.data(),
                    refMat.stride_1(),
                    (ordinal_type*)pivVec.data(),
                    ortMat.data(),
                    ortMat.stride_1(),
                    &info);

        if (info) {
          std::stringstream ss;
          ss << ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HGRAD): "
             << "LAPACK return with error code: "
             << info;
          INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str().c_str() );
        }

        //Transpose matrix and clean up numerical noise
        {
          const double eps = threshold();
          for (ordinal_type i=0;i<ndofSubcell;++i)
            for (ordinal_type j=i;j<ndofSubcell;++j) {
              auto intOrtMat = std::round(ortMat(i,j));
              ortMat(i,j) = (std::abs(ortMat(i,j) - std::round(ortMat(i,j))) < eps) ? intOrtMat : ortMat(i,j);
            }
        }
      }

      {
        // move the data to original device memory
        const Kokkos::pair<ordinal_type,ordinal_type> range(0, ndofSubcell);
        Kokkos::deep_copy(Kokkos::subview(output, range, range), ortMat);
      }
    }
  }

}
#endif
