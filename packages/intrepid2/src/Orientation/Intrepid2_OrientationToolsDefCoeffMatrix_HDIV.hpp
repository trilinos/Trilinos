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


/** \file  Intrepid2_OrientationToolsDefCoeffMatrix_HDIV.hpp
    \brief  Creation of orientation matrix A of a face or edge for HDIV elements

    \li     \f$\sum_k A_ik \psi_k(F_s (\eta_o (\xi_j))) \cdot (n_s det(J_\eta)) = \phi_i (\xi_j)\f$, where
    \li     \f$\psi_k\f$ are the basis functions of the reference cell,
    \li     \f$\phi_i\f$ are the basis function of the subcell,
    \li     \f$n_j\f$ is the side normal of subcell s (s == subcellId)
    \li     \f$F_s\f$ is the Jacobian of the map from the reference subcell to the subcell s of the reference cell (s == subcellId),
    \li     \f$\eta_o\f$ is the orientation map o associated to the subcell s (o == subcellOrt),
    \li     \f$J_{\eta}\f$ is the Jacobian of the map \f$\eta_o\f$,
    \li     \f$\xi_j\f$ are dof points of the subcell basis so that \f$\phi_i (\xi_j) \dot t_j = \delta_ij\f$

    \author Created by Kyungjoo Kim
*/
#ifndef __INTREPID2_ORIENTATIONTOOLS_DEF_COEFF_MATRIX_HDIV_HPP__
#define __INTREPID2_ORIENTATIONTOOLS_DEF_COEFF_MATRIX_HDIV_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {

  template<typename ExecSpaceType,
           typename outputValueType,
           typename pointValueType>
  class Basis_HDIV_TET_In_FEM;
  
  namespace Impl {

    template<typename OutputViewType,
             typename subcellBasisType,
             typename cellBasisType>
    inline
    void
    OrientationTools::
    getCoeffMatrix_HDIV(OutputViewType &output,
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
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): " \
                                    "cellDim must be greater than subcellDim.");

      const auto subcellBaseKey = subcellTopo.getBaseKey();
      const auto cellBaseKey = cellTopo.getBaseKey();

      INTREPID2_TEST_FOR_EXCEPTION( subcellBaseKey != shards::Line<>::key &&
                                    subcellBaseKey != shards::Quadrilateral<>::key &&
                                    subcellBaseKey != shards::Triangle<>::key,
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): " \
                                    "subcellBasis must have line, quad, or triangle topology.");

      INTREPID2_TEST_FOR_EXCEPTION( cellBaseKey != shards::Quadrilateral<>::key &&
                                    cellBaseKey != shards::Triangle<>::key &&
                                    cellBaseKey != shards::Hexahedron<>::key &&
                                    cellBaseKey != shards::Tetrahedron<>::key,
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): " \
                                    "cellBasis must have quad, triangle, hexhedron or tetrahedron topology.");

      //
      // Function space
      //
      {
        const std::string cellBasisName(cellBasis.getName());
        INTREPID2_TEST_FOR_EXCEPTION( cellBasisName.find("HDIV") == std::string::npos,
                                      std::logic_error,
                                      ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): "
                                      "cellBasis is not HDIV.");
        {
          const std::string subcellBasisName(subcellBasis.getName());
          switch (subcellDim) {
          case 1: {
            //TODO: Hex, QUAD, TET and TRI element should have the same 1d basis
            if ((cellBasisName.find("HEX") != std::string::npos) || (cellBasisName.find("QUAD") != std::string::npos)) {
              INTREPID2_TEST_FOR_EXCEPTION( subcellBasisName.find("HGRAD") == std::string::npos,
                                          std::logic_error,
                                          ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_DIV): " 
                                          "subcellBasis function space (1d) is not consistent to cellBasis, which should be open line hgrad, order -1.");
            } else if ((cellBasisName.find("TET") != std::string::npos) || (cellBasisName.find("TRI") != std::string::npos)) {
              INTREPID2_TEST_FOR_EXCEPTION( subcellBasisName.find("HVOL") == std::string::npos,
                                          std::logic_error,
                                          ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_DIV): " 
                                          "subcellBasis function space (1d) is not consistent to cellBasis, which should be HVOL line, order -1.");
            }
            break;
          }
          case 2: {
            if        (subcellBaseKey == shards::Quadrilateral<>::key) {
              // quad face basis is tensor product of open line basis functions
              INTREPID2_TEST_FOR_EXCEPTION( subcellBasisName.find("HGRAD") == std::string::npos,
                                            std::logic_error,
                                            ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): "
                                            "subcellBasis function space is not compatible, which should be open line hgrad, order -1.");
            } else if (subcellBaseKey == shards::Triangle<>::key) {
              // triangle face basis comes from HVOL basis
              INTREPID2_TEST_FOR_EXCEPTION( subcellBasisName.find("HVOL") == std::string::npos,
                                            std::logic_error,
                                            ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): "
                                            "subcellBasis function space is not compatible, which should HVOL, order-1.");
            }
            break;
          }
          }
        }
      }

      //
      // Collocation points
      //
      const ordinal_type numCellBasis = cellBasis.getCardinality();
      const ordinal_type numSubcellBasis = subcellBasis.getCardinality();

      const ordinal_type ordSubcell = cellBasis.getDofOrdinal(subcellDim, subcellId, 0);
      INTREPID2_TEST_FOR_EXCEPTION( ordSubcell == -1,
                                    std::logic_error,
                                    ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): " \
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
      // evaluate values on the reference cell
      DynRankViewHostType refValues("refValues", numCellBasis, ndofSubcell, cellDim);
      cellBasis.getValues(refValues, refPtsCell, OPERATOR_VALUE);


      DynRankViewHostType ortJacobian("ortJacobian", subcellDim, subcellDim);
      Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobian, subcellTopo, subcellOrt);
      auto ortJacobianDet = (subcellDim == 1) ? ortJacobian(0,0) :
                                 ortJacobian(0,0)*ortJacobian(1,1)-ortJacobian(1,0)*ortJacobian(0,1);

      //
      // Restrict vector valued basis functions on the subcell dimensions
      //
      DynRankViewHostType sideNormal("sideNormal", cellDim);
      CellTools<host_space_type>::getReferenceSideNormal(sideNormal, subcellId, cellTopo);

      if((cellBaseKey == shards::Triangle<>::key))
        for(ordinal_type i=0; i<cellDim; ++i)
          sideNormal(i) *= 2;  //scale by reference edge length

      Kokkos::View<value_type**,Kokkos::LayoutLeft,host_space_type>
        refMat("refMat", ndofSubcell, ndofSubcell),
        ortMat("ortMat", ndofSubcell, ndofSubcell);
      DynRankViewHostType tmpValues("tmpValues", numCellBasis, ndofSubcell);
      for (ordinal_type i=0;i<ndofSubcell;++i) {
        const ordinal_type ii = cellBasis.getDofOrdinal(subcellDim, subcellId, i);
        for (ordinal_type j=0;j<ndofSubcell;++j)
          for (ordinal_type k=0;k<cellDim;++k) {
            refMat(i,j) += refValues(ii,j,k)*sideNormal(k)* ortJacobianDet;
            ortMat(i,j) = (i==j);  //ortMat is the identity thanks to the kronecker property of the bases
          }
      }

       // solve the system
      {
        Teuchos::LAPACK<ordinal_type,value_type> lapack;
        ordinal_type info = 0;
        Kokkos::View<value_type**,Kokkos::LayoutLeft,host_space_type> pivVec("pivVec", ndofSubcell, 1);

        lapack.GESV(ndofSubcell, ndofSubcell,
                    refMat.data(),
                    refMat.stride_1(),
                    (ordinal_type*)pivVec.data(),
                    ortMat.data(),
                    ortMat.stride_1(),
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
          for (ordinal_type i=0;i<ndofSubcell;++i)
            for (ordinal_type j=0;j<ndofSubcell;++j) {
              auto intOrtMat = std::round(ortMat(i,j));
              ortMat(i,j) = (std::abs(ortMat(i,j) - std::round(ortMat(i,j))) < eps) ? intOrtMat : ortMat(i,j);
            }
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
