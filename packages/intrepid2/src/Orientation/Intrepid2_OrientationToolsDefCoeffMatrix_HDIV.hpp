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
    \li     \f$\xi_j\f$ are points in the reference subcell such that subcell bases are uniquely
    determined by the values they take on these points.

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

namespace Debug {

#ifdef HAVE_INTREPID2_DEBUG
template<typename subcellBasisType,
typename cellBasisType>
inline
void
check_getCoeffMatrix_HDIV(const subcellBasisType& subcellBasis,
    const cellBasisType& cellBasis,
    const ordinal_type subcellId,
    const ordinal_type subcellOrt) {
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
      cellBaseKey != shards::Wedge<>::key &&
      cellBaseKey != shards::Tetrahedron<>::key,
      std::logic_error,
      ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): " \
      "cellBasis must have quad, triangle, hexhedron or tetrahedron topology.");

  //
  // Function space
  //
  {
    const bool isHDIV = cellBasis.getFunctionSpace() == FUNCTION_SPACE_HDIV;
    INTREPID2_TEST_FOR_EXCEPTION( !isHDIV,
        std::logic_error,
        ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): "
        "cellBasis is not HDIV.");
    {
      const bool subcellBasisIsHGRAD = subcellBasis.getFunctionSpace() == FUNCTION_SPACE_HGRAD;
      const bool subcellBasisIsHVOL  = subcellBasis.getFunctionSpace() == FUNCTION_SPACE_HVOL;
      const bool cellIsTri  = cellBaseKey == shards::Triangle<>::key;
      const bool cellIsTet  = cellBaseKey == shards::Tetrahedron<>::key;
      const bool cellIsHex  = cellBaseKey == shards::Hexahedron<>::key;
      const bool cellIsQuad = cellBaseKey == shards::Quadrilateral<>::key;

      switch (subcellDim) {
      case 1: {
        //TODO: Hex, QUAD, TET and TRI element should have the same 1d basis
        if (cellIsHex || cellIsQuad) {
          INTREPID2_TEST_FOR_EXCEPTION( !(subcellBasisIsHGRAD||subcellBasisIsHVOL),
              std::logic_error,
              ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_DIV): "
              "subcellBasis function space (1d) is not consistent to cellBasis, which should be open line hgrad, order -1.");
        } else if (cellIsTet || cellIsTri) {
          INTREPID2_TEST_FOR_EXCEPTION( !subcellBasisIsHVOL,
              std::logic_error,
              ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_DIV): "
              "subcellBasis function space (1d) is not consistent to cellBasis, which should be HVOL line, order -1.");
        }
        break;
      }
      case 2: {
        if        (subcellBaseKey == shards::Quadrilateral<>::key) {
          // quad face basis is tensor product of open line basis functions
          INTREPID2_TEST_FOR_EXCEPTION( !(subcellBasisIsHGRAD||subcellBasisIsHVOL),
              std::logic_error,
              ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): "
              "subcellBasis function space is not compatible, which should be open line hgrad, order -1.");
        } else if (subcellBaseKey == shards::Triangle<>::key) {
          // triangle face basis comes from HVOL basis
          INTREPID2_TEST_FOR_EXCEPTION( !subcellBasisIsHVOL,
              std::logic_error,
              ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): "
              "subcellBasis function space is not compatible, which should HVOL, order-1.");
        }
        break;
      }
      }
    }
  }
}
#endif
} //Debug Namespace

template<typename OutputViewType,
typename subcellBasisType,
typename cellBasisType>
inline
void
OrientationTools::
getCoeffMatrix_HDIV(OutputViewType &output,
    const subcellBasisType& subcellBasis,
    const cellBasisType& cellBasis,
    const ordinal_type subcellId,
    const ordinal_type subcellOrt) {

#ifdef HAVE_INTREPID2_DEBUG
  Debug::check_getCoeffMatrix_HDIV(subcellBasis,cellBasis,subcellId,subcellOrt);
#endif

  using ScalarType = typename cellBasisType::scalarType;
  using ExecutionSpace = typename cellBasisType::ExecutionSpace;
  using HostExecutionSpace =
      typename Kokkos::Impl::is_space<ExecutionSpace>::host_mirror_space::execution_space;
  using OutputValueType = typename cellBasisType::OutputValueType;
  using PointValueType = typename cellBasisType::PointValueType;
  using BasisViewType = Kokkos::DynRankView<OutputValueType,ExecutionSpace>;
  using PointViewType = Kokkos::DynRankView<PointValueType,ExecutionSpace>;
  using ScalarViewType = Kokkos::DynRankView<ScalarType,ExecutionSpace>;


  //
  // Collocation points
  //
  const shards::CellTopology cellTopo = cellBasis.getBaseCellTopology();
  const shards::CellTopology subcellTopo = subcellBasis.getBaseCellTopology();
  const ordinal_type cellDim = cellTopo.getDimension();
  const ordinal_type subcellDim = subcellTopo.getDimension();
  const auto subcellBaseKey = subcellTopo.getBaseKey();
  const ordinal_type numCellBasis = cellBasis.getCardinality();
  const ordinal_type numSubcellBasis = subcellBasis.getCardinality();
  const ordinal_type ndofSubcell = cellBasis.getDofCount(subcellDim,subcellId);

  //
  // Reference points
  //

  //use lattice to compute reference subcell points \xi_j
  auto latticeDegree = (subcellBaseKey == shards::Triangle<>::key) ?
      cellBasis.getDegree()+2 : cellBasis.getDegree()+1;

  // Reference points \xi_j on the subcell
  PointViewType refPtsSubcell("refPtsSubcell", ndofSubcell, subcellDim);
  auto latticeSize=PointTools::getLatticeSize(subcellTopo, latticeDegree, 1);
  INTREPID2_TEST_FOR_EXCEPTION( latticeSize != ndofSubcell,
      std::logic_error,
      ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): " \
      "Lattice size should be equal to the onber of subcell internal DoFs");
  PointTools::getLattice(refPtsSubcell, subcellTopo, latticeDegree, 1);//, POINTTYPE_WARPBLEND);

  // evaluate values on the modified cell
  typename CellTools<ExecutionSpace>::subcellParamViewType subcellParam;
  CellTools<ExecutionSpace>::getSubcellParametrization(subcellParam, subcellDim, cellTopo);

  // refPtsCell = F_s (\eta_o (refPtsSubcell))
  PointViewType refPtsCell("refPtsCell", ndofSubcell, cellDim);
  // map points from the subcell manifold into the cell one
  mapSubcellCoordsToRefCell(refPtsCell,refPtsSubcell, subcellParam, subcellBaseKey, subcellId, subcellOrt);

  //computing normal to the subcell accounting for orientation
  ScalarViewType tangentsAndNormal("trJacobianF", cellDim, cellDim );
  OrientationTools::getRefSideTangentsAndNormal(tangentsAndNormal, subcellParam, subcellBaseKey, subcellId, subcellOrt);
  auto sideNormal = Kokkos::subview(tangentsAndNormal, cellDim-1, Kokkos::ALL());


  //
  // Basis evaluation on the collocation points
  //

  // cellBasisValues = \psi_k(F_s (\eta_o (\xi_j)))
  PointViewType cellBasisValues("cellBasisValues", numCellBasis, ndofSubcell, cellDim);
  cellBasis.getValues(cellBasisValues, refPtsCell, OPERATOR_VALUE);

  // subcellBasisValues = \phi_i (\xi_j)
  BasisViewType subCellValues("subCellValues", numSubcellBasis, ndofSubcell);
  subcellBasis.getValues(subCellValues, refPtsSubcell, OPERATOR_VALUE);
  ExecutionSpace().fence();

  //
  // Compute Psi_jk = \psi_k(F_s (\eta_o (\xi_j))) \cdot (n_s det(J_\eta))
  // and Phi_ji = \phi_i (\xi_j), and solve
  // Psi A^T = Phi
  //

  // construct Psi and Phi  matrices.  LAPACK wants left layout
  Kokkos::View<ScalarType**,Kokkos::LayoutLeft,ExecutionSpace>
  PsiMat("PsiMat", ndofSubcell, ndofSubcell),
  PhiMat("PhiMat", ndofSubcell, ndofSubcell);

  auto cellTagToOrdinal = Kokkos::create_mirror_view_and_copy(typename ExecutionSpace::memory_space(), cellBasis.getAllDofOrdinal());
  auto subcellTagToOrdinal = Kokkos::create_mirror_view_and_copy(typename ExecutionSpace::memory_space(), subcellBasis.getAllDofOrdinal());

  for (ordinal_type i=0;i<ndofSubcell;++i) {
    const ordinal_type ic = cellTagToOrdinal(subcellDim, subcellId, i);
    const ordinal_type isc = subcellTagToOrdinal(subcellDim, 0, i);
    for (ordinal_type j=0;j<ndofSubcell;++j) {
      PhiMat(j,i) = get_scalar_value(subCellValues(isc,j));
      for (ordinal_type k=0;k<cellDim;++k)
        PsiMat(j,i) += get_scalar_value(cellBasisValues(ic,j,k))*sideNormal(k);
    }
  }

  auto hostRefMat = Kokkos::create_mirror_view_and_copy(typename ExecutionSpace::memory_space(), PsiMat);
  auto hostOrtMat = Kokkos::create_mirror_view_and_copy(typename ExecutionSpace::memory_space(), PhiMat);

  // solve the system
  {
    Teuchos::LAPACK<ordinal_type, ScalarType> lapack;
    ordinal_type info = 0;
    Kokkos::View<ordinal_type*,Kokkos::LayoutLeft,HostExecutionSpace> pivVec("pivVec", ndofSubcell);

    lapack.GESV(ndofSubcell, ndofSubcell,
        hostRefMat.data(),
        hostRefMat.stride_1(),
        pivVec.data(),
        hostOrtMat.data(),
        hostOrtMat.stride_1(),
        &info);

    if (info) {
      std::stringstream ss;
      ss << ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): "
          << "LAPACK return with error code: "
          << info;
      INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str().c_str() );
    }

    //After solving the system w/ LAPACK, Phi contains A^T

    // transpose and clean up numerical noise (for permutation matrices)
    const double eps = tolerence();
    for (ordinal_type i=0;i<ndofSubcell;++i) {
      auto intmatii = std::round(hostOrtMat(i,i));
      hostOrtMat(i,i) = (std::abs(hostOrtMat(i,i) - intmatii) < eps) ? intmatii : hostOrtMat(i,i);
      for (ordinal_type j=i+1;j<ndofSubcell;++j) {
        auto matij = hostOrtMat(i,j);

        auto intmatji = std::round(hostOrtMat(j,i));
        hostOrtMat(i,j) = (std::abs(hostOrtMat(j,i) - intmatji) < eps) ? intmatji : hostOrtMat(j,i);

        auto intmatij = std::round(matij);
        hostOrtMat(j,i) = (std::abs(matij - intmatij) < eps) ? intmatij : matij;
      }
    }

    // Print Matrix A
    /*
    {
      std::cout  << "|";
      for (ordinal_type i=0;i<ndofSubcell;++i) {
        for (ordinal_type j=0;j<ndofSubcell;++j) {
          std::cout << hostOrtMat(i,j) << " ";
        }
        std::cout  << "| ";
      }
      std::cout <<std::endl;
    }
    */

  }

  {
    // move the data to original device memory
    const Kokkos::pair<ordinal_type,ordinal_type> range(0, ndofSubcell);
    Kokkos::deep_copy(Kokkos::subview(output, range, range), hostOrtMat);
  }
}
}

}
#endif
