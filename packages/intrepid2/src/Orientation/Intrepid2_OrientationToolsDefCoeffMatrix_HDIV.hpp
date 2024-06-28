// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
      cellBaseKey != shards::Tetrahedron<>::key &&
      cellBaseKey != shards::Pyramid<>::key,
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
typename subcellBasisHostType,
typename cellBasisHostType>
inline
void
OrientationTools::
getCoeffMatrix_HDIV(OutputViewType &output,
                    const subcellBasisHostType& subcellBasis,
                    const cellBasisHostType& cellBasis,
                    const ordinal_type subcellId,
                    const ordinal_type subcellOrt,
                    const bool inverse) {
  
#ifdef HAVE_INTREPID2_DEBUG
  Debug::check_getCoeffMatrix_HDIV(subcellBasis,cellBasis,subcellId,subcellOrt);
#endif

  using value_type = typename OutputViewType::non_const_value_type;
  using host_device_type = Kokkos::HostSpace::device_type;

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
  Kokkos::DynRankView<value_type,host_device_type> refPtsSubcell("refPtsSubcell", ndofSubcell, subcellDim);
  auto latticeSize=PointTools::getLatticeSize(subcellTopo, latticeDegree, 1);
  INTREPID2_TEST_FOR_EXCEPTION( latticeSize != ndofSubcell,
      std::logic_error,
      ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HDIV): " \
      "Lattice size should be equal to the onber of subcell internal DoFs");
  PointTools::getLattice(refPtsSubcell, subcellTopo, latticeDegree, 1);//, POINTTYPE_WARPBLEND);

  // evaluate values on the modified cell
  auto subcellParam = Intrepid2::RefSubcellParametrization<host_device_type>::get(subcellDim, cellTopo.getKey());

  // refPtsCell = F_s (\eta_o (refPtsSubcell))
  Kokkos::DynRankView<value_type,host_device_type> refPtsCell("refPtsCell", ndofSubcell, cellDim);
  // map points from the subcell manifold into the cell one
  mapSubcellCoordsToRefCell(refPtsCell,refPtsSubcell, subcellParam, subcellBaseKey, subcellId, subcellOrt);

  //computing normal to the subcell accounting for orientation
  Kokkos::DynRankView<value_type,host_device_type> tangentsAndNormal("trJacobianF", cellDim, cellDim );
  OrientationTools::getRefSideTangentsAndNormal(tangentsAndNormal, subcellParam, subcellBaseKey, subcellId, subcellOrt);
  auto sideNormal = Kokkos::subview(tangentsAndNormal, cellDim-1, Kokkos::ALL());


  //
  // Basis evaluation on the collocation points
  //

  // cellBasisValues = \psi_k(F_s (\eta_o (\xi_j)))
  Kokkos::DynRankView<value_type,host_device_type> cellBasisValues("cellBasisValues", numCellBasis, ndofSubcell, cellDim);
  cellBasis.getValues(cellBasisValues, refPtsCell, OPERATOR_VALUE);

  // subcellBasisValues = \phi_i (\xi_j)
  Kokkos::DynRankView<value_type,host_device_type> subCellValues("subCellValues", numSubcellBasis, ndofSubcell);
  subcellBasis.getValues(subCellValues, refPtsSubcell, OPERATOR_VALUE);

  //
  // Compute Psi_jk = \psi_k(F_s (\eta_o (\xi_j))) \cdot (n_s det(J_\eta))
  // and Phi_ji = \phi_i (\xi_j), and solve
  // Psi A^T = Phi
  //

  // construct Psi and Phi  matrices.  LAPACK wants left layout
  Kokkos::DynRankView<value_type,Kokkos::LayoutLeft,host_device_type>
    PsiMat("PsiMat", ndofSubcell, ndofSubcell),
    PhiMat("PhiMat", ndofSubcell, ndofSubcell),
    RefMat,
    OrtMat;

  auto cellTagToOrdinal = cellBasis.getAllDofOrdinal();
  auto subcellTagToOrdinal = subcellBasis.getAllDofOrdinal();

  for (ordinal_type i=0;i<ndofSubcell;++i) {
    const ordinal_type ic = cellTagToOrdinal(subcellDim, subcellId, i);
    const ordinal_type isc = subcellTagToOrdinal(subcellDim, 0, i);
    for (ordinal_type j=0;j<ndofSubcell;++j) {
      PhiMat(j,i) = get_scalar_value(subCellValues(isc,j));
      for (ordinal_type k=0;k<cellDim;++k)
        PsiMat(j,i) += get_scalar_value(cellBasisValues(ic,j,k))*sideNormal(k);
    }
  }

  RefMat = inverse ? PhiMat : PsiMat;
  OrtMat = inverse ? PsiMat : PhiMat;

  // solve the system
  {
    Teuchos::LAPACK<ordinal_type,value_type> lapack;
    ordinal_type info = 0;
    Kokkos::DynRankView<ordinal_type,host_device_type> pivVec("pivVec", ndofSubcell);

    lapack.GESV(ndofSubcell, ndofSubcell,
        RefMat.data(),
        RefMat.stride_1(),
        pivVec.data(),
        OrtMat.data(),
        OrtMat.stride_1(),
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
      auto intmatii = std::round(OrtMat(i,i));
      OrtMat(i,i) = (std::abs(OrtMat(i,i) - intmatii) < eps) ? intmatii : OrtMat(i,i);
      for (ordinal_type j=i+1;j<ndofSubcell;++j) {
        auto matij = OrtMat(i,j);

        auto intmatji = std::round(OrtMat(j,i));
        OrtMat(i,j) = (std::abs(OrtMat(j,i) - intmatji) < eps) ? intmatji : OrtMat(j,i);

        auto intmatij = std::round(matij);
        OrtMat(j,i) = (std::abs(matij - intmatij) < eps) ? intmatij : matij;
      }
    }

    // Print Matrix A
    /*
    {
      std::cout  << "|";
      for (ordinal_type i=0;i<ndofSubcell;++i) {
        for (ordinal_type j=0;j<ndofSubcell;++j) {
          std::cout << OrtMat(i,j) << " ";
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
    auto suboutput = Kokkos::subview(output, range, range);
    auto tmp = Kokkos::create_mirror_view_and_copy(typename OutputViewType::device_type::memory_space(), OrtMat);
    Kokkos::deep_copy(suboutput, tmp);
  }
}
}

}
#endif
