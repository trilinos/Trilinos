// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file   Intrepid2_OrientationToolsDefCoeffMatrix_HVOL.hpp
    \brief  Creation of orientation matrix A of a face or edge for HGRAD elements
    \li     \f$\sum_k A_ik \psi_k(\eta_o (\xi_j)) det(J_\eta) = \phi_i (\xi_j) \f$ where
    \li     \f$\psi_k\f$ are the basis functions of the reference cell ,
    \li     \f$\phi_i\f$ are the basis function of the subcell
    \li     \f$\eta_o\f$ is the orientation map o associated to the subcell s (o == subcellOrt)
    \li     \f$J_{\eta}\f$ is the Jacobian of the map \f$\eta_o\f$,
    \li     \f$\xi_j\f$ are points in the reference subcell such that subcell bases are uniquely
    determined by the values they take on these points.

    \author Created by Kyungjoo Kim
 */
#ifndef __INTREPID2_ORIENTATIONTOOLS_DEF_COEFF_MATRIX_HVOL_HPP__
#define __INTREPID2_ORIENTATIONTOOLS_DEF_COEFF_MATRIX_HVOL_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {

namespace Impl {
namespace Debug {

#ifdef HAVE_INTREPID2_DEBUG
template<typename cellBasisType>
inline
void
check_getCoeffMatrix_HVOL(const cellBasisType& cellBasis,
    const ordinal_type cellOrt) {

  // populate points on a subcell and map to subcell
  const shards::CellTopology cellTopo = cellBasis.getBaseCellTopology();
  const ordinal_type cellDim = cellTopo.getDimension();

  INTREPID2_TEST_FOR_EXCEPTION( cellDim > 2,
      std::logic_error,
      ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HVOL): " \
      "HVOL orientation supported only for (side) cells with dimension less than 3.");

  const auto cellBaseKey = cellTopo.getBaseKey();

  INTREPID2_TEST_FOR_EXCEPTION( cellBaseKey != shards::Line<>::key &&
      cellBaseKey != shards::Quadrilateral<>::key &&
      cellBaseKey != shards::Triangle<>::key,
      std::logic_error,
      ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HVOL): " \
      "cellBasis must have line, quad, or triangle topology.");


  //
  // Function space
  //

  {
    INTREPID2_TEST_FOR_EXCEPTION( cellBasis.getFunctionSpace() != FUNCTION_SPACE_HVOL,
        std::logic_error,
        ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HVOL): " \
        "cellBasis function space is not HVOL.");
    }
  }
#endif
} // Debug namespace

template<typename OutputViewType,
typename cellBasisHostType>
inline
void
OrientationTools::
getCoeffMatrix_HVOL(OutputViewType &output, /// this is device view
                    const cellBasisHostType& cellBasis, /// this also must be host basis object
                    const ordinal_type cellOrt,
                    const bool inverse) {

#ifdef HAVE_INTREPID2_DEBUG
  Debug::check_getCoeffMatrix_HVOL(cellBasis,cellOrt);
#endif

  using host_device_type = typename Kokkos::HostSpace::device_type;
  using value_type = typename OutputViewType::non_const_value_type;

  //
  // Topology
  //

  const shards::CellTopology cellTopo = cellBasis.getBaseCellTopology();
  const ordinal_type cellDim = cellTopo.getDimension();
  const auto cellBaseKey = cellTopo.getBaseKey();
  const ordinal_type cardinality = cellBasis.getCardinality();

  //
  // Reference points
  //

  // Reference points \xi_j on the subcell
  Kokkos::DynRankView<value_type,host_device_type> refPtsCell("refPtsCell", cardinality, cellDim),refPtsCellNotOriented("refPtsCellNotOriented", cardinality, cellDim);

  ordinal_type latticeOffset(1);

  // this work for line and quadrilateral topologies
  ordinal_type latticeOrder = (cellTopo.getBaseKey() == shards::Triangle<>::key) ?
      cellBasis.getDegree() + 3 * latticeOffset :  // triangle
      cellBasis.getDegree() + 2 * latticeOffset;   // line and quad

  PointTools::getLattice(refPtsCellNotOriented, cellTopo, latticeOrder, 1, POINTTYPE_WARPBLEND);

  // map the points into the parent, cell accounting for orientation
  mapToModifiedReference(refPtsCell,refPtsCellNotOriented,cellBaseKey,cellOrt);

  //
  // Bases evaluation on the reference points
  //

  // cellBasisValues = \psi_k(\eta_o (\xi_j))
  Kokkos::DynRankView<value_type,host_device_type> cellBasisValues("cellBasisValues", cardinality, cardinality);

  // basisValues = \phi_i (\xi_j)
  Kokkos::DynRankView<value_type,host_device_type> nonOrientedBasisValues("subcellBasisValues", cardinality, cardinality);

  cellBasis.getValues(cellBasisValues, refPtsCell, OPERATOR_VALUE);
  cellBasis.getValues(nonOrientedBasisValues, refPtsCellNotOriented, OPERATOR_VALUE);

  //
  // Compute Psi_jk = \psi_k(\eta_o (\xi_j)) det(J_\eta) and Phi_ji = \phi_i (\xi_j),
  // and solve
  // Psi A^T = Phi
  //

  // construct Psi and Phi  matrices.  LAPACK wants left layout
  Kokkos::DynRankView<value_type,Kokkos::LayoutLeft,host_device_type>
    PsiMat("PsiMat", cardinality, cardinality),
    PhiMat("PhiMat", cardinality, cardinality),
    RefMat,
    OrtMat;
  
  auto cellTagToOrdinal = cellBasis.getAllDofOrdinal();

  Kokkos::DynRankView<value_type,host_device_type> jac("jacobian",cellDim,cellDim);
  Intrepid2::Impl::OrientationTools::getJacobianOfOrientationMap(jac,cellBaseKey,cellOrt);
  value_type jacDet(0);
  if(cellDim == 2) {
    jacDet = jac(0,0)*jac(1,1)-jac(0,1)*jac(1,0);
  } else { //celldim == 1
    jacDet = jac(0,0);
  }

  for (ordinal_type i=0;i<cardinality;++i) {
    const ordinal_type ic = cellTagToOrdinal(cellDim, 0, i);
    for (ordinal_type j=0;j<cardinality;++j) {
      PsiMat(j, i) = cellBasisValues(ic,j)*jacDet;
      PhiMat(j, i) = nonOrientedBasisValues(ic,j);
    }
  }

  RefMat = inverse ? PhiMat : PsiMat;
  OrtMat = inverse ? PsiMat : PhiMat;

  // Solve the system
  {
    Teuchos::LAPACK<ordinal_type,value_type> lapack;
    ordinal_type info = 0;

    Kokkos::DynRankView<ordinal_type,Kokkos::LayoutLeft,host_device_type> pivVec("pivVec", cardinality);
    lapack.GESV(cardinality, cardinality,
                RefMat.data(),
                RefMat.stride_1(),
                pivVec.data(),
                OrtMat.data(),
                OrtMat.stride_1(),
                &info);
    
    if (info) {
      std::stringstream ss;
      ss << ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HVOL): "
         << "LAPACK return with error code: "
         << info;
      INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str().c_str() );
    }
    
    //After solving the system w/ LAPACK, Phi contains A^T
    
    // transpose B and clean up numerical noise (for permutation matrices)
    const double eps = tolerence();
    for (ordinal_type i=0;i<cardinality;++i) {
      auto intmatii = std::round(OrtMat(i,i));
      OrtMat(i,i) = (std::abs(OrtMat(i,i) - intmatii) < eps) ? intmatii : OrtMat(i,i);
      for (ordinal_type j=i+1;j<cardinality;++j) {
        auto matij = OrtMat(i,j);

        auto intmatji = std::round(OrtMat(j,i));
        OrtMat(i,j) = (std::abs(OrtMat(j,i) - intmatji) < eps) ? intmatji : OrtMat(j,i);

        auto intmatij = std::round(matij);
        OrtMat(j,i) = (std::abs(matij - intmatij) < eps) ? intmatij : matij;
      }
    }

  }

  // Print A Matrix
  /*
  {
    std::cout  << "Ort: " << cellOrt << ": |";
    for (ordinal_type i=0;i<cardinality;++i) {
      for (ordinal_type j=0;j<cardinality;++j) {
        std::cout << OrtMat(i,j) << " ";
      }
      std::cout  << "| ";
    }
    std::cout <<std::endl;
  }
  */

  {
    // move the data to original device memory
    const Kokkos::pair<ordinal_type,ordinal_type> range(0, cardinality);
    auto suboutput = Kokkos::subview(output, range, range);
    auto tmp = Kokkos::create_mirror_view_and_copy(typename OutputViewType::device_type::memory_space(), OrtMat);
    Kokkos::deep_copy(suboutput, tmp);
  }
}
}

}
#endif
