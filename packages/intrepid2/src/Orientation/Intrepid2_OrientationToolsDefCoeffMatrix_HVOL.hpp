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
                    const ordinal_type cellOrt) {

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
    PhiMat("PhiMat", cardinality, cardinality);
  
  auto cellTagToOrdinal = cellBasis.getAllDofOrdinal();

  double jacDet;
  Intrepid2::Impl::OrientationTools::getJacobianDetOfOrientationMap(&jacDet,cellBaseKey,cellOrt);

  for (ordinal_type i=0;i<cardinality;++i) {
    const ordinal_type ic = cellTagToOrdinal(cellDim, 0, i);
    for (ordinal_type j=0;j<cardinality;++j) {
      PsiMat(j, i) = cellBasisValues(ic,j)*jacDet;
      PhiMat(j, i) = nonOrientedBasisValues(ic,j);
    }
  }

  // Solve the system
  {
    Teuchos::LAPACK<ordinal_type,value_type> lapack;
    ordinal_type info = 0;

    Kokkos::DynRankView<ordinal_type,Kokkos::LayoutLeft,host_device_type> pivVec("pivVec", cardinality);
    lapack.GESV(cardinality, cardinality,
                PsiMat.data(),
                PsiMat.stride_1(),
                pivVec.data(),
                PhiMat.data(),
                PhiMat.stride_1(),
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
      auto intmatii = std::round(PhiMat(i,i));
      PhiMat(i,i) = (std::abs(PhiMat(i,i) - intmatii) < eps) ? intmatii : PhiMat(i,i);
      for (ordinal_type j=i+1;j<cardinality;++j) {
        auto matij = PhiMat(i,j);

        auto intmatji = std::round(PhiMat(j,i));
        PhiMat(i,j) = (std::abs(PhiMat(j,i) - intmatji) < eps) ? intmatji : PhiMat(j,i);

        auto intmatij = std::round(matij);
        PhiMat(j,i) = (std::abs(matij - intmatij) < eps) ? intmatij : matij;
      }
    }

  }

  // Print A Matrix
  /*
  {
    std::cout  << "Ort: " << cellOrt << ": |";
    for (ordinal_type i=0;i<cardinality;++i) {
      for (ordinal_type j=0;j<cardinality;++j) {
        std::cout << PhiMat(i,j) << " ";
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
    auto tmp = Kokkos::create_mirror_view_and_copy(typename OutputViewType::device_type::memory_space(), PhiMat);
    Kokkos::deep_copy(suboutput, tmp);
  }
}
}

}
#endif
