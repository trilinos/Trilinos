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
    \li     \f$\xi_j\f$ are points in the reference subcell such that subcell bases are uniquely
    determined by the values they take on these points.

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
namespace Debug {

#ifdef HAVE_INTREPID2_DEBUG
template<typename subcellBasisType,
typename cellBasisType>
inline
void
check_getCoeffMatrix_HGRAD(const subcellBasisType& subcellBasis,
    const cellBasisType& cellBasis,
    const ordinal_type subcellId,
    const ordinal_type subcellOrt) {

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
    const bool cellBasisIsHGRAD = cellBasis.getFunctionSpace() == FUNCTION_SPACE_HGRAD;
    const bool subcellBasisIsHGRAD = subcellBasis.getFunctionSpace() == FUNCTION_SPACE_HGRAD;
    if (cellBasisIsHGRAD) {
      INTREPID2_TEST_FOR_EXCEPTION( !subcellBasisIsHGRAD,
          std::logic_error,
          ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HGRAD): " \
          "subcellBasis function space is not consistent to cellBasis.");
    }

    INTREPID2_TEST_FOR_EXCEPTION( subcellBasis.getDegree() != cellBasis.getDegree(),
        std::logic_error,
        ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HGRAD): " \
        "subcellBasis has a different polynomial degree from cellBasis' degree.");
  }
}
#endif
} // Debug namespace

template<typename OutputViewType,
typename subcellBasisType,
typename cellBasisType>
inline
void
OrientationTools::
getCoeffMatrix_HGRAD(OutputViewType &output,
    const subcellBasisType& subcellBasis,
    const cellBasisType& cellBasis,
    const ordinal_type subcellId,
    const ordinal_type subcellOrt) {

#ifdef HAVE_INTREPID2_DEBUG
  Debug::check_getCoeffMatrix_HGRAD(subcellBasis,cellBasis,subcellId,subcellOrt);
#endif

  using ScalarType = typename cellBasisType::scalarType;
  using ExecutionSpace = typename cellBasisType::ExecutionSpace;
  using HostExecutionSpace =
      typename Kokkos::Impl::is_space<ExecutionSpace>::host_mirror_space::execution_space;
  using OutputValueType = typename cellBasisType::OutputValueType;
  using PointValueType = typename cellBasisType::PointValueType;
  using BasisViewType = Kokkos::DynRankView<OutputValueType,ExecutionSpace>;
  using PointViewType = Kokkos::DynRankView<PointValueType,ExecutionSpace>;

  //
  // Topology
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

  // Reference points \xi_j on the subcell
  PointViewType refPtsSubcell("refPtsSubcell", ndofSubcell, subcellDim);
  auto latticeSize=PointTools::getLatticeSize(subcellTopo, subcellBasis.getDegree(), 1);

  INTREPID2_TEST_FOR_EXCEPTION( latticeSize != ndofSubcell,
      std::logic_error,
      ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HGRAD): " \
      "Lattice size should be equal to the onber of subcell internal DoFs");
  PointTools::getLattice(refPtsSubcell, subcellTopo, subcellBasis.getDegree(), 1, POINTTYPE_WARPBLEND);

  // map the points into the parent, cell accounting for orientation
  typename CellTools<ExecutionSpace>::subcellParamViewType subcellParam;
  CellTools<ExecutionSpace>::getSubcellParametrization(subcellParam, subcellDim, cellTopo);
  PointViewType refPtsCell("refPtsCell", ndofSubcell, cellDim);
  // refPtsCell = F_s (\eta_o (refPtsSubcell))
  mapSubcellCoordsToRefCell(refPtsCell,refPtsSubcell, subcellParam, subcellBaseKey, subcellId, subcellOrt);

  //
  // Bases evaluation on the reference points
  //

  // cellBasisValues = \psi_k(F_s (\eta_o (\xi_j)))
  BasisViewType cellBasisValues("cellBasisValues", numCellBasis, ndofSubcell);

  // subcellBasisValues = \phi_i (\xi_j)
  BasisViewType subcellBasisValues("subcellBasisValues", numSubcellBasis, ndofSubcell);

  cellBasis.getValues(cellBasisValues, refPtsCell, OPERATOR_VALUE);
  subcellBasis.getValues(subcellBasisValues, refPtsSubcell, OPERATOR_VALUE);
  ExecutionSpace().fence();

  //
  // Compute Psi_jk = \psi_k(F_s (\eta_o (\xi_j))) and Phi_ji = \phi_i (\xi_j),
  // and solve
  // Psi A^T = Phi
  //

  // construct Psi and Phi  matrices.  LAPACK wants left layout
  Kokkos::View<ScalarType**,Kokkos::LayoutLeft,ExecutionSpace>
  PsiMat("PsiMat", ndofSubcell, ndofSubcell),
  PhiMat("PhiMat", ndofSubcell, ndofSubcell);

  auto cellTagToOrdinal = Kokkos::create_mirror_view_and_copy(typename ExecutionSpace::memory_space(), cellBasis.getAllDofOrdinal());
  auto subcellTagToOrdinal = Kokkos::create_mirror_view_and_copy(typename ExecutionSpace::memory_space(), subcellBasis.getAllDofOrdinal());

  Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>(0, ndofSubcell),
        KOKKOS_LAMBDA (const int i){
    const ordinal_type ic = cellTagToOrdinal(subcellDim, subcellId, i);
    const ordinal_type isc = subcellTagToOrdinal(subcellDim, 0, i);
    for (ordinal_type j=0;j<ndofSubcell;++j) {
      PsiMat(j, i) = get_scalar_value(cellBasisValues(ic,j));
      PhiMat(j, i) = get_scalar_value(subcellBasisValues(isc,j));
    }
  });

  auto hostPsiMat = Kokkos::create_mirror_view_and_copy(typename ExecutionSpace::memory_space(), PsiMat);
  auto hostPhiMat = Kokkos::create_mirror_view_and_copy(typename ExecutionSpace::memory_space(), PhiMat);

  // Solve the system
  {
    Teuchos::LAPACK<ordinal_type,ScalarType> lapack;
    ordinal_type info = 0;

    Kokkos::View<ordinal_type*,Kokkos::LayoutLeft,HostExecutionSpace> pivVec("pivVec", ndofSubcell);
    lapack.GESV(ndofSubcell, ndofSubcell,
        hostPsiMat.data(),
        hostPsiMat.stride_1(),
        pivVec.data(),
        hostPhiMat.data(),
        hostPhiMat.stride_1(),
        &info);

    if (info) {
      std::stringstream ss;
      ss << ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HGRAD): "
          << "LAPACK return with error code: "
          << info;
      INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error, ss.str().c_str() );
    }

    //After solving the system w/ LAPACK, Phi contains A^T


    // transpose B and clean up numerical noise (for permutation matrices)
    const double eps = tolerence();
    for (ordinal_type i=0;i<ndofSubcell;++i) {
      auto intmatii = std::round(hostPhiMat(i,i));
      hostPhiMat(i,i) = (std::abs(hostPhiMat(i,i) - intmatii) < eps) ? intmatii : hostPhiMat(i,i);
      for (ordinal_type j=i+1;j<ndofSubcell;++j) {
        auto matij = hostPhiMat(i,j);

        auto intmatji = std::round(hostPhiMat(j,i));
        hostPhiMat(i,j) = (std::abs(hostPhiMat(j,i) - intmatji) < eps) ? intmatji : hostPhiMat(j,i);

        auto intmatij = std::round(matij);
        hostPhiMat(j,i) = (std::abs(matij - intmatij) < eps) ? intmatij : matij;
      }
    }

  }

  // Print A Matrix
  /*
  {
    std::cout  << "|";
    for (ordinal_type i=0;i<ndofSubcell;++i) {
      for (ordinal_type j=0;j<ndofSubcell;++j) {
        std::cout << PhiMat(i,j) << " ";
      }
      std::cout  << "| ";
    }
    std::cout <<std::endl;
  }
  */

  {
    // move the data to original device memory
    const Kokkos::pair<ordinal_type,ordinal_type> range(0, ndofSubcell);
    Kokkos::deep_copy(Kokkos::subview(output, range, range), hostPhiMat);
  }
}
}

}
#endif
