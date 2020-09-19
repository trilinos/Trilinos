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


/** \file   Intrepid2_OrientationToolsDefCoeffMatrix_HCURL.hpp
    \brief  Creation of orientation matrix A of a face or edge for HCURL elements

    \li     \f$\sum_k A_ik \psi_k(F_s (\eta_o (\xi_j))) \cdot (J_F J_\eta t_j) = \phi_i (\xi_j) \dot t_j\f$ where
    \li     \f$\psi_k\f$ are the basis functions of the reference cell,
    \li     \f$\phi_i\f$ are the basis function of the subcell,
    \li     \f$t_j\f$ is the tangent on the reference subcell associated to dof j
    \li     \f$F_s\f$ is the Jacobian of the map from the reference subcell to the subcell s of the reference cell (s == subcellId),
    \li     \f$J_F\f$ is the Jacobian of the map \f$F_s\f$,
    \li     \f$\eta_o\f$ is the orientation map o associated to the subcell s (o == subcellOrt),
    \li     \f$J_{\eta}\f$ is the Jacobian of the map \eta,
    \li     \f$\xi_j\f$ are points of the subcell manifold

    Note: the points \xi_j and tangent vectors t_j are chosen such that the bases \phi_i are
    uniquely identified by the values \phi_i(\xi_j) \dot t_j.

    \author Created by Kyungjoo Kim
 */
#ifndef __INTREPID2_ORIENTATIONTOOLS_DEF_COEFF_MATRIX_HCURL_HPP__
#define __INTREPID2_ORIENTATIONTOOLS_DEF_COEFF_MATRIX_HCURL_HPP__

#include "Intrepid2_HCURL_QUAD_In_FEM.hpp"
#include "Intrepid2_HCURL_TRI_In_FEM.hpp"
#include "Intrepid2_HVOL_LINE_Cn_FEM.hpp"

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
check_getCoeffMatrix_HCURL(const subcellBasisType& subcellBasis,
    const cellBasisType& cellBasis,
    const ordinal_type subcellId,
    const ordinal_type subcellOrt) {
  const shards::CellTopology cellTopo = cellBasis.getBaseCellTopology();
  const shards::CellTopology subcellTopo = subcellBasis.getBaseCellTopology();

  const ordinal_type cellDim = cellTopo.getDimension();
  const ordinal_type subcellDim = subcellTopo.getDimension();



  INTREPID2_TEST_FOR_EXCEPTION( subcellDim >= cellDim,
      std::logic_error,
      ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HCURL): " \
      "cellDim must be greater than subcellDim.");

  const auto subcellBaseKey = subcellTopo.getBaseKey();
  const auto cellBaseKey = cellTopo.getBaseKey();

  INTREPID2_TEST_FOR_EXCEPTION( subcellBaseKey != shards::Line<>::key &&
      subcellBaseKey != shards::Quadrilateral<>::key &&
      subcellBaseKey != shards::Triangle<>::key,
      std::logic_error,
      ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HCURL): " \
      "subcellBasis must have line, quad, or triangle topology.");

  INTREPID2_TEST_FOR_EXCEPTION( cellBaseKey != shards::Quadrilateral<>::key &&
      cellBaseKey != shards::Triangle<>::key &&
      cellBaseKey != shards::Hexahedron<>::key &&
      cellBaseKey != shards::Wedge<>::key &&
      cellBaseKey != shards::Tetrahedron<>::key,
      std::logic_error,
      ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HCURL): " \
      "cellBasis must have quad, triangle, hexhedron or tetrahedron topology.");

  //
  // Function space
  //
  {
    const bool cellBasisIsHCURL = cellBasis.getFunctionSpace() == FUNCTION_SPACE_HCURL;


    if (cellBasisIsHCURL) {
      const bool subcellBasisIsHGRAD = subcellBasis.getFunctionSpace() == FUNCTION_SPACE_HGRAD;
      const bool subcellBasisIsHVOL  = subcellBasis.getFunctionSpace() == FUNCTION_SPACE_HVOL;
      const bool subcellBasisIsHCURL = subcellBasis.getFunctionSpace() == FUNCTION_SPACE_HCURL;
      const bool cellIsTri  = cellBaseKey == shards::Triangle<>::key;
      const bool cellIsTet  = cellBaseKey == shards::Tetrahedron<>::key;
      const bool cellIsHex  = cellBaseKey == shards::Hexahedron<>::key;
      const bool cellIsQuad = cellBaseKey == shards::Quadrilateral<>::key;


      // edge hcurl is hgrad with gauss legendre points
      switch (subcellDim) {
      case 1: {
        //TODO: Hex, QUAD, TET and TRI element should have the same 1d basis
        if (cellIsHex || cellIsQuad) {
          INTREPID2_TEST_FOR_EXCEPTION( !(subcellBasisIsHGRAD||subcellBasisIsHVOL),
              std::logic_error,
              ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HCURL): "
              "subcellBasis function space (1d) is not consistent to cellBasis.");
        } else if (cellIsTet || cellIsTri) {
          INTREPID2_TEST_FOR_EXCEPTION( !subcellBasisIsHVOL,
              std::logic_error,
              ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HCURL): "
              "subcellBasis function space (1d) is not consistent to cellBasis.");
        }
        break;
      }
      case 2: {
        INTREPID2_TEST_FOR_EXCEPTION( !subcellBasisIsHCURL,
            std::logic_error,
            ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HCURL): "
            "subcellBasis function space (2d) is not consistent to cellBasis.");
        break;
      }
      }
    }
  }
}
#endif
}

template<typename OutputViewType,
typename subcellBasisType,
typename cellBasisType>
inline
void
OrientationTools::
getCoeffMatrix_HCURL(OutputViewType &output,
    const subcellBasisType& subcellBasis,
    const cellBasisType& cellBasis,
    const ordinal_type subcellId,
    const ordinal_type subcellOrt) {

#ifdef HAVE_INTREPID2_DEBUG
  Debug::check_getCoeffMatrix_HCURL(subcellBasis,cellBasis,subcellId,subcellOrt);
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
  // Topology
  //
  // populate points on a subcell and map to subcell
  const shards::CellTopology cellTopo = cellBasis.getBaseCellTopology();
  const shards::CellTopology subcellTopo = subcellBasis.getBaseCellTopology();
  const ordinal_type cellDim = cellTopo.getDimension();
  const ordinal_type subcellDim = subcellTopo.getDimension();
  const auto subcellBaseKey = subcellTopo.getBaseKey();
  const ordinal_type numCellBasis = cellBasis.getCardinality();
  const ordinal_type numSubcellBasis = subcellBasis.getCardinality();
  const ordinal_type ndofSubcell = cellBasis.getDofCount(subcellDim,subcellId);

  // Compute reference points \xi_j and tangents t_j on the subcell
  // To do so we use the DoF coordinates and DoF coefficients of a Lagrangian HCURL basis
  // on the subcell spanning the same space as the bases \phi_j

  // Tangents t_j
  ScalarViewType subcellTangents("subcellTangents", numSubcellBasis, subcellDim);
  auto degree = subcellBasis.getDegree();
  BasisPtr<ExecutionSpace, ScalarType, ScalarType> basisPtr;
  if(subcellBaseKey == shards::Line<>::key) {
    basisPtr = Teuchos::rcp(new Intrepid2::Basis_HVOL_LINE_Cn_FEM<ExecutionSpace, ScalarType, ScalarType>(degree));
    basisPtr->getDofCoeffs(Kokkos::subview(subcellTangents, Kokkos::ALL(),0));
  } else if (subcellBaseKey == shards::Triangle<>::key) {
    basisPtr = Teuchos::rcp(new Intrepid2::Basis_HCURL_TRI_In_FEM<ExecutionSpace, ScalarType, ScalarType>(degree));
    basisPtr->getDofCoeffs(subcellTangents);
  } else if (subcellBaseKey == shards::Quadrilateral<>::key) {
    basisPtr =  Teuchos::rcp(new Intrepid2::Basis_HCURL_QUAD_In_FEM<ExecutionSpace, ScalarType, ScalarType>(degree));
    basisPtr->getDofCoeffs(subcellTangents);
  }

  // coordinates \xi_j
  ScalarViewType subcellDofCoords("subcellDofCoords", basisPtr->getCardinality(), subcellDim);
  basisPtr->getDofCoords(subcellDofCoords);
  INTREPID2_TEST_FOR_EXCEPTION( basisPtr->getDofCount(subcellDim,0) != ndofSubcell,
      std::logic_error,
      ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HCURL): " \
      "the number of basisPtr internal DoFs should equate those of the subcell");

  // restrict \xi_j (and corresponding t_j) to points internal to the HCURL basis
  PointViewType refPtsSubcell("refPtsSubcell", ndofSubcell, subcellDim);
  ScalarViewType refSubcellTangents("subcellTangents", ndofSubcell, subcellDim);
  auto tagToOrdinal = Kokkos::create_mirror_view_and_copy(typename ExecutionSpace::memory_space(), basisPtr->getAllDofOrdinal());
  Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>(0, ndofSubcell),
      KOKKOS_LAMBDA (const int i){
    ordinal_type isc = tagToOrdinal(subcellDim, 0, i);
    for(ordinal_type d=0; d <subcellDim; ++d){
      refPtsSubcell(i,d) = subcellDofCoords(isc,d);
      for(ordinal_type k=0; k <subcellDim; ++k)
        refSubcellTangents(i,d) = subcellTangents(isc,d);
    }
  });

  //
  // Bases evaluation on the reference points
  //

  // subcellBasisValues = \phi_i (\xi_j)
  BasisViewType subCellValues("subCellValues", numSubcellBasis, ndofSubcell, subcellDim);
  if(subcellDim==1) {
    auto lineValues = Kokkos::subview(subCellValues, Kokkos::ALL(), Kokkos::ALL(), 0);
    subcellBasis.getValues(lineValues, refPtsSubcell, OPERATOR_VALUE);
  } else {
    subcellBasis.getValues(subCellValues, refPtsSubcell, OPERATOR_VALUE);
  }

  //
  // Basis evaluation on the reference points
  //

  typename CellTools<ExecutionSpace>::subcellParamViewType subcellParam;
  CellTools<ExecutionSpace>::getSubcellParametrization(subcellParam, subcellDim, cellTopo);

  // refPtsCell = F_s (\eta_o (refPtsSubcell))
  PointViewType refPtsCell("refPtsCell", ndofSubcell, cellDim);
  mapSubcellCoordsToRefCell(refPtsCell,refPtsSubcell, subcellParam, subcellBaseKey, subcellId, subcellOrt);


  //mapping tangents t_j into parent cell, i.e. computing J_F J_\eta t_j
  ScalarViewType trJacobianF("trJacobianF", subcellDim, cellDim );
  OrientationTools::getRefSubcellTangents(trJacobianF, subcellParam, subcellBaseKey, subcellId, subcellOrt);



  // cellBasisValues = \psi_k(F_s (\eta_o (\xi_j)))
  BasisViewType cellBasisValues("cellBasisValues", numCellBasis, ndofSubcell, cellDim);
  cellBasis.getValues(cellBasisValues, refPtsCell, OPERATOR_VALUE);
  ExecutionSpace().fence();

  //
  // Compute Psi_jk = \psi_k(F_s (\eta_o (\xi_j))) \cdot (J_F J_\eta t_j)
  // and Phi_ji = \phi_i (\xi_j) \ctot t_j, and solve
  // Psi A^T = Phi
  //

  // construct Psi and Phi  matrices.  LAPACK wants left layout
  Kokkos::View<ScalarType**,Kokkos::LayoutLeft,ExecutionSpace> // left layout for lapack
  PsiMat("PsiMat", ndofSubcell, ndofSubcell),
  PhiMat("PhiMat", ndofSubcell, ndofSubcell);

  auto cellTagToOrdinal = Kokkos::create_mirror_view_and_copy(typename ExecutionSpace::memory_space(), cellBasis.getAllDofOrdinal());
  auto subcellTagToOrdinal = Kokkos::create_mirror_view_and_copy(typename ExecutionSpace::memory_space(), subcellBasis.getAllDofOrdinal());

  Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>(0, ndofSubcell),
      KOKKOS_LAMBDA (const int i){
    const ordinal_type ic = cellTagToOrdinal(subcellDim, subcellId, i);
    for (ordinal_type j=0;j<ndofSubcell;++j) {
      const ordinal_type isc = subcellTagToOrdinal(subcellDim, 0, i);
      ScalarType refEntry = 0, ortEntry =0;
      for (ordinal_type k=0;k<subcellDim;++k) {
        ortEntry += get_scalar_value(subCellValues(isc,j,k))*refSubcellTangents(j,k);
        for (ordinal_type d=0; d<cellDim; ++d)
          refEntry +=  get_scalar_value(cellBasisValues(ic,j,d))*trJacobianF(k,d)*refSubcellTangents(j,k);
      }
      PsiMat(j,i) = refEntry;
      PhiMat(j,i) = ortEntry;
    }
  });

  auto hostRefMat = Kokkos::create_mirror_view_and_copy(typename ExecutionSpace::memory_space(), PsiMat);
  auto hostOrtMat = Kokkos::create_mirror_view_and_copy(typename ExecutionSpace::memory_space(), PhiMat);


  // Solve the system using Lapack
  {
    Teuchos::LAPACK<ordinal_type,ScalarType> lapack;
    ordinal_type info = 0;


    /*
        Kokkos::View<value_type**,Kokkos::LayoutLeft,host_space_type> work("pivVec", 2*ndofSubcell, 1);
        lapack.GELS('N', ndofSubcell*card, ndofSubcell, ndofSubcell,
            PsiMat.data(),
            PsiMat.stride_1(),
            PhiMat.data(),
            PhiMat.stride_1(),
            work.data(), work.extent(0),
            &info);

        */
    Kokkos::View<ordinal_type*,Kokkos::LayoutLeft,HostExecutionSpace> pivVec("pivVec", ndofSubcell);
    lapack.GESV(ndofSubcell, ndofSubcell,
        hostRefMat.data(),
        hostOrtMat.stride_1(),
        pivVec.data(),
        PhiMat.data(),
        PhiMat.stride_1(),
        &info);
    //*/
    if (info) {
      std::stringstream ss;
      ss << ">>> ERROR (Intrepid::OrientationTools::getCoeffMatrix_HCURL): "
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



    // Print A matrix
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
