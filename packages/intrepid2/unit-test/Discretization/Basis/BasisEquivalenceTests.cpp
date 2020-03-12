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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   BasisEquivalenceTests.cpp
    \brief  Tests to verify that bases that span the same space are equivalent.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_NodalBasisFamily.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_TestUtils.hpp"

#include "Intrepid2_HGRAD_LINE_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"

#include <Kokkos_Core.hpp>
//#include "KokkosBlas3.hpp"

#include <Teuchos_LAPACK.hpp>
#include <Teuchos_BLAS.hpp>

namespace
{
  //! solve a system Ax=b using LAPACK Cholesky factorization on host
  template<class AViewType, class BXViewType>
  int solveSystemUsingHostLapack(const AViewType &A_device, const BXViewType &bx_device)
  {
    using Scalar = typename AViewType::value_type;
    
    auto A_host  = getHostCopy(A_device);
    auto bx_host = getHostCopy(bx_device);
    
//    std::cout << "\nA = [";
//    for (int i=0; i<A_host.extent_int(0); i++)
//    {
//      std::cout << "[ ";
//      for (int j=0; j<A_host.extent(1); j++)
//      {
//        std::cout << A_host(i,j) << " ";
//      }
//      std::cout << "]; ";
//    }
//    std::cout << "]\n";
//
//    std::cout << "b = [";
//    for (int i=0; i<bx_host.extent_int(0); i++)
//    {
//      std::cout << "[ ";
//      for (int j=0; j<bx_host.extent(1); j++)
//      {
//        std::cout << bx_host(i,j) << " ";
//      }
//      std::cout << "]; ";
//    }
//    std::cout << "]\n";
    
    int N = A_host.extent_int(0);
    TEUCHOS_TEST_FOR_EXCEPTION(N != A_host.extent_int(1), std::invalid_argument, "A must be square!");
    
    int M = (bx_host.rank() == 1) ? 1 : bx_host.extent_int(1);
    
    // since A is SPD, col/row major has no effect on the data
    // but B's data may be transposed relative to what LAPACK expects (column-major order)
    // so we allocate our own storage for B to make sure of the ordering
    double B[N*M];
    
    for (int j=0; j<M; j++)
    {
      for (int i=0; i<N; i++)
      {
        B[i+j*N] = bx_host(i,j);
      }
    }
    
    char UPLO = 'L'; // lower-triangular
    
    int result = 0;
    int INFO;
    
    Teuchos::LAPACK<int, Scalar> lapack;
    
    // factor
    lapack.POTRF(UPLO, N, A_host.data(), N, &INFO);
    
    if (INFO != 0) result = INFO;
    if (INFO != 0) std::cout << "ERROR: got " << INFO << " from POTRF.\n";
    
    // back-substitute
    lapack.POTRS(UPLO, N, M, A_host.data(), N, B, N, &INFO);
    
    if (INFO != 0) result = INFO;
    
    // copy from our B container back to bx_host
    for (int j=0; j<M; j++)
    {
      for (int i=0; i<N; i++)
      {
        bx_host(i,j) = B[i+j*N];
      }
    }
    
//    std::cout << "x = [";
//    for (int i=0; i<bx_host.extent_int(0); i++)
//    {
//      std::cout << "[ ";
//      for (int j=0; j<bx_host.extent(1); j++)
//      {
//        std::cout << bx_host(i,j) << " ";
//      }
//      std::cout << "]; ";
//    }
//    std::cout << "]\n";
    
    // solution should be in bx_host; copy to device
    Kokkos::deep_copy(bx_device, bx_host);
    
    return result;
  }
  
  //! Computes C := A*B or C := A^T*B
  //! B is allowed to be either a (rank-2) matrix, or a higher-rank object.  Either way, B's first (row) index should match A's column count.
  template<class Rank2View>
  void deviceGeneralizedMatrixMultiply(Rank2View &A, bool transposeA, Rank2View &B, Rank2View &C)
  {
    using Scalar = typename Rank2View::value_type;
    const int N0 = transposeA ? A.extent_int(1) : A.extent_int(0);
    const int N1 = transposeA ? A.extent_int(0) : A.extent_int(1);
    TEUCHOS_TEST_FOR_EXCEPTION(N1 != B.extent_int(0), std::invalid_argument, "A column count should match B row count");
    const int N2 = B.extent_int(1);
    
    TEUCHOS_TEST_FOR_EXCEPTION(N0 != C.extent_int(0), std::invalid_argument, "C row count should match A row count");
    TEUCHOS_TEST_FOR_EXCEPTION(N2 != C.extent_int(1), std::invalid_argument, "C column count should match B column count");
    
    TEUCHOS_TEST_FOR_EXCEPTION(B.rank() != C.rank(), std::invalid_argument, "B's rank must match C's rank");
    for (unsigned dim=2; dim<B.rank(); dim++)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(B.extent_int(dim) != C.extent_int(dim), std::invalid_argument, "B and C must agree in all dimensions beyond the first two");
    }
    using ViewIteratorScalar = Intrepid2::ViewIterator<ViewType<Scalar>, Scalar>;
    Kokkos::parallel_for(N0, KOKKOS_LAMBDA(const int A_row_ordinal)
    {
      ViewIteratorScalar B_viewIterator(B);
      ViewIteratorScalar C_viewIterator(C);
      
      for (int B_col_ordinal=0; B_col_ordinal<N2; B_col_ordinal++)
      {
        C_viewIterator.setLocation({A_row_ordinal,B_col_ordinal,0,0,0,0,0});
        do
        {
          Scalar value = 0.0;
          auto B_location = C_viewIterator.getLocation();
          B_location[0] = 0;
          B_location[1] = B_col_ordinal;
          B_viewIterator.setLocation(B_location);
          for (int k=0; k<N1; k++)
          {
            B_viewIterator.getLocation()[0] = k;
            if (transposeA)
            {
              value += A(k,A_row_ordinal) * B_viewIterator.get();
            }
            else
            {
              value += A(A_row_ordinal,k) * B_viewIterator.get();
            }
          }
          C_viewIterator.set(value);
        }
        while (C_viewIterator.increment() > 1); // increment returns the rank of the leftmost index that was changed; when it changes B_col_ordinal, it's time to stop.
      }
    });
    Kokkos::fence();
  }

  using namespace Intrepid2;

  //! tests that two bases are equivalent; computes a conversion from one to the other and then uses that to confirm that the equivalence
  //! holds for OPERATOR_VALUE (as it should by construction), as well as the operators passed in in opsToTest.
  template<class Basis1, class Basis2>
  void testBasisEquivalence(Basis1 &basis1, Basis2 &basis2, const std::vector<EOperator> &opsToTest,
                            const double relTol, const double absTol, Teuchos::FancyOStream &out, bool &success)
  {
    // first, check that the bases agree on cardinality
    TEST_ASSERT(basis1.getCardinality() == basis2.getCardinality());
    
    const int basisCardinality = basis1.getCardinality();
    
    // get quadrature points for integrating up to 2*polyOrder
    const int quadratureDegree = 2*basis1.getDegree();
    using ExecutionSpace = typename Basis1::ExecutionSpace;
    using PointScalar = typename Basis1::PointValueType;
    using WeightScalar = typename Basis1::OutputValueType;
    using Scalar = WeightScalar;
    Intrepid2::DefaultCubatureFactory cub_factory;
    auto cellTopoKey = basis1.getBaseCellTopology().getKey();
    auto quadrature = cub_factory.create<ExecutionSpace, PointScalar, WeightScalar>(cellTopoKey, quadratureDegree);
    Intrepid2::ordinal_type numRefPoints = quadrature->getNumPoints();
    const int spaceDim = basis1.getBaseCellTopology().getDimension();
    ViewType<PointScalar> points = ViewType<PointScalar>("quadrature points 1D ref cell", numRefPoints, spaceDim);
    ViewType<WeightScalar> weights = ViewType<WeightScalar>("quadrature weights 1D ref cell", numRefPoints);
    quadrature->getCubature(points, weights);
    
    auto functionSpace = basis1.getFunctionSpace();
    
    // set up a projection of basis2 onto basis1
    using Scalar = typename Basis1::OutputValueType;
    ViewType<Scalar> basis1Values = getOutputView<Scalar>(functionSpace, OPERATOR_VALUE, basisCardinality, numRefPoints, spaceDim);
    ViewType<Scalar> basis2Values = getOutputView<Scalar>(functionSpace, OPERATOR_VALUE, basisCardinality, numRefPoints, spaceDim);
    
    basis1.getValues(basis1Values, points, OPERATOR_VALUE);
    basis2.getValues(basis2Values, points, OPERATOR_VALUE);
    
    // integrate basis1 against itself to compute the SPD matrix A that we'll use to set up the basis conversion system
    ViewType<Scalar> basis1_vs_basis1 = ViewType<Scalar>("basis 1 vs basis 1", basisCardinality, basisCardinality);
    // integrate basis1 against basis2 to compute the RHS b for the basis conversion system
    ViewType<Scalar> basis1_vs_basis2 = ViewType<Scalar>("basis 1 vs basis 2", basisCardinality, basisCardinality);
    
    Kokkos::parallel_for(basisCardinality, KOKKOS_LAMBDA(const int basisOrdinal1)
    {
      // we could use hierarchical parallelism to speed this up
      for (int basisOrdinal2=0; basisOrdinal2<basisCardinality; basisOrdinal2++)
      {
        Scalar integral1v1 = 0.0, integral1v2 = 0.0;
        for (int pointOrdinal=0; pointOrdinal<numRefPoints; pointOrdinal++)
        {
          const auto quadratureWeight = weights(pointOrdinal);
          integral1v1 += quadratureWeight * basis1Values(basisOrdinal1,pointOrdinal) * basis1Values(basisOrdinal2,pointOrdinal);
          integral1v2 += quadratureWeight * basis1Values(basisOrdinal1,pointOrdinal) * basis2Values(basisOrdinal2,pointOrdinal);
        }
        basis1_vs_basis1(basisOrdinal1,basisOrdinal2) = integral1v1;
        basis1_vs_basis2(basisOrdinal1,basisOrdinal2) = integral1v2;
      }
    });
    
    // each column in the following matrix will represent the corresponding member of basis 2 in terms of members of basis 1
    ViewType<Scalar> basis1Coefficients("basis 1 coefficients to represent basis 2", basisCardinality, basisCardinality);
    Kokkos::deep_copy(basis1Coefficients, basis1_vs_basis2);
    solveSystemUsingHostLapack(basis1_vs_basis1, basis1Coefficients);
    
    // compute the values for basis2 using basis1Coefficients, basis1Values, and confirm that these agree with basisValues2
    ViewType<Scalar> basis2ValuesFromBasis1 = getOutputView<Scalar>(functionSpace, OPERATOR_VALUE, basisCardinality, numRefPoints, spaceDim);
    deviceGeneralizedMatrixMultiply(basis1Coefficients, true, basis1Values, basis2ValuesFromBasis1); // true: transpose
    
    testViewFloatingEquality(basis2Values, basis2ValuesFromBasis1, relTol, absTol, out, success, "expected", "actual");
    
    for (auto op : opsToTest)
    {
      out << "** Testing operator " << EOperatorToString(op) << " **\n";
      ViewType<Scalar> basis1OpValues = getOutputView<Scalar>(functionSpace, op, basisCardinality, numRefPoints, spaceDim);
      ViewType<Scalar> basis2OpValues = getOutputView<Scalar>(functionSpace, op, basisCardinality, numRefPoints, spaceDim);
      
      basis1.getValues(basis1OpValues, points, op);
      basis2.getValues(basis2OpValues, points, op);
      
      // compute the values for basis2 using basis1Coefficients, basis1Values, and confirm that these agree with basisValues2
      ViewType<Scalar> basis2OpValuesFromBasis1 = getOutputView<Scalar>(functionSpace, op, basisCardinality, numRefPoints, spaceDim);
      deviceGeneralizedMatrixMultiply(basis1Coefficients, true, basis1OpValues, basis2OpValuesFromBasis1); // true: transpose
      
      testViewFloatingEquality(basis2OpValues, basis2OpValuesFromBasis1, relTol, absTol, out, success, "expected", "actual");
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, LineNodalVersusHierarchicalCG_HGRAD )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<>::HGRAD_LINE;
    using NodalBasis        = NodalBasisFamily<>::HGRAD_LINE;
    
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-13; // 2e-14 is sharp on development setup for polyOrder=4; relaxing for potential architectural differences
    const double absTol=1e-13; // 4e-14 is sharp on development setup for polyOrder=4; relaxing for potential architectural differences
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      testBasisEquivalence(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, LineNodalCnVersusNodalC1_HGRAD )
  {
    using CnBasis = Intrepid2::Basis_HGRAD_LINE_Cn_FEM<Kokkos::DefaultExecutionSpace>;
    using C1Basis = Intrepid2::Basis_HGRAD_LINE_C1_FEM<Kokkos::DefaultExecutionSpace>;
    
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-15; // 3e-16 is sharp on development setup for polyOrder=1; relaxing for potential architectural differences
    const double absTol=0.0;   // since there are no quadrature points for polyOrder=1 for which the (analytical) c1Basis evaluates to 0, absTol turns out not to enter into it.
    
    CnBasis cnBasis(1);
    C1Basis c1Basis;
    testBasisEquivalence(cnBasis, c1Basis, opsToTest, relTol, absTol, out, success);
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, LineHierarchicalDGVersusHierarchicalCG_HGRAD )
  {
    using CGBasis = HierarchicalBasisFamily<>::HGRAD_LINE;
    using DGBasis = DGHierarchicalBasisFamily<>::HGRAD_LINE;
    
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-13; // 6e-15 is sharp on development setup for polyOrder=4; relaxing for potential architectural differences
    const double absTol=1e-12; // 1e-13 is sharp on development setup for polyOrder=4; relaxing for potential architectural differences
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      CGBasis cgBasis(polyOrder);
      DGBasis dgBasis(polyOrder);
      testBasisEquivalence(cgBasis, dgBasis, opsToTest, relTol, absTol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, QuadrilateralHierarchicalDGVersusHierarchicalCG_HGRAD )
  {
    using CGBasis = HierarchicalBasisFamily<>::HGRAD_QUAD;
    using DGBasis = DGHierarchicalBasisFamily<>::HGRAD_QUAD;
    
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-12; // 2e-13 is sharp on development setup for polyOrder=3; relaxing for potential architectural differences
    const double absTol=1e-12; // 5e-13 is sharp on development setup for polyOrder=3; relaxing for potential architectural differences
    
    for (int polyOrder=1; polyOrder<4; polyOrder++)
    {
      CGBasis cgBasis(polyOrder);
      DGBasis dgBasis(polyOrder);
      testBasisEquivalence(cgBasis, dgBasis, opsToTest, relTol, absTol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, QuadrilateralNodalVersusHierarchicalCG_HGRAD )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<>::HGRAD_QUAD;
    using NodalBasis        = NodalBasisFamily<>::HGRAD_QUAD;
    
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-13; // 8e-14 is sharp on development setup for polyOrder=3; relaxing for potential architectural differences
    const double absTol=1e-12; // 2e-13 is sharp on development setup for polyOrder=3; relaxing for potential architectural differences
    
    for (int polyOrder=1; polyOrder<4; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      testBasisEquivalence(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, HexahedronHierarchicalDGVersusHierarchicalCG_HGRAD )
  {
    using CGBasis = HierarchicalBasisFamily<>::HGRAD_HEX;
    using DGBasis = DGHierarchicalBasisFamily<>::HGRAD_HEX;
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-11; // 2e-12 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences
    const double absTol=1e-12; // 9e-13 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences
    
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    for (int polyOrder=1; polyOrder<3; polyOrder++)
    {
      CGBasis cgBasis(polyOrder);
      DGBasis dgBasis(polyOrder);
      testBasisEquivalence(cgBasis, dgBasis, opsToTest, relTol, absTol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, HexahedronNodalVersusHierarchicalCG_HGRAD )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<>::HGRAD_HEX;
    using NodalBasis        = NodalBasisFamily<>::HGRAD_HEX;
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-12; // 2e-13 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences
    const double absTol=1e-13; // 2e-14 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences
    
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    for (int polyOrder=1; polyOrder<3; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      testBasisEquivalence(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, HexahedronNodalCnVersusNodalC1_HGRAD )
  {
    using CnBasis = Intrepid2::Basis_HGRAD_HEX_Cn_FEM<Kokkos::DefaultExecutionSpace>;
    using C1Basis = Intrepid2::Basis_HGRAD_HEX_C1_FEM<Kokkos::DefaultExecutionSpace>;

    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-13; // ____ is sharp on development setup for polyOrder=1; relaxing for potential architectural differences
    const double absTol=1e-13; // ____ is sharp on development setup for polyOrder=1; relaxing for potential architectural differences

    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    const int polyOrder = 1;
    CnBasis cnBasis(polyOrder);
    C1Basis c1Basis;
    testBasisEquivalence(cnBasis, c1Basis, opsToTest, relTol, absTol, out, success);
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, HexahedronNodalCnVersusNodalC2_HGRAD )
  {
    using CnBasis = Intrepid2::Basis_HGRAD_HEX_Cn_FEM<Kokkos::DefaultExecutionSpace>;
    using C2Basis = Intrepid2::Basis_HGRAD_HEX_C2_FEM<Kokkos::DefaultExecutionSpace>;

    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-13; // 4e-14 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences
    const double absTol=1e-14; // 2e-15 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences

    // C2 throws an exception for OPERATOR_D5 and OPERATOR_D6, with a message that these are unsupported.
    // I'm not sure why that is, but for that reason we don't test with OPERATOR_D5 here, as we do in other tests
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2, OPERATOR_D3, OPERATOR_D4};
    const int polyOrder = 2;
    CnBasis cnBasis(polyOrder);
    C2Basis c2Basis;
    testBasisEquivalence(cnBasis, c2Basis, opsToTest, relTol, absTol, out, success);
  }
  
} // namespace
