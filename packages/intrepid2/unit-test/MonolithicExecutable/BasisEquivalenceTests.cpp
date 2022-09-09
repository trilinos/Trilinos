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

using namespace Intrepid2;

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
    std::vector<double> B(N*M);
    
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
    lapack.POTRS(UPLO, N, M, A_host.data(), N, B.data(), N, &INFO);
    
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
  template<class Rank2View, class BFunctor, int rankB = 2>
  void deviceGeneralizedMatrixMultiply(Rank2View &A, bool transposeA, BFunctor &B, Rank2View &C)
  {
    using DeviceType = DefaultTestDeviceType;
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
    TEUCHOS_TEST_FOR_EXCEPTION(getFunctorRank(B) != rankB, std::invalid_argument, "run-time rank of B does not match compile-time rank");
    using  ViewIteratorScalar = ViewIterator<Rank2View, Scalar>;
    using     BIteratorScalar = FunctorIterator<BFunctor, Scalar, rankB>;
    Kokkos::parallel_for(N0, KOKKOS_LAMBDA(const int A_row_ordinal)
    {
      BIteratorScalar     B_viewIterator(B);
      ViewIteratorScalar  C_viewIterator(C);
      
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
  //! If ordinalMap is non-empty, it should contain a mapping from the ordinals in basis1 to those ordinals in basis2 that form a sub-basis spanning the same space as basis1.
  template<class DeviceType, class Basis1, class Basis2>
  void testBasisEquivalence(Basis1 &basis1, Basis2 &basis2, typename Basis1::OrdinalTypeArray1D ordinalMap, const std::vector<EOperator> &opsToTest,
                            const double relTol, const double absTol, Teuchos::FancyOStream &out, bool &success)
  {
    if (ordinalMap.size() == 0)
    {
      // first, check that the bases agree on cardinality
      TEST_EQUALITY(basis1.getCardinality(), basis2.getCardinality());
    }
    else
    {
      TEST_EQUALITY(basis1.getCardinality(), ordinalMap.extent_int(0));
    }
    
    const int basisCardinality  = basis1.getCardinality();
    const int basisCardinality2 = basis2.getCardinality();
    
    typename Basis1::OrdinalTypeArray1D reverseOrdinalMap;
    reverseOrdinalMap = typename Basis1::OrdinalTypeArray1D("reverseOrdinalMap - empty", 0);
    auto reverseOrdinalMapHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), reverseOrdinalMap);
    
    if (ordinalMap.size() > 0)
    {
      // then set up reverse ordinal map
      reverseOrdinalMap = typename Basis1::OrdinalTypeArray1D("reverseOrdinalMap", basisCardinality2);
      reverseOrdinalMapHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), reverseOrdinalMap);
      
      auto ordinalMapHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ordinalMap);
      
      // initialize with -1; these are ordinals for members in basis2 that are not included in basis1
      Kokkos::deep_copy(reverseOrdinalMapHost, -1);
      
      for (int i=0; i<basisCardinality; i++)
      {
        int i_mapped = ordinalMapHost(i);
        reverseOrdinalMapHost(i_mapped) = i;
      }
      Kokkos::deep_copy(reverseOrdinalMap, reverseOrdinalMapHost);
      printFunctor1(ordinalMap,        out, "ordinalMap");
      printFunctor1(reverseOrdinalMap, out, "reverseOrdinalMap");
    }
    
    // get quadrature points for integrating up to 2*polyOrder
    const int quadratureDegree = 2*basis1.getDegree();
    using PointScalar = typename Basis1::PointValueType;
    using WeightScalar = typename Basis1::OutputValueType;
    using Scalar = WeightScalar;
    DefaultCubatureFactory cub_factory;
    auto cellTopoKey = basis1.getBaseCellTopology().getKey();
    
    using Cubature       = Intrepid2::Cubature<DeviceType, PointScalar, WeightScalar>;
    using CubatureTensor = Intrepid2::CubatureTensor<DeviceType, PointScalar, WeightScalar>;
    using CubatureDirect = Intrepid2::CubatureDirect<DeviceType, PointScalar, WeightScalar>;
    
    Teuchos::RCP<Cubature> lineTopoQuadrature = cub_factory.create<DeviceType, PointScalar, WeightScalar>(shards::Line<>::key, quadratureDegree);
    Teuchos::RCP<Cubature> baseTopoQuadrature = cub_factory.create<DeviceType, PointScalar, WeightScalar>(cellTopoKey, quadratureDegree);
    
    Teuchos::RCP<Intrepid2::Cubature<DeviceType, PointScalar, WeightScalar>> quadrature;
    
    const int numTensorialExtrusions = basis1.getNumTensorialExtrusions();
    if (numTensorialExtrusions == 0)
    {
      quadrature = baseTopoQuadrature;
    }
    else
    {
      const CubatureDirect* baseTopoQuadratureDirect = dynamic_cast<CubatureDirect*>(baseTopoQuadrature.get());
      const CubatureDirect* lineTopoQuadratureDirect = dynamic_cast<CubatureDirect*>(lineTopoQuadrature.get());
      
      Teuchos::RCP<CubatureTensor> tensorCubature = Teuchos::rcp( new CubatureTensor(*baseTopoQuadratureDirect, *lineTopoQuadratureDirect));
      
      for (int d=1; d<numTensorialExtrusions; d++)
      {
        tensorCubature = Teuchos::rcp( new CubatureTensor(*tensorCubature, *lineTopoQuadratureDirect) );
      }
      
      quadrature = tensorCubature;
    }
    
    ordinal_type numRefPoints = quadrature->getNumPoints();
    const int spaceDim = basis1.getBaseCellTopology().getDimension() + basis1.getNumTensorialExtrusions();
    
    auto points  = quadrature->allocateCubaturePoints();
    auto weights = quadrature->allocateCubatureWeights();
    
    quadrature->getCubature(points, weights);
    
    using HostExecSpace = Kokkos::HostSpace::execution_space;
    TensorPoints<PointScalar,HostExecSpace> pointsHost(points);
    
    out << "Points being tested:\n";
    for (int pointOrdinal=0; pointOrdinal<numRefPoints; pointOrdinal++)
    {
      out << pointOrdinal << ": " << "(";
      for (int d=0; d<spaceDim; d++)
      {
        out << pointsHost(pointOrdinal,d);
        if (d < spaceDim-1) out << ",";
      }
      out << ")\n";
    }
    
    auto functionSpace = basis1.getFunctionSpace();
    
    // set up a projection of basis2 onto basis1
    using Scalar = typename Basis1::OutputValueType;
    auto basis1Values = basis1.allocateBasisValues(points, OPERATOR_VALUE);
    auto basis2Values = basis2.allocateBasisValues(points, OPERATOR_VALUE);
    
    basis1.getValues(basis1Values, points, OPERATOR_VALUE);
    basis2.getValues(basis2Values, points, OPERATOR_VALUE);
    
    basis2Values.setOrdinalFilter(ordinalMap);
    
    // integrate basis1 against itself to compute the SPD matrix A that we'll use to set up the basis conversion system
    ViewType<Scalar,DeviceType> basis1_vs_basis1 = getView<Scalar, DeviceType>("basis 1 vs basis 1", basisCardinality, basisCardinality);
    // integrate basis1 against basis2 to compute the RHS b for the basis conversion system
    ViewType<Scalar,DeviceType> basis1_vs_basis2 = getView<Scalar, DeviceType>("basis 1 vs basis 2", basisCardinality, basisCardinality);
    
    // HDIV and HCURL under OPERATOR_VALUE are vector-valued; "scalarValued" tests for this case.
    const bool scalarValued = (functionSpace == FUNCTION_SPACE_HGRAD) || (functionSpace == FUNCTION_SPACE_HVOL);
    
    if (scalarValued)
    {
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
    }
    else
    {
      int spaceDim = basis1Values.extent_int(2);
      Kokkos::parallel_for(basisCardinality, KOKKOS_LAMBDA(const int basisOrdinal1)
      {
        // we could use hierarchical parallelism to speed this up
        for (int basisOrdinal2=0; basisOrdinal2<basisCardinality; basisOrdinal2++)
        {
          Scalar integral1v1 = 0.0, integral1v2 = 0.0;
          for (int pointOrdinal=0; pointOrdinal<numRefPoints; pointOrdinal++)
          {
            const auto quadratureWeight = weights(pointOrdinal);
            for (int d=0; d<spaceDim; d++)
            {
              integral1v1 += quadratureWeight * basis1Values(basisOrdinal1,pointOrdinal,d) * basis1Values(basisOrdinal2,pointOrdinal,d);
              integral1v2 += quadratureWeight * basis1Values(basisOrdinal1,pointOrdinal,d) * basis2Values(basisOrdinal2,pointOrdinal,d);
            }
          }
          basis1_vs_basis1(basisOrdinal1,basisOrdinal2) = integral1v1;
          basis1_vs_basis2(basisOrdinal1,basisOrdinal2) = integral1v2;
        }
      });
    }
    
    // each column in the following matrix will represent the corresponding member of basis 2 in terms of members of basis 1
    ViewType<Scalar,DeviceType> basis1Coefficients = getView<Scalar, DeviceType>("basis 1 coefficients", basisCardinality, basisCardinality);
    Kokkos::deep_copy(basis1Coefficients, basis1_vs_basis2);
    solveSystemUsingHostLapack(basis1_vs_basis1, basis1Coefficients);
    
    // for non-"DG" bases, check that the topological associations match up
    // to check for DG-ness of basis, examine how many dofs are associated with the interior
    auto cellTopo = ::Intrepid2::CellTopology::cellTopology(basis1.getBaseCellTopology(), basis1.getNumTensorialExtrusions());
    const int interiorDim = cellTopo->getDimension();
    const int interiorSubcellOrdinal = 0;
    const int firstDofOrdinalForSubcell = 0;
    
    int basis1InteriorCardinality = 0, basis2InteriorCardinality = 0;
    int basis1FirstInterior = -1, basis2FirstInterior = -1;
    if (basis1.getAllDofOrdinal().extent_int(0) > interiorDim)
    {
      basis1FirstInterior = basis1.getAllDofOrdinal()(interiorDim, interiorSubcellOrdinal, firstDofOrdinalForSubcell);
    }
    if (basis2.getAllDofOrdinal().extent_int(0) > interiorDim)
    {
      basis2FirstInterior = basis2.getAllDofOrdinal()(interiorDim, interiorSubcellOrdinal, firstDofOrdinalForSubcell);
    }
    
    // if there are no interior dofs, we'll get a -1 value
    if (basis1FirstInterior != -1)
    {
      // at index 3, dof tag stores the total dof count associated with the subcell (here, the interior)
      basis1InteriorCardinality = basis1.getDofTag(basis1FirstInterior)(3);
    }
    if (basis2FirstInterior != -1)
    {
      // at index 3, dof tag stores the total dof count associated with the subcell (here, the interior)
      basis2InteriorCardinality = basis2.getDofTag(basis2FirstInterior)(3);
    }
    const bool basis1IsDG = (basis1InteriorCardinality == basisCardinality);
    const bool basis2IsDG = (basis2InteriorCardinality == basisCardinality2);
    const bool neitherBasisIsDG = !basis1IsDG && !basis2IsDG;
    
    // Hypercube bases can be defined on a strictly tensorial topology, as opposed to the shards topology.  These have different subcell numbering; and the
    // block below relies on the subcell numbering being the same between the two.  To do the equivalent thing with differing topologies, we would need a
    // mapping from one cell topology to the other.  We have something like this in TensorTopologyMap, so it might not be too hard to add this…
    const bool basesAgreeOnCellTopo = (basis1.getBaseCellTopology() == basis2.getBaseCellTopology()) && (basis1.getNumTensorialExtrusions() == basis2.getNumTensorialExtrusions());
    if (neitherBasisIsDG && basesAgreeOnCellTopo)
    {
      auto basis1CoefficientsHost = getHostCopy(basis1Coefficients);
      // if neither basis is DG, then we can expect them to have matching counts on basis members
      // associated with a given subcell.  We can also expect the representation of of members of basis1
      // on a given subcell in terms of basis2 to involve at least some members of basis2 associated with the same
      // subcell.  (We can check this by examining the weights in basis1Coefficients.)
      for (int subcellDim=0; subcellDim<=interiorDim; subcellDim++)
      {
        out << "checking subcells of dimension " << subcellDim << std::endl;
        // if there are no dofs for this subcell dimension and greater, then getAllDofOrdinal() won't have
        // an entry for subcellDim.  The following guards against that:
        if (basis1.getAllDofOrdinal().extent_int(0) <= subcellDim)
        {
          // basis1 has no dofs for this subcell dim.  Check that basis2 also does not, at least among the members that form the sub-basis corresponding to basis1:
          if (reverseOrdinalMapHost.size() == 0) // basis2 and basis1 are supposed to be equivalent
          {
            TEST_ASSERT(basis2.getAllDofOrdinal().extent_int(0) <= subcellDim);
          }
          else // basis1 corresponds to a non-trivial sub-basis of basis2
          {
            if (basis2.getAllDofOrdinal().extent_int(0) > subcellDim)
            {
              const int subcellCount = cellTopo->getSubcellCount(subcellDim);
              for (int subcellOrdinal=0; subcellOrdinal<subcellCount; subcellOrdinal++)
              {
                const int basis2FirstDofOrdinal = basis2.getAllDofOrdinal()(subcellDim, subcellOrdinal, firstDofOrdinalForSubcell);
                const int basis2SubcellCardinality = basis2.getDofTag(basis2FirstDofOrdinal)(3);
                for (int subcellDofOrdinal=0; subcellDofOrdinal<basis2SubcellCardinality; subcellDofOrdinal++)
                {
                  const int basis2DofOrdinal = basis2.getAllDofOrdinal()(subcellDim, subcellOrdinal, subcellDofOrdinal);
                  TEST_ASSERT(reverseOrdinalMapHost(basis2DofOrdinal) == -1); // indicates that basis2DofOrdinal corresponds to basis 2 member that lies outside the span of basis 1
                }
              }
            }
          }
          break; // we've already checked all subcell dimensions that have any dofs associated with them
        }
        
        const int subcellCount = cellTopo->getSubcellCount(subcellDim);
        auto basis1AllDofOrdinal = basis1.getAllDofOrdinal();
        auto basis2AllDofOrdinal = basis2.getAllDofOrdinal();
        for (int subcellOrdinal=0; subcellOrdinal<subcellCount; subcellOrdinal++)
        {
          // need to find the first dof ordinal for the subcell to get the basis cardinality on the subcell
          int basis1FirstDofOrdinal, basis2FirstDofOrdinal;
          if ((subcellDim < basis1AllDofOrdinal.extent_int(0)) && (subcellOrdinal < basis1AllDofOrdinal.extent_int(1)))
          {
            basis1FirstDofOrdinal = basis1AllDofOrdinal(subcellDim, subcellOrdinal, firstDofOrdinalForSubcell);
          }
          else
          {
            basis1FirstDofOrdinal = -1;
          }
          if ((subcellDim < basis2AllDofOrdinal.extent_int(0)) && (subcellOrdinal < basis2AllDofOrdinal.extent_int(1)))
          {
            basis2FirstDofOrdinal = basis2AllDofOrdinal(subcellDim, subcellOrdinal, firstDofOrdinalForSubcell);
          }
          else
          {
            basis2FirstDofOrdinal = -1;
          }
          // if there are no dofs on the subcell, we'll get a -1 value
          if ((basis1FirstDofOrdinal == -1) || (basis2FirstDofOrdinal == -1))
          {
            // if one of the bases has no dofs on the subcell, then both should:
            TEST_ASSERT((basis1FirstDofOrdinal == -1) && (basis2FirstDofOrdinal == -1));
            // if either of the bases has no dofs on the subcell, we should continue with the next subcell
            continue;
          }
          // at index 3, dof tag stores the total dof count associated with the subcell
          const int basis1SubcellCardinality = basis1.getDofTag(basis1FirstDofOrdinal)(3);
          const int basis2SubcellCardinality = basis2.getDofTag(basis2FirstDofOrdinal)(3);
          
          out << "subcell " << subcellOrdinal << " has " << basis1SubcellCardinality << " dofs.\n";
          
          std::vector<ordinal_type> basis1SubcellDofOrdinals(basis1SubcellCardinality);
          std::vector<ordinal_type> basis2SubcellDofOrdinals(basis2SubcellCardinality);
          int basis2SubcellCardinalityFiltered = 0;
          for (int subcellDofOrdinal=0; subcellDofOrdinal<basis1SubcellCardinality; subcellDofOrdinal++)
          {
            basis1SubcellDofOrdinals[subcellDofOrdinal] = basis1.getAllDofOrdinal()(subcellDim, subcellOrdinal, subcellDofOrdinal);
          }
          for (int subcellDofOrdinal=0; subcellDofOrdinal<basis2SubcellCardinality; subcellDofOrdinal++)
          {
            int basis2SubcellDofOrdinal                 = basis2.getAllDofOrdinal()(subcellDim, subcellOrdinal, subcellDofOrdinal);
            basis2SubcellDofOrdinals[subcellDofOrdinal] = basis2SubcellDofOrdinal;
            
            if ((reverseOrdinalMapHost.size() == 0) || (reverseOrdinalMapHost(basis2SubcellDofOrdinal) != -1))
            {
              basis2SubcellCardinalityFiltered++;
            }
          }
          
          TEST_EQUALITY(basis1SubcellCardinality, basis2SubcellCardinalityFiltered);
          // if we fail the cardinality check, not much point in checking the coefficients
          if (basis1SubcellCardinality != basis2SubcellCardinalityFiltered) continue;
          
          for (int subcellDofOrdinal=0; subcellDofOrdinal<basis2SubcellCardinality; subcellDofOrdinal++)
          {
            // each column in basis1Coefficients represents the corresponding member of basis 2 in terms of members of basis 1
            const int basis2DofOrdinal = basis2SubcellDofOrdinals[subcellDofOrdinal];
            out << "checking representation of basis2 dof ordinal " << basis2DofOrdinal << std::endl;
            
            int basis2DofOrdinalFiltered;
            if (reverseOrdinalMapHost.size() == 0)
            {
              basis2DofOrdinalFiltered = basis2DofOrdinal;
            }
            else
            {
              basis2DofOrdinalFiltered = reverseOrdinalMapHost(basis2DofOrdinal);
              if (basis2DofOrdinalFiltered == -1)
              {
                // indicates a basis 2 member outside the span of basis 1.
                continue;
              }
            }
            // look for at least one member of basis1's representation on the subcell
            bool hasNonzeroCoefficient = false;
            for (const int basis1DofOrdinal : basis1SubcellDofOrdinals)
            {
              Scalar basisCoefficient = basis1CoefficientsHost(basis1DofOrdinal,basis2DofOrdinalFiltered);
              bool nonzeroCoefficient = (std::abs(basisCoefficient) > absTol);
              if (nonzeroCoefficient)
              {
                out << "basis coefficient " << basisCoefficient << " is nonzero (absTol = " << absTol << ").\n";
                hasNonzeroCoefficient = true;
              }
            }
            TEST_ASSERT(hasNonzeroCoefficient);
          }
        }
      }
    }
    
    // compute the values for basis2 using basis1Coefficients, basis1Values, and confirm that these agree with basisValues2
    auto basis2ValuesFromBasis1 = getOutputView<Scalar,DeviceType>(functionSpace, OPERATOR_VALUE, basisCardinality, numRefPoints, spaceDim);
    {
      const int rankB = basis1Values.rank();
      using ViewType = decltype(basis1Coefficients);
      using BType    = decltype(basis1Values);
      switch(rankB)
      {
        case 2: deviceGeneralizedMatrixMultiply<ViewType,BType,2>(basis1Coefficients, true, basis1Values, basis2ValuesFromBasis1); break;
        case 3: deviceGeneralizedMatrixMultiply<ViewType,BType,3>(basis1Coefficients, true, basis1Values, basis2ValuesFromBasis1); break;
        default:
          TEUCHOS_TEST_FOR_EXCEPTION((rankB < 2) || (rankB > 3), std::invalid_argument, "Unhandled rank for basis1Values");
      }
    }
    
    testViewFloatingEquality(basis2ValuesFromBasis1, basis2Values, relTol, absTol, out, success, "actual", "expected");
    
    for (auto op : opsToTest)
    {
      out << "** Testing operator " << EOperatorToString(op) << " **\n";
      auto basis1OpValues = basis1.allocateBasisValues(points, op);
      auto basis2OpValues = basis2.allocateBasisValues(points, op);
      
      basis1.getValues(basis1OpValues, points, op);
      basis2.getValues(basis2OpValues, points, op);
      
      // compute the values for basis2 using basis1Coefficients, basis1Values, and confirm that these agree with basisValues2
      auto basis2OpValuesFromBasis1 = getOutputView<Scalar,DeviceType>(functionSpace, op, basisCardinality, numRefPoints, spaceDim);
      int rankB = basis1OpValues.rank();
      using ViewType = decltype(basis1Coefficients);
      using BType    = decltype(basis1OpValues);
      const bool transpose = true;
      switch(rankB)
      {
        case 2: deviceGeneralizedMatrixMultiply<ViewType,BType,2>(basis1Coefficients, transpose, basis1OpValues, basis2OpValuesFromBasis1); break;
        case 3: deviceGeneralizedMatrixMultiply<ViewType,BType,3>(basis1Coefficients, transpose, basis1OpValues, basis2OpValuesFromBasis1); break;
        default:
          TEUCHOS_TEST_FOR_EXCEPTION((rankB < 2) || (rankB > 3), std::invalid_argument, "Unhandled rank for basis1OpValues");
      }
      
      basis2OpValues.setOrdinalFilter(ordinalMap);
      testViewFloatingEquality(basis2OpValuesFromBasis1, basis2OpValues, relTol, absTol, out, success, "actual", "expected");
    }
  }

  //! tests that two bases are equivalent; computes a conversion from one to the other and then uses that to confirm that the equivalence
  //! holds for OPERATOR_VALUE (as it should by construction), as well as the operators passed in in opsToTest.
  template<class DeviceType, class Basis1, class Basis2>
  void testBasisEquivalence(Basis1 &basis1, Basis2 &basis2, const std::vector<EOperator> &opsToTest,
                            const double relTol, const double absTol, Teuchos::FancyOStream &out, bool &success)
  {
    typename Basis1::OrdinalTypeArray1D ordinalMap; // empty map
    testBasisEquivalence<DeviceType,Basis1,Basis2>(basis1, basis2, ordinalMap, opsToTest, relTol, absTol, out, success);
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, LineNodalVersusHierarchicalCG_HGRAD )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_LINE;
    using NodalBasis        = NodalBasisFamily<DefaultTestDeviceType>::HGRAD_LINE;
    
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-13; // 2e-14 is sharp on development setup for polyOrder=4; relaxing for potential architectural differences
    const double absTol=1e-13; // 4e-14 is sharp on development setup for polyOrder=4; relaxing for potential architectural differences
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, LineNodalCnVersusNodalC1_HGRAD )
  {
    using CnBasis = Intrepid2::Basis_HGRAD_LINE_Cn_FEM<DefaultTestDeviceType>;
    using C1Basis = Intrepid2::Basis_HGRAD_LINE_C1_FEM<DefaultTestDeviceType>;
    
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-15; // 3e-16 is sharp on development setup for polyOrder=1; relaxing for potential architectural differences
    const double absTol=0.0;   // since there are no quadrature points for polyOrder=1 for which the (analytical) c1Basis evaluates to 0, absTol turns out not to enter into it.
    
    CnBasis cnBasis(1);
    C1Basis c1Basis;
    testBasisEquivalence<DefaultTestDeviceType>(cnBasis, c1Basis, opsToTest, relTol, absTol, out, success);
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, LineHierarchicalDGVersusHierarchicalCG_HGRAD )
  {
    using CGBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_LINE;
    using DGBasis = DGHierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_LINE;
    
    // NOTE: for the moment, OPERATOR_Dn for n > 2 not supported by DerivedBasis.  We can support more by either increasing
    //       Parameters::MaxVectorComponents (which is 7 right now), or by changing VectorData to allow a dynamic number of
    //       components.  (We were doing the latter using Kokkos::vector, but have switched to a Kokkos::Array instead to
    //       avoid using UVM.)
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2}; //, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-13; // 6e-15 is sharp on development setup for polyOrder=4; relaxing for potential architectural differences
    const double absTol=1e-12; // 1e-13 is sharp on development setup for polyOrder=4; relaxing for potential architectural differences
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      CGBasis cgBasis(polyOrder);
      DGBasis dgBasis(polyOrder);
      testBasisEquivalence<DefaultTestDeviceType>(cgBasis, dgBasis, opsToTest, relTol, absTol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, QuadrilateralHierarchicalDGVersusHierarchicalCG_HGRAD )
  {
    using CGBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_QUAD;
    using DGBasis = DGHierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_QUAD;
    
    // NOTE: for the moment, OPERATOR_Dn for n > 2 not supported by DerivedBasis.  We can support more by either increasing
    //       Parameters::MaxVectorComponents (which is 7 right now), or by changing VectorData to allow a dynamic number of
    //       components.  (We were doing the latter using Kokkos::vector, but have switched to a Kokkos::Array instead to
    //       avoid using UVM.)
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2}; //, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-12; // 2e-13 is sharp on development setup for polyOrder=3; relaxing for potential architectural differences
    const double absTol=1e-11; // 5e-13 is sharp on development setup for polyOrder=3; relaxing for potential architectural differences
    
    for (int polyOrder=1; polyOrder<4; polyOrder++)
    {
      CGBasis cgBasis(polyOrder);
      DGBasis dgBasis(polyOrder);
      testBasisEquivalence<DefaultTestDeviceType>(cgBasis, dgBasis, opsToTest, relTol, absTol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, QuadrilateralNodalVersusHierarchicalCG_HGRAD )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_QUAD;
    using NodalBasis        = NodalBasisFamily<DefaultTestDeviceType>::HGRAD_QUAD;
    
    // NOTE: for the moment, OPERATOR_Dn for n > 2 not supported by DerivedBasis.  We can support more by either increasing
    //       Parameters::MaxVectorComponents (which is 7 right now), or by changing VectorData to allow a dynamic number of
    //       components.  (We were doing the latter using Kokkos::vector, but have switched to a Kokkos::Array instead to
    //       avoid using UVM.)
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2}; //, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-13; // 8e-14 is sharp on development setup for polyOrder=3; relaxing for potential architectural differences
    const double absTol=1e-12; // 2e-13 is sharp on development setup for polyOrder=3; relaxing for potential architectural differences
    
    for (int polyOrder=1; polyOrder<4; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, HexahedronHierarchicalDGVersusHierarchicalCG_HGRAD )
  {
    using CGBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_HEX;
    using DGBasis = DGHierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_HEX;
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-11; // 2e-12 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences
    const double absTol=1e-11; // 9e-13 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences
    
    // NOTE: for the moment, OPERATOR_Dn for n > 2 not supported by DerivedBasis.  We can support more by either increasing
    //       Parameters::MaxVectorComponents (which is 7 right now), or by changing VectorData to allow a dynamic number of
    //       components.  (We were doing the latter using Kokkos::vector, but have switched to a Kokkos::Array instead to
    //       avoid using UVM.)
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2}; //, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    for (int polyOrder=1; polyOrder<3; polyOrder++)
    {
      CGBasis cgBasis(polyOrder);
      DGBasis dgBasis(polyOrder);
      testBasisEquivalence<DefaultTestDeviceType>(cgBasis, dgBasis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, HexahedronHierarchicalCGVersusHypercube3D_HGRAD )
  {
    using CGBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_HEX;
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-11;
    const double absTol=1e-11;
    
    // NOTE: for the moment, OPERATOR_Dn for n > 2 not supported by DerivedBasis.  We can support more by either increasing
    //       Parameters::MaxVectorComponents (which is 7 right now), or by changing VectorData to allow a dynamic number of
    //       components.  (We were doing the latter using Kokkos::vector, but have switched to a Kokkos::Array instead to
    //       avoid using UVM.)
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2}; //, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    for (int polyOrder=1; polyOrder<3; polyOrder++)
    {
      CGBasis cgBasis(polyOrder);
      const int spaceDim = 3;
      auto hypercubeBasis = getHypercubeBasis_HGRAD<HierarchicalBasisFamily<DefaultTestDeviceType>>(polyOrder, spaceDim);
      testBasisEquivalence<DefaultTestDeviceType>(cgBasis, *hypercubeBasis, opsToTest, relTol, absTol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, HexahedronNodalVersusHierarchicalCG_HGRAD )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_HEX;
    using NodalBasis        = NodalBasisFamily<DefaultTestDeviceType>::HGRAD_HEX;
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-12; // 2e-13 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences
    const double absTol=1e-13; // 2e-14 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences
    
    // NOTE: for the moment, OPERATOR_Dn for n > 2 not supported by DerivedBasis.  We can support more by either increasing
    //       Parameters::MaxVectorComponents (which is 7 right now), or by changing VectorData to allow a dynamic number of
    //       components.  (We were doing the latter using Kokkos::vector, but have switched to a Kokkos::Array instead to
    //       avoid using UVM.)
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2}; //, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    for (int polyOrder=1; polyOrder<3; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, HexahedronNodalCnVersusNodalC1_HGRAD )
  {
    using CnBasis = Intrepid2::Basis_HGRAD_HEX_Cn_FEM<DefaultTestDeviceType>;
    using C1Basis = Intrepid2::Basis_HGRAD_HEX_C1_FEM<DefaultTestDeviceType>;

    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-13; // ____ is sharp on development setup for polyOrder=1; relaxing for potential architectural differences
    const double absTol=1e-13; // ____ is sharp on development setup for polyOrder=1; relaxing for potential architectural differences

    // NOTE: for the moment, OPERATOR_Dn for n > 2 on Hexahedron not supported by BasisValues.  We can support more by either increasing
    //       Parameters::MaxVectorComponents (which is 7 right now), or by changing VectorData to allow a dynamic number of
    //       components.  (We were doing the latter using Kokkos::vector, but have switched to a Kokkos::Array instead to
    //       avoid using UVM.)
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2}; //, OPERATOR_D3, OPERATOR_D4, OPERATOR_D5};
    const int polyOrder = 1;
    CnBasis cnBasis(polyOrder);
    C1Basis c1Basis;
    testBasisEquivalence<DefaultTestDeviceType>(cnBasis, c1Basis, opsToTest, relTol, absTol, out, success);
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, HexahedronNodalCnVersusNodalC2_HGRAD )
  {
    using CnBasis = Intrepid2::Basis_HGRAD_HEX_Cn_FEM<DefaultTestDeviceType>;
    using C2Basis = Intrepid2::Basis_HGRAD_HEX_C2_FEM<DefaultTestDeviceType>;

    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-13; // 4e-14 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences
    const double absTol=1e-14; // 2e-15 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences

    // C2 throws an exception for OPERATOR_D5 and OPERATOR_D6, with a message that these are unsupported.
    // I'm not sure why that is, but for that reason we don't test with OPERATOR_D5 here, as we do in other tests
    // NOTE: for the moment, OPERATOR_Dn for n > 2 on Hexahedron not supported by BasisValues.  We can support more by either increasing
    //       Parameters::MaxVectorComponents (which is 7 right now), or by changing VectorData to allow a dynamic number of
    //       components.  (We were doing the latter using Kokkos::vector, but have switched to a Kokkos::Array instead to
    //       avoid using UVM.)
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1, OPERATOR_D2}; //, OPERATOR_D3, OPERATOR_D4};
    const int polyOrder = 2;
    CnBasis cnBasis(polyOrder);
    C2Basis c2Basis;
    testBasisEquivalence<DefaultTestDeviceType>(cnBasis, c2Basis, opsToTest, relTol, absTol, out, success);
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, HypercubeNodalVersusHypercubeHierarchical_HGRAD )
  {
    const int maxSpaceDim = 7; // we only test polyDegree = 1 for spaceDim > 5, due to test performance considerations
    const int maxDegree = 2;
    
    using HierarchicalBasisFamily = HierarchicalBasisFamily<DefaultTestDeviceType>;
    using NodalBasisFamily = NodalBasisFamily<DefaultTestDeviceType>;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_GRAD};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-14;
    const double absTol=1e-14;
    
    Intrepid2::EFunctionSpace fs = FUNCTION_SPACE_HGRAD;
    
    const int minDegree = (fs == FUNCTION_SPACE_HVOL) ? 0 : 1;
    for (int polyDegree = minDegree; polyDegree <= maxDegree; polyDegree++)
    {
      out << "** polyDegree " << polyDegree << " **\n";
      for (int spaceDim = 1; spaceDim <= maxSpaceDim; spaceDim++)
      {
        if ((polyDegree > 1) && (spaceDim > 5)) continue; // skip this case in the interest of test performance
        out << "** spaceDim " << spaceDim << " **\n";
        auto hierarchicalBasis = getHypercubeBasis_HGRAD<HierarchicalBasisFamily>(polyDegree, spaceDim);
        auto nodalBasis        = getHypercubeBasis_HGRAD<NodalBasisFamily>(polyDegree, spaceDim);
        testBasisEquivalence<DefaultTestDeviceType>(*nodalBasis, *hierarchicalBasis, opsToTest, relTol, absTol, out, success);
      }
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, HypercubeLowestOrderVersusSerendipity_HGRAD )
  {
    const int maxSpaceDim = 7; // lowest-order hypercube basis should be identitical to its serendipity basis
    const int minDegree = 1;
    const int maxDegree = 1;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_GRAD};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-14;
    const double absTol=1e-14;
    
    using BasisFamily = HierarchicalBasisFamily<DefaultTestDeviceType>;
    using BasisBase = typename BasisFamily::HGRAD_LINE::BasisBase;
    
    for (int polyDegree = minDegree; polyDegree <= maxDegree; polyDegree++)
    {
      out << "** polyDegree " << polyDegree << " **\n";
      for (int spaceDim = 1; spaceDim <= maxSpaceDim; spaceDim++)
      {
        out << "** spaceDim " << spaceDim << " **\n";
        auto hierarchicalBasis = getHypercubeBasis_HGRAD<BasisFamily>(polyDegree, spaceDim);
        auto serendipityBasis = Teuchos::rcp(new SerendipityBasis<BasisBase>(hierarchicalBasis));
        
        testBasisEquivalence<DefaultTestDeviceType>(*hierarchicalBasis, *serendipityBasis, opsToTest, relTol, absTol, out, success);
      }
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, HypercubeHigherOrderVersusSerendipity_HGRAD )
  {
    const int maxSpaceDim = 4;
    const int minDegree = 2;
    const int maxDegree = 4;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_GRAD};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    using BasisFamily = HierarchicalBasisFamily<DefaultTestDeviceType>;
    using BasisBase = typename BasisFamily::HGRAD_LINE::BasisBase;
    
    for (int polyDegree = minDegree; polyDegree <= maxDegree; polyDegree++)
    {
      out << "** polyDegree " << polyDegree << " **\n";
      for (int spaceDim = 1; spaceDim <= maxSpaceDim; spaceDim++)
      {
        out << "** spaceDim " << spaceDim << " **\n";
        auto hierarchicalBasis = getHypercubeBasis_HGRAD<BasisFamily>(polyDegree, spaceDim);
        auto serendipityBasis = Teuchos::rcp(new SerendipityBasis<BasisBase>(hierarchicalBasis));
        
        testBasisEquivalence<DefaultTestDeviceType>(*serendipityBasis, *hierarchicalBasis, serendipityBasis->ordinalMap(), opsToTest, relTol, absTol, out, success);
      }
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, HypercubeHigherOrderVersusSerendipity_HVOL )
  {
    const int maxSpaceDim = 4;
    const int minDegree = 2;
    const int maxDegree = 4;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE};
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    using BasisFamily = HierarchicalBasisFamily<DefaultTestDeviceType>;
    using BasisBase = typename BasisFamily::HGRAD_LINE::BasisBase;
    
    for (int polyDegree = minDegree; polyDegree <= maxDegree; polyDegree++)
    {
      out << "** polyDegree " << polyDegree << " **\n";
      for (int spaceDim = 1; spaceDim <= maxSpaceDim; spaceDim++)
      {
        out << "** spaceDim " << spaceDim << " **\n";
        auto hierarchicalBasis = getHypercubeBasis_HVOL<BasisFamily>(polyDegree, spaceDim);
        auto serendipityBasis = Teuchos::rcp(new SerendipityBasis<BasisBase>(hierarchicalBasis));
        
        testBasisEquivalence<DefaultTestDeviceType>(*serendipityBasis, *hierarchicalBasis, serendipityBasis->ordinalMap(), opsToTest, relTol, absTol, out, success);
      }
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, HypercubeNodalVersusHypercubeHierarchical_HVOL )
  {
    const int maxSpaceDim = 7; // we only test polyDegree = 0,1 for spaceDim > 5, due to test performance considerations
    const int maxDegreeForCardinalityTests = 2;
    
    using HierarchicalBasisFamily = HierarchicalBasisFamily<DefaultTestDeviceType>;
    using NodalBasisFamily = NodalBasisFamily<DefaultTestDeviceType>;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher-order bases in higher dimensions)
    const double relTol=1e-10;
    const double absTol=1e-10;
    
    Intrepid2::EFunctionSpace fs = FUNCTION_SPACE_HVOL;
    
    const int minDegree = (fs == FUNCTION_SPACE_HVOL) ? 0 : 1;
    for (int polyDegree = minDegree; polyDegree <= maxDegreeForCardinalityTests; polyDegree++)
    {
      out << "** polyDegree " << polyDegree << " **\n";
      for (int spaceDim = 1; spaceDim <= maxSpaceDim; spaceDim++)
      {
        if ((polyDegree > 1) && (spaceDim > 5)) continue; // skip this case in the interest of test performance
        out << "** spaceDim " << spaceDim << " **\n";
        auto hierarchicalBasis = getHypercubeBasis_HVOL<HierarchicalBasisFamily>(polyDegree, spaceDim);
        auto nodalBasis        = getHypercubeBasis_HVOL<NodalBasisFamily>(polyDegree, spaceDim);
        testBasisEquivalence<DefaultTestDeviceType>(*nodalBasis, *hierarchicalBasis, opsToTest, relTol, absTol, out, success);
      }
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, TetrahedronNodalCnVersusNodalC2_HGRAD )
  {
    using CnBasis = Intrepid2::Basis_HGRAD_TET_Cn_FEM<DefaultTestDeviceType>;
    using C2Basis = Intrepid2::Basis_HGRAD_TET_C2_FEM<DefaultTestDeviceType>;
    
    // OPERATOR_D2 and above are not supported by either the nodal or the hierarchical basis at present...
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-12; // 2e-14 is sharp on development setup; relaxing for potential architectural differences
    const double absTol=1e-13; // 3e-15 is sharp on development setup; relaxing for potential architectural differences
    
    CnBasis cnBasis(2);
    C2Basis c2Basis;
    testBasisEquivalence<DefaultTestDeviceType>(cnBasis, c2Basis, opsToTest, relTol, absTol, out, success);
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, TetrahedronHierarchicalDGVersusHierarchicalCG_HGRAD )
  {
    using CGBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_TET;
    using DGBasis = DGHierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_TET;
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-6; // 3e-8 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences
    const double absTol=1e-9; // 5e-11 is sharp on development setup for polyOrder=2; relaxing for potential architectural differences
    
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1};
    for (int polyOrder=1; polyOrder<7; polyOrder++)
    {
      CGBasis cgBasis(polyOrder);
      DGBasis dgBasis(polyOrder);
      testBasisEquivalence<DefaultTestDeviceType>(cgBasis, dgBasis, opsToTest, relTol, absTol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, TetrahedronNodalVersusHierarchicalCG_HGRAD )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_TET;
    using NodalBasis        = NodalBasisFamily<DefaultTestDeviceType>::HGRAD_TET;
    
    // OPERATOR_D2 and above are not supported by either the nodal or the hierarchical basis at present...
    std::vector<EOperator> opsToTest {OPERATOR_GRAD, OPERATOR_D1};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-6;  // 3e-08 is sharp on development setup for polyOrder=6; relaxing for potential architectural differences
    const double absTol=1e-10; // 3e-12 is sharp on development setup for polyOrder=6; relaxing for potential architectural differences
    
    for (int polyOrder=1; polyOrder<7; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      
      testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, TetrahedronNodalVersusHierarchical_HCURL )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HCURL_TET;
    using NodalBasis        = NodalBasisFamily<DefaultTestDeviceType>::HCURL_TET;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_CURL};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-6;
    const double absTol=1e-10;
    
    for (int polyOrder=1; polyOrder<7; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      
      testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, TetrahedronNodalVersusHierarchical_HDIV )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HDIV_TET;
    using NodalBasis        = NodalBasisFamily<DefaultTestDeviceType>::HDIV_TET;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_DIV};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-6;
    const double absTol=1e-9;
    
    for (int polyOrder=1; polyOrder<7; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      
      testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, TetrahedronNodalVersusHierarchical_HVOL )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HVOL_TET;
    using NodalBasis        = NodalBasisFamily<DefaultTestDeviceType>::HVOL_TET;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-6;
    const double absTol=1e-10;
    
    for (int polyOrder=1; polyOrder<7; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      
      testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, TriangleNodalVersusHierarchicalCG_HGRAD )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_TRI;
    using NodalBasis        = NodalBasisFamily<DefaultTestDeviceType>::HGRAD_TRI;
    
    // OPERATOR_D2 and above are not supported by either the nodal or the hierarchical basis at present...
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_GRAD, OPERATOR_D1};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-11; // 7e-13 is sharp on development setup for polyOrder=4; relaxing for potential architectural differences
    const double absTol=1e-12; // 2e-14 is sharp on development setup for polyOrder=4; relaxing for potential architectural differences
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisEquivalence, TriangleHierarchicalCGVersusHierarchicalDG_HGRAD )
  {
    using CGBasis =   HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_TRI;
    using DGBasis = DGHierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_TRI;
    
    // OPERATOR_D2 and above are not supported by either the nodal or the hierarchical basis at present...
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_GRAD, OPERATOR_D1};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-11; // 7e-13 is sharp on development setup for polyOrder=4; relaxing for potential architectural differences
    const double absTol=1e-12; // 2e-14 is sharp on development setup for polyOrder=4; relaxing for potential architectural differences
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      CGBasis cgBasis(polyOrder);
      DGBasis dgBasis(polyOrder);
      testBasisEquivalence<DefaultTestDeviceType>(cgBasis, dgBasis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, TriangleNodalVersusHierarchical_HCURL )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HCURL_TRI;
    using NodalBasis        = NodalBasisFamily<DefaultTestDeviceType>::HCURL_TRI;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_CURL};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-10;
    const double absTol=1e-10;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      testBasisEquivalence<DefaultTestDeviceType>(hierarchicalBasis, nodalBasis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, TriangleNodalVersusHierarchical_HDIV )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HDIV_TRI;
    using NodalBasis        = NodalBasisFamily<DefaultTestDeviceType>::HDIV_TRI;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE, OPERATOR_DIV};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-10;
    const double absTol=1e-10;
    
    for (int polyOrder=1; polyOrder<5; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      testBasisEquivalence<DefaultTestDeviceType>(hierarchicalBasis, nodalBasis, opsToTest, relTol, absTol, out, success);
    }
  }

  TEUCHOS_UNIT_TEST( BasisEquivalence, TriangleNodalVersusHierarchical_HVOL )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HVOL_TRI;
    using NodalBasis        = NodalBasisFamily<DefaultTestDeviceType>::HVOL_TRI;
    
    std::vector<EOperator> opsToTest {OPERATOR_VALUE};
    
    // these tolerances are selected such that we have a little leeway for architectural differences
    // (It is true, though, that we incur a fair amount of floating point error for higher order bases in higher dimensions)
    const double relTol=1e-11;
    const double absTol=1e-12;
    
    for (int polyOrder=0; polyOrder<5; polyOrder++)
    {
      HierarchicalBasis hierarchicalBasis(polyOrder);
      NodalBasis        nodalBasis(polyOrder);
      testBasisEquivalence<DefaultTestDeviceType>(nodalBasis, hierarchicalBasis, opsToTest, relTol, absTol, out, success);
    }
  }
  
} // namespace
