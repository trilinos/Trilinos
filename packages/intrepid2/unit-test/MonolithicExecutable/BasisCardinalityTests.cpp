// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   BasisCardinalityTests.cpp
    \brief  Tests to verify that basis implementations return the expected cardinality.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_NodalBasisFamily.hpp"
#include "Intrepid2_SerendipityBasisFamily.hpp"
#include "Intrepid2_Types.hpp"

#include "Intrepid2_TestUtils.hpp"

#include "Kokkos_Core.hpp"

namespace
{
  using namespace Intrepid2;

  void testQuadBasisCardinality(Intrepid2::EFunctionSpace fs, int polyOrder_x, int polyOrder_y, Teuchos::FancyOStream &out, bool &success)
  {
    shards::CellTopology cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<> >() );
    const bool defineVertexFunctions = true;
    int expectedCardinality = -1;
    switch (fs)
    {
      case Intrepid2::FUNCTION_SPACE_HVOL:
      case Intrepid2::FUNCTION_SPACE_HGRAD:
        expectedCardinality = (polyOrder_x + 1) * (polyOrder_y + 1);
        break;
      case Intrepid2::FUNCTION_SPACE_HCURL:
      case Intrepid2::FUNCTION_SPACE_HDIV:
        expectedCardinality = (polyOrder_x    ) * (polyOrder_y + 1)
                            + (polyOrder_x + 1) * (polyOrder_y    );
        break;
      default:
        out << "Unsupported function space.\n";
        success = false;
        return;
    }
    auto basis = getHierarchicalBasis<defineVertexFunctions>(cellTopo, fs, polyOrder_x, polyOrder_y, -1);
    if (basis->getCardinality() != expectedCardinality)
    {
      out << "FAILURE: expected cardinality of " << expectedCardinality << " but got " << basis->getCardinality() << std::endl;
      success = false;
    }
  }

  void testSerendipityQuadBasisCardinality(Intrepid2::EFunctionSpace fs, int polyOrder_x, int polyOrder_y, Teuchos::FancyOStream &out, bool &success)
  {
    using BasisFamily = SerendipityBasisFamily<DefaultTestDeviceType>;
    
    int expectedCardinality = 0;
    int maxDegree = std::max(polyOrder_x,polyOrder_y);
    int maxH1Degree = (fs == FUNCTION_SPACE_HVOL) ? maxDegree + 1 : maxDegree;
    
    switch (fs)
    {
      case Intrepid2::FUNCTION_SPACE_HVOL:  out << "Testing HVOL"; break;
      case Intrepid2::FUNCTION_SPACE_HGRAD: out << "Testing HGRAD"; break;
      case Intrepid2::FUNCTION_SPACE_HDIV:  out << "Testing HDIV"; break;
      case Intrepid2::FUNCTION_SPACE_HCURL: out << "Testing HCURL"; break;
      default:
        out << "Unhandled function space\n";
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unhandled function space");
    }
    out << " for polyOrder_x = " << polyOrder_x << ", polyOrder_y = " << polyOrder_y << "\n";
    
    switch (fs)
    {
      case Intrepid2::FUNCTION_SPACE_HVOL:
      case Intrepid2::FUNCTION_SPACE_HGRAD:
        for (int i=0; i<=polyOrder_x; i++)
        {
          const int i_H1 = (fs == FUNCTION_SPACE_HVOL) ? i + 1 : i;
          const int i_sl = (i_H1 > 1) ? i_H1 : 0; // superlinear part of i_H1
          for (int j=0; j<=polyOrder_y; j++)
          {
            const int j_H1 = (fs == FUNCTION_SPACE_HVOL) ? j + 1 : j;
            const int j_sl = (j_H1 > 1) ? j_H1 : 0;
            const int superlinearDegree = i_sl + j_sl;
            if (superlinearDegree <= maxH1Degree)
            {
              expectedCardinality++;
            }
          }
        }
        break;
      case Intrepid2::FUNCTION_SPACE_HCURL:
      case Intrepid2::FUNCTION_SPACE_HDIV:
      {
        /*
         In the H(div) and H(curl) bases, we want to include any members whose H^1 orders satisfy the Serendipity criterion.
         
         Have not worked out what the closed form for this is; let's just count in a brute-force manner.
         */
        
        // family H(vol) x H(grad)
        for (int i_vol=0; i_vol<polyOrder_x; i_vol++)
        {
          int i_grad = i_vol + 1;
          const int i_sl = (i_grad > 1) ? i_grad : 0; // superlinear part of i_grad
          for (int j=0; j<=polyOrder_y; j++)
          {
            const int j_sl = (j > 1) ? j : 0;
            const int superlinearDegree = i_sl + j_sl;
            if (superlinearDegree <= maxH1Degree)
            {
              expectedCardinality++;
            }
          }
        }
        // family H(grad) x H(vol)
        for (int i=0; i<=polyOrder_x; i++)
        {
          const int i_sl = (i > 1) ? i : 0; // superlinear part of i_grad
          for (int j_vol=0; j_vol<polyOrder_y; j_vol++)
          {
            const int j_grad = j_vol + 1;
            const int j_sl = (j_grad > 1) ? j_grad : 0;
            const int superlinearDegree = i_sl + j_sl;
            if (superlinearDegree <= maxDegree)
            {
              expectedCardinality++;
            }
          }
        }
      }
        break;
      default:
        out << "Unsupported function space.\n";
        success = false;
        return;
    }
    auto basis = getQuadrilateralBasis<BasisFamily>(fs, polyOrder_x, polyOrder_y);
    if (basis->getCardinality() != expectedCardinality)
    {
      out << "FAILURE: expected cardinality of " << expectedCardinality << " but got " << basis->getCardinality() << std::endl;
      success = false;
    }
  }
  
  void testHexBasisCardinality(Intrepid2::EFunctionSpace fs, int polyOrder_x, int polyOrder_y, int polyOrder_z, Teuchos::FancyOStream &out, bool &success)
  {
    shards::CellTopology cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<> >() );
    const bool defineVertexFunctions = true;
    int expectedCardinality = -1;
    switch (fs)
    {
      case Intrepid2::FUNCTION_SPACE_HVOL:
      case Intrepid2::FUNCTION_SPACE_HGRAD:
        expectedCardinality = (polyOrder_x + 1) * (polyOrder_y + 1) * (polyOrder_z + 1);
        break;
      case Intrepid2::FUNCTION_SPACE_HCURL:
        expectedCardinality = (polyOrder_x    ) * (polyOrder_y + 1) * (polyOrder_z + 1)
                            + (polyOrder_x + 1) * (polyOrder_y    ) * (polyOrder_z + 1)
                            + (polyOrder_x + 1) * (polyOrder_y + 1) * (polyOrder_z    );
        break;
      case Intrepid2::FUNCTION_SPACE_HDIV:
        expectedCardinality = (polyOrder_x + 1) * (polyOrder_y    ) * (polyOrder_z    )
                            + (polyOrder_x    ) * (polyOrder_y + 1) * (polyOrder_z    )
                            + (polyOrder_x    ) * (polyOrder_y    ) * (polyOrder_z + 1);
        break;
      default:
        out << "Unsupported function space.\n";
        success = false;
        return;
    }
    auto basis = getHierarchicalBasis<defineVertexFunctions>(cellTopo, fs, polyOrder_x, polyOrder_y, polyOrder_z);
    if (basis->getCardinality() != expectedCardinality)
    {
      out << "FAILURE: expected cardinality of " << expectedCardinality << " but got " << basis->getCardinality() << std::endl;
      success = false;
    }
  }

  void testSerendipityHexBasisCardinality(Intrepid2::EFunctionSpace fs, int polyOrder_x, int polyOrder_y, int polyOrder_z, Teuchos::FancyOStream &out, bool &success)
  {
    using BasisFamily = SerendipityBasisFamily<DefaultTestDeviceType>;
    
    int expectedCardinality = 0;
    using std::max;
    int maxDegree = max(max(polyOrder_x,polyOrder_y),polyOrder_z);
    int maxH1Degree = (fs == FUNCTION_SPACE_HVOL) ? maxDegree + 1 : maxDegree;
    
    switch (fs)
    {
      case Intrepid2::FUNCTION_SPACE_HVOL:  out << "Testing HVOL"; break;
      case Intrepid2::FUNCTION_SPACE_HGRAD: out << "Testing HGRAD"; break;
      case Intrepid2::FUNCTION_SPACE_HDIV:  out << "Testing HDIV"; break;
      case Intrepid2::FUNCTION_SPACE_HCURL: out << "Testing HCURL"; break;
      default:
        out << "Unhandled function space\n";
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unhandled function space");
    }
    out << " for polyOrder_x = " << polyOrder_x << ", polyOrder_y = " << polyOrder_y << ", polyOrder_z = " << polyOrder_z << "\n";
    
    switch (fs)
    {
      case Intrepid2::FUNCTION_SPACE_HVOL:
      case Intrepid2::FUNCTION_SPACE_HGRAD:
        for (int i=0; i<=polyOrder_x; i++)
        {
          const int i_H1 = (fs == FUNCTION_SPACE_HVOL) ? i + 1 : i;
          const int i_sl = (i_H1 > 1) ? i_H1 : 0; // superlinear part of i_H1
          for (int j=0; j<=polyOrder_y; j++)
          {
            const int j_H1 = (fs == FUNCTION_SPACE_HVOL) ? j + 1 : j;
            const int j_sl = (j_H1 > 1) ? j_H1 : 0;
            for (int k=0; k<=polyOrder_z; k++)
            {
              const int k_H1 = (fs == FUNCTION_SPACE_HVOL) ? k + 1 : k;
              const int k_sl = (k_H1 > 1) ? k_H1 : 0;
              const int superlinearDegree = i_sl + j_sl + k_sl;
              if (superlinearDegree <= maxH1Degree)
              {
                expectedCardinality++;
              }
            }
          }
        }
        break;
      case Intrepid2::FUNCTION_SPACE_HCURL:
      {
        /*
         In the H(div) and H(curl) bases, we want to include any members whose H^1 orders satisfy the Serendipity criterion.
         
         Have not worked out what the closed form for this is; let's just count in a brute-force manner.
         */
        
        // H(curl) is H(vol) x H(grad) x H(grad), and permutations thereof
        
        // family H(vol) x H(grad) x H(grad)
        for (int i_vol=0; i_vol<polyOrder_x; i_vol++)
        {
          int i_grad = i_vol + 1;
          const int i_sl = (i_grad > 1) ? i_grad : 0; // superlinear part of i_grad
          for (int j=0; j<=polyOrder_y; j++)
          {
            const int j_sl = (j > 1) ? j : 0;
            for (int k=0; k<=polyOrder_z; k++)
            {
              const int k_sl = (k > 1) ? k : 0;
              const int superlinearDegree = i_sl + j_sl + k_sl;
              if (superlinearDegree <= maxH1Degree)
              {
                expectedCardinality++;
              }
            }
          }
        }
        // family H(grad) x H(vol) x H(grad)
        for (int i=0; i<=polyOrder_x; i++)
        {
          const int i_sl = (i > 1) ? i : 0;
          for (int j_vol=0; j_vol<polyOrder_y; j_vol++)
          {
            int j_grad = j_vol + 1;
            const int j_sl = (j_grad > 1) ? j_grad : 0; // superlinear part of j_grad
            for (int k=0; k<=polyOrder_z; k++)
            {
              const int k_sl = (k > 1) ? k : 0;
              const int superlinearDegree = i_sl + j_sl + k_sl;
              if (superlinearDegree <= maxH1Degree)
              {
                expectedCardinality++;
              }
            }
          }
        }
        
        // family H(grad) x H(grad) x H(vol)
        for (int i=0; i<=polyOrder_x; i++)
        {
          const int i_sl = (i > 1) ? i : 0;
          for (int j=0; j<=polyOrder_y; j++)
          {
            const int j_sl = (j > 1) ? j : 0;
            for (int k_vol=0; k_vol<polyOrder_z; k_vol++)
            {
              int k_grad = k_vol + 1;
              const int k_sl = (k_grad > 1) ? k_grad : 0; // superlinear part of k_grad
              const int superlinearDegree = i_sl + j_sl + k_sl;
              if (superlinearDegree <= maxH1Degree)
              {
                expectedCardinality++;
              }
            }
          }
        }
      }
      break;
      case Intrepid2::FUNCTION_SPACE_HDIV:
      {
        /*
         In the H(div) and H(curl) bases, we want to include any members whose H^1 orders satisfy the Serendipity criterion.
         
         Have not worked out what the closed form for this is; let's just count in a brute-force manner.
         */
        
        // H(vol) is H(vol) x H(vol) x H(grad), and permutations thereof
        
        // family H(vol) x H(vol) x H(grad)
        for (int i_vol=0; i_vol<polyOrder_x; i_vol++)
        {
          const int i_grad = i_vol + 1;
          const int i_sl = (i_grad > 1) ? i_grad : 0; // superlinear part of i_grad
          for (int j_vol=0; j_vol<polyOrder_y; j_vol++)
          {
            const int j_grad = j_vol + 1;
            const int j_sl = (j_grad > 1) ? j_grad : 0; // superlinear part of j_grad
            for (int k=0; k<=polyOrder_z; k++)
            {
              const int k_sl = (k > 1) ? k : 0;
              const int superlinearDegree = i_sl + j_sl + k_sl;
              if (superlinearDegree <= maxH1Degree)
              {
                expectedCardinality++;
              }
            }
          }
        }
        // family H(vol) x H(grad) x H(vol)
        for (int i_vol=0; i_vol<polyOrder_x; i_vol++)
        {
          const int i_grad = i_vol + 1;
          const int i_sl = (i_grad > 1) ? i_grad : 0; // superlinear part of i_grad
          for (int j=0; j<=polyOrder_y; j++)
          {
            const int j_sl = (j > 1) ? j : 0;
            for (int k_vol=0; k_vol<polyOrder_z; k_vol++)
            {
              int k_grad = k_vol + 1;
              const int k_sl = (k_grad > 1) ? k_grad : 0; // superlinear part of k_grad
              const int superlinearDegree = i_sl + j_sl + k_sl;
              if (superlinearDegree <= maxH1Degree)
              {
                expectedCardinality++;
              }
            }
          }
        }
        
        // family H(grad) x H(vol) x H(vol)
        for (int i=0; i<=polyOrder_x; i++)
        {
          const int i_sl = (i > 1) ? i : 0;
          for (int j_vol=0; j_vol<polyOrder_y; j_vol++)
          {
            const int j_grad = j_vol + 1;
            const int j_sl = (j_grad > 1) ? j_grad : 0; // superlinear part of j_grad
            for (int k_vol=0; k_vol<polyOrder_z; k_vol++)
            {
              int k_grad = k_vol + 1;
              const int k_sl = (k_grad > 1) ? k_grad : 0; // superlinear part of k_grad
              const int superlinearDegree = i_sl + j_sl + k_sl;
              if (superlinearDegree <= maxH1Degree)
              {
                expectedCardinality++;
              }
            }
          }
        }
      }
        break;
      default:
        out << "Unsupported function space.\n";
        success = false;
        return;
    }
    auto basis = getHexahedronBasis<BasisFamily>(fs, polyOrder_x, polyOrder_y, polyOrder_z);
    TEST_EQUALITY(basis->getCardinality(), expectedCardinality);
    if (basis->getCardinality() != expectedCardinality)
    {
      out << "FAILURE: expected cardinality of " << expectedCardinality << " but got " << basis->getCardinality() << std::endl;
      success = false;
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisCardinality, Quadrilateral )
  {
    std::vector<Intrepid2::EFunctionSpace> functionSpaces_2D = {FUNCTION_SPACE_HGRAD,FUNCTION_SPACE_HCURL,FUNCTION_SPACE_HDIV,FUNCTION_SPACE_HVOL};
    
    const int maxDegreeForCardinalityTests = 3;

    for (auto fs : functionSpaces_2D)
    {
      const int minDegree = (fs == FUNCTION_SPACE_HVOL) ? 0 : 1;
      int spaceDim = 2;
      auto cardinalityTestCases = getBasisTestCasesUpToDegree(spaceDim, minDegree, maxDegreeForCardinalityTests, maxDegreeForCardinalityTests);
      for (auto testCase : cardinalityTestCases)
      {
        testQuadBasisCardinality(fs, testCase[0], testCase[1], out, success);
      }
    }
  }

  TEUCHOS_UNIT_TEST( BasisCardinality, Quadrilateral_Serendipity )
  {
    std::vector<Intrepid2::EFunctionSpace> functionSpaces_2D = {FUNCTION_SPACE_HGRAD,FUNCTION_SPACE_HCURL,FUNCTION_SPACE_HDIV,FUNCTION_SPACE_HVOL};
    
    const int maxDegreeForCardinalityTests = 3;

    for (auto fs : functionSpaces_2D)
    {
      const int minDegree = (fs == FUNCTION_SPACE_HVOL) ? 0 : 1;
      int spaceDim = 2;
      auto cardinalityTestCases = getBasisTestCasesUpToDegree(spaceDim, minDegree, maxDegreeForCardinalityTests, maxDegreeForCardinalityTests);
      for (auto testCase : cardinalityTestCases)
      {
        testSerendipityQuadBasisCardinality(fs, testCase[0], testCase[1], out, success);
      }
    }
  }
  
  TEUCHOS_UNIT_TEST( BasisCardinality, Hexahedron )
  {
    std::vector<Intrepid2::EFunctionSpace> functionSpaces_3D = {FUNCTION_SPACE_HGRAD,FUNCTION_SPACE_HCURL,FUNCTION_SPACE_HDIV,FUNCTION_SPACE_HVOL};

    const int maxDegreeForCardinalityTests = 3;
    
    for (auto fs : functionSpaces_3D)
    {
      const int minDegree = (fs == FUNCTION_SPACE_HVOL) ? 0 : 1;
      const int spaceDim = 3;
      auto cardinalityTestCases = getBasisTestCasesUpToDegree(spaceDim, minDegree, maxDegreeForCardinalityTests, maxDegreeForCardinalityTests, maxDegreeForCardinalityTests);
      for (auto testCase : cardinalityTestCases)
      {
        testHexBasisCardinality(fs, testCase[0], testCase[1], testCase[2], out, success);
      }
    }
  }

  TEUCHOS_UNIT_TEST( BasisCardinality, Hexahedron_Serendipity )
  {
    std::vector<Intrepid2::EFunctionSpace> functionSpaces_3D = {FUNCTION_SPACE_HGRAD,FUNCTION_SPACE_HCURL,FUNCTION_SPACE_HDIV,FUNCTION_SPACE_HVOL};
    
    const int maxDegreeForCardinalityTests = 3;

    for (auto fs : functionSpaces_3D)
    {
      const int minDegree = (fs == FUNCTION_SPACE_HVOL) ? 0 : 1;
      int spaceDim = 3;
      auto cardinalityTestCases = getBasisTestCasesUpToDegree(spaceDim, minDegree, maxDegreeForCardinalityTests, maxDegreeForCardinalityTests, maxDegreeForCardinalityTests);
      for (auto testCase : cardinalityTestCases)
      {
        testSerendipityHexBasisCardinality(fs, testCase[0], testCase[1], testCase[2], out, success);
      }
    }
  }

  TEUCHOS_UNIT_TEST( BasisCardinality, Hypercube )
  {
    std::vector<Intrepid2::EFunctionSpace> functionSpaces = {FUNCTION_SPACE_HGRAD,FUNCTION_SPACE_HVOL};

    const int maxSpaceDim = 7;
    const int maxDegreeForCardinalityTests = 2;
    
    using HierarchicalBasisFamily = HierarchicalBasisFamily<DefaultTestDeviceType>;
    using NodalBasisFamily = NodalBasisFamily<DefaultTestDeviceType>;
    
    for (auto fs : functionSpaces)
    {
      if (fs == FUNCTION_SPACE_HGRAD)
      {
        out << "Testing HGRAD.\n";
      }
      else if (fs == FUNCTION_SPACE_HVOL)
      {
        out << "Testing HVOL.\n";
      }
      
      const int minDegree = (fs == FUNCTION_SPACE_HVOL) ? 0 : 1;
      for (int polyDegree = minDegree; polyDegree <= maxDegreeForCardinalityTests; polyDegree++)
      {
        out << "** polyDegree " << polyDegree << " **\n";
        for (int spaceDim = 1; spaceDim <= maxSpaceDim; spaceDim++)
        {
          out << "** spaceDim " << spaceDim << " **\n";
          int expectedCardinality = 1;
          for (int d=0; d<spaceDim; d++)
          {
            expectedCardinality *= (polyDegree + 1);
          }
          int expectedExtrusionCount = spaceDim - 1;
          if (fs == FUNCTION_SPACE_HGRAD)
          {
            auto hierarchicalBasis = getHypercubeBasis_HGRAD<HierarchicalBasisFamily>(polyDegree, spaceDim);
            auto nodalBasis        = getHypercubeBasis_HGRAD<NodalBasisFamily>(polyDegree, spaceDim);
            
            TEST_EQUALITY(expectedCardinality, hierarchicalBasis->getCardinality());
            TEST_EQUALITY(expectedCardinality, nodalBasis->getCardinality());
            
            TEST_EQUALITY(expectedExtrusionCount, hierarchicalBasis->getNumTensorialExtrusions());
            TEST_EQUALITY(expectedExtrusionCount, nodalBasis->getNumTensorialExtrusions());
          }
          else if (fs == FUNCTION_SPACE_HVOL)
          {
            auto hierarchicalBasis = getHypercubeBasis_HVOL<HierarchicalBasisFamily>(polyDegree, spaceDim);
            auto nodalBasis        = getHypercubeBasis_HVOL<NodalBasisFamily>(polyDegree, spaceDim);
            
            TEST_EQUALITY(expectedCardinality, hierarchicalBasis->getCardinality());
            TEST_EQUALITY(expectedCardinality, nodalBasis->getCardinality());
            
            TEST_EQUALITY(expectedExtrusionCount, hierarchicalBasis->getNumTensorialExtrusions());
            TEST_EQUALITY(expectedExtrusionCount, nodalBasis->getNumTensorialExtrusions());
          }
        }
      }
    }
  }

  //! Compute (n choose k)
  int choose(const int &n, const int &k)
  {
    int result = 1;
    if (n - k < k)
    {
      return choose(n, n-k);
    }
    else
    {
      for (int i=1; i<=k; i++)
      {
        result *= (n - k + i);
      }
      for (int i=2; i<=k; i++)
      {
        result /= i;
      }
    }
    return result;
  }

  TEUCHOS_UNIT_TEST( BasisCardinality, Pyramid_HGRAD )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_PYR;
    
    for (ordinal_type p=1; p<10; p++)
    {
      // expected cardinality is p^3 + 3p + 1
      const ordinal_type expectedCardinality = p * p * p + 3 * p + 1;
      
      HierarchicalBasis hierarchicalBasis(p);
      const ordinal_type actualCardinality = hierarchicalBasis.getCardinality();
      TEST_EQUALITY(expectedCardinality, actualCardinality);
    }
  }

  TEUCHOS_UNIT_TEST( BasisCardinality, Pyramid_HDIV )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HDIV_PYR;
    
    for (ordinal_type p=1; p<10; p++)
    {
      // expected cardinality is p^2 + 2p(p+1) + 3p^2(p-1):
      // - quadrilateral face: p^2
      // - triangular faces:   2p(p+1)
      // - interior bubbles:   3p^2(p-1)
      const ordinal_type expectedCardinality = p * p + 2 * p * (p+1) + 3 * p * p * (p-1);
      
      HierarchicalBasis hierarchicalBasis(p);
      const ordinal_type actualCardinality = hierarchicalBasis.getCardinality();
      TEST_EQUALITY(expectedCardinality, actualCardinality);
    }
  }

  TEUCHOS_UNIT_TEST( BasisCardinality, Pyramid_HVOL )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HVOL_PYR;
    
    for (ordinal_type p=1; p<10; p++)
    {
      // expected cardinality is (p+1)^3
      const ordinal_type expectedCardinality = (p+1) * (p+1) * (p+1);
      
      HierarchicalBasis hierarchicalBasis(p);
      const ordinal_type actualCardinality = hierarchicalBasis.getCardinality();
      TEST_EQUALITY(expectedCardinality, actualCardinality);
    }
  }

  TEUCHOS_UNIT_TEST( BasisCardinality, Serendipity_HDIV_HEX )
  {
    // TODO: finish this test.  (Can we do something templated on the BasisFamily??)
    
    std::vector<Intrepid2::EFunctionSpace> functionSpaces = {FUNCTION_SPACE_HDIV};

    const int maxSpaceDim = 7;
    const int maxDegreeForCardinalityTests = 7;
    
    using BasisFamily = SerendipityBasisFamily<DefaultTestDeviceType>;
    using BasisBase = typename BasisFamily::HGRAD_LINE::BasisBase;
    
    for (auto fs : functionSpaces)
    {
      if (fs == FUNCTION_SPACE_HGRAD)
      {
        out << "Testing HGRAD.\n";
      }
      else if (fs == FUNCTION_SPACE_HVOL)
      {
        out << "Testing HVOL.\n";
      }
      
      const int minDegree = (fs == FUNCTION_SPACE_HVOL) ? 0 : 1;
      for (int polyDegree = minDegree; polyDegree <= maxDegreeForCardinalityTests; polyDegree++)
      {
        out << "** polyDegree " << polyDegree << " **\n";
        for (int spaceDim = 1; spaceDim <= maxSpaceDim; spaceDim++)
        {
          out << "** spaceDim " << spaceDim << " **\n";
          int expectedCardinality = 0; // we'll sum into this
          int i_max = std::min(spaceDim,polyDegree/2);
          
          if (polyDegree == 0)
          {
            expectedCardinality = 1; // serendipity of constant basis is a constant basis
          }
          else
          {
            for (int i = 0; i <= i_max; i++)
            {
              int d_choose_i = choose(spaceDim, i);
              int p_minus_i_choose_i = choose(polyDegree - i, i);
              int two_to_the_d_minus_i = 1 << (spaceDim-i);
              expectedCardinality += two_to_the_d_minus_i * d_choose_i * p_minus_i_choose_i;
            }
          }
          if (fs == FUNCTION_SPACE_HGRAD)
          {
            auto fullBasis        = getHypercubeBasis_HGRAD<BasisFamily>(polyDegree, spaceDim);
            auto serendipityBasis = Teuchos::rcp(new SerendipityBasis<BasisBase>(fullBasis));
            
            TEST_EQUALITY(expectedCardinality, serendipityBasis->getCardinality());
          }
          else if (fs == FUNCTION_SPACE_HVOL)
          {
            auto fullBasis        = getHypercubeBasis_HVOL<BasisFamily>(polyDegree, spaceDim);
            auto serendipityBasis = Teuchos::rcp(new SerendipityBasis<BasisBase>(fullBasis));
            
            TEST_EQUALITY(expectedCardinality, serendipityBasis->getCardinality());
          }
        }
      }
    }
  }

  TEUCHOS_UNIT_TEST( BasisCardinality, Serendipity_Hypercube_HGRAD )
  {
    std::vector<Intrepid2::EFunctionSpace> functionSpaces = {FUNCTION_SPACE_HGRAD};

    const int maxSpaceDim = 7;
    const int maxDegreeForCardinalityTests = 7;
    
    // Serendipity requires a hierarchical basis
    using BasisFamily = HierarchicalBasisFamily<DefaultTestDeviceType>;
    using BasisBase = typename BasisFamily::HGRAD_LINE::BasisBase;
    
    for (auto fs : functionSpaces)
    {
      const int minDegree = (fs == FUNCTION_SPACE_HVOL) ? 0 : 1;
      for (int polyDegree = minDegree; polyDegree <= maxDegreeForCardinalityTests; polyDegree++)
      {
        out << "** polyDegree " << polyDegree << " **\n";
        for (int spaceDim = 1; spaceDim <= maxSpaceDim; spaceDim++)
        {
          out << "** spaceDim " << spaceDim << " **\n";
          int expectedCardinality = 0; // we'll sum into this
          int i_max = std::min(spaceDim,polyDegree/2);
          
          if (polyDegree == 0)
          {
            expectedCardinality = 1; // serendipity of constant basis is a constant basis
          }
          else
          {
            // closed form for H^1 -- see Juno et al. doi:10.1016/j.jcp.2017.10.009, Eqn. (83).
            for (int i = 0; i <= i_max; i++)
            {
              int d_choose_i = choose(spaceDim, i);
              int p_minus_i_choose_i = choose(polyDegree - i, i);
              int two_to_the_d_minus_i = 1 << (spaceDim-i);
              expectedCardinality += two_to_the_d_minus_i * d_choose_i * p_minus_i_choose_i;
            }
          }
          auto fullBasis        = getHypercubeBasis_HGRAD<BasisFamily>(polyDegree, spaceDim);
          auto serendipityBasis = Teuchos::rcp(new SerendipityBasis<BasisBase>(fullBasis));
          
          TEST_EQUALITY(expectedCardinality, serendipityBasis->getCardinality());
        }
      }
    }
  }
  
} // namespace
