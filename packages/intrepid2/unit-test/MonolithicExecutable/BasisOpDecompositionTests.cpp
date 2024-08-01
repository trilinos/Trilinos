// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   BasisOpDecompositionTests.cpp
    \brief  Tests against operator decompositions implemented by various TensorBasis subclasses.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_NodalBasisFamily.hpp"
#include "Intrepid2_ScalarView.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_TestUtils.hpp"

#include <Kokkos_Core.hpp>

using namespace Intrepid2;

namespace
{
  using namespace Intrepid2;

  // OPERATOR_D1 and OPERATOR_GRAD mean the same thing.  This method allows us to replace D1 with GRAD in testing that we match the expected decomposition.
  EOperator canonicalOp(const EOperator &op)
  {
    if (op == OPERATOR_D1)
    {
      return OPERATOR_GRAD;
    }
    else
    {
      return op;
    }
  }

  //! tests that hierarchical basis has the expected decomposition into component bases
  template<class Basis>
  void testBasisOpDecomposition(Basis &basis, const EOperator &opToTest, std::vector< std::vector<EOperator> > &expectedOps,
                                std::vector<double> &expectedWeights, Teuchos::FancyOStream &out, bool &success)
  {
    auto decomp = basis.getOperatorDecomposition(opToTest);
    
    const ordinal_type expectedOpsSize=expectedOps.size();
    TEST_EQUALITY(decomp.numVectorComponents(), expectedOpsSize);
    
    for (ordinal_type vectorEntry=0; vectorEntry<decomp.numVectorComponents(); vectorEntry++)
    {
      const ordinal_type size=expectedOps[vectorEntry].size();
      if (decomp.identicallyZeroComponent(vectorEntry))
      {
        TEST_ASSERT(size == 0);
        continue;
      }
      else
      {
        TEST_EQUALITY(size, decomp.numBasisComponents());
      }
      if (size == decomp.numBasisComponents())
      {
        for (ordinal_type basisEntry=0; basisEntry<decomp.numBasisComponents(); basisEntry++)
        {
          // for easier reading of test failures, compare as strings
          auto expectedOp       = canonicalOp(expectedOps[vectorEntry][basisEntry]);
          auto expectedOpString = EOperatorToString(expectedOp);
          auto actualOp         = canonicalOp(decomp.op(vectorEntry, basisEntry));
          auto actualOpString   = EOperatorToString(actualOp);
          TEST_EQUALITY( expectedOpString, actualOpString );
        }
      }
      TEST_EQUALITY( expectedWeights[vectorEntry], decomp.weight(vectorEntry) );
    }
  }
  
TEUCHOS_UNIT_TEST( BasisOpDecomposition, HierarchicalHGRAD_QUAD )
{
  using Basis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_QUAD;
  
  auto GRAD  = OPERATOR_GRAD;
  auto VALUE = OPERATOR_VALUE;
  
  EOperator opToTest = OPERATOR_GRAD;
  std::vector< std::vector<EOperator> > expectedDecomposition;
  expectedDecomposition.push_back( std::vector<EOperator>{GRAD,  VALUE} );
  expectedDecomposition.push_back( std::vector<EOperator>{VALUE,  GRAD} );
  
  std::vector<double> expectedWeights(2, 1.0);
  
  const int polyOrder = 3; // arbitrary; should not matter
  Basis basis(polyOrder);
  
  testBasisOpDecomposition(basis, opToTest, expectedDecomposition, expectedWeights, out, success);
}

  TEUCHOS_UNIT_TEST( BasisOpDecomposition, HierarchicalHGRAD_HEX )
  {
    using Basis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_HEX;
    
    auto GRAD  = OPERATOR_GRAD;
    auto VALUE = OPERATOR_VALUE;
    
    EOperator opToTest = OPERATOR_GRAD;
    std::vector< std::vector<EOperator> > expectedDecomposition;
    expectedDecomposition.push_back( std::vector<EOperator>{GRAD,  VALUE, VALUE} );
    expectedDecomposition.push_back( std::vector<EOperator>{VALUE,  GRAD, VALUE} );
    expectedDecomposition.push_back( std::vector<EOperator>{VALUE, VALUE,  GRAD} );
    
    std::vector<double> expectedWeights(3, 1.0);
    
    const int polyOrder = 3; // arbitrary; should not matter
    Basis basis(polyOrder);
    
    testBasisOpDecomposition(basis, opToTest, expectedDecomposition, expectedWeights, out, success);
  }

  TEUCHOS_UNIT_TEST( BasisOpDecomposition, HierarchicalHDIV_HEX_Family1 )
  {
    using BasisGRAD = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_LINE;
    using BasisVOL  = HierarchicalBasisFamily<DefaultTestDeviceType>::HVOL_LINE;
    
    using Basis = Basis_Derived_HDIV_Family1_HEX<BasisGRAD,BasisVOL>;
    const int familyOrdinal = 0;
    
    auto DIV   = OPERATOR_DIV;
    auto GRAD  = OPERATOR_GRAD;
    auto VALUE = OPERATOR_VALUE;
    
    const int polyOrder = 3; // arbitrary; should not matter
    Basis basis(polyOrder,polyOrder,polyOrder);
    
    // test VALUE:
    EOperator opToTest = VALUE;
    out << "Testing " << EOperatorToString(opToTest) << std::endl;
    std::vector< std::vector<EOperator> > expectedDecomposition(3);
    
    // VALUE goes in each slot
    std::vector<EOperator> ops(3, VALUE);
    
    expectedDecomposition[familyOrdinal] = ops; // nonzero in the diagonal slot

    std::vector<double> expectedWeights(3, 0.0);
    expectedWeights[familyOrdinal] = 1.0;
    
    testBasisOpDecomposition(basis, opToTest, expectedDecomposition, expectedWeights, out, success);
    
    // test DIV:
    opToTest = DIV;
    out << "Testing " << EOperatorToString(opToTest) << std::endl;
    expectedDecomposition = std::vector< std::vector<EOperator> >(1); // scalar-valued
    
    // GRAD goes in the "diagonal" slot
    ops = std::vector<EOperator>(3, VALUE);
    ops[familyOrdinal] = GRAD;
    
    expectedDecomposition[0] = ops;

    expectedWeights = std::vector<double>(1, 1.0);
    
    testBasisOpDecomposition(basis, opToTest, expectedDecomposition, expectedWeights, out, success);
  }

  TEUCHOS_UNIT_TEST( BasisOpDecomposition, HierarchicalHDIV_HEX_Family2 )
  {
    
    using BasisGRAD = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_LINE;
    using BasisVOL  = HierarchicalBasisFamily<DefaultTestDeviceType>::HVOL_LINE;
    
    using Basis = Basis_Derived_HDIV_Family2_HEX<BasisGRAD,BasisVOL>;
    const int familyOrdinal = 1;
    
    auto DIV   = OPERATOR_DIV;
    auto GRAD  = OPERATOR_GRAD;
    auto VALUE = OPERATOR_VALUE;
    
    const int polyOrder = 3; // arbitrary; should not matter
    Basis basis(polyOrder,polyOrder,polyOrder);
    
    // test VALUE:
    EOperator opToTest = VALUE;
    out << "Testing " << EOperatorToString(opToTest) << std::endl;
    std::vector< std::vector<EOperator> > expectedDecomposition(3);
    
    // VALUE goes in each slot
    std::vector<EOperator> ops(3, VALUE);
    
    expectedDecomposition[familyOrdinal] = ops; // nonzero in the diagonal slot

    std::vector<double> expectedWeights(3, 0.0);
    expectedWeights[familyOrdinal] = 1.0;
    
    testBasisOpDecomposition(basis, opToTest, expectedDecomposition, expectedWeights, out, success);
    
    // test DIV:
    opToTest = DIV;
    out << "Testing " << EOperatorToString(opToTest) << std::endl;
    expectedDecomposition = std::vector< std::vector<EOperator> >(1); // scalar-valued
    
    // GRAD goes in the "diagonal" slot
    ops = std::vector<EOperator>(3, VALUE);
    ops[familyOrdinal] = GRAD;
    
    expectedDecomposition[0] = ops;

    expectedWeights = std::vector<double>(1, 1.0);
    
    testBasisOpDecomposition(basis, opToTest, expectedDecomposition, expectedWeights, out, success);
  }

  TEUCHOS_UNIT_TEST( BasisOpDecomposition, HierarchicalHDIV_HEX_Family3 )
  {
    
    using BasisGRAD = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_LINE;
    using BasisVOL  = HierarchicalBasisFamily<DefaultTestDeviceType>::HVOL_LINE;
    
    using Basis = Basis_Derived_HDIV_Family3_HEX<BasisGRAD,BasisVOL>;
    const int familyOrdinal = 2;
    
    auto DIV   = OPERATOR_DIV;
    auto GRAD  = OPERATOR_GRAD;
    auto VALUE = OPERATOR_VALUE;
    
    const int polyOrder = 3; // arbitrary; should not matter
    Basis basis(polyOrder,polyOrder,polyOrder);
    
    // test VALUE:
    EOperator opToTest = VALUE;
    out << "Testing " << EOperatorToString(opToTest) << std::endl;
    std::vector< std::vector<EOperator> > expectedDecomposition(3);
    
    // VALUE goes in each slot
    std::vector<EOperator> ops(3, VALUE);
    
    expectedDecomposition[familyOrdinal] = ops; // nonzero in the diagonal slot

    std::vector<double> expectedWeights(3, 0.0);
    expectedWeights[familyOrdinal] = 1.0;
    
    testBasisOpDecomposition(basis, opToTest, expectedDecomposition, expectedWeights, out, success);
    
    // test DIV:
    opToTest = DIV;
    out << "Testing " << EOperatorToString(opToTest) << std::endl;
    expectedDecomposition = std::vector< std::vector<EOperator> >(1); // scalar-valued
    
    // GRAD goes in the "diagonal" slot
    ops = std::vector<EOperator>(3, VALUE);
    ops[familyOrdinal] = GRAD;
    
    expectedDecomposition[0] = ops;

    expectedWeights = std::vector<double>(1, 1.0);
    
    testBasisOpDecomposition(basis, opToTest, expectedDecomposition, expectedWeights, out, success);
  }

  TEUCHOS_UNIT_TEST( BasisOpDecomposition, HierarchicalHCURL_QUAD_Family1 )
  {
    using BasisGRAD = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_LINE;
    using BasisVOL  = HierarchicalBasisFamily<DefaultTestDeviceType>::HVOL_LINE;
    
    using Basis = Basis_Derived_HCURL_Family1_QUAD<BasisGRAD,BasisVOL>;
    const int familyOrdinal = 0;
    
    auto CURL  = OPERATOR_CURL;
    auto GRAD  = OPERATOR_GRAD;
    auto VALUE = OPERATOR_VALUE;
    
    const int polyOrder = 3; // arbitrary; should not matter
    Basis basis(polyOrder,polyOrder);
    
    // test VALUE:
    EOperator opToTest = VALUE;
    out << "Testing " << EOperatorToString(opToTest) << std::endl;
    std::vector< std::vector<EOperator> > expectedDecomposition(2);
    
    std::vector<EOperator> ops(2, VALUE);
    expectedDecomposition[0] = std::vector<EOperator>{VALUE,VALUE};
    expectedDecomposition[1] = std::vector<EOperator>{};

    std::vector<double> expectedWeights(2, 0.0);
    expectedWeights[familyOrdinal] = 1.0;
    
    testBasisOpDecomposition(basis, opToTest, expectedDecomposition, expectedWeights, out, success);
    
    // test CURL:
    opToTest = CURL;
    out << "Testing " << EOperatorToString(opToTest) << std::endl;
    expectedDecomposition = std::vector< std::vector<EOperator> >(1);
    
    // Family 1:
    // 2D curl: scalar-valued
    // x component of curl is -d/dy: (VALUE,GRAD), weight = -1.0
    
    expectedDecomposition[0] = std::vector<EOperator>{VALUE,GRAD};

    expectedWeights = std::vector<double> {-1.0};
    
    testBasisOpDecomposition(basis, opToTest, expectedDecomposition, expectedWeights, out, success);
  }

  TEUCHOS_UNIT_TEST( BasisOpDecomposition, HierarchicalHCURL_QUAD_Family2 )
  {
    
    using BasisGRAD = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_LINE;
    using BasisVOL  = HierarchicalBasisFamily<DefaultTestDeviceType>::HVOL_LINE;
    
    using Basis = Basis_Derived_HCURL_Family2_QUAD<BasisGRAD,BasisVOL>;
    const int familyOrdinal = 1;
    
    auto CURL  = OPERATOR_CURL;
    auto GRAD  = OPERATOR_GRAD;
    auto VALUE = OPERATOR_VALUE;
    
    const int polyOrder = 3; // arbitrary; should not matter
    Basis basis(polyOrder,polyOrder);
    
    // test VALUE:
    EOperator opToTest = VALUE;
    out << "Testing " << EOperatorToString(opToTest) << std::endl;
    std::vector< std::vector<EOperator> > expectedDecomposition(2);
    
    std::vector<EOperator> ops(2, VALUE);
    expectedDecomposition[0] = std::vector<EOperator>{};
    expectedDecomposition[1] = std::vector<EOperator>{VALUE,VALUE};

    std::vector<double> expectedWeights(2, 0.0);
    expectedWeights[familyOrdinal] = 1.0;
    
    testBasisOpDecomposition(basis, opToTest, expectedDecomposition, expectedWeights, out, success);
    
    // test CURL:
    opToTest = CURL;
    out << "Testing " << EOperatorToString(opToTest) << std::endl;
    expectedDecomposition = std::vector< std::vector<EOperator> >(1);
    
    // Family 2:
    // 2D curl: scalar-valued
    // x component iz zero
    // y component is d/dx: (GRAD,VALUE), weight = 1.0
    
    expectedDecomposition[0] = std::vector<EOperator>{GRAD,VALUE};

    expectedWeights = std::vector<double> {1.0};
    
    testBasisOpDecomposition(basis, opToTest, expectedDecomposition, expectedWeights, out, success);
  }

  TEUCHOS_UNIT_TEST( BasisOpDecomposition, HierarchicalHCURL_HEX_Family1 )
  {
    using BasisGRAD = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_LINE;
    using BasisVOL  = HierarchicalBasisFamily<DefaultTestDeviceType>::HVOL_LINE;
    
    using Basis = Basis_Derived_HCURL_Family1_HEX<BasisGRAD,BasisVOL>;
    const int familyOrdinal = 0;
    
    auto CURL  = OPERATOR_CURL;
    auto GRAD  = OPERATOR_GRAD;
    auto VALUE = OPERATOR_VALUE;
    
    const int polyOrder = 3; // arbitrary; should not matter
    Basis basis(polyOrder,polyOrder,polyOrder);
    
    // test VALUE:
    EOperator opToTest = VALUE;
    out << "Testing " << EOperatorToString(opToTest) << std::endl;
    std::vector< std::vector<EOperator> > expectedDecomposition(3);
    
    std::vector<EOperator> ops(3, VALUE);
    expectedDecomposition[0] = std::vector<EOperator>{VALUE,VALUE,VALUE};
    expectedDecomposition[1] = std::vector<EOperator>{};
    expectedDecomposition[2] = std::vector<EOperator>{};

    std::vector<double> expectedWeights(3, 0.0);
    expectedWeights[familyOrdinal] = 1.0;
    
    testBasisOpDecomposition(basis, opToTest, expectedDecomposition, expectedWeights, out, success);
    
    // test CURL:
    opToTest = CURL;
    out << "Testing " << EOperatorToString(opToTest) << std::endl;
    expectedDecomposition = std::vector< std::vector<EOperator> >(3);
    
    // Family 1:
    // x component is zero
    // y component is  d/dz: (VALUE,VALUE,GRAD), weight =  1.0
    // z component is -d/dy: (VALUE,GRAD,VALUE), weight = -1.0
    
    expectedDecomposition[0] = std::vector<EOperator>{};
    expectedDecomposition[1] = std::vector<EOperator>{VALUE,VALUE,GRAD};
    expectedDecomposition[2] = std::vector<EOperator>{VALUE,GRAD,VALUE};

    expectedWeights = std::vector<double> {0.0, 1.0, -1.0};
    
    testBasisOpDecomposition(basis, opToTest, expectedDecomposition, expectedWeights, out, success);
  }

  TEUCHOS_UNIT_TEST( BasisOpDecomposition, HierarchicalHCURL_HEX_Family2 )
  {
    
    using BasisGRAD = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_LINE;
    using BasisVOL  = HierarchicalBasisFamily<DefaultTestDeviceType>::HVOL_LINE;
    
    using Basis = Basis_Derived_HCURL_Family2_HEX<BasisGRAD,BasisVOL>;
    const int familyOrdinal = 1;
    
    auto CURL  = OPERATOR_CURL;
    auto GRAD  = OPERATOR_GRAD;
    auto VALUE = OPERATOR_VALUE;
    
    const int polyOrder = 3; // arbitrary; should not matter
    Basis basis(polyOrder,polyOrder,polyOrder);
    
    // test VALUE:
    EOperator opToTest = VALUE;
    out << "Testing " << EOperatorToString(opToTest) << std::endl;
    std::vector< std::vector<EOperator> > expectedDecomposition(3);
    
    std::vector<EOperator> ops(3, VALUE);
    expectedDecomposition[0] = std::vector<EOperator>{};
    expectedDecomposition[1] = std::vector<EOperator>{VALUE,VALUE,VALUE};
    expectedDecomposition[2] = std::vector<EOperator>{};

    std::vector<double> expectedWeights(3, 0.0);
    expectedWeights[familyOrdinal] = 1.0;
    
    testBasisOpDecomposition(basis, opToTest, expectedDecomposition, expectedWeights, out, success);
    
    // test CURL:
    opToTest = CURL;
    out << "Testing " << EOperatorToString(opToTest) << std::endl;
    expectedDecomposition = std::vector< std::vector<EOperator> >(3);
    
    // Family 2:
    // x component is -d/dz: (VALUE,VALUE,GRAD), weight = -1.0
    // y component iz zero
    // z component is  d/dx: (VALUE,GRAD,VALUE), weight =  1.0
    
    expectedDecomposition[0] = std::vector<EOperator>{VALUE,VALUE,GRAD};
    expectedDecomposition[1] = std::vector<EOperator>{};
    expectedDecomposition[2] = std::vector<EOperator>{GRAD,VALUE,VALUE};

    expectedWeights = std::vector<double> {-1.0, 0.0, 1.0};
    
    testBasisOpDecomposition(basis, opToTest, expectedDecomposition, expectedWeights, out, success);
  }

  TEUCHOS_UNIT_TEST( BasisOpDecomposition, HierarchicalHCURL_HEX_Family3 )
  {
    
    using BasisGRAD = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_LINE;
    using BasisVOL  = HierarchicalBasisFamily<DefaultTestDeviceType>::HVOL_LINE;
    
    using Basis = Basis_Derived_HCURL_Family3_HEX<BasisGRAD,BasisVOL>;
    const int familyOrdinal = 2;
    
    auto CURL   = OPERATOR_CURL;
    auto GRAD  = OPERATOR_GRAD;
    auto VALUE = OPERATOR_VALUE;
    
    const int polyOrder = 3; // arbitrary; should not matter
    Basis basis(polyOrder,polyOrder,polyOrder);
    
    // test VALUE:
    EOperator opToTest = VALUE;
    out << "Testing " << EOperatorToString(opToTest) << std::endl;
    std::vector< std::vector<EOperator> > expectedDecomposition(3);
    
    std::vector<EOperator> ops(3, VALUE);
    expectedDecomposition[0] = std::vector<EOperator>{};
    expectedDecomposition[1] = std::vector<EOperator>{};
    expectedDecomposition[2] = std::vector<EOperator>{VALUE,VALUE,VALUE};

    std::vector<double> expectedWeights(3, 0.0);
    expectedWeights[familyOrdinal] = 1.0;
    
    testBasisOpDecomposition(basis, opToTest, expectedDecomposition, expectedWeights, out, success);
    
    // test CURL:
    opToTest = CURL;
    out << "Testing " << EOperatorToString(opToTest) << std::endl;
    expectedDecomposition = std::vector< std::vector<EOperator> >(3);
    
    // Family 3:
    // x component is  d/dy: (VALUE,GRAD,VALUE), weight =  1.0
    // y component is -d/dx: (GRAD,VALUE,VALUE), weight = -1.0
    // z component is zero
    
    expectedDecomposition[0] = std::vector<EOperator>{VALUE,GRAD,VALUE};
    expectedDecomposition[1] = std::vector<EOperator>{GRAD,VALUE,VALUE};
    expectedDecomposition[2] = std::vector<EOperator>{};

    expectedWeights = std::vector<double> {1.0, -1.0, 0.0};
    
    testBasisOpDecomposition(basis, opToTest, expectedDecomposition, expectedWeights, out, success);
  }
} // namespace
