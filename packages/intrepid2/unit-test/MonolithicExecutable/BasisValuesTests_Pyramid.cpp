// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   BasisValuesTests_Pyramid.cpp
    \brief  Some tests of the values of the H^1 hierarchical basis on the pyramid.
    \author Created by N.V. Roberts.
 
 For now, these tests are somewhat limited: we have one test that verifies that the lowest-order basis is nodal at the vertices, and one test, primarily intended for debugging, which exercises computation of the p=9 basis at a single point.
 
 We have done considerably more testing offline, verifying that values and gradients agree with the ESEAS implementation of the same basis on a 1240-point lattice, for p up to 9.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "BasisEquivalenceHelpers.hpp"

#include "Intrepid2_HGRAD_PYR_C1_FEM.hpp"

using namespace Intrepid2;

namespace
{
  TEUCHOS_UNIT_TEST( BasisValues, PyramidHierarchical_HGRAD_VertexValues )
  {
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_PYR;
    using DeviceType = DefaultTestDeviceType;
    
    const double relTol=1e-13;
    const double absTol=1e-13;
    
    const ordinal_type p = 1;
    
    HierarchicalBasis hierarchicalBasis(p);
    
    using std::vector;
    const int numPoints = 5;
    vector< vector<double> > ESEAS_points(5, vector<double>(3)); // these have a pyramid base of (0,1)^2
    
    ESEAS_points[0] = {0,0,0};
    ESEAS_points[1] = {1,0,0};
    ESEAS_points[2] = {1,1,0};
    ESEAS_points[3] = {0,1,0};
    ESEAS_points[4] = {0,0,1};
    
    ViewType<double,DeviceType> intrepid2_points = getView<double, DeviceType>("points", numPoints, 3);
    
    auto intrepid2_points_host = getHostCopy(intrepid2_points);
    
    for (int i=0; i<numPoints; i++)
    {
      const double z = ESEAS_points[i][2];
      const double x = ESEAS_points[i][0] * 2 - 1 + z;
      const double y = ESEAS_points[i][1] * 2 - 1 + z;
      
      intrepid2_points_host(i,0) = x;
      intrepid2_points_host(i,1) = y;
      intrepid2_points_host(i,2) = z;
      
      // use Intrepid2's transformation function to convert back (a check on the transformation we just did above).
      double x_eseas, y_eseas, z_eseas;
      transformToESEASPyramid(x_eseas, y_eseas, z_eseas, x, y, z);
      TEST_FLOATING_EQUALITY(ESEAS_points[i][0], x_eseas, relTol);
      TEST_FLOATING_EQUALITY(ESEAS_points[i][1], y_eseas, relTol);
      TEST_FLOATING_EQUALITY(ESEAS_points[i][2], z_eseas, relTol);
    }
    Kokkos::deep_copy(intrepid2_points, intrepid2_points_host);
    
    const int basisCardinality = hierarchicalBasis.getCardinality();
    ViewType<double,DeviceType> expectedValues = getView<double, DeviceType>("expected values", basisCardinality, numPoints);
    auto expectedValuesHost = getHostCopy(expectedValues);
    
    for (int basisOrdinal=0; basisOrdinal<basisCardinality; basisOrdinal++)
    {
      for (int pointOrdinal=0; pointOrdinal<numPoints; pointOrdinal++)
      {
        // expect basis to be nodal at vertices
        expectedValuesHost(basisOrdinal,pointOrdinal) = (basisOrdinal == pointOrdinal) ? 1 : 0;
      }
    }
    Kokkos::deep_copy(expectedValues, expectedValuesHost);
    
    ViewType<double,DeviceType> actualValues   = getView<double, DeviceType>("expected values", basisCardinality, numPoints);
    hierarchicalBasis.getValues(actualValues, intrepid2_points);
    
    testViewFloatingEquality(actualValues, expectedValues, relTol, absTol, out, success, "actual", "expected");
  }

  TEUCHOS_UNIT_TEST( BasisValues, PyramidHierarchical_HGRAD_PointTest )
  {
    // Test to compute values at a given point, for debugging purposes
    using HierarchicalBasis = HierarchicalBasisFamily<DefaultTestDeviceType>::HGRAD_PYR;
    using DeviceType = DefaultTestDeviceType;
    
    const ordinal_type p = 9;
    
    HierarchicalBasis hierarchicalBasis(p);
    
    using std::vector;
    const int numPoints = 1;
    vector< vector<double> > ESEAS_points(numPoints, vector<double>(3));
    
    ESEAS_points[0] = {0,0,0};
    
    ViewType<double,DeviceType> intrepid2_points = getView<double, DeviceType>("points", numPoints, 3);
    
    auto intrepid2_points_host = getHostCopy(intrepid2_points);
    
    for (int i=0; i<numPoints; i++)
    {
      const double z = ESEAS_points[i][2];
      const double x = ESEAS_points[i][0] * 2 - 1 + z;
      const double y = ESEAS_points[i][1] * 2 - 1 + z;
      intrepid2_points_host(i,0) = x;
      intrepid2_points_host(i,1) = y;
      intrepid2_points_host(i,2) = z;
    }
    Kokkos::deep_copy(intrepid2_points, intrepid2_points_host);
    
    const int basisCardinality = hierarchicalBasis.getCardinality();
    
    ViewType<double,DeviceType> actualValues   = getView<double, DeviceType>("actual values", basisCardinality, numPoints);
    hierarchicalBasis.getValues(actualValues, intrepid2_points, OPERATOR_VALUE);
    
    ViewType<double,DeviceType> actualGradients   = getView<double, DeviceType>("actual gradients", basisCardinality, numPoints, 3);
    hierarchicalBasis.getValues(actualGradients, intrepid2_points, OPERATOR_GRAD);
  }
} // namespace
