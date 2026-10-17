// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_ParametricElement.hpp>

#include <stk_topology/topology.hpp>

#include <array>
#include <cmath>
#include <vector>

namespace krino {

// Pure master-element math for ParametricElementHybrid (the public
// ParametricElement typedef).  Everything here is process-local, so the tests
// carry no parallel guards and never touch BulkData.  Kokkos/Intrepid2 is
// initialized once for the test executable by stk::initialize in
// krino/unit_main/Akri_Unit_main.cpp.
//
// Gold values are exact analytic invariants of the reference elements:
//   parametric volumes TET_4 = 1/6, TRI_3 = 1/2, QUAD_4 = 4;
//   shape functions partition unity and satisfy the Kronecker-delta property;
//   the Jacobian determinant of the reference element is 1 and scales as
//   s^dim under a uniform scaling; a linear field's recovered gradient is the
//   field's (constant) coefficients.

namespace {
const double exactTol = 1.e-9;
const double cubatureTol = 1.e-6;

// Pack per-node physical coordinates into the (numCoordDims, npe, nelem=1)
// layout expected by determinant()/gradient_operator():
//   coords[d + numCoordDims*(node + npe*elem)].
std::vector<double> pack_single_element_coords(const std::vector<std::array<double,3>> & nodeLocs, const int numCoordDims)
{
  const int npe = static_cast<int>(nodeLocs.size());
  std::vector<double> coords(numCoordDims * npe, 0.);
  for (int n = 0; n < npe; ++n)
    for (int d = 0; d < numCoordDims; ++d)
      coords[d + numCoordDims * n] = nodeLocs[n][d];
  return coords;
}

// Reference nodal coordinates of a parametric element as a vector of points,
// so they can be reused to build a "reference-shaped" physical element.
std::vector<std::array<double,3>> reference_node_locs(const ParametricElement & pe)
{
  const int dim = static_cast<int>(pe.topology_dimension());
  const int npe = static_cast<int>(pe.num_nodes());
  const double * refCoords = pe.nodal_parametric_coordinates();
  std::vector<std::array<double,3>> locs(npe, {{0.,0.,0.}});
  for (int n = 0; n < npe; ++n)
    for (int d = 0; d < dim; ++d)
      locs[n][d] = refCoords[n*dim + d];
  return locs;
}
}

//--------------------------------------------------------------------------
// Basic per-topology invariants
//--------------------------------------------------------------------------

TEST(ParametricElement, tet4_basicInvariants)
{
  const ParametricElement pe(stk::topology::TET_4);

  EXPECT_EQ(stk::topology::TET_4, pe.get_topology());
  EXPECT_EQ(3u, pe.topology_dimension());
  EXPECT_EQ(4u, pe.num_nodes());
  EXPECT_NEAR(1./6., pe.parametric_volume(), exactTol);

  const double gold[4][3] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1}};
  const double * refCoords = pe.nodal_parametric_coordinates();
  for (int n = 0; n < 4; ++n)
    for (int d = 0; d < 3; ++d)
      EXPECT_NEAR(gold[n][d], refCoords[n*3 + d], exactTol);
}

TEST(ParametricElement, tri3_basicInvariants)
{
  const ParametricElement pe(stk::topology::TRI_3_2D);

  EXPECT_EQ(2u, pe.topology_dimension());
  EXPECT_EQ(3u, pe.num_nodes());
  EXPECT_NEAR(0.5, pe.parametric_volume(), exactTol);

  const double gold[3][2] = {{0,0},{1,0},{0,1}};
  const double * refCoords = pe.nodal_parametric_coordinates();
  for (int n = 0; n < 3; ++n)
    for (int d = 0; d < 2; ++d)
      EXPECT_NEAR(gold[n][d], refCoords[n*2 + d], exactTol);
}

TEST(ParametricElement, quad4_basicInvariants)
{
  const ParametricElement pe(stk::topology::QUAD_4_2D);

  EXPECT_EQ(2u, pe.topology_dimension());
  EXPECT_EQ(4u, pe.num_nodes());
  EXPECT_NEAR(4.0, pe.parametric_volume(), exactTol);

  const double gold[4][2] = {{-1,-1},{1,-1},{1,1},{-1,1}};
  const double * refCoords = pe.nodal_parametric_coordinates();
  for (int n = 0; n < 4; ++n)
    for (int d = 0; d < 2; ++d)
      EXPECT_NEAR(gold[n][d], refCoords[n*2 + d], exactTol);
}

//--------------------------------------------------------------------------
// Integration weights integrate the reference cell exactly.
//--------------------------------------------------------------------------

TEST(ParametricElement, integrationWeightsSumToParametricVolume)
{
  for (const stk::topology topo : {stk::topology::TET_4, stk::topology::TRI_3_2D, stk::topology::QUAD_4_2D})
  {
    const ParametricElement pe(topo);
    double weightSum = 0.;
    const double * weights = pe.intg_weights();
    for (unsigned q = 0; q < pe.num_intg_pts(); ++q)
      weightSum += weights[q];
    EXPECT_NEAR(pe.parametric_volume(), weightSum, cubatureTol) << "topology " << topo.name();
  }
}

//--------------------------------------------------------------------------
// centroid_parametric_coordinates is the sum (not average) of nodal coords.
//--------------------------------------------------------------------------

TEST(ParametricElement, centroidIsSumOfNodalParametricCoordinates)
{
  const ParametricElement pe(stk::topology::TET_4);
  const int dim = static_cast<int>(pe.topology_dimension());
  const int npe = static_cast<int>(pe.num_nodes());
  const double * refCoords = pe.nodal_parametric_coordinates();
  const double * centroid = pe.centroid_parametric_coordinates();

  for (int d = 0; d < dim; ++d)
  {
    double sum = 0.;
    for (int n = 0; n < npe; ++n)
      sum += refCoords[n*dim + d];
    EXPECT_NEAR(sum, centroid[d], exactTol);
  }
}

//--------------------------------------------------------------------------
// Shape functions: partition of unity and the Kronecker-delta property.
//--------------------------------------------------------------------------

TEST(ParametricElement, shapeFunctionsPartitionUnityAndAreKroneckerDeltaAtNodes)
{
  for (const stk::topology topo : {stk::topology::TET_4, stk::topology::TRI_3_2D, stk::topology::QUAD_4_2D})
  {
    const ParametricElement pe(topo);
    const int dim = static_cast<int>(pe.topology_dimension());
    const int npe = static_cast<int>(pe.num_nodes());
    const double * refCoords = pe.nodal_parametric_coordinates();

    // Evaluate the npe shape functions at all npe nodal parametric coords in
    // one call; result is laid out result[node_of_eval*npe + basis].
    std::vector<double> shape(npe * npe, 0.);
    pe.shape_fcn(npe, refCoords, shape.data());

    for (int evalNode = 0; evalNode < npe; ++evalNode)
    {
      double sum = 0.;
      for (int basis = 0; basis < npe; ++basis)
      {
        const double val = shape[evalNode*npe + basis];
        sum += val;
        const double expected = (basis == evalNode) ? 1.0 : 0.0;
        EXPECT_NEAR(expected, val, exactTol) << "topology " << topo.name()
            << " N_" << basis << " at node " << evalNode;
      }
      EXPECT_NEAR(1.0, sum, exactTol) << "topology " << topo.name() << " partition of unity at node " << evalNode;
    }

    // Partition of unity also holds at an arbitrary interior point.
    std::vector<double> interior(dim, 0.);
    const double * centroid = pe.centroid_parametric_coordinates();
    for (int d = 0; d < dim; ++d)
      interior[d] = centroid[d] / npe; // average of nodal coords: safely interior
    std::vector<double> shapeInterior(npe, 0.);
    pe.shape_fcn(1, interior.data(), shapeInterior.data());
    double interiorSum = 0.;
    for (int basis = 0; basis < npe; ++basis)
      interiorSum += shapeInterior[basis];
    EXPECT_NEAR(1.0, interiorSum, exactTol) << "topology " << topo.name() << " partition of unity at interior point";
  }
}

//--------------------------------------------------------------------------
// determinant: reference element -> 1, uniform scaling -> s^dim, inverted -> error.
//--------------------------------------------------------------------------

TEST(ParametricElement, determinantOfReferenceElementIsOne)
{
  for (const stk::topology topo : {stk::topology::TET_4, stk::topology::TRI_3_2D, stk::topology::QUAD_4_2D})
  {
    const ParametricElement pe(topo);
    const int dim = static_cast<int>(pe.topology_dimension());
    const int nint = static_cast<int>(pe.num_intg_pts());

    const std::vector<double> coords = pack_single_element_coords(reference_node_locs(pe), dim);
    std::vector<double> detJ(nint, 0.);
    double error = -1.;
    pe.determinant(dim, 1, coords.data(), detJ.data(), &error);

    EXPECT_NEAR(0., error, exactTol) << "topology " << topo.name();
    for (int q = 0; q < nint; ++q)
      EXPECT_NEAR(1.0, detJ[q], exactTol) << "topology " << topo.name() << " ip " << q;
  }
}

TEST(ParametricElement, determinantScalesAsSToTheDim)
{
  const double s = 2.5;
  for (const stk::topology topo : {stk::topology::TET_4, stk::topology::TRI_3_2D, stk::topology::QUAD_4_2D})
  {
    const ParametricElement pe(topo);
    const int dim = static_cast<int>(pe.topology_dimension());
    const int nint = static_cast<int>(pe.num_intg_pts());

    std::vector<std::array<double,3>> locs = reference_node_locs(pe);
    for (auto & loc : locs)
      for (int d = 0; d < dim; ++d)
        loc[d] *= s;

    const std::vector<double> coords = pack_single_element_coords(locs, dim);
    std::vector<double> detJ(nint, 0.);
    double error = -1.;
    pe.determinant(dim, 1, coords.data(), detJ.data(), &error);

    const double gold = std::pow(s, dim);
    EXPECT_NEAR(0., error, exactTol) << "topology " << topo.name();
    for (int q = 0; q < nint; ++q)
      EXPECT_NEAR(gold, detJ[q], cubatureTol) << "topology " << topo.name() << " ip " << q;
  }
}

TEST(ParametricElement, determinantFlagsInvertedElement)
{
  const ParametricElement pe(stk::topology::TET_4);
  const int dim = static_cast<int>(pe.topology_dimension());
  const int nint = static_cast<int>(pe.num_intg_pts());

  // Swap two nodes of the reference tet to invert its orientation (det < 0).
  std::vector<std::array<double,3>> locs = reference_node_locs(pe);
  std::swap(locs[1], locs[2]);

  const std::vector<double> coords = pack_single_element_coords(locs, dim);
  std::vector<double> detJ(nint, 0.);
  double error = -1.;
  pe.determinant(dim, 1, coords.data(), detJ.data(), &error);

  EXPECT_NEAR(1.0, error, exactTol);
  for (int q = 0; q < nint; ++q)
    EXPECT_LT(detJ[q], 0.);
}

//--------------------------------------------------------------------------
// interpolate_point reproduces an affine map of the reference element.
//--------------------------------------------------------------------------

TEST(ParametricElement, interpolatePointReproducesAffineMap)
{
  const ParametricElement pe(stk::topology::TET_4);
  const int npe = static_cast<int>(pe.num_nodes());
  const double * refCoords = pe.nodal_parametric_coordinates();

  // A generic affine map phys = A*xi + t.
  auto affine = [](double x, double y, double z, double out[3]) {
    out[0] = 2.0*x + 0.5*y - 1.0*z + 1.0;
    out[1] = -1.0*x + 0.0*y + 3.0*z + 2.0;
    out[2] = 1.0*x + 1.0*y + 1.0*z - 0.5;
  };

  // Nodal physical coords as an (ncomp=3, num_nodes) field: field[3*node + comp].
  std::vector<double> field(3 * npe, 0.);
  for (int n = 0; n < npe; ++n)
    affine(refCoords[n*3+0], refCoords[n*3+1], refCoords[n*3+2], &field[3*n]);

  const double queryPar[3] = {0.25, 0.5, 0.1};
  double gold[3];
  affine(queryPar[0], queryPar[1], queryPar[2], gold);

  double result[3] = {0.,0.,0.};
  pe.interpolate_point(3, queryPar, 3, field.data(), result);

  for (int c = 0; c < 3; ++c)
    EXPECT_NEAR(gold[c], result[c], exactTol) << "component " << c;
}

//--------------------------------------------------------------------------
// LINE_3 shape-function derivatives match the analytic quadratic derivatives.
//
// Regression: the three end/mid-node derivatives always sum to zero (gradient
// partition of unity), so a bug that drops the constant terms is invisible to
// a sum check.  Pin the individual values instead.  For the LINE_3 shape
// functions N0 = -x(1-x)/2, N1 = x(1+x)/2, N2 = (1-x)(1+x), the derivatives are
// dN0/dx = x - 1/2, dN1/dx = x + 1/2, dN2/dx = -2x.
//--------------------------------------------------------------------------

TEST(ParametricElement, line3ShapeFunctionDerivativesMatchAnalytic)
{
  const ParametricElement pe(stk::topology::LINE_3);
  ASSERT_EQ(3u, pe.num_nodes());
  ASSERT_EQ(1u, pe.topology_dimension());

  // Evaluate at the two end nodes, the mid node, and an interior point.
  const std::vector<double> parCoords = {-1.0, 1.0, 0.0, 0.3};
  const int nint = static_cast<int>(parCoords.size());

  std::vector<double> deriv(nint * 3, 0.);
  pe.shape_fcn_deriv(nint, parCoords.data(), deriv.data());

  for (int ip = 0; ip < nint; ++ip)
  {
    const double x = parCoords[ip];
    EXPECT_NEAR(x - 0.5, deriv[ip*3 + 0], exactTol) << "N0' at x=" << x;
    EXPECT_NEAR(x + 0.5, deriv[ip*3 + 1], exactTol) << "N1' at x=" << x;
    EXPECT_NEAR(-2.0*x,  deriv[ip*3 + 2], exactTol) << "N2' at x=" << x;
  }
}

//--------------------------------------------------------------------------
// HEX_27 cubature is exact for the quartic polynomials a triquadratic basis
// produces.  The cubature order is 2*degree, so a basis mis-declared as degree
// 1 requests only order-2 cubature and under-integrates.  Integrating x^4 over
// the reference cube [-1,1]^3 (exact value 8/5) fails at order 2 and passes at
// order 4, so this pins the declared degree.
//--------------------------------------------------------------------------

TEST(ParametricElement, hex27CubatureIntegratesQuarticExactly)
{
  const ParametricElement pe(stk::topology::HEX_27);
  const int nint = static_cast<int>(pe.num_intg_pts());
  const double * weights = pe.intg_weights();
  const double * points = pe.intg_pt_locations();

  double integral = 0.;
  for (int q = 0; q < nint; ++q)
  {
    const double x = points[q*3 + 0];
    integral += weights[q] * x*x*x*x;
  }

  // integral of x^4 over [-1,1]^3 = (2/5)*2*2 = 8/5.
  EXPECT_NEAR(8./5., integral, cubatureTol);
}

//--------------------------------------------------------------------------
// gradient_operator + scalar_gradient recover the gradient of a linear field.
//--------------------------------------------------------------------------

TEST(ParametricElement, scalarGradientRecoversLinearFieldGradient)
{
  const ParametricElement pe(stk::topology::TET_4);
  const int dim = static_cast<int>(pe.topology_dimension());
  const int npe = static_cast<int>(pe.num_nodes());
  const int nint = static_cast<int>(pe.num_intg_pts());

  // Physical element == reference tet.
  const std::vector<std::array<double,3>> locs = reference_node_locs(pe);
  const std::vector<double> coords = pack_single_element_coords(locs, dim);

  // Linear field f(x,y,z) = a*x + b*y + c*z + d sampled at the nodes.
  const double a = 3.0, b = -2.0, c = 5.0, d = 1.5;
  std::vector<double> sfield(npe, 0.);
  for (int n = 0; n < npe; ++n)
    sfield[n] = a*locs[n][0] + b*locs[n][1] + c*locs[n][2] + d;

  // Reference-point shape-function gradients serve as both geometry and field
  // derivative operators (isoparametric).
  const double * deriv = pe.shape_fcn_deriv();

  std::vector<double> gradop(dim * npe * 1 * nint, 0.);
  std::vector<double> detJ(nint, 0.);
  double error = -1.;
  pe.gradient_operator(dim, nint, npe, deriv, npe, deriv, 1, coords.data(), gradop.data(), detJ.data(), &error);
  EXPECT_NEAR(0., error, exactTol);

  std::vector<double> grad(dim * nint, 0.);
  pe.scalar_gradient(nint, 1, gradop.data(), detJ.data(), sfield.data(), grad.data());

  const double gold[3] = {a, b, c};
  for (int q = 0; q < nint; ++q)
    for (int dd = 0; dd < dim; ++dd)
      EXPECT_NEAR(gold[dd], grad[q*dim + dd], cubatureTol) << "ip " << q << " dim " << dd;
}

}
