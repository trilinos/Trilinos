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

/** \file   Intrepid2_TestUtils.hpp
    \brief  Utility methods for Intrepid2 unit tests.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_TestUtils_h
#define Intrepid2_TestUtils_h

#include "Kokkos_Core.hpp"
#include "Kokkos_DynRankView.hpp"

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_DerivedBasisFamily.hpp"
#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_Sacado.hpp" // Sacado includes, guarded by the appropriate preprocessor variable
#include "Intrepid2_Utils.hpp"

#include "Teuchos_UnitTestHarness.hpp"

static const double TEST_TOLERANCE_TIGHT = 1.e2 * std::numeric_limits<double>::epsilon();

// we use DynRankView for both input points and values
template<typename ScalarType>
using ViewType = Kokkos::DynRankView<ScalarType,Kokkos::DefaultExecutionSpace>;

template<typename ScalarType>
inline bool valuesAreSmall(ScalarType a, ScalarType b, double epsilon)
{
  using std::abs;
  return (abs(a) < epsilon) && (abs(b) < epsilon);
}

inline bool approximatelyEqual(double a, double b, double epsilon)
{
  const double larger_magnitude = (std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a));
  return std::abs(a - b) <= larger_magnitude * epsilon;
}

inline bool essentiallyEqual(double a, double b, double epsilon)
{
  const double smaller_magnitude = (std::abs(a) > std::abs(b) ? std::abs(b) : std::abs(a));
  return std::abs(a - b) <= smaller_magnitude * epsilon;
}

// conversion from the ref element [0,1] (used by ESEAS) to the ref element [-1,1] (used by Intrepid2)
KOKKOS_INLINE_FUNCTION double fromZeroOne(double x_zero_one)
{
  return x_zero_one * 2.0 - 1.0;
}

// conversion from the ref element [-1,1] (used by Intrepid2) to the ref element [0,1] (used by ESEAS)
KOKKOS_INLINE_FUNCTION double toZeroOne(double x_minus_one_one)
{
  return (x_minus_one_one + 1.0) / 2.0;
}

// conversion from the ref element [0,1] (used by ESEAS) to the ref element [-1,1] (used by Intrepid2)
KOKKOS_INLINE_FUNCTION double fromZeroOne_dx(double dx_zero_one)
{
  return dx_zero_one / 2.0;
}

// conversion from the ref element [-1,1] (used by Intrepid2) to the ref element [0,1] (used by ESEAS)
KOKKOS_INLINE_FUNCTION double toZeroOne_dx(double dx_minus_one_one)
{
  return dx_minus_one_one * 2.0;
}

template<class DeviceViewType>
typename DeviceViewType::HostMirror getHostCopy(const DeviceViewType &deviceView)
{
  typename DeviceViewType::HostMirror hostView = Kokkos::create_mirror(deviceView);
  Kokkos::deep_copy(hostView, deviceView);
  return hostView;
}

template<class BasisFamily>
inline Teuchos::RCP< Intrepid2::Basis<Kokkos::DefaultExecutionSpace,double,double> > getBasisUsingFamily(shards::CellTopology cellTopo, Intrepid2::EFunctionSpace fs,
                                                                                                         int polyOrder_x, int polyOrder_y=-1, int polyOrder_z = -1)
{
  using BasisPtr = typename BasisFamily::BasisPtr;
  
  BasisPtr basis;
  using namespace Intrepid2;
  
  if (cellTopo.getBaseKey() == shards::Line<>::key)
  {
    basis = getLineBasis<BasisFamily>(fs,polyOrder_x);
  }
  else if (cellTopo.getBaseKey() == shards::Quadrilateral<>::key)
  {
    INTREPID2_TEST_FOR_EXCEPTION(polyOrder_y < 0, std::invalid_argument, "polyOrder_y must be specified");
    basis = getQuadrilateralBasis<BasisFamily>(fs,polyOrder_x,polyOrder_y);
  }
  else if (cellTopo.getBaseKey() == shards::Triangle<>::key)
  {
    basis = getTriangleBasis<BasisFamily>(fs,polyOrder_x);
  }
  else if (cellTopo.getBaseKey() == shards::Hexahedron<>::key)
  {
    INTREPID2_TEST_FOR_EXCEPTION(polyOrder_y < 0, std::invalid_argument, "polyOrder_y must be specified");
    INTREPID2_TEST_FOR_EXCEPTION(polyOrder_z < 0, std::invalid_argument, "polyOrder_z must be specified");
    basis = getHexahedronBasis<BasisFamily>(fs,polyOrder_x,polyOrder_y,polyOrder_z);
  }
  else if (cellTopo.getBaseKey() == shards::Tetrahedron<>::key)
  {
    basis = getTetrahedronBasis<BasisFamily>(fs, polyOrder_x);
  }
  else
  {
    INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported cell topology");
  }
  return basis;
}

template<bool defineVertexFunctions>
inline Teuchos::RCP< Intrepid2::Basis<Kokkos::DefaultExecutionSpace,double,double> > getHierarchicalBasis(shards::CellTopology cellTopo, Intrepid2::EFunctionSpace fs,
                                                                                                          int polyOrder_x, int polyOrder_y=-1, int polyOrder_z = -1)
{
  using ExecSpace = Kokkos::DefaultExecutionSpace;
  using Scalar = double;
  using namespace Intrepid2;
  
  using LineBasisGrad = Intrepid2::IntegratedLegendreBasis_HGRAD_LINE<ExecSpace, Scalar, Scalar, defineVertexFunctions, true>;
  using LineBasisVol  = Intrepid2::LegendreBasis_HVOL_LINE< ExecSpace, Scalar, Scalar>;
  using TriangleBasisFamily = Intrepid2::HierarchicalTriangleBasisFamily<ExecSpace, Scalar, Scalar, defineVertexFunctions>;
  using TetrahedronBasisFamily = Intrepid2::HierarchicalTetrahedronBasisFamily<ExecSpace, Scalar, Scalar, defineVertexFunctions>;
  
  using BasisFamily = DerivedBasisFamily<LineBasisGrad, LineBasisVol, TriangleBasisFamily, TetrahedronBasisFamily>;
  
  return getBasisUsingFamily<BasisFamily>(cellTopo, fs, polyOrder_x, polyOrder_y, polyOrder_z);
}

template<typename ValueType, class ... DimArgs>
inline ViewType<ValueType> getView(const std::string &label, DimArgs... dims)
{
  const bool allocateFadStorage = !std::is_pod<ValueType>::value;
  constexpr int maxDerivatives = 3; // maximum derivatives used for fad types in unit tests
  if (!allocateFadStorage)
  {
    return ViewType<ValueType>(label,dims...);
  }
  else
  {
    return ViewType<ValueType>(label,dims...,maxDerivatives+1);
  }
}

/** \brief Returns a View containing equispaced points on the line.
 \param [in] numPointsBase - the number of points that will be defined along each edge.
 */
template <typename PointValueType>
inline ViewType<PointValueType> lineInputPointsView(int numPoints)
{
  ViewType<PointValueType> inputPoints = getView<PointValueType>("line input points",numPoints,1);
  Kokkos::parallel_for(numPoints, KOKKOS_LAMBDA(const int i)
  {
    double x_zero_one = i * 1.0 / (numPoints-1); // the x value on [0,1]
    inputPoints(i,0) = PointValueType(fromZeroOne(x_zero_one));   // standard Intrepid2 element
  });
  return inputPoints;
}

/** \brief Returns a View containing equispaced points on the hexahedron.
 \param [in] numPointsBase - the number of points that will be defined along each edge.
 
 The total number of points defined will be a cubic number; if n=numPointsBase, then the point count is n^3.
 */
template <typename PointValueType>
inline ViewType<PointValueType> hexInputPointsView(int numPoints_1D)
{
  ViewType<PointValueType> inputPoints = getView<PointValueType>("hex input points",numPoints_1D*numPoints_1D*numPoints_1D,3);
  Kokkos::parallel_for(numPoints_1D, KOKKOS_LAMBDA(const int i)
  {
    double x_zero_one = i * 1.0 / (numPoints_1D-1); // the x value on [0,1]
    for (int j=0; j<numPoints_1D; j++)
    {
      double y_zero_one = j * 1.0 / (numPoints_1D-1); // the y value on [0,1]
      for (int k=0; k<numPoints_1D; k++)
      {
        double z_zero_one = k * 1.0 / (numPoints_1D-1); // the z value on [0,1]
        int pointOrdinal = (i*numPoints_1D+j)*numPoints_1D+k;
        inputPoints(pointOrdinal,0) = PointValueType(fromZeroOne(x_zero_one));   // standard Intrepid2 element
        inputPoints(pointOrdinal,1) = PointValueType(fromZeroOne(y_zero_one));   // standard Intrepid2 element
        inputPoints(pointOrdinal,2) = PointValueType(fromZeroOne(z_zero_one));   // standard Intrepid2 element
      }
    }
  });
  return inputPoints;
}

/** \brief Returns a View containing equispaced points on the quadrilateral.
 \param [in] numPointsBase - the number of points that will be defined along each edge.
 
 The total number of points defined will be a square number; if n=numPointsBase, then the point count is n^2.
 */
template <typename PointValueType>
inline ViewType<PointValueType> quadInputPointsView(int numPoints_1D)
{
  ViewType<PointValueType> inputPoints = getView<PointValueType>("quad input points",numPoints_1D*numPoints_1D,2);
  
  Kokkos::parallel_for(numPoints_1D, KOKKOS_LAMBDA(const int i)
  {
    double x_zero_one = i * 1.0 / (numPoints_1D-1); // the x value on [0,1]
    for (int j=0; j<numPoints_1D; j++)
    {
      double y_zero_one = j * 1.0 / (numPoints_1D-1); // the y value on [0,1]
      int pointOrdinal = i*numPoints_1D+j;
      inputPoints(pointOrdinal,0) = PointValueType(fromZeroOne(x_zero_one));   // standard Intrepid2 element
      inputPoints(pointOrdinal,1) = PointValueType(fromZeroOne(y_zero_one));   // standard Intrepid2 element
    }
  });
  return inputPoints;
}

/** \brief Returns a View containing regularly-spaced points on the tetrahedron.
 \param [in] numPointsBase - the number of points that will be defined along each edge.
 
 The total number of points defined will be a tetrahedral number; if n=numPointsBase, then the point count is the nth tetrahedral number, given by n*(n+1)*(n+2)/6.
 */
template <typename PointValueType>
inline ViewType<PointValueType> tetInputPointsView(int numPointsBase)
{
  const int numPoints = numPointsBase*(numPointsBase+1)*(numPointsBase+2)/6;
  ViewType<PointValueType> inputPoints = getView<PointValueType>("tetrahedron input points",numPoints,3);
  Kokkos::parallel_for(numPointsBase, KOKKOS_LAMBDA(const int d0) // d0 generalizes row
  {
    // the following might be the formula for the nested for loop below, but we need to check this
    // for now, we comment this out and do the clearer thing that requires more computation
//    const int n = numPointsBase-d0;
//    const int pointOrdinalOffset = n*(n+1)*(n+2)/6;
    int pointOrdinalOffset = 0;
    for (int i=0; i<d0; i++)
    {
      for (int j=0; j<numPointsBase-i; j++)
      {
        pointOrdinalOffset += (numPointsBase - i - j);
      }
    }
    
    double z_zero_one = d0 * 1.0 / (numPointsBase-1); // z value on [0,1]
    int pointOrdinal = pointOrdinalOffset;
    for (int d1=0; d1<numPointsBase-d0; d1++) // d1 generalizes column
    {
      double y_zero_one = d1 * 1.0 / (numPointsBase-1); // the y value on [0,1]
      for (int d2=0; d2<numPointsBase-d0-d1; d2++)
      {
        double x_zero_one = d2 * 1.0 / (numPointsBase-1); // the x value on [0,1]
        
        inputPoints(pointOrdinal,0) = PointValueType(x_zero_one);
        inputPoints(pointOrdinal,1) = PointValueType(y_zero_one);
        inputPoints(pointOrdinal,2) = PointValueType(z_zero_one);
        
        pointOrdinal++;
      }
    }
  });
  return inputPoints;
}

/** \brief Returns a View containing regularly-spaced points on the triangle.
 \param [in] numPointsBase - the number of points that will be defined along each edge.
 
 The total number of points defined will be a triangular number; if n=numPointsBase, then the point count is the nth triangular number, given by n*(n+1)/2.
 */
template <typename PointValueType>
inline ViewType<PointValueType> triInputPointsView(int numPointsBase)
{
  const int numPoints = numPointsBase*(numPointsBase+1)/2;
  ViewType<PointValueType> inputPoints = getView<PointValueType>("triangle input points",numPoints,2);
  Kokkos::parallel_for(numPointsBase, KOKKOS_LAMBDA(const int row)
  {
    int rowPointOrdinalOffset = 0;
    for (int i=0; i<row; i++)
    {
      rowPointOrdinalOffset += (numPointsBase - i);
    }
    double y_zero_one = row * 1.0 / (numPointsBase-1); // y value on [0,1]
    for (int col=0; col<numPointsBase-row; col++)
    {
      const int pointOrdinal = rowPointOrdinalOffset + col;
      double x_zero_one = col * 1.0 / (numPointsBase-1); // the x value on [0,1]
      inputPoints(pointOrdinal,0) = PointValueType(x_zero_one);
      inputPoints(pointOrdinal,1) = PointValueType(y_zero_one);
    }
  });
  return inputPoints;
}

template <typename PointValueType>
inline ViewType<PointValueType> getInputPointsView(shards::CellTopology &cellTopo, int numPoints_1D)
{
  if (cellTopo.getBaseKey() == shards::Line<>::key)
  {
    return lineInputPointsView<PointValueType>(numPoints_1D);
  }
  else if (cellTopo.getBaseKey() == shards::Quadrilateral<>::key)
  {
    return quadInputPointsView<PointValueType>(numPoints_1D);
  }
  else if (cellTopo.getBaseKey() == shards::Hexahedron<>::key)
  {
    return hexInputPointsView<PointValueType>(numPoints_1D);
  }
  else if (cellTopo.getBaseKey() == shards::Triangle<>::key)
  {
    return triInputPointsView<PointValueType>(numPoints_1D);
  }
  else if (cellTopo.getBaseKey() == shards::Tetrahedron<>::key)
  {
    return tetInputPointsView<PointValueType>(numPoints_1D);
  }
  else
  {
    INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported Cell Topology");
  }
}

template<typename OutputValueType>
inline ViewType<OutputValueType> getOutputView(Intrepid2::EFunctionSpace fs, Intrepid2::EOperator op, int basisCardinality, int numPoints, int spaceDim)
{
  switch (fs) {
    case Intrepid2::FUNCTION_SPACE_HGRAD:
      switch (op) {
        case Intrepid2::OPERATOR_VALUE:
          return getView<OutputValueType>("H^1 value output",basisCardinality,numPoints);
        case Intrepid2::OPERATOR_GRAD:
          return getView<OutputValueType>("H^1 derivative output",basisCardinality,numPoints,spaceDim);
        case Intrepid2::OPERATOR_D1:
        case Intrepid2::OPERATOR_D2:
        case Intrepid2::OPERATOR_D3:
        case Intrepid2::OPERATOR_D4:
        case Intrepid2::OPERATOR_D5:
        case Intrepid2::OPERATOR_D6:
        case Intrepid2::OPERATOR_D7:
        case Intrepid2::OPERATOR_D8:
        case Intrepid2::OPERATOR_D9:
        case Intrepid2::OPERATOR_D10:
        {
          const auto dkcard = getDkCardinality(op, spaceDim);
          return getView<OutputValueType>("H^1 derivative output",basisCardinality,numPoints,dkcard);
        }
        default:
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported op/fs combination");
      }
    case Intrepid2::FUNCTION_SPACE_HCURL:
      switch (op) {
        case Intrepid2::OPERATOR_VALUE:
          return getView<OutputValueType>("H(curl) value output",basisCardinality,numPoints,spaceDim);
        case Intrepid2::OPERATOR_CURL:
          if (spaceDim == 2)
            return getView<OutputValueType>("H(curl) derivative output",basisCardinality,numPoints);
          else if (spaceDim == 3)
            return getView<OutputValueType>("H(curl) derivative output",basisCardinality,numPoints,spaceDim);
        default:
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported op/fs combination");
      }
    case Intrepid2::FUNCTION_SPACE_HDIV:
      switch (op) {
        case Intrepid2::OPERATOR_VALUE:
          return getView<OutputValueType>("H(div) value output",basisCardinality,numPoints,spaceDim);
        case Intrepid2::OPERATOR_DIV:
          return getView<OutputValueType>("H(div) derivative output",basisCardinality,numPoints);
        default:
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported op/fs combination");
      }
      
    case Intrepid2::FUNCTION_SPACE_HVOL:
      switch (op) {
        case Intrepid2::OPERATOR_VALUE:
          return getView<OutputValueType>("H(vol) value output",basisCardinality,numPoints);
        case Intrepid2::OPERATOR_D1:
        case Intrepid2::OPERATOR_D2:
        case Intrepid2::OPERATOR_D3:
        case Intrepid2::OPERATOR_D4:
        case Intrepid2::OPERATOR_D5:
        case Intrepid2::OPERATOR_D6:
        case Intrepid2::OPERATOR_D7:
        case Intrepid2::OPERATOR_D8:
        case Intrepid2::OPERATOR_D9:
        case Intrepid2::OPERATOR_D10:
        {
          const auto dkcard = getDkCardinality(op, spaceDim);
          return getView<OutputValueType>("H(vol) derivative output",basisCardinality,numPoints,dkcard);
        }
        default:
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported op/fs combination");
      }
    default:
      INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported op/fs combination");
  }
}

// ! This returns a vector whose entries are vector<int>s containing 1-3 polynomial orders from 1 up to and including those specified
// ! Intended for testing bases that support anisotropic polynomial degree, such as the hierarchical bases
inline std::vector< std::vector<int> > getBasisTestCasesUpToDegree(int spaceDim, int minDegree, int polyOrder_x, int polyOrder_y=-1, int polyOrder_z = -1)
{
  std::vector<int> degrees(spaceDim);
  degrees[0] = polyOrder_x;
  if (spaceDim > 1) degrees[1] = polyOrder_y;
  if (spaceDim > 2) degrees[2] = polyOrder_z;
  
  int numCases = degrees[0];
  for (unsigned d=1; d<degrees.size(); d++)
  {
    numCases = numCases * (degrees[d] + 1 - minDegree);
  }
  std::vector< std::vector<int> > subBasisDegreeTestCases(numCases);
  for (int caseOrdinal=0; caseOrdinal<numCases; caseOrdinal++)
  {
    std::vector<int> subBasisDegrees(degrees.size());
    int caseRemainder = caseOrdinal;
    for (int d=degrees.size()-1; d>=0; d--)
    {
      int subBasisDegree = caseRemainder % (degrees[d] + 1 - minDegree);
      caseRemainder = caseRemainder / (degrees[d] + 1 - minDegree);
      subBasisDegrees[d] = subBasisDegree + minDegree;
    }
    subBasisDegreeTestCases[caseOrdinal] = subBasisDegrees;
  }
  return subBasisDegreeTestCases;
}

template <class ViewType>
void testViewFloatingEquality(ViewType &view1, ViewType &view2, double relTol, double absTol, Teuchos::FancyOStream &out, bool &success,
                              std::string view1Name = "View 1", std::string view2Name = "View 2")
{
  using Scalar   = typename ViewType::value_type;
  using HostView = typename ViewType::HostMirror;
  using HostViewIteratorScalar = Intrepid2::ViewIterator<HostView, Scalar>;
  
  auto hostView1 = getHostCopy(view1);
  auto hostView2 = getHostCopy(view2);
  
  // check that rank/size match
  TEUCHOS_TEST_FOR_EXCEPTION(view1.rank() != view2.rank(), std::invalid_argument, "views must agree in rank");
  for (unsigned i=0; i<view1.rank(); i++)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(view1.extent_int(i) != view2.extent_int(i), std::invalid_argument, "views must agree in size in each dimension");
  }
  
  if (view1.size() == 0) return; // nothing to test
  
  HostViewIteratorScalar vi1(hostView1);
  HostViewIteratorScalar vi2(hostView2);
  
  double maxDiff = 0.0;
  double maxRelativeDiff = 0.0;
  bool differencesFound = false;
  
  bool moreEntries = (vi1.nextIncrementRank() != -1) && (vi2.nextIncrementRank() != -1);
  while (moreEntries)
  {
    Scalar value1 = vi1.get();
    Scalar value2 = vi2.get();
    
    const bool small = valuesAreSmall(value1, value2, absTol);
    const bool valuesMatch = small || essentiallyEqual(value1, value2, relTol);
    
    if (!valuesMatch)
    {
      differencesFound = true;
      success = false;
      auto location = vi1.getLocation();
      out << "At location (";
      for (unsigned i=0; i<view1.rank(); i++)
      {
        out << location[i];
        if (i<view1.rank()-1)
        {
          out << ",";
        }
      }
      out << ") " << view1Name << " value != " << view2Name << " value (";
      out << value1 << " != " << value2 << "); diff is " << std::abs(value1-value2) << std::endl;
      
      maxDiff = std::max(maxDiff, std::abs(value1-value2));
      double smallerMagnitude = (std::abs(value1) > std::abs(value2) ? std::abs(value2) : std::abs(value1));
      if (smallerMagnitude > 0)
      {
        const double relativeDiff = std::abs(value1 - value2) / smallerMagnitude;
        maxRelativeDiff = std::max(maxRelativeDiff, relativeDiff);
      }
    }
    
    moreEntries =                (vi1.increment() != -1);
    moreEntries = moreEntries && (vi2.increment() != -1);
  }
  
  if (differencesFound)
  {
    out << "max difference = " << maxDiff << std::endl;
    out << "max relative difference = " << maxRelativeDiff << std::endl;
  }
}

#ifdef HAVE_INTREPID2_SACADO
#include <Sacado.hpp>
using Sacado_Fad_DFadType = Sacado::Fad::DFad<double>;
#define INTREPID2_OUTPUTSCALAR_POINTSCALAR_TEST_INSTANT(GROUP_NAME, TEST_NAME) \
\
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( GROUP_NAME, TEST_NAME, double, double ) \
\
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( GROUP_NAME, TEST_NAME, Sacado_Fad_DFadType, Sacado_Fad_DFadType ) \

#else
#define INTREPID2_OUTPUTSCALAR_POINTSCALAR_TEST_INSTANT(GROUP_NAME, TEST_NAME) \
\
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( GROUP_NAME, TEST_NAME, double, double ) \

#endif

#endif /* Intrepid2_TestUtils_h */
