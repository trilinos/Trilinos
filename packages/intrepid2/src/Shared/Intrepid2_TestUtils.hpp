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
  
  if (cellTopo.getBaseKey() == shards::Line<2>::key)
  {
    basis = getLineBasis<BasisFamily>(fs,polyOrder_x);
  }
  else if (cellTopo.getBaseKey() == shards::Quadrilateral<4>::key)
  {
    basis = getQuadrilateralBasis<BasisFamily>(fs,polyOrder_x,polyOrder_y);
  }
  else if (cellTopo.getBaseKey() == shards::Hexahedron<>::key)
  {
    INTREPID2_TEST_FOR_EXCEPTION(polyOrder_z < 0, std::invalid_argument, "polyOrder_z must be specified");
    basis = getHexahedralBasis<BasisFamily>(fs,polyOrder_x,polyOrder_y,polyOrder_z);
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
  
  using BasisFamily = DerivedBasisFamily<LineBasisGrad, LineBasisVol>;
  
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

template <typename PointValueType>
inline ViewType<PointValueType> lineInputPointsView(int numPoints)
{
  ViewType<PointValueType> inputPoints = getView<PointValueType>("line input points",numPoints,1);
  Kokkos::parallel_for(numPoints, KOKKOS_LAMBDA(const int i)
  {
    double x_zero_one = i * 1.0 / (numPoints-1); // the x value on [0,1]
    inputPoints(i,0) = fromZeroOne(x_zero_one);   // standard Intrepid2 element
  });
  return inputPoints;
}

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
        inputPoints(pointOrdinal,0) = fromZeroOne(x_zero_one);   // standard Intrepid2 element
        inputPoints(pointOrdinal,1) = fromZeroOne(y_zero_one);   // standard Intrepid2 element
        inputPoints(pointOrdinal,2) = fromZeroOne(z_zero_one);   // standard Intrepid2 element
      }
    }
  });
  return inputPoints;
}

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
      inputPoints(pointOrdinal,0) = fromZeroOne(x_zero_one);   // standard Intrepid2 element
      inputPoints(pointOrdinal,1) = fromZeroOne(y_zero_one);   // standard Intrepid2 element
    }
  });
  return inputPoints;
}

template <typename PointValueType>
inline ViewType<PointValueType> getInputPointsView(shards::CellTopology &cellTopo, int numPoints_1D)
{
  if (cellTopo.getBaseKey() == shards::Line<2>::key)
  {
    return lineInputPointsView<PointValueType>(numPoints_1D);
  }
  else if (cellTopo.getBaseKey() == shards::Quadrilateral<4>::key)
  {
    return quadInputPointsView<PointValueType>(numPoints_1D);
  }
  else if (cellTopo.getBaseKey() == shards::Hexahedron<8>::key)
  {
    return hexInputPointsView<PointValueType>(numPoints_1D);
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
