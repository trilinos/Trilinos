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
// Questions? Contact Mauro Perego  (mperego@sandia.gov) or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file   StructuredIntegrationTests_Utils.hpp
    \brief  Helper methods for use in StructuredIntegrationTests.
    \author Nathan V. Roberts
*/

#include <Intrepid2_Types.hpp>
#include <Intrepid2_CellGeometry.hpp>

#include "GRADGRADStandardAssembly.hpp"
#include "GRADGRADStructuredAssembly.hpp"
#include "H1StandardAssembly.hpp"
#include "H1StructuredAssembly.hpp"
#include "HDIVStandardAssembly.hpp"
#include "HDIVStructuredAssembly.hpp"
#include "HCURLStandardAssembly.hpp"
#include "HCURLStructuredAssembly.hpp"
#include "HVOLStandardAssembly.hpp"
#include "HVOLStructuredAssembly.hpp"

template< typename PointScalar, int spaceDim, typename DeviceType >
inline
CellGeometry<PointScalar, spaceDim, DeviceType> getMesh(AlgorithmChoice algorithmChoice, const Kokkos::Array<int,spaceDim> &gridCellCounts)
{
  Kokkos::Array<PointScalar,spaceDim> domainExtents;
  for (int d=0; d<spaceDim; d++)
  {
    domainExtents[d]  = 1.0;
  }
  auto uniformTensorGeometry = uniformCartesianMesh<PointScalar,spaceDim,DeviceType>(domainExtents, gridCellCounts);
  
  switch (algorithmChoice)
  {
    case Standard:
    case NonAffineTensor:
    {
      // Standard and non-affine tensor use the same geometry; the difference is how this is used in assembly
      const bool copyAffineness = false;
      auto genericGeometry = getNodalCellGeometry(uniformTensorGeometry, copyAffineness);
      return genericGeometry;
    }
    case Uniform:
      return uniformTensorGeometry;
    case AffineNonTensor:
    case AffineTensor:
    {
      const bool copyAffineness = true;
      auto affineNonTensorGeometry = getNodalCellGeometry(uniformTensorGeometry, copyAffineness);
      return affineNonTensorGeometry;
    }
    case DiagonalJacobian:
    {
      INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "DiagonalJacobian case not yet implemented");
    }
  }
  return uniformTensorGeometry; // this line should be unreachable; included to avoid compiler warnings from nvcc
}

template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType>
Intrepid2::ScalarView<Scalar,DeviceType> performStandardQuadrature(FormulationChoice formulation,
                                        Intrepid2::CellGeometry<PointScalar, spaceDim, DeviceType> &geometry, const int &polyOrder, const int &worksetSize,
                                        double &transformIntegrateFlopCount, double &jacobianCellMeasureFlopCount)
{
  switch (formulation)
  {
    case Poisson:
      return performStandardQuadratureGRADGRAD<Scalar,BasisFamily>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
    case Hgrad:
      return performStandardQuadratureH1<Scalar, BasisFamily>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
    case Hdiv:
      return performStandardQuadratureHDIV<Scalar, BasisFamily>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
    case Hcurl:
      return performStandardQuadratureHCURL<Scalar, BasisFamily>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
    case L2:
      return performStandardQuadratureHVOL<Scalar, BasisFamily>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
    default:
      INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported formulation");
  }
}

template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType>
Intrepid2::ScalarView<Scalar,DeviceType> performStructuredQuadrature(FormulationChoice formulation,
                                          Intrepid2::CellGeometry<PointScalar, spaceDim, DeviceType> &geometry, const int &polyOrder, const int &worksetSize,
                                          double &transformIntegrateFlopCount, double &jacobianCellMeasureFlopCount)
{
  switch (formulation)
  {
    case Poisson:
      return performStructuredQuadratureGRADGRAD<Scalar,BasisFamily>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
    case Hgrad:
      return performStructuredQuadratureH1<Scalar, BasisFamily>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
    case Hdiv:
      return performStructuredQuadratureHDIV<Scalar, BasisFamily>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
    case Hcurl:
      return performStructuredQuadratureHCURL<Scalar, BasisFamily>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
    case L2:
      return performStructuredQuadratureHVOL<Scalar, BasisFamily>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
    default:
      INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported formulation");
  }
}

//! version of integrate that performs a standard integration for affine meshes; does not take advantage of the tensor product structure at all
//! this version can be used to verify correctness of other versions
template<class Scalar, typename DeviceType>
void integrate_baseline(Data<Scalar,DeviceType> integrals, const TransformedBasisValues<Scalar,DeviceType> vectorDataLeft,
                      const TensorData<Scalar,DeviceType> cellMeasures, const TransformedBasisValues<Scalar,DeviceType> vectorDataRight)
{
const int spaceDim       = vectorDataLeft.spaceDim();

// use the CFPD operator() provided by the vector data objects; don't take advantage of tensor product structure at all
const int numPoints = vectorDataLeft.numPoints();

INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(vectorDataLeft.rank()  != 4, std::invalid_argument, "vectorDataLeft must be of shape (C,F,P,D)");
INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(vectorDataRight.rank() != 4, std::invalid_argument, "vectorDataRight must be of shape (C,F,P,D)");
INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(cellMeasures.rank()    != 2, std::invalid_argument, "cellMeasures must be of shape (C,P)");

INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(vectorDataLeft.spaceDim() != vectorDataRight.spaceDim(), std::invalid_argument, "vectorDataLeft and vectorDataRight must agree on the spatial dimension");

INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(vectorDataRight.extent_int(2) != vectorDataLeft.extent_int(2), std::invalid_argument, "vectorData point dimensions must match");
INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(cellMeasures.extent_int(1) != vectorDataLeft.extent_int(2), std::invalid_argument,
                                         std::string("cellMeasures point dimension (" + std::to_string(cellMeasures.extent_int(1)) +
                                                     ") must match vectorData point dimension (" + std::to_string(vectorDataLeft.extent_int(2)) + ")").c_str());

//  printFunctor4(vectorDataLeft, std::cout, "vectorDataLeft");
//  printFunctor2(cellMeasures, std::cout, "cellMeasures");

// integral data may have shape (C,F1,F2) or (if the variation type is CONSTANT in the cell dimension) shape (F1,F2)
const int integralViewRank = integrals.getUnderlyingViewRank();

using ExecutionSpace = typename DeviceType::execution_space;
auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{integrals.getDataExtent(0),vectorDataLeft.numFields(),vectorDataRight.numFields()});
Kokkos::parallel_for("fill expanded cell nodes", policy,
KOKKOS_LAMBDA (const int &cellDataOrdinal, const int &fieldOrdinalLeft, const int &fieldOrdinalRight)
{
  Scalar integral=0;
  for (int pointOrdinal=0; pointOrdinal<numPoints; pointOrdinal++)
  {
    Scalar pointSum = 0.0;
    for (int d=0; d<spaceDim; d++)
    {
      pointSum += vectorDataLeft(cellDataOrdinal,fieldOrdinalLeft,pointOrdinal,d) * vectorDataRight(cellDataOrdinal,fieldOrdinalRight,pointOrdinal,d);
    }
    integral += pointSum * cellMeasures(cellDataOrdinal,pointOrdinal);
  }
  if (integralViewRank == 3)
  {
    // shape (C,F1,F2)
    auto integralView = integrals.getUnderlyingView3();
    integralView(cellDataOrdinal,fieldOrdinalLeft,fieldOrdinalRight) = integral;
  }
  else
  {
    // shape (F1,F2)
    auto integralView = integrals.getUnderlyingView2();
    integralView(fieldOrdinalLeft,fieldOrdinalRight) = integral;
  }
});
}

template<class Scalar, typename DeviceType>
void testIntegrateMatchesBaseline(const TransformedBasisValues<Scalar,DeviceType> vectorDataLeft,
                                  const TensorData<Scalar,DeviceType> cellMeasures, const TransformedBasisValues<Scalar,DeviceType> vectorDataRight,
                                  Teuchos::FancyOStream &out, bool &success)
{
  const double relTol = 1e-12;
  const double absTol = 1e-12;
  
  using IntegrationTools = Intrepid2::IntegrationTools<DeviceType>;
  
  auto integralsBaseline  = IntegrationTools::allocateIntegralData(vectorDataLeft, cellMeasures, vectorDataRight);
  auto integralsIntegrate = IntegrationTools::allocateIntegralData(vectorDataLeft, cellMeasures, vectorDataRight);
  
  integrate_baseline(integralsBaseline, vectorDataLeft, cellMeasures, vectorDataRight);
  IntegrationTools::integrate(integralsIntegrate, vectorDataLeft, cellMeasures, vectorDataRight);
  
  const int integralsBaselineViewRank = integralsBaseline.getUnderlyingViewRank();
  
  printFunctor3( integralsBaseline, out, "integralsBaseline");
  printFunctor3(integralsIntegrate, out, "integralsIntegrate");
//    printView(integralsBaselineView, out);
  
  if (integralsBaselineViewRank == 3)
  {
    auto integralsBaselineView  = integralsBaseline.getUnderlyingView3();
    auto integralsIntegrateView = integralsIntegrate.getUnderlyingView3();
    
    testFloatingEquality3(integralsBaselineView, integralsIntegrateView, relTol, absTol, out, success, "baseline integral", "sum factorized integral");
  }
  else
  {
    auto integralsBaselineView  = integralsBaseline.getUnderlyingView2();
    auto integralsIntegrateView = integralsIntegrate.getUnderlyingView2();
    
    testFloatingEquality2(integralsBaselineView, integralsIntegrateView, relTol, absTol, out, success, "baseline integral", "sum factorized integral");
  }
}
