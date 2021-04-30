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

/** \file   StructuredIntegrationPerformance.cpp
    \brief  Main for performance tests comparing structured integration performance to standard.
 */

#include "Teuchos_GlobalMPISession.hpp"

#include "Teuchos_StackedTimer.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_DefaultComm.hpp"

#include "Kokkos_Core.hpp"

#include "Intrepid2_CellGeometry.hpp"
#include "Intrepid2_CellGeometryTestUtils.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_IntegrationTools.hpp"

enum AlgorithmChoice
{
  Standard,
  AffineNonTensor,
  NonAffineTensor,
  AffineTensor,
  DiagonalJacobian,
  Uniform
};

std::string to_string(AlgorithmChoice choice)
{
  switch (choice) {
    case Standard:         return "Standard";
    case AffineNonTensor:  return "AffineNonTensor";
    case NonAffineTensor:  return "NonAffineTensor";
    case AffineTensor:     return "AffineTensor";
    case DiagonalJacobian: return "DiagonalJacobian";
    case Uniform:          return "Uniform";
    
    default:               return "Unknown AlgorithmChoice";
  }
}

using namespace Intrepid2;

template< typename PointScalar, int spaceDim, typename ExecutionSpace >
inline
CellGeometry<PointScalar, spaceDim, ExecutionSpace> getMesh(AlgorithmChoice algorithmChoice, const int &meshWidth)
{
  Kokkos::Array<PointScalar,spaceDim> domainExtents;
  Kokkos::Array<int,spaceDim> gridCellCounts;
  for (int d=0; d<spaceDim; d++)
  {
    domainExtents[d]  = 0.0;
    gridCellCounts[d] = meshWidth;
  }
  auto uniformTensorGeometry = uniformCartesianMesh<PointScalar,spaceDim,ExecutionSpace>(domainExtents, gridCellCounts);
  
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

double flopsPerJacobian(const int &spaceDim, const int &numPoints, const int &numGeometryNodes)
{
  // implementation looks like:
//  for (ordinal_type i=0;i<dim;++i)
//  for (ordinal_type j=0;j<dim;++j) {
//    _jacobian(cell, point, i, j) = 0;
//    for (ordinal_type bf=0;bf<cardinality;++bf)
//      _jacobian(cell, point, i, j) += _worksetCells(cell+_startCell, bf, i) * _basisGrads(bf, point, j); // 2 flops: one multiply, one add
//  }
  return 2.0 * spaceDim * spaceDim * numPoints * numGeometryNodes;
}

double flopsPerJacobianDet(const int &spaceDim, const int &numPoints)
{
  //implementation in RealSpaceTools:
  /*value_type r_val = 0.0;
  switch (dim) {
  case 3:
    r_val = ( inMat(0,0) * inMat(1,1) * inMat(2,2) + // 3 flops: 2 mults, 1 add
              inMat(1,0) * inMat(2,1) * inMat(0,2) + // 3 flops: 2 mults, 1 add
              inMat(2,0) * inMat(0,1) * inMat(1,2) - // 3 flops: 2 mults, 1 subtract
              inMat(2,0) * inMat(1,1) * inMat(0,2) - // 3 flops: 2 mults, 1 subtract
              inMat(0,0) * inMat(2,1) * inMat(1,2) - // 3 flops: 2 mults, 1 subtract
              inMat(1,0) * inMat(0,1) * inMat(2,2) ); // 2 flops: 2 mults
    break;
  case 2:
    r_val = ( inMat(0,0) * inMat(1,1) -
              inMat(0,1) * inMat(1,0) );
    break;
  case 1:
    r_val = ( inMat(0,0) );
    break;
  }
  return r_val;*/
  int r_val;
  switch (spaceDim) {
    case 3: r_val = 17.0 * numPoints; break;
    case 2: r_val = 3.0 * numPoints; break;
    case 1: r_val = 0.0; break;
    default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unhandled spaceDim");
  }
  return r_val;
}

double flopsPerJacobianInverse(const int &spaceDim, const int &numPoints)
{
  // implementation looks like:
  // const value_type val = RealSpaceTools<>::Serial::det(mat);
  double totalFlops = flopsPerJacobianDet(spaceDim, numPoints);
  
  if (spaceDim == 3)
  {
  //  val0 =   mat(1,1)*mat(2,2) - mat(2,1)*mat(1,2); // 3 flops: 2 mults, 1 subtraction
  //  val1 = - mat(1,0)*mat(2,2) + mat(2,0)*mat(1,2); // 4 flops: 2 mults, 1 negation, 1 add
  //  val2 =   mat(1,0)*mat(2,1) - mat(2,0)*mat(1,1); // 3 flops: 2 mults, 1 subtraction
  //
  //  inv(0,0) = val0/val; // 1 flop
  //  inv(1,0) = val1/val; // 1 flop
  //  inv(2,0) = val2/val; // 1 flop
  //
  //  val0 =   mat(2,1)*mat(0,2) - mat(0,1)*mat(2,2); // 3
  //  val1 =   mat(0,0)*mat(2,2) - mat(2,0)*mat(0,2); // 3
  //  val2 = - mat(0,0)*mat(2,1) + mat(2,0)*mat(0,1); // 4
  //
  //  inv(0,1) = val0/val; // 1
  //  inv(1,1) = val1/val; // 1
  //  inv(2,1) = val2/val; // 1
  //
  //  val0 =   mat(0,1)*mat(1,2) - mat(1,1)*mat(0,2); // 3
  //  val1 = - mat(0,0)*mat(1,2) + mat(1,0)*mat(0,2); // 4
  //  val2 =   mat(0,0)*mat(1,1) - mat(1,0)*mat(0,1); // 3
  //
  //  inv(0,2) = val0/val; // 1
  //  inv(1,2) = val1/val; // 1
  //  inv(2,2) = val2/val; // 1
    totalFlops += 36.0 * numPoints;
  }
  else
  {
    INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unhandled spaceDim");
  }
  
  return totalFlops;
}

//! version that uses the classic, generic Intrepid2 paths
template<class Scalar, class PointScalar, int spaceDim, typename ExecSpaceType>
ScalarView<Scalar,ExecSpaceType> performStandardQuadratureHypercubeGRADGRAD(CellGeometry<PointScalar, spaceDim, ExecSpaceType> &geometry,
                                                                            const int &polyOrder, int worksetSize,
                                                                            double &transformIntegrateFlopCount, double &jacobianCellMeasureFlopCount)
{
  int numVertices = 1;
  for (int d=0; d<spaceDim; d++)
  {
    numVertices *= 2;
  }
  
  auto jacobianAndCellMeasureTimer = Teuchos::TimeMonitor::getNewTimer("Jacobians");
  auto fstIntegrateCall = Teuchos::TimeMonitor::getNewTimer("transform + integrate()");
  auto initialSetupTimer = Teuchos::TimeMonitor::getNewTimer("Initial Setup");
  initialSetupTimer->start();
  
  using CellTools = Intrepid2::CellTools<Kokkos::DefaultExecutionSpace>;
  using FunctionSpaceTools = Intrepid2::FunctionSpaceTools<Kokkos::DefaultExecutionSpace>;
  
  using namespace std;
  // dimensions of the returned view are (C,F,F)
  auto fs = Intrepid2::FUNCTION_SPACE_HGRAD;

  shards::CellTopology cellTopo = geometry.cellTopology();
  
  auto basis = Intrepid2::getBasis< Intrepid2::NodalBasisFamily<Kokkos::DefaultExecutionSpace> >(cellTopo, fs, polyOrder);
  
  int numFields = basis->getCardinality();
  int numCells = geometry.numCells();
  
  if (worksetSize > numCells) worksetSize = numCells;
  
  // local stiffness matrices:
  ScalarView<Scalar,ExecSpaceType> cellStiffness("cell stiffness matrices",numCells,numFields,numFields);
  
  using Kokkos::DefaultExecutionSpace;
  auto cubature = Intrepid2::DefaultCubatureFactory::create<ExecSpaceType>(cellTopo,polyOrder*2);
  int numPoints = cubature->getNumPoints();
  ScalarView<PointScalar,ExecSpaceType> cubaturePoints("cubature points",numPoints,spaceDim);
  ScalarView<double,ExecSpaceType> cubatureWeights("cubature weights", numPoints);
  
  cubature->getCubature(cubaturePoints, cubatureWeights);
  
  const double flopsPerJacobianPerCell    = flopsPerJacobian(spaceDim, numPoints, numVertices);
  const double flopsPerJacobianDetPerCell = flopsPerJacobianDet(spaceDim, numPoints);
  const double flopsPerJacobianInvPerCell = flopsPerJacobianInverse(spaceDim, numPoints);
  
  // Allocate some intermediate containers
  ScalarView<Scalar,ExecSpaceType> basisValues    ("basis values", numFields, numPoints );
  ScalarView<Scalar,ExecSpaceType> basisGradValues("basis grad values", numFields, numPoints, spaceDim);

  ScalarView<Scalar,ExecSpaceType> transformedGradValues("transformed grad values", worksetSize, numFields, numPoints, spaceDim);
  ScalarView<Scalar,ExecSpaceType> transformedWeightedGradValues("transformed weighted grad values", worksetSize, numFields, numPoints, spaceDim);
  
  basis->getValues(basisValues,     cubaturePoints, Intrepid2::OPERATOR_VALUE );
  basis->getValues(basisGradValues, cubaturePoints, Intrepid2::OPERATOR_GRAD  );
  
  const int numNodesPerCell = geometry.numNodesPerCell();
  ScalarView<PointScalar,ExecSpaceType> expandedCellNodes("expanded cell nodes",numCells,numNodesPerCell,spaceDim);
  Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numCells),
  KOKKOS_LAMBDA (const int &cellOrdinal) {
    for (int nodeOrdinal=0; nodeOrdinal<numNodesPerCell; nodeOrdinal++)
    {
      for (int d=0; d<spaceDim; d++)
      {
        expandedCellNodes(cellOrdinal,nodeOrdinal,d) = geometry(cellOrdinal,nodeOrdinal,d);
      }
    }
  });
  
  ScalarView<Scalar,ExecSpaceType> cellMeasures("cell measures", worksetSize, numPoints);
  ScalarView<Scalar,ExecSpaceType> jacobianDeterminant("jacobian determinant", worksetSize, numPoints);
  ScalarView<Scalar,ExecSpaceType> jacobian("jacobian", worksetSize, numPoints, spaceDim, spaceDim);
  ScalarView<Scalar,ExecSpaceType> jacobianInverse("jacobian inverse", worksetSize, numPoints, spaceDim, spaceDim);

  initialSetupTimer->stop();
  
  transformIntegrateFlopCount  = 0;
  jacobianCellMeasureFlopCount  = numCells * flopsPerJacobianPerCell;    // jacobian itself
  jacobianCellMeasureFlopCount += numCells * flopsPerJacobianInvPerCell; // inverse
  jacobianCellMeasureFlopCount += numCells * flopsPerJacobianDetPerCell; // determinant
  jacobianCellMeasureFlopCount += numCells * numPoints; // cell measure: (C,P) gets weighted with cubature weights of shape (P)
  
  int cellOffset = 0;
  while (cellOffset < numCells)
  {
    int startCell         = cellOffset;
    int numCellsInWorkset = (cellOffset + worksetSize - 1 < numCells) ? worksetSize : numCells - startCell;
    
    std::pair<int,int> cellRange = {startCell, startCell+numCellsInWorkset};
    auto cellWorkset = Kokkos::subview(expandedCellNodes, cellRange, Kokkos::ALL(), Kokkos::ALL());
    
    if (numCellsInWorkset != worksetSize)
    {
      Kokkos::resize(jacobian,                      numCellsInWorkset, numPoints, spaceDim, spaceDim);
      Kokkos::resize(jacobianInverse,               numCellsInWorkset, numPoints, spaceDim, spaceDim);
      Kokkos::resize(jacobianDeterminant,           numCellsInWorkset, numPoints);
      Kokkos::resize(cellMeasures,                  numCellsInWorkset, numPoints);
      Kokkos::resize(transformedGradValues,         numCellsInWorkset, numFields, numPoints, spaceDim);
      Kokkos::resize(transformedWeightedGradValues, numCellsInWorkset, numFields, numPoints, spaceDim);
    }
    jacobianAndCellMeasureTimer->start();
    CellTools::setJacobian(jacobian, cubaturePoints, cellWorkset, cellTopo); // accounted for outside loop, as numCells * flopsPerJacobianPerCell.
    CellTools::setJacobianInv(jacobianInverse, jacobian);
    CellTools::setJacobianDet(jacobianDeterminant, jacobian);
    
    FunctionSpaceTools::computeCellMeasure(cellMeasures, jacobianDeterminant, cubatureWeights);
    ExecSpaceType().fence();
    jacobianAndCellMeasureTimer->stop();
    
    // because structured integration performs transformations within integrate(), to get a fairer comparison here we include the transformation calls.
    fstIntegrateCall->start();
    FunctionSpaceTools::HGRADtransformGRAD(transformedGradValues, jacobianInverse, basisGradValues);
    transformIntegrateFlopCount += double(numCellsInWorkset) * double(numFields) * double(numPoints) * double(spaceDim) * (spaceDim - 1) * 2.0; // 2: one multiply, one add per (P,D) entry in the contraction.
    FunctionSpaceTools::multiplyMeasure(transformedWeightedGradValues, cellMeasures, transformedGradValues);
    transformIntegrateFlopCount += double(numCellsInWorkset) * double(numFields) * double(numPoints) * double(spaceDim); // multiply each entry of transformedGradValues: one flop for each.
        
    auto cellStiffnessSubview = Kokkos::subview(cellStiffness, cellRange, Kokkos::ALL(), Kokkos::ALL());
    
    FunctionSpaceTools::integrate(cellStiffnessSubview, transformedGradValues, transformedWeightedGradValues);
    ExecSpaceType().fence();
    fstIntegrateCall->stop();
    
    transformIntegrateFlopCount += double(numCellsInWorkset) * double(numFields) * double(numFields) * double(numPoints * spaceDim * 2 - 1); // 2: one multiply, one add per (P,D) entry in the contraction.
    
    cellOffset += worksetSize;
  }
//  std::cout << "standard integration, approximateFlopCount: " << approximateFlopCount << std::endl;
  return cellStiffness;
}

//! returns an estimated count of the floating point operations performed.
template<class Scalar, class PointScalar, int spaceDim, typename ExecSpaceType>
void performStructuredQuadratureHypercubeGRADGRAD(CellGeometry<PointScalar, spaceDim, ExecSpaceType> &geometry, const int &polyOrder, const int &worksetSize,
                                                  double &transformIntegrateFlopCount, double &jacobianCellMeasureFlopCount)
{
  int numVertices = 1;
  for (int d=0; d<spaceDim; d++)
  {
    numVertices *= 2;
  }
  
  auto initialSetupTimer = Teuchos::TimeMonitor::getNewTimer("Initial Setup");
  initialSetupTimer->start();
  using namespace std;
  using FunctionSpaceTools = Intrepid2::FunctionSpaceTools<ExecSpaceType>;
  using IntegrationTools   = Intrepid2::IntegrationTools<ExecSpaceType>;
  // dimensions of the returned view are (C,F,F)
  auto fs = Intrepid2::FUNCTION_SPACE_HGRAD;
  
  shards::CellTopology cellTopo = geometry.cellTopology();
  
  auto basis = Intrepid2::getBasis< Intrepid2::DerivedNodalBasisFamily<Kokkos::DefaultExecutionSpace> >(cellTopo, fs, polyOrder);
  
  int numFields = basis->getCardinality();
  int numCells = geometry.numCells();
    
  // local stiffness matrix:
  ScalarView<Scalar,ExecSpaceType> cellStiffness("cell stiffness matrices",numCells,numFields,numFields);
  
  auto cubature = Intrepid2::DefaultCubatureFactory::create<ExecSpaceType>(cellTopo,polyOrder*2);
  auto tensorCubatureWeights = cubature->allocateCubatureWeights();
  TensorPoints<PointScalar,ExecSpaceType> tensorCubaturePoints  = cubature->allocateCubaturePoints();
  
  cubature->getCubature(tensorCubaturePoints, tensorCubatureWeights);
  
  EOperator op = OPERATOR_GRAD;
  BasisValues<Scalar,ExecSpaceType> gradientValues = basis->allocateBasisValues(tensorCubaturePoints, op);
  basis->getValues(gradientValues, tensorCubaturePoints, op);
  
  // goal here is to do a weighted Poisson; i.e. (f grad u, grad v) on each cell
    
  int cellOffset = 0;
  
  auto jacobianAndCellMeasureTimer = Teuchos::TimeMonitor::getNewTimer("Jacobians");
  auto fstIntegrateCall = Teuchos::TimeMonitor::getNewTimer("transform + integrate()");
  
  Data<PointScalar,ExecSpaceType> jacobian = geometry.allocateJacobianData(tensorCubaturePoints, 0, worksetSize);
  Data<PointScalar,ExecSpaceType> jacobianDet = CellTools<ExecSpaceType>::allocateJacobianDet(jacobian);
  Data<PointScalar,ExecSpaceType> jacobianInv = CellTools<ExecSpaceType>::allocateJacobianInv(jacobian);
  TensorData<PointScalar,ExecSpaceType> cellMeasures = geometry.allocateCellMeasure(jacobianDet, tensorCubatureWeights);
  
  // lazily-evaluated transformed gradient values (temporary to allow integralData allocation)
  auto transformedGradientValuesTemp = FunctionSpaceTools::getHGRADtransformGRAD(jacobianInv, gradientValues.vectorData());
  auto integralData = IntegrationTools::allocateIntegralData(transformedGradientValuesTemp, cellMeasures, transformedGradientValuesTemp);
  
  const int numPoints = jacobian.getDataExtent(1); // data extent will be 1 for affine, numPoints for other cases
  
  // TODO: make the below determination accurate for diagonal/block-diagonal cases… (right now, will overcount)
  const double flopsPerJacobianPerCell    = flopsPerJacobian(spaceDim, numPoints, numVertices);
  const double flopsPerJacobianDetPerCell = flopsPerJacobianDet(spaceDim, numPoints);
  const double flopsPerJacobianInvPerCell = flopsPerJacobianInverse(spaceDim, numPoints);
  
  transformIntegrateFlopCount = 0;
  jacobianCellMeasureFlopCount  = numCells * flopsPerJacobianPerCell;    // jacobian itself
  jacobianCellMeasureFlopCount += numCells * flopsPerJacobianInvPerCell; // inverse
  jacobianCellMeasureFlopCount += numCells * flopsPerJacobianDetPerCell; // determinant
  jacobianCellMeasureFlopCount += numCells * numPoints; // cell measure: (C,P) gets weighted with cubature weights of shape (P)
  
  auto refData = geometry.getJacobianRefData(tensorCubaturePoints);
  
  initialSetupTimer->stop();
  while (cellOffset < numCells)
  {
    int startCell         = cellOffset;
    int numCellsInWorkset = (cellOffset + worksetSize - 1 < numCells) ? worksetSize : numCells - startCell;
    int endCell           = numCellsInWorkset + startCell;
    
    jacobianAndCellMeasureTimer->start();
    if (numCellsInWorkset != worksetSize)
    {
      const int CELL_DIM = 0; // first dimension corresponds to cell
      jacobian.setExtent(    CELL_DIM, numCellsInWorkset);
      jacobianDet.setExtent( CELL_DIM, numCellsInWorkset);
      jacobianInv.setExtent( CELL_DIM, numCellsInWorkset);
      integralData.setExtent(CELL_DIM, numCellsInWorkset);
      
      // cellMeasures is a TensorData object with separateFirstComponent_ = true; the below sets the cell dimension…
      cellMeasures.setFirstComponentExtentInDimension0(numCellsInWorkset);
    }
    
    geometry.setJacobian(jacobian, tensorCubaturePoints, refData, startCell, endCell);
    CellTools<ExecSpaceType>::setJacobianDet(jacobianDet, jacobian);
    CellTools<ExecSpaceType>::setJacobianInv(jacobianInv, jacobian);
    
    // lazily-evaluated transformed gradient values:
    auto transformedGradientValues = FunctionSpaceTools::getHGRADtransformGRAD(jacobianInv, gradientValues.vectorData());
    
    geometry.computeCellMeasure(cellMeasures, jacobianDet, tensorCubatureWeights);
    ExecSpaceType().fence();
    jacobianAndCellMeasureTimer->stop();
    
    bool sumInto = false;
    double approximateFlopCountIntegrateWorkset = 0;
    fstIntegrateCall->start();
    IntegrationTools::integrate(integralData, transformedGradientValues, cellMeasures, transformedGradientValues, sumInto, &approximateFlopCountIntegrateWorkset);
    ExecSpaceType().fence();
    fstIntegrateCall->stop();
    
    // copy into cellStiffness container.  (Alternately, do something like allocateIntegralData, but outside this loop, and take a subview to construct the workset integralData.)
    if (integralData.getUnderlyingViewRank() == 3)
    {
      std::pair<int,int> cellRange = {startCell, endCell};
      auto cellStiffnessSubview = Kokkos::subview(cellStiffness, cellRange, Kokkos::ALL(), Kokkos::ALL());
      Kokkos::deep_copy(cellStiffnessSubview, integralData.getUnderlyingView3());
    }
    else // underlying view rank is 2; copy to each cell in destination stiffness matrix
    {
      auto integralView2 = integralData.getUnderlyingView2();
      auto policy = Kokkos::MDRangePolicy<ExecSpaceType,Kokkos::Rank<3>>({0,0,0},{numCellsInWorkset,numFields,numFields});
      Kokkos::parallel_for("copy uniform data to expanded container", policy,
                       KOKKOS_LAMBDA (const int &cellOrdinal, const int &leftFieldOrdinal, const int &rightFieldOrdinal) {
        cellStiffness(startCell + cellOrdinal, leftFieldOrdinal, rightFieldOrdinal) = integralView2(leftFieldOrdinal,rightFieldOrdinal);
      });
    }
    
    transformIntegrateFlopCount  += approximateFlopCountIntegrateWorkset;
    
    cellOffset += worksetSize;
  }
}

int main( int argc, char* argv[] )
{
  // Note that the dtor for GlobalMPISession will call Kokkos::finalize_all() but does not call Kokkos::initialize()...
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize(argc,argv);
  
  using std::cout;
  using std::endl;
  using std::string;
  using std::vector;
  
  {
    // For now, we focus on 3D Poisson test.
    // The non-affine tensor case is one that, like the Standard case, does *not* make any algorithmic assumptions about the geometry.
    // This should be a good proxy for Poisson with unstructured material data on a curvilinear hexahedral mesh.
    const int spaceDim = 3;
    
    enum Mode
    {
      Calibration,
      Test,
      BestSerial,
      BestCuda
    };
    Mode mode;
    
    vector<AlgorithmChoice> allAlgorithmChoices {Standard, NonAffineTensor, AffineTensor, Uniform};
    
    Teuchos::CommandLineProcessor cmdp(false,true); // false: don't throw exceptions; true: do return errors for unrecognized options
    
    string algorithmChoiceString = "All"; // alternatives: Standard, NonAffineTensor, AffineTensor, Uniform
    
    int polyOrderFixed = -1;
    int polyOrderMin = 1;
    int polyOrderMax = 8;
    
    string modeChoiceString = "Test"; // alternatives: Calibration, BestSerial, BestCuda
    
    cmdp.setOption("algorithm", &algorithmChoiceString, "Options: All, Standard, NonAffineTensor, AffineTensor, Uniform");
    cmdp.setOption("polyOrder", &polyOrderFixed, "Single polynomial degree to run at");
    cmdp.setOption("minPolyOrder", &polyOrderMin, "Starting polynomial degree to run at");
    cmdp.setOption("maxPolyOrder", &polyOrderMax, "Maximum polynomial degree to run at");
    cmdp.setOption("mode", &modeChoiceString);
    
    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
    {
  #ifdef HAVE_MPI
      MPI_Finalize();
  #endif
      return -1;
    }

    vector<AlgorithmChoice> algorithmChoices;
    if (algorithmChoiceString == "All")
    {
      algorithmChoices = allAlgorithmChoices;
    }
    else if (algorithmChoiceString == "Standard")
    {
      algorithmChoices = vector<AlgorithmChoice>{Standard};
    }
    else if (algorithmChoiceString == "NonAffineTensor")
    {
      algorithmChoices = vector<AlgorithmChoice>{NonAffineTensor};
    }
    else if (algorithmChoiceString == "AffineTensor")
    {
      algorithmChoices = vector<AlgorithmChoice>{AffineTensor};
    }
    else if (algorithmChoiceString == "Uniform")
    {
      algorithmChoices = vector<AlgorithmChoice>{Uniform};
    }
    else
    {
      cout << "Unrecognized algorithm choice: " << algorithmChoiceString << endl;
#ifdef HAVE_MPI
      MPI_Finalize();
#endif
      return -1;
    }
    
    if (polyOrderFixed > 0)
    {
      polyOrderMin = polyOrderFixed;
      polyOrderMax = polyOrderFixed;
    }
    
    if (modeChoiceString == "Calibration")
    {
      mode = Calibration;
    }
    else if (modeChoiceString == "BestCuda")
    {
      mode = BestCuda;
    }
    else if (modeChoiceString == "BestSerial")
    {
      mode = BestSerial;
    }
    else if (modeChoiceString == "Test")
    {
      mode = Test;
    }
    else
    {
      cout << "Unrecognized mode choice: " << modeChoiceString << endl;
#ifdef HAVE_MPI
      MPI_Finalize();
#endif
      return -1;
    }
    
    using Scalar = double;
    using ExecutionSpace = Kokkos::DefaultExecutionSpace;
    
    using std::vector;
    using std::map;
    using std::pair;
    using std::make_pair;
    using std::tuple;
    using std::cout;
    using std::endl;
    using std::setw;
    
    using WorksetForAlgorithmChoice = map<AlgorithmChoice, int>;
    
    vector< tuple<int,int,WorksetForAlgorithmChoice> > polyOrderMeshWidthWorksetTestCases;
    
    const int meshWidth = 16;
    vector<int> worksetSizes {1,2,4,8,16,32,64,128,256,512,1024,2048,4096};
    
    // due to memory constraints, restrict the workset size for higher orders
    map<int,int> maxWorksetSizeForPolyOrder;
    maxWorksetSizeForPolyOrder[1] = 4096;
    maxWorksetSizeForPolyOrder[2] = 4096;
    maxWorksetSizeForPolyOrder[3] = 4096;
    maxWorksetSizeForPolyOrder[4] = 4096;
    maxWorksetSizeForPolyOrder[5] = 4096;
    maxWorksetSizeForPolyOrder[6] = 2048;
    maxWorksetSizeForPolyOrder[7] = 1024;
    maxWorksetSizeForPolyOrder[8] = 512;
    
    map<int,int> minWorksetSizeForPolyOrder;
    minWorksetSizeForPolyOrder[1] = 1;
    minWorksetSizeForPolyOrder[2] = 1;
    minWorksetSizeForPolyOrder[3] = 1;
    minWorksetSizeForPolyOrder[4] = 1;
    minWorksetSizeForPolyOrder[5] = 1;
    minWorksetSizeForPolyOrder[6] = 1;
    minWorksetSizeForPolyOrder[7] = 1;
    minWorksetSizeForPolyOrder[8] = 1;
    
    switch (mode) {
      case Calibration:
      {
        for (int polyOrder=polyOrderMin; polyOrder<=polyOrderMax; polyOrder++)
        {
          for (int worksetSize : worksetSizes)
          {
            if (worksetSize > maxWorksetSizeForPolyOrder[polyOrder])
            {
              continue;
            }
            if (worksetSize < minWorksetSizeForPolyOrder[polyOrder])
            {
              continue;
            }
            // use the same worksetSize for all AlgorithmChoice's
            WorksetForAlgorithmChoice worksetForAlgorithmChoice;
            for (auto algorithmChoice : algorithmChoices)
            {
              worksetForAlgorithmChoice[algorithmChoice] = worksetSize;
            }
            polyOrderMeshWidthWorksetTestCases.push_back(tuple<int,int,WorksetForAlgorithmChoice>{polyOrder,meshWidth,worksetForAlgorithmChoice} );
          }
        }
      }
      break;
      case Test:
      {
        // for test run, use the same modestly-sized tuples for each AlgorithmChoice
        // (note that meshWidth varies here)
        vector< tuple<int,int,int> > testCases { tuple<int,int,int> {1,8,512},
                                                 tuple<int,int,int> {2,8,256},
                                                 tuple<int,int,int> {3,4,64},
                                                 tuple<int,int,int> {3,4,9} // test case whose workset size does not evenly divide the cell count
        };
        
        for (auto testCase : testCases )
        {
          int polyOrder   = std::get<0>(testCase);
          int meshWidth   = std::get<1>(testCase);
          
          int numCells = 1;
          for (int d=0; d<spaceDim; d++)
          {
            numCells *= meshWidth;
          }
          
          WorksetForAlgorithmChoice worksetForAlgorithmChoice;
          for (auto algorithmChoice : algorithmChoices)
          {
            worksetForAlgorithmChoice[algorithmChoice] = std::get<2>(testCase);
          }
          worksetForAlgorithmChoice[Uniform] = numCells;
          polyOrderMeshWidthWorksetTestCases.push_back(tuple<int,int,WorksetForAlgorithmChoice>{polyOrder,meshWidth,worksetForAlgorithmChoice} );
        }
      }
      break;
      case BestSerial:
      {
        // manually calibrated workset sizes on iMac Pro (2.3 GHz Xeon W, 18-core, running in serial)
        // (these were calibrated without much tuning for the affine tensor case; if/when that happens, will want to recalibrate.)
        
        map<int,int> standardWorksetForPolyOrder;
        standardWorksetForPolyOrder[1] = 64;
        standardWorksetForPolyOrder[2] = 64;
        standardWorksetForPolyOrder[3] = 128;
        standardWorksetForPolyOrder[4] = 64;
        standardWorksetForPolyOrder[5] = 16;
        standardWorksetForPolyOrder[6] = 2;
        standardWorksetForPolyOrder[7] = 2;
        standardWorksetForPolyOrder[8] = 1;
        
        // Non-Affine Tensor
        // it turns out 4096 is the best choice for the PointValueCache algorithm for any polyOrder from 1 to 5
        // this likely means we're not exposing enough parallelism within the cell.
        map<int,int> nonAffineTensorWorksetForPolyOrder;
        nonAffineTensorWorksetForPolyOrder[1] = 256;
        nonAffineTensorWorksetForPolyOrder[2] = 256;
        nonAffineTensorWorksetForPolyOrder[3] = 128;
        nonAffineTensorWorksetForPolyOrder[4] = 64;
        nonAffineTensorWorksetForPolyOrder[5] = 16;
        nonAffineTensorWorksetForPolyOrder[6] = 8;
        nonAffineTensorWorksetForPolyOrder[7] = 2;
        nonAffineTensorWorksetForPolyOrder[8] = 2;
        
        map<int,int> affineTensorWorksetForPolyOrder;
        affineTensorWorksetForPolyOrder[1] = 256;
        affineTensorWorksetForPolyOrder[2] = 128;
        affineTensorWorksetForPolyOrder[3] = 16;
        affineTensorWorksetForPolyOrder[4] = 4;
        affineTensorWorksetForPolyOrder[5] = 2;
        affineTensorWorksetForPolyOrder[6] = 1;
        affineTensorWorksetForPolyOrder[7] = 1;
        affineTensorWorksetForPolyOrder[8] = 1;
        
        // for the cases that we have not tried yet, we try to choose sensible guesses for workset size:
        // 1 is best, we think, for polyOrder 8, so it'll be the best for the rest.
        int worksetSize = 1;
        for (int polyOrder=9; polyOrder <= polyOrderMax; polyOrder++)
        {
          nonAffineTensorWorksetForPolyOrder[polyOrder] = worksetSize;
          affineTensorWorksetForPolyOrder[polyOrder]    = worksetSize;
          standardWorksetForPolyOrder[polyOrder]        = worksetSize;
        }
        
        int numCells = 1;
        for (int d=0; d<spaceDim; d++)
        {
          numCells *= meshWidth;
        }
        
        for (int polyOrder=polyOrderMin; polyOrder<=polyOrderMax; polyOrder++)
        {
          WorksetForAlgorithmChoice worksetForAlgorithmChoice;
          worksetForAlgorithmChoice[Standard]        = standardWorksetForPolyOrder       [polyOrder];
          worksetForAlgorithmChoice[NonAffineTensor] = nonAffineTensorWorksetForPolyOrder[polyOrder];
          worksetForAlgorithmChoice[AffineTensor]    = nonAffineTensorWorksetForPolyOrder[polyOrder];
          worksetForAlgorithmChoice[Uniform]         = numCells;
          
          polyOrderMeshWidthWorksetTestCases.push_back(tuple<int,int,WorksetForAlgorithmChoice>{polyOrder,meshWidth,worksetForAlgorithmChoice} );
        }
      }
        break;
      case BestCuda:
      {
        {
          // STANDARD
          // manually calibrated workset size on P100 (ride) - best for Standard
          // these are for 4096-element meshes
          map<int,int> standardWorksetForPolyOrder;
          standardWorksetForPolyOrder[1] = 4096;
          standardWorksetForPolyOrder[2] = 1024;
          standardWorksetForPolyOrder[3] = 128;
          standardWorksetForPolyOrder[4] = 64;
          standardWorksetForPolyOrder[5] = 8;
          standardWorksetForPolyOrder[6] = 4;
          standardWorksetForPolyOrder[7] = 2;
          standardWorksetForPolyOrder[8] = 1;
          
          // Non-Affine Tensor
          // it turns out 4096 is the best choice for the PointValueCache algorithm for any polyOrder from 1 to 5
          // this likely means we're not exposing enough parallelism within the cell.
          map<int,int> nonAffineTensorWorksetForPolyOrder;
          nonAffineTensorWorksetForPolyOrder[1] = 4096;
          nonAffineTensorWorksetForPolyOrder[2] = 4096;
          nonAffineTensorWorksetForPolyOrder[3] = 4096;
          nonAffineTensorWorksetForPolyOrder[4] = 4096;
          nonAffineTensorWorksetForPolyOrder[5] = 4096;
          nonAffineTensorWorksetForPolyOrder[6] = 2048;
          nonAffineTensorWorksetForPolyOrder[7] = 512;
          nonAffineTensorWorksetForPolyOrder[8] = 512;
          
          // for the cases that we have not tried yet, we try to choose sensible guesses for workset size:
          int nonAffineWorksetSize = 256; // divide by 2 for each polyOrder beyond 8
          int standardWorksetSize  = 1;  // 1 is best, we think, for polyOrder 8, so it'll be the best for the rest.
          for (int polyOrder=9; polyOrder <= polyOrderMax; polyOrder++)
          {
            nonAffineTensorWorksetForPolyOrder[polyOrder] = nonAffineWorksetSize;
            nonAffineWorksetSize = (nonAffineWorksetSize > 1) ? nonAffineWorksetSize / 2 : 1;
            standardWorksetForPolyOrder[polyOrder] = standardWorksetSize;
          }
          
          int numCells = 1;
          for (int d=0; d<spaceDim; d++)
          {
            numCells *= meshWidth;
          }
          
          for (int polyOrder=polyOrderMin; polyOrder<=polyOrderMax; polyOrder++)
          {
            WorksetForAlgorithmChoice worksetForAlgorithmChoice;
            worksetForAlgorithmChoice[Standard]        = standardWorksetForPolyOrder       [polyOrder];
            worksetForAlgorithmChoice[NonAffineTensor] = nonAffineTensorWorksetForPolyOrder[polyOrder];
            worksetForAlgorithmChoice[AffineTensor]    = nonAffineTensorWorksetForPolyOrder[polyOrder];
            worksetForAlgorithmChoice[Uniform]         = numCells;
            
            polyOrderMeshWidthWorksetTestCases.push_back(tuple<int,int,WorksetForAlgorithmChoice>{polyOrder,meshWidth,worksetForAlgorithmChoice} );
          }
        }
        break;
        
      default:
        break;
    }
    }
    
    cout << std::setprecision(2) << std::scientific;
    
    map< AlgorithmChoice, map<int, pair<double,int> > > maxAlgorithmThroughputForPolyOrder; // values are (throughput in GFlops/sec, worksetSize)
    
    const int charWidth = 15;
    
    for (auto & testCase : polyOrderMeshWidthWorksetTestCases)
    {
      int polyOrder       = std::get<0>(testCase);
      int meshWidth       = std::get<1>(testCase);
      auto worksetSizeMap = std::get<2>(testCase);
      std::cout << "\n\n";
      std::cout << "Running with polyOrder = " << polyOrder << ", meshWidth = " << meshWidth << std::endl;
      for (auto algorithmChoice : algorithmChoices)
      {
        int worksetSize = worksetSizeMap[algorithmChoice];
        auto geometry = getMesh<Scalar, spaceDim, ExecutionSpace>(algorithmChoice, meshWidth);
        
        // timers recorded in performStructuredQuadratureHypercubeGRADGRAD, performStandardQuadratureHypercubeGRADGRAD
        auto jacobianAndCellMeasureTimer = Teuchos::TimeMonitor::getNewTimer("Jacobians");
        auto fstIntegrateCall = Teuchos::TimeMonitor::getNewTimer("transform + integrate()");
        auto initialSetupTimer = Teuchos::TimeMonitor::getNewTimer("Initial Setup");

        jacobianAndCellMeasureTimer->reset();
        fstIntegrateCall->reset();
        initialSetupTimer->reset();
        
        double elapsedTimeSeconds = 0;
        double jacobianCellMeasureFlopCount = 0;
        double transformIntegrateFlopCount = 0;
        
        if (algorithmChoice == Standard)
        {
          // each cell needs on the order of polyOrder^N quadrature points, each of which has a Jacobian of size N * N.
          auto timer = Teuchos::TimeMonitor::getNewTimer("Standard Integration");
          timer->start();
          performStandardQuadratureHypercubeGRADGRAD<Scalar,Scalar,spaceDim,ExecutionSpace>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
          timer->stop();
          elapsedTimeSeconds = timer->totalElapsedTime();
          
          cout << "Standard, workset size:          " << setw(charWidth) << worksetSize << endl;
          
          timer->reset();
        }
        else if (algorithmChoice == AffineTensor)
        {
          auto timer = Teuchos::TimeMonitor::getNewTimer("Affine tensor Integration");
          timer->start();
          performStructuredQuadratureHypercubeGRADGRAD<Scalar>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
          timer->stop();
          
          elapsedTimeSeconds = timer->totalElapsedTime();
          
          cout << "Affine Tensor, workset size:     " << setw(charWidth) << worksetSize << endl;
                    
          timer->reset();
        }
        else if (algorithmChoice == NonAffineTensor)
        {
          auto timer = Teuchos::TimeMonitor::getNewTimer("Non-affine tensor Integration");
          timer->start();
          performStructuredQuadratureHypercubeGRADGRAD<Scalar>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
          timer->stop();
          
          elapsedTimeSeconds = timer->totalElapsedTime();
          
          cout << "Non-Affine Tensor, workset size: " << setw(charWidth) << worksetSize << endl;
          
          timer->reset();
        }
        else if (algorithmChoice == Uniform)
        {
          // for uniform, override worksetSize: no loss in taking maximal worksetSize
          int numCells = 1;
          for (int d=0; d<spaceDim; d++)
          {
            numCells *= meshWidth;
          }
          auto timer = Teuchos::TimeMonitor::getNewTimer("Uniform Integration");
          timer->start();
          performStructuredQuadratureHypercubeGRADGRAD<Scalar>(geometry, polyOrder, numCells, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
          timer->stop();
          
          elapsedTimeSeconds = timer->totalElapsedTime();
          
          cout << "Uniform, workset size:           " << setw(charWidth) << worksetSize << endl;
          
          timer->reset();
        }
        
        const double approximateFlopCountTotal = transformIntegrateFlopCount + jacobianCellMeasureFlopCount;
        const double overallThroughputInGFlops = approximateFlopCountTotal / elapsedTimeSeconds / 1.0e9;
        
        const double previousMaxThroughput = maxAlgorithmThroughputForPolyOrder[algorithmChoice][polyOrder].first;
        if (overallThroughputInGFlops > previousMaxThroughput)
        {
          maxAlgorithmThroughputForPolyOrder[algorithmChoice][polyOrder] = make_pair(overallThroughputInGFlops,worksetSize);
        }
        
        // timing details
        double integrateCallTime       = fstIntegrateCall->totalElapsedTime();
        double integrateCallPercentage = integrateCallTime / elapsedTimeSeconds * 100.0;
        double jacobianTime            = jacobianAndCellMeasureTimer->totalElapsedTime();
        double jacobianPercentage      = jacobianTime / elapsedTimeSeconds * 100.0;
        double initialSetupTime        = initialSetupTimer->totalElapsedTime();
        double initialSetupPercentage  = initialSetupTime / elapsedTimeSeconds * 100.0;
        double remainingTime           = elapsedTimeSeconds - (integrateCallTime + jacobianTime + initialSetupTime);
        double remainingPercentage     = remainingTime / elapsedTimeSeconds * 100.0;
        
        const double transformIntegrateThroughputInGFlops = transformIntegrateFlopCount  / integrateCallTime / 1.0e9;
        const double jacobiansThroughputInGFlops          = jacobianCellMeasureFlopCount / jacobianTime      / 1.0e9;
        cout << "Time (core integration)          " << setw(charWidth) << std::scientific << integrateCallTime << " seconds (" << std::fixed << integrateCallPercentage << "%)." << endl;
        cout << "flop estimate (core):            " << setw(charWidth) << std::scientific << transformIntegrateFlopCount << endl;
        cout << "estimated throughput (core):     " << setw(charWidth) << std::scientific << transformIntegrateThroughputInGFlops << " GFlops" << endl;
        cout << std::fixed;
        cout << "Time (Jacobians)                 " << setw(charWidth) << std::scientific << jacobianTime      << " seconds (" << std::fixed << jacobianPercentage      << "%)." << endl;
        cout << "flop estimate (Jacobians):       " << setw(charWidth) << std::scientific << jacobianCellMeasureFlopCount << endl;
        cout << "estimated throughput (Jac.):     " << setw(charWidth) << std::scientific << jacobiansThroughputInGFlops << " GFlops" << endl;
        cout << "Time (initial setup)             " << setw(charWidth) << std::scientific << initialSetupTime  << " seconds (" << std::fixed << initialSetupPercentage  << "%)." << endl;
        cout << "Time (other)                     " << setw(charWidth) << std::scientific << remainingTime     << " seconds (" << std::fixed << remainingPercentage     << "%)." << endl;
        cout << "Time (total):                    " << setw(charWidth) << std::scientific << elapsedTimeSeconds   << " seconds.\n";
        cout << "flop estimate (total):           " << setw(charWidth) << std::scientific << approximateFlopCountTotal << endl;
        cout << "estimated throughput (total):    " << setw(charWidth) << std::scientific << overallThroughputInGFlops << " GFlops" << endl;
        
        cout << endl;
      }
    }
    
    if (mode == Calibration)
    {
      for (auto & algorithmChoice : algorithmChoices)
      {
        cout << "Best workset sizes for " << to_string(algorithmChoice) << ":" << endl;
        
        for (auto & maxThroughputEntry : maxAlgorithmThroughputForPolyOrder[algorithmChoice])
        {
          int polyOrder   = maxThroughputEntry.first;
          int worksetSize = maxThroughputEntry.second.second;
          double throughput = maxThroughputEntry.second.first;
          cout << "p = " << polyOrder << ":" << setw(5) << worksetSize << " (" << throughput << " GFlops/sec)\n";
        }
      }
    }
  }
  
  return 0;
}
