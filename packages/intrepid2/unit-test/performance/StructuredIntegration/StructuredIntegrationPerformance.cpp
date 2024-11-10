// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   StructuredIntegrationPerformance.cpp
    \brief  Driver for performance tests comparing structured integration performance to standard.
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
#include "Intrepid2_TestUtils.hpp"

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
#include "VectorWeightedGRADGRADStandardAssembly.hpp"
#include "VectorWeightedGRADGRADStructuredAssembly.hpp"

enum FormulationChoice
{
  Poisson, // (grad, grad)
  Hgrad,   // (grad, grad) + (value, value)
  Hdiv,    // (div, div)   + (value, value)
  Hcurl,   // (curl, curl) + (value, value)
  L2,      // (value, value)
  VectorWeightedPoisson,
  UnknownFormulation
};

enum AlgorithmChoice
{
  Standard,
  AffineNonTensor,
  NonAffineTensor,
  AffineTensor,
  DiagonalJacobian,
  Uniform
};

enum BasisFamilyChoice
{
  Nodal,
  Hierarchical,
  Serendipity
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

std::string to_string(FormulationChoice choice)
{
  switch (choice) {
    case Poisson:               return "Poisson";
    case Hgrad:                 return "Hgrad";
    case Hdiv:                  return "Hdiv";
    case Hcurl:                 return "Hcurl";
    case L2:                    return "L2";
    case VectorWeightedPoisson: return "VectorWeightedPoisson";
    
    default:      return "Unknown FormulationChoice";
  }
}

AlgorithmChoice algorithm_from_string(const std::string &algorithmString)
{
  if (algorithmString == "Standard")
  {
    return Standard;
  }
  else if (algorithmString == "AffineNonTensor")
  {
    return AffineNonTensor;
  }
  else if (algorithmString == "NonAffineTensor")
  {
    return NonAffineTensor;
  }
  else if (algorithmString == "AffineTensor")
  {
    return AffineTensor;
  }
  else if (algorithmString == "DiagonalJacobian")
  {
    return DiagonalJacobian;
  }
  else if (algorithmString == "Uniform")
  {
    return Uniform;
  }
  else
  {
    INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unrecognized algorithm");
  }
}

FormulationChoice formulation_from_string(const std::string &formulationString)
{
  if (formulationString == "Poisson")
  {
    return Poisson;
  }
  else if (formulationString == "Hgrad")
  {
    return Hgrad;
  }
  else if (formulationString == "Hcurl")
  {
    return Hcurl;
  }
  else if (formulationString == "Hdiv")
  {
    return Hdiv;
  }
  else if (formulationString == "L2")
  {
    return L2;
  }
  else
  {
    INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unrecognized formulation");
  }
}

using namespace Intrepid2;

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

template< typename PointScalar, int spaceDim, typename DeviceType >
inline
CellGeometry<PointScalar, spaceDim, DeviceType> getMesh(AlgorithmChoice algorithmChoice, const int &meshWidth)
{
  Kokkos::Array<int,spaceDim> gridCellCounts;
  for (int d=0; d<spaceDim; d++)
  {
    gridCellCounts[d] = meshWidth;
  }
  return getMesh<PointScalar, spaceDim, DeviceType>(algorithmChoice, gridCellCounts);
}

//! Returns an Array of roughly isotropic grid dimensions for which (C,F,F) stifness matrix will have at most maxStiffnessEntryCount entries, and at most maxElements total cells.
template<int spaceDim>
Kokkos::Array<int,spaceDim>
getMeshWidths(int basisCardinality, int maxStiffnessEntryCount, int maxElements)
{
  Kokkos::Array<int,spaceDim> meshWidths;
  const int entriesPerElement = basisCardinality * basisCardinality;
  const int maxElementCount   = std::min(maxStiffnessEntryCount / entriesPerElement, maxElements);
  
  // initialize meshWidths:
  for (int d=0; d<spaceDim; d++)
  {
    meshWidths[d] = 1;
  }
  
  // double in each dimension until doing so would make the number of elements would exceed maxElementCount:
  int numElements = 1;
  int d = 0;
  while (numElements * 2 <= maxElementCount)
  {
    meshWidths[d] *= 2;
    d = (d + 1) % spaceDim;
    numElements *= 2;
  }
  return meshWidths;
}

template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType, unsigned long spaceDim2=spaceDim>
Intrepid2::ScalarView<Scalar,DeviceType> performStandardQuadrature(FormulationChoice formulation,
                                                                   Intrepid2::CellGeometry<PointScalar, spaceDim, DeviceType> &geometry, const int &polyOrder, const int &worksetSize,
                                                                   double &transformIntegrateFlopCount, double &jacobianCellMeasureFlopCount,
                                                                   Teuchos::RCP<Kokkos::Array<Scalar,spaceDim2>> vectorWeight1 = Teuchos::null,
                                                                   Teuchos::RCP<Kokkos::Array<Scalar,spaceDim2>> vectorWeight2 = Teuchos::null)
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
    case VectorWeightedPoisson:
      return performStandardQuadratureVectorWeightedGRADGRAD<Scalar, BasisFamily>(geometry, polyOrder, worksetSize, vectorWeight1, vectorWeight2, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
    default:
      return Intrepid2::ScalarView<Scalar,DeviceType>();
  }
}

template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType, unsigned long spaceDim2=spaceDim>
Intrepid2::ScalarView<Scalar,DeviceType> performStructuredQuadrature(FormulationChoice formulation,
                                                                     Intrepid2::CellGeometry<PointScalar, spaceDim, DeviceType> &geometry, const int &polyOrder, const int &worksetSize,
                                                                     double &transformIntegrateFlopCount, double &jacobianCellMeasureFlopCount,
                                                                     Teuchos::RCP<Kokkos::Array<Scalar,spaceDim2>> vectorWeight1 = Teuchos::null,
                                                                     Teuchos::RCP<Kokkos::Array<Scalar,spaceDim2>> vectorWeight2 = Teuchos::null)
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
    case VectorWeightedPoisson:
      return performStructuredQuadratureVectorWeightedGRADGRAD<Scalar, BasisFamily>(geometry, polyOrder, worksetSize, vectorWeight1, vectorWeight2, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
    default:
      return Intrepid2::ScalarView<Scalar,DeviceType>();
  }
}

template<class Scalar, class BasisFamily, class PointScalar, int spaceDim>
typename BasisFamily::BasisPtr getBasisForFormulation(FormulationChoice formulation, shards::CellTopology &cellTopo, const int polyOrder)
{
  Intrepid2::EFunctionSpace fs;
  switch (formulation)
  {
    case Poisson:               fs = FUNCTION_SPACE_HGRAD; break;
    case Hgrad:                 fs = FUNCTION_SPACE_HGRAD; break;
    case Hdiv:                  fs = FUNCTION_SPACE_HDIV;  break;
    case Hcurl:                 fs = FUNCTION_SPACE_HCURL; break;
    case L2:                    fs = FUNCTION_SPACE_HVOL;  break;
    case VectorWeightedPoisson: fs = FUNCTION_SPACE_HGRAD; break;
    case UnknownFormulation:    INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unknown formulation");
  }
  
  auto basis = getBasis< BasisFamily >(cellTopo, fs, polyOrder);
  return basis;
}

template<class Scalar, class PointScalar, class DeviceType, int spaceDim>
BasisPtr<DeviceType,Scalar,Scalar> getHypercubeBasisForFormulation(FormulationChoice formulation, BasisFamilyChoice basisFamilyChoice, const int polyOrder)
{
  shards::CellTopology cellTopo;
  switch (spaceDim)
  {
    case 1: cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Line<2> >());          break;
    case 2: cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >()); break;
    case 3: cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >());    break;
    default:
      INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported spaceDim");
  }
  
  switch (basisFamilyChoice)
  {
    case Nodal:
    {
      using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
      return getBasisForFormulation<Scalar, BasisFamily, Scalar, spaceDim>(formulation, cellTopo, polyOrder);
    }
      break;
    case Hierarchical:
    {
      using BasisFamily = HierarchicalBasisFamily<DeviceType>;
      return getBasisForFormulation<Scalar, BasisFamily, Scalar, spaceDim>(formulation, cellTopo, polyOrder);
    }
      break;
    case Serendipity:
      INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "basis family choice not yet implemented");
  }
  
  return Teuchos::null;
}

using std::map;
using std::tuple;
using std::vector;

enum Mode
{
  Calibration,
  Test,
  BestSerial,
  BestOpenMP_16,
  BestCuda,
  Precalibrated
};

map<tuple<Mode,FormulationChoice,AlgorithmChoice>,map<int,int> > getWorksetSizeMap(const std::string &precalibrationFile,
                                                                                   const map<int,int> &cellCountForPolyOrder,
                                                                                   const int &polyOrderMin,
                                                                                   const int &polyOrderMax)
{
  const int spaceDim = 3;
  
  map<tuple<Mode,FormulationChoice,AlgorithmChoice>,map<int,int> > worksetSizeMap; // keys are maps p -> worksetSize
  
  vector<AlgorithmChoice> allAlgorithmChoices {Standard, NonAffineTensor, AffineTensor, Uniform};
  vector<FormulationChoice> allFormulationChoices {Poisson, Hgrad, Hdiv, Hcurl, L2, VectorWeightedPoisson};
  
  // skip calibration case; want that to span workset sizes in a particular way…
  vector<Mode> allModes {Test,BestSerial,BestOpenMP_16,BestCuda,Precalibrated};
  
  for (auto mode : allModes)
  {
    // for the cases that we have not tried yet (polyOrder > 8), we try to choose sensible guesses for workset size:
    // 1 is best for polyOrder 8, so it'll be the best for the rest.
    for (int polyOrder=9; polyOrder <= polyOrderMax; polyOrder++)
    {
      for (auto formulation : allFormulationChoices)
      {
        for (auto algorithm : allAlgorithmChoices)
        {
          tuple<Mode,FormulationChoice,AlgorithmChoice> key {mode,formulation,algorithm};
          worksetSizeMap[key][polyOrder] = 1;
        }
      }
    }
    
    // best choice for Uniform: whole mesh (cellCount)
    for (auto formulation : allFormulationChoices)
    {
      tuple<Mode,FormulationChoice,AlgorithmChoice> key {mode,formulation,Uniform};
      worksetSizeMap[key] = cellCountForPolyOrder;
    }
    
    switch (mode) {
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
          
          for (auto formulation : allFormulationChoices)
          {
            for (auto algorithm : allAlgorithmChoices)
            {
              tuple<Mode,FormulationChoice,AlgorithmChoice> key {mode,formulation,algorithm};
              
              worksetSizeMap[key][polyOrder] = (algorithm != Uniform) ? std::get<2>(testCase) : numCells;
            }
          }
        }
      } // Test mode case
      break;
      case BestSerial:
      {
        // manually calibrated workset sizes on Mac Pro (2.5 GHz Xeon W, 28-core, running in serial)
        {
          // Poisson formulation
          FormulationChoice formulation = Poisson;
          tuple<Mode,FormulationChoice,AlgorithmChoice> standardKey {mode,formulation,Standard};
          tuple<Mode,FormulationChoice,AlgorithmChoice> nonAffineTensorKey {mode,formulation,NonAffineTensor};
          tuple<Mode,FormulationChoice,AlgorithmChoice> affineTensorKey {mode,formulation,AffineTensor};
          
          // best for Poisson - these are for meshes that range from 32768 for p=1 to 256 for p=8
          worksetSizeMap[standardKey][1] = 8192;
          worksetSizeMap[standardKey][2] = 4096;
          worksetSizeMap[standardKey][3] =   64;
          worksetSizeMap[standardKey][4] =   16;
          worksetSizeMap[standardKey][5] =   16;
          worksetSizeMap[standardKey][6] =    1;
          worksetSizeMap[standardKey][7] =    1;
          worksetSizeMap[standardKey][8] =    1;
          
          worksetSizeMap[nonAffineTensorKey][1] = 2048;
          worksetSizeMap[nonAffineTensorKey][2] =  256;
          worksetSizeMap[nonAffineTensorKey][3] =  128;
          worksetSizeMap[nonAffineTensorKey][4] =   16;
          worksetSizeMap[nonAffineTensorKey][5] =    2;
          worksetSizeMap[nonAffineTensorKey][6] =    1;
          worksetSizeMap[nonAffineTensorKey][7] =    1;
          worksetSizeMap[nonAffineTensorKey][8] =    1;
          
          worksetSizeMap[affineTensorKey][1] = 4096;
          worksetSizeMap[affineTensorKey][2] =   64;
          worksetSizeMap[affineTensorKey][3] =   32;
          worksetSizeMap[affineTensorKey][4] =    4;
          worksetSizeMap[affineTensorKey][5] =    2;
          worksetSizeMap[affineTensorKey][6] =    1;
          worksetSizeMap[affineTensorKey][7] =    1;
          worksetSizeMap[affineTensorKey][8] =    1;
        }
        {
          // Hgrad formulation
          FormulationChoice formulation = Hgrad;
          tuple<Mode,FormulationChoice,AlgorithmChoice> standardKey {mode,formulation,Standard};
          tuple<Mode,FormulationChoice,AlgorithmChoice> nonAffineTensorKey {mode,formulation,NonAffineTensor};
          tuple<Mode,FormulationChoice,AlgorithmChoice> affineTensorKey {mode,formulation,AffineTensor};
          
          // best for Hgrad - these are for meshes that range from 32768 for p=1 to 256 for p=8
          worksetSizeMap[standardKey][1] = 32768;
          worksetSizeMap[standardKey][2] = 16384;
          worksetSizeMap[standardKey][3] =   512;
          worksetSizeMap[standardKey][4] =   512;
          worksetSizeMap[standardKey][5] =   512;
          worksetSizeMap[standardKey][6] =     2;
          worksetSizeMap[standardKey][7] =     1;
          worksetSizeMap[standardKey][8] =     1;
          
          worksetSizeMap[nonAffineTensorKey][1] = 4096;
          worksetSizeMap[nonAffineTensorKey][2] =  512;
          worksetSizeMap[nonAffineTensorKey][3] =  128;
          worksetSizeMap[nonAffineTensorKey][4] =   32;
          worksetSizeMap[nonAffineTensorKey][5] =   16;
          worksetSizeMap[nonAffineTensorKey][6] =    1;
          worksetSizeMap[nonAffineTensorKey][7] =    1;
          worksetSizeMap[nonAffineTensorKey][8] =    1;
          
          worksetSizeMap[affineTensorKey][1] = 8192;
          worksetSizeMap[affineTensorKey][2] =  512;
          worksetSizeMap[affineTensorKey][3] =  128;
          worksetSizeMap[affineTensorKey][4] =   64;
          worksetSizeMap[affineTensorKey][5] =   16;
          worksetSizeMap[affineTensorKey][6] =    1;
          worksetSizeMap[affineTensorKey][7] =    1;
          worksetSizeMap[affineTensorKey][8] =    1;
        }
        {
          // Hdiv formulation
          FormulationChoice formulation = Hdiv;
          tuple<Mode,FormulationChoice,AlgorithmChoice> standardKey {mode,formulation,Standard};
          tuple<Mode,FormulationChoice,AlgorithmChoice> nonAffineTensorKey {mode,formulation,NonAffineTensor};
          tuple<Mode,FormulationChoice,AlgorithmChoice> affineTensorKey {mode,formulation,AffineTensor};
          
          // best for Hdiv - these are for meshes that range from 32768 for p=1 to 64 for p=8
          worksetSizeMap[standardKey][1] = 256;
          worksetSizeMap[standardKey][2] =  64;
          worksetSizeMap[standardKey][3] =  64;
          worksetSizeMap[standardKey][4] =  16;
          worksetSizeMap[standardKey][5] =   4;
          worksetSizeMap[standardKey][6] =   1;
          worksetSizeMap[standardKey][7] =   1;
          worksetSizeMap[standardKey][8] =   1;
          
          worksetSizeMap[nonAffineTensorKey][1] = 4096;
          worksetSizeMap[nonAffineTensorKey][2] =  256;
          worksetSizeMap[nonAffineTensorKey][3] =   64;
          worksetSizeMap[nonAffineTensorKey][4] =   16;
          worksetSizeMap[nonAffineTensorKey][5] =    4;
          worksetSizeMap[nonAffineTensorKey][6] =    1;
          worksetSizeMap[nonAffineTensorKey][7] =    1;
          worksetSizeMap[nonAffineTensorKey][8] =    1;
          
          worksetSizeMap[affineTensorKey][1] = 8192;
          worksetSizeMap[affineTensorKey][2] =  512;
          worksetSizeMap[affineTensorKey][3] =   64;
          worksetSizeMap[affineTensorKey][4] =   16;
          worksetSizeMap[affineTensorKey][5] =    8;
          worksetSizeMap[affineTensorKey][6] =    1;
          worksetSizeMap[affineTensorKey][7] =    1;
          worksetSizeMap[affineTensorKey][8] =    1;
        }
        {
          // Hcurl formulation
          FormulationChoice formulation = Hcurl;
          tuple<Mode,FormulationChoice,AlgorithmChoice> standardKey {mode,formulation,Standard};
          tuple<Mode,FormulationChoice,AlgorithmChoice> nonAffineTensorKey {mode,formulation,NonAffineTensor};
          tuple<Mode,FormulationChoice,AlgorithmChoice> affineTensorKey {mode,formulation,AffineTensor};
          
          // best for Hcurl - these are for meshes that range from 32768 for p=1 to 64 for p=8
          worksetSizeMap[standardKey][1] = 1024;
          worksetSizeMap[standardKey][2] =  512;
          worksetSizeMap[standardKey][3] =  256;
          worksetSizeMap[standardKey][4] =    4;
          worksetSizeMap[standardKey][5] =    1;
          worksetSizeMap[standardKey][6] =    1;
          worksetSizeMap[standardKey][7] =    1;
          worksetSizeMap[standardKey][8] =    1;
          
          worksetSizeMap[nonAffineTensorKey][1] = 512;
          worksetSizeMap[nonAffineTensorKey][2] =  64;
          worksetSizeMap[nonAffineTensorKey][3] =  16;
          worksetSizeMap[nonAffineTensorKey][4] =   4;
          worksetSizeMap[nonAffineTensorKey][5] =   1;
          worksetSizeMap[nonAffineTensorKey][6] =   1;
          worksetSizeMap[nonAffineTensorKey][7] =   1;
          worksetSizeMap[nonAffineTensorKey][8] =   1;
          
          worksetSizeMap[affineTensorKey][1] = 1024;
          worksetSizeMap[affineTensorKey][2] =  128;
          worksetSizeMap[affineTensorKey][3] =   16;
          worksetSizeMap[affineTensorKey][4] =    4;
          worksetSizeMap[affineTensorKey][5] =    1;
          worksetSizeMap[affineTensorKey][6] =    1;
          worksetSizeMap[affineTensorKey][7] =    1;
          worksetSizeMap[affineTensorKey][8] =    1;
        }
        {
          // L^2 formulation
          FormulationChoice formulation = L2;
          tuple<Mode,FormulationChoice,AlgorithmChoice> standardKey {mode,formulation,Standard};
          tuple<Mode,FormulationChoice,AlgorithmChoice> nonAffineTensorKey {mode,formulation,NonAffineTensor};
          tuple<Mode,FormulationChoice,AlgorithmChoice> affineTensorKey {mode,formulation,AffineTensor};
          
          // best for L^2 - these are for meshes that range from 32768 for p=1 to 256 for p=8
          worksetSizeMap[standardKey][1] = 1024;
          worksetSizeMap[standardKey][2] =  256;
          worksetSizeMap[standardKey][3] =   64;
          worksetSizeMap[standardKey][4] =   16;
          worksetSizeMap[standardKey][5] =   16;
          worksetSizeMap[standardKey][6] =   16;
          worksetSizeMap[standardKey][7] =    1;
          worksetSizeMap[standardKey][8] =    1;
          
          worksetSizeMap[nonAffineTensorKey][1] = 16384;
          worksetSizeMap[nonAffineTensorKey][2] =   512;
          worksetSizeMap[nonAffineTensorKey][3] =   256;
          worksetSizeMap[nonAffineTensorKey][4] =    64;
          worksetSizeMap[nonAffineTensorKey][5] =    16;
          worksetSizeMap[nonAffineTensorKey][6] =     8;
          worksetSizeMap[nonAffineTensorKey][7] =     2;
          worksetSizeMap[nonAffineTensorKey][8] =     1;
          
          worksetSizeMap[affineTensorKey][1] = 32768;
          worksetSizeMap[affineTensorKey][2] =  1024;
          worksetSizeMap[affineTensorKey][3] =   256;
          worksetSizeMap[affineTensorKey][4] =   128;
          worksetSizeMap[affineTensorKey][5] =    16;
          worksetSizeMap[affineTensorKey][6] =     8;
          worksetSizeMap[affineTensorKey][7] =     1;
          worksetSizeMap[affineTensorKey][8] =     1;
        }
        {
          // VectorWeightedPoisson
          // These calibrations were run 5-25-24 on an M2 Ultra, on a fork expected to be merged into Trilinos develop soon.
          FormulationChoice formulation = VectorWeightedPoisson;
          tuple<Mode,FormulationChoice,AlgorithmChoice> standardKey {mode,formulation,Standard};
          tuple<Mode,FormulationChoice,AlgorithmChoice> nonAffineTensorKey {mode,formulation,NonAffineTensor};
          tuple<Mode,FormulationChoice,AlgorithmChoice> affineTensorKey {mode,formulation,AffineTensor};
          
          // best for VectorWeightedPoisson - these are for meshes that range from 32,768 for p=1 to 128 for p=10
          worksetSizeMap[standardKey][1]  = 4096;
          worksetSizeMap[standardKey][2]  = 1024;
          worksetSizeMap[standardKey][3]  =   32;
          worksetSizeMap[standardKey][4]  =    4;
          worksetSizeMap[standardKey][5]  =    1;
          worksetSizeMap[standardKey][6]  =    1;
          worksetSizeMap[standardKey][7]  =    1;
          worksetSizeMap[standardKey][8]  =    1;
          worksetSizeMap[standardKey][9]  =    1;
          worksetSizeMap[standardKey][10] =    1;
          
          worksetSizeMap[nonAffineTensorKey][1]  = 2048;
          worksetSizeMap[nonAffineTensorKey][2]  = 2048;
          worksetSizeMap[nonAffineTensorKey][3]  = 128;
          worksetSizeMap[nonAffineTensorKey][4]  = 16;
          worksetSizeMap[nonAffineTensorKey][5]  = 2;
          worksetSizeMap[nonAffineTensorKey][6]  = 1;
          worksetSizeMap[nonAffineTensorKey][7]  = 1;
          worksetSizeMap[nonAffineTensorKey][8]  = 1;
          worksetSizeMap[nonAffineTensorKey][9]  = 1;
          worksetSizeMap[nonAffineTensorKey][10] = 1;
           
          worksetSizeMap[affineTensorKey][1]  = 32768;
          worksetSizeMap[affineTensorKey][2]  =  8192;
          worksetSizeMap[affineTensorKey][3]  =   128;
          worksetSizeMap[affineTensorKey][4]  =     8;
          worksetSizeMap[affineTensorKey][5]  =     2;
          worksetSizeMap[affineTensorKey][6]  =     1;
          worksetSizeMap[affineTensorKey][7]  =     1;
          worksetSizeMap[affineTensorKey][8]  =     1;
          worksetSizeMap[affineTensorKey][9]  =     1;
          worksetSizeMap[affineTensorKey][10] =     1;
        }
      } // BestSerial case
        break;
      case BestOpenMP_16:
      {
        // manually calibrated workset sizes on Mac Pro (2.5 GHz Xeon W, 28-core, running with OpenMP, OMP_NUM_THREADS=16)
        // Calibration for sum factorization cases was run while usePointCacheForRank3Tensor = true.
        
        
        // manually calibrated workset sizes on Mac Pro (2.5 GHz Xeon W, 28-core, running in serial)
        {
          // Poisson formulation
          FormulationChoice formulation = Poisson;
          tuple<Mode,FormulationChoice,AlgorithmChoice> standardKey {mode,formulation,Standard};
          tuple<Mode,FormulationChoice,AlgorithmChoice> nonAffineTensorKey {mode,formulation,NonAffineTensor};
          tuple<Mode,FormulationChoice,AlgorithmChoice> affineTensorKey {mode,formulation,AffineTensor};
          
          // best for Poisson - these are for meshes that range from 32768 for p=1 to 256 for p=8
          worksetSizeMap[standardKey][1] = 4096;
          worksetSizeMap[standardKey][2] = 2048;
          worksetSizeMap[standardKey][3] = 2048;
          worksetSizeMap[standardKey][4] = 2048;
          worksetSizeMap[standardKey][5] = 2048;
          worksetSizeMap[standardKey][6] = 2048;
          worksetSizeMap[standardKey][7] =    4;
          worksetSizeMap[standardKey][8] =    2;
          
          worksetSizeMap[nonAffineTensorKey][1] = 2048;
          worksetSizeMap[nonAffineTensorKey][2] =  512;
          worksetSizeMap[nonAffineTensorKey][3] =  256;
          worksetSizeMap[nonAffineTensorKey][4] =  128;
          worksetSizeMap[nonAffineTensorKey][5] =   64;
          worksetSizeMap[nonAffineTensorKey][6] =   32;
          worksetSizeMap[nonAffineTensorKey][7] =   16;
          worksetSizeMap[nonAffineTensorKey][8] =   16;
          
          worksetSizeMap[affineTensorKey][1] = 8192;
          worksetSizeMap[affineTensorKey][2] = 4096;
          worksetSizeMap[affineTensorKey][3] = 1024;
          worksetSizeMap[affineTensorKey][4] =  256;
          worksetSizeMap[affineTensorKey][5] =   64;
          worksetSizeMap[affineTensorKey][6] =   32;
          worksetSizeMap[affineTensorKey][7] =   16;
          worksetSizeMap[affineTensorKey][8] =   16;
        }
        {
          // Hgrad formulation
          FormulationChoice formulation = Hgrad;
          tuple<Mode,FormulationChoice,AlgorithmChoice> standardKey {mode,formulation,Standard};
          tuple<Mode,FormulationChoice,AlgorithmChoice> nonAffineTensorKey {mode,formulation,NonAffineTensor};
          tuple<Mode,FormulationChoice,AlgorithmChoice> affineTensorKey {mode,formulation,AffineTensor};
          
          // best for Hgrad - these are for meshes that range from 32768 for p=1 to 256 for p=8
          worksetSizeMap[standardKey][1] = 16384;
          worksetSizeMap[standardKey][2] =  8192;
          worksetSizeMap[standardKey][3] =  8192;
          worksetSizeMap[standardKey][4] =  2048;
          worksetSizeMap[standardKey][5] =   512;
          worksetSizeMap[standardKey][6] =   512;
          worksetSizeMap[standardKey][7] =   512;
          worksetSizeMap[standardKey][8] =     1;
          
          worksetSizeMap[nonAffineTensorKey][1] = 16384;
          worksetSizeMap[nonAffineTensorKey][2] =  8192;
          worksetSizeMap[nonAffineTensorKey][3] =   256;
          worksetSizeMap[nonAffineTensorKey][4] =   256;
          worksetSizeMap[nonAffineTensorKey][5] =    64;
          worksetSizeMap[nonAffineTensorKey][6] =    32;
          worksetSizeMap[nonAffineTensorKey][7] =    16;
          worksetSizeMap[nonAffineTensorKey][8] =    16;
          
          worksetSizeMap[affineTensorKey][1] = 8192;
          worksetSizeMap[affineTensorKey][2] = 4096;
          worksetSizeMap[affineTensorKey][3] = 1024;
          worksetSizeMap[affineTensorKey][4] =  256;
          worksetSizeMap[affineTensorKey][5] =   64;
          worksetSizeMap[affineTensorKey][6] =   32;
          worksetSizeMap[affineTensorKey][7] =   16;
          worksetSizeMap[affineTensorKey][8] =   16;
        }
        {
          // Hdiv formulation
          FormulationChoice formulation = Hdiv;
          tuple<Mode,FormulationChoice,AlgorithmChoice> standardKey {mode,formulation,Standard};
          tuple<Mode,FormulationChoice,AlgorithmChoice> nonAffineTensorKey {mode,formulation,NonAffineTensor};
          tuple<Mode,FormulationChoice,AlgorithmChoice> affineTensorKey {mode,formulation,AffineTensor};
          
          // best for Hdiv - these are for meshes that range from 32768 for p=1 to 64 for p=8
          worksetSizeMap[standardKey][1] = 32768;
          worksetSizeMap[standardKey][2] = 32768;
          worksetSizeMap[standardKey][3] =   512;
          worksetSizeMap[standardKey][4] =   256;
          worksetSizeMap[standardKey][5] =    64;
          worksetSizeMap[standardKey][6] =     2;
          worksetSizeMap[standardKey][7] =     2;
          worksetSizeMap[standardKey][8] =     1;
          
          worksetSizeMap[nonAffineTensorKey][1] = 32768;
          worksetSizeMap[nonAffineTensorKey][2] = 16384;
          worksetSizeMap[nonAffineTensorKey][3] =  8192;
          worksetSizeMap[nonAffineTensorKey][4] =    64;
          worksetSizeMap[nonAffineTensorKey][5] =    16;
          worksetSizeMap[nonAffineTensorKey][6] =    16;
          worksetSizeMap[nonAffineTensorKey][7] =    16;
          worksetSizeMap[nonAffineTensorKey][8] =    16;
          
          worksetSizeMap[affineTensorKey][1] = 16384;
          worksetSizeMap[affineTensorKey][2] =  4096;
          worksetSizeMap[affineTensorKey][3] =   256;
          worksetSizeMap[affineTensorKey][4] =   128;
          worksetSizeMap[affineTensorKey][5] =    64;
          worksetSizeMap[affineTensorKey][6] =    16;
          worksetSizeMap[affineTensorKey][7] =    16;
          worksetSizeMap[affineTensorKey][8] =    16;
        }
        {
          // Hcurl formulation
          FormulationChoice formulation = Hcurl;
          tuple<Mode,FormulationChoice,AlgorithmChoice> standardKey {mode,formulation,Standard};
          tuple<Mode,FormulationChoice,AlgorithmChoice> nonAffineTensorKey {mode,formulation,NonAffineTensor};
          tuple<Mode,FormulationChoice,AlgorithmChoice> affineTensorKey {mode,formulation,AffineTensor};
          
          // best for Hcurl - these are for meshes that range from 32768 for p=1 to 64 for p=8
          worksetSizeMap[standardKey][1] = 4096;
          worksetSizeMap[standardKey][2] =  128;
          worksetSizeMap[standardKey][3] =  128;
          worksetSizeMap[standardKey][4] =   32;
          worksetSizeMap[standardKey][5] =    4;
          worksetSizeMap[standardKey][6] =    1;
          worksetSizeMap[standardKey][7] =    1;
          worksetSizeMap[standardKey][8] =    1;
          
          worksetSizeMap[nonAffineTensorKey][1] = 16384;
          worksetSizeMap[nonAffineTensorKey][2] =   512;
          worksetSizeMap[nonAffineTensorKey][3] =   128;
          worksetSizeMap[nonAffineTensorKey][4] =    64;
          worksetSizeMap[nonAffineTensorKey][5] =    32;
          worksetSizeMap[nonAffineTensorKey][6] =    16;
          worksetSizeMap[nonAffineTensorKey][7] =    16;
          worksetSizeMap[nonAffineTensorKey][8] =    16;
          
          worksetSizeMap[affineTensorKey][1] = 32768;
          worksetSizeMap[affineTensorKey][2] =  4096;
          worksetSizeMap[affineTensorKey][3] =   128;
          worksetSizeMap[affineTensorKey][4] =    64;
          worksetSizeMap[affineTensorKey][5] =    16;
          worksetSizeMap[affineTensorKey][6] =    16;
          worksetSizeMap[affineTensorKey][7] =    16;
          worksetSizeMap[affineTensorKey][8] =    16;
        }
        {
          // L^2 formulation
          FormulationChoice formulation = L2;
          tuple<Mode,FormulationChoice,AlgorithmChoice> standardKey {mode,formulation,Standard};
          tuple<Mode,FormulationChoice,AlgorithmChoice> nonAffineTensorKey {mode,formulation,NonAffineTensor};
          tuple<Mode,FormulationChoice,AlgorithmChoice> affineTensorKey {mode,formulation,AffineTensor};
          
          // best for L^2 - these are for meshes that range from 32768 for p=1 to 256 for p=8
          worksetSizeMap[standardKey][1] = 8192;
          worksetSizeMap[standardKey][2] =  512;
          worksetSizeMap[standardKey][3] =   32;
          worksetSizeMap[standardKey][4] =   32;
          worksetSizeMap[standardKey][5] =   32;
          worksetSizeMap[standardKey][6] =    1;
          worksetSizeMap[standardKey][7] =    1;
          worksetSizeMap[standardKey][8] =    1;
          
          worksetSizeMap[nonAffineTensorKey][1] = 16384;
          worksetSizeMap[nonAffineTensorKey][2] =  4096;
          worksetSizeMap[nonAffineTensorKey][3] =  1024;
          worksetSizeMap[nonAffineTensorKey][4] =   256;
          worksetSizeMap[nonAffineTensorKey][5] =    64;
          worksetSizeMap[nonAffineTensorKey][6] =    32;
          worksetSizeMap[nonAffineTensorKey][7] =    16;
          worksetSizeMap[nonAffineTensorKey][8] =    16;
          
          worksetSizeMap[affineTensorKey][1] = 32768;
          worksetSizeMap[affineTensorKey][2] =  4096;
          worksetSizeMap[affineTensorKey][3] =  1024;
          worksetSizeMap[affineTensorKey][4] =   256;
          worksetSizeMap[affineTensorKey][5] =   128;
          worksetSizeMap[affineTensorKey][6] =    32;
          worksetSizeMap[affineTensorKey][7] =    16;
          worksetSizeMap[affineTensorKey][8] =    16;
        }
        {
          // VectorWeightedPoisson
          // These calibrations were run 5-25-24 on an M2 Ultra, on a fork expected to be merged into Trilinos develop soon.
          FormulationChoice formulation = VectorWeightedPoisson;
          tuple<Mode,FormulationChoice,AlgorithmChoice> standardKey {mode,formulation,Standard};
          tuple<Mode,FormulationChoice,AlgorithmChoice> nonAffineTensorKey {mode,formulation,NonAffineTensor};
          tuple<Mode,FormulationChoice,AlgorithmChoice> affineTensorKey {mode,formulation,AffineTensor};
          
          // best for VectorWeightedPoisson - these are for meshes that range from 32,768 for p=1 to 128 for p=10
          worksetSizeMap[standardKey][1]  = 16384;
          worksetSizeMap[standardKey][2]  = 16384;
          worksetSizeMap[standardKey][3]  =  8192;
          worksetSizeMap[standardKey][4]  =  1024;
          worksetSizeMap[standardKey][5]  =  1024;
          worksetSizeMap[standardKey][6]  =  1024;
          worksetSizeMap[standardKey][7]  =   512;
          worksetSizeMap[standardKey][8]  =   256;
          worksetSizeMap[standardKey][9]  =   128;
          worksetSizeMap[standardKey][10] =    32;
          
          worksetSizeMap[nonAffineTensorKey][1]  = 32768;
          worksetSizeMap[nonAffineTensorKey][2]  =  8192;
          worksetSizeMap[nonAffineTensorKey][3]  =  8192;
          worksetSizeMap[nonAffineTensorKey][4]  =  4096;
          worksetSizeMap[nonAffineTensorKey][5]  =  4096;
          worksetSizeMap[nonAffineTensorKey][6]  =    64;
          worksetSizeMap[nonAffineTensorKey][7]  =    32;
          worksetSizeMap[nonAffineTensorKey][8]  =    32;
          worksetSizeMap[nonAffineTensorKey][9]  =    16;
          worksetSizeMap[nonAffineTensorKey][10] =    16;
           
          worksetSizeMap[affineTensorKey][1]  = 32768;
          worksetSizeMap[affineTensorKey][2]  = 16384;
          worksetSizeMap[affineTensorKey][3]  =  8192;
          worksetSizeMap[affineTensorKey][4]  =  4096;
          worksetSizeMap[affineTensorKey][5]  =  4096;
          worksetSizeMap[affineTensorKey][6]  =  2048;
          worksetSizeMap[affineTensorKey][7]  =    32;
          worksetSizeMap[affineTensorKey][8]  =    16;
          worksetSizeMap[affineTensorKey][9]  =    16;
          worksetSizeMap[affineTensorKey][10] =    16;
        }
      } // BestOpenMP_16 case
        break;
      case BestCuda:
      {
        {
          // Poisson formulation
          FormulationChoice formulation = Poisson;
          tuple<Mode,FormulationChoice,AlgorithmChoice> standardKey {mode,formulation,Standard};
          tuple<Mode,FormulationChoice,AlgorithmChoice> nonAffineTensorKey {mode,formulation,NonAffineTensor};
          tuple<Mode,FormulationChoice,AlgorithmChoice> affineTensorKey {mode,formulation,AffineTensor};
          
          // best for Poisson - these are for meshes that range from 32768 for p=1 to 256 for p=8
          worksetSizeMap[standardKey][1] = 16384;
          worksetSizeMap[standardKey][2] =   512;
          worksetSizeMap[standardKey][3] =   128;
          worksetSizeMap[standardKey][4] =     8;
          worksetSizeMap[standardKey][5] =     4;
          worksetSizeMap[standardKey][6] =     1;
          worksetSizeMap[standardKey][7] =     1;
          worksetSizeMap[standardKey][8] =     1;
          
          worksetSizeMap[nonAffineTensorKey][1] = 32768;
          worksetSizeMap[nonAffineTensorKey][2] = 32768;
          worksetSizeMap[nonAffineTensorKey][3] = 16384;
          worksetSizeMap[nonAffineTensorKey][4] =  8192;
          worksetSizeMap[nonAffineTensorKey][5] =  4096;
          worksetSizeMap[nonAffineTensorKey][6] =  2048;
          worksetSizeMap[nonAffineTensorKey][7] =   256;
          worksetSizeMap[nonAffineTensorKey][8] =   256;
          
          worksetSizeMap[affineTensorKey][1] = 8192;
          worksetSizeMap[affineTensorKey][2] = 8192;
          worksetSizeMap[affineTensorKey][3] = 8192;
          worksetSizeMap[affineTensorKey][4] = 8192;
          worksetSizeMap[affineTensorKey][5] = 4096;
          worksetSizeMap[affineTensorKey][6] = 2048;
          worksetSizeMap[affineTensorKey][7] =  256;
          worksetSizeMap[affineTensorKey][8] =  128;
        }
        {
          // Hgrad formulation
          FormulationChoice formulation = Hgrad;
          tuple<Mode,FormulationChoice,AlgorithmChoice> standardKey {mode,formulation,Standard};
          tuple<Mode,FormulationChoice,AlgorithmChoice> nonAffineTensorKey {mode,formulation,NonAffineTensor};
          tuple<Mode,FormulationChoice,AlgorithmChoice> affineTensorKey {mode,formulation,AffineTensor};
          
          // best for Hgrad - these are for meshes that range from 32768 for p=1 to 256 for p=8
          worksetSizeMap[standardKey][1] = 32768;
          worksetSizeMap[standardKey][2] =   512;
          worksetSizeMap[standardKey][3] =   128;
          worksetSizeMap[standardKey][4] =    16;
          worksetSizeMap[standardKey][5] =     4;
          worksetSizeMap[standardKey][6] =     1;
          worksetSizeMap[standardKey][7] =     1;
          worksetSizeMap[standardKey][8] =     1;
          
          worksetSizeMap[nonAffineTensorKey][1] = 32768;
          worksetSizeMap[nonAffineTensorKey][2] = 32768;
          worksetSizeMap[nonAffineTensorKey][3] = 16384;
          worksetSizeMap[nonAffineTensorKey][4] =  8192;
          worksetSizeMap[nonAffineTensorKey][5] =  4096;
          worksetSizeMap[nonAffineTensorKey][6] =  2048;
          worksetSizeMap[nonAffineTensorKey][7] =   256;
          worksetSizeMap[nonAffineTensorKey][8] =   256;
          
          worksetSizeMap[affineTensorKey][1] = 32768;
          worksetSizeMap[affineTensorKey][2] = 32768;
          worksetSizeMap[affineTensorKey][3] =  8192;
          worksetSizeMap[affineTensorKey][4] =  8192;
          worksetSizeMap[affineTensorKey][5] =  4096;
          worksetSizeMap[affineTensorKey][6] =  2048;
          worksetSizeMap[affineTensorKey][7] =   256;
          worksetSizeMap[affineTensorKey][8] =   256;
        }
        {
          // Hdiv formulation
          FormulationChoice formulation = Hdiv;
          tuple<Mode,FormulationChoice,AlgorithmChoice> standardKey {mode,formulation,Standard};
          tuple<Mode,FormulationChoice,AlgorithmChoice> nonAffineTensorKey {mode,formulation,NonAffineTensor};
          tuple<Mode,FormulationChoice,AlgorithmChoice> affineTensorKey {mode,formulation,AffineTensor};
          
          // best for Hdiv - these are for meshes that range from 32768 for p=1 to 64 for p=8
          worksetSizeMap[standardKey][1] = 32768;
          worksetSizeMap[standardKey][2] =   512;
          worksetSizeMap[standardKey][3] =    32;
          worksetSizeMap[standardKey][4] =     4;
          worksetSizeMap[standardKey][5] =     1;
          worksetSizeMap[standardKey][6] =     1;
          worksetSizeMap[standardKey][7] =     1;
          worksetSizeMap[standardKey][8] =     1;
          
          worksetSizeMap[nonAffineTensorKey][1] = 32768;
          worksetSizeMap[nonAffineTensorKey][2] = 32768;
          worksetSizeMap[nonAffineTensorKey][3] = 16384;
          worksetSizeMap[nonAffineTensorKey][4] =  4096;
          worksetSizeMap[nonAffineTensorKey][5] =  1024;
          worksetSizeMap[nonAffineTensorKey][6] =   256;
          worksetSizeMap[nonAffineTensorKey][7] =   128;
          worksetSizeMap[nonAffineTensorKey][8] =    64;
          
          worksetSizeMap[affineTensorKey][1] = 32768;
          worksetSizeMap[affineTensorKey][2] = 32768;
          worksetSizeMap[affineTensorKey][3] = 16384;
          worksetSizeMap[affineTensorKey][4] =  4096;
          worksetSizeMap[affineTensorKey][5] =  1024;
          worksetSizeMap[affineTensorKey][6] =   256;
          worksetSizeMap[affineTensorKey][7] =   128;
          worksetSizeMap[affineTensorKey][8] =    64;
        }
        {
          // Hcurl formulation
          FormulationChoice formulation = Hcurl;
          tuple<Mode,FormulationChoice,AlgorithmChoice> standardKey {mode,formulation,Standard};
          tuple<Mode,FormulationChoice,AlgorithmChoice> nonAffineTensorKey {mode,formulation,NonAffineTensor};
          tuple<Mode,FormulationChoice,AlgorithmChoice> affineTensorKey {mode,formulation,AffineTensor};
          
          // best for Hcurl - these are for meshes that range from 32768 for p=1 to 64 for p=8
          worksetSizeMap[standardKey][1] = 1024;
          worksetSizeMap[standardKey][2] =  128;
          worksetSizeMap[standardKey][3] =   16;
          worksetSizeMap[standardKey][4] =    4;
          worksetSizeMap[standardKey][5] =    1;
          worksetSizeMap[standardKey][6] =    1;
          worksetSizeMap[standardKey][7] =    1;
          worksetSizeMap[standardKey][8] =    1;
          
          worksetSizeMap[nonAffineTensorKey][1] = 32768;
          worksetSizeMap[nonAffineTensorKey][2] = 32768;
          worksetSizeMap[nonAffineTensorKey][3] =  8192;
          worksetSizeMap[nonAffineTensorKey][4] =  2048;
          worksetSizeMap[nonAffineTensorKey][5] =   512;
          worksetSizeMap[nonAffineTensorKey][6] =   256;
          worksetSizeMap[nonAffineTensorKey][7] =   128;
          worksetSizeMap[nonAffineTensorKey][8] =    64;
          
          worksetSizeMap[affineTensorKey][1] = 32768;
          worksetSizeMap[affineTensorKey][2] = 32768;
          worksetSizeMap[affineTensorKey][3] =  8192;
          worksetSizeMap[affineTensorKey][4] =  2048;
          worksetSizeMap[affineTensorKey][5] =   512;
          worksetSizeMap[affineTensorKey][6] =   256;
          worksetSizeMap[affineTensorKey][7] =   128;
          worksetSizeMap[affineTensorKey][8] =    64;
        }
        {
          // L^2 formulation
          FormulationChoice formulation = L2;
          tuple<Mode,FormulationChoice,AlgorithmChoice> standardKey {mode,formulation,Standard};
          tuple<Mode,FormulationChoice,AlgorithmChoice> nonAffineTensorKey {mode,formulation,NonAffineTensor};
          tuple<Mode,FormulationChoice,AlgorithmChoice> affineTensorKey {mode,formulation,AffineTensor};
          
          // best for L^2 - these are for meshes that range from 32768 for p=1 to 256 for p=8
          worksetSizeMap[standardKey][1] = 32768;
          worksetSizeMap[standardKey][2] =  1024;
          worksetSizeMap[standardKey][3] =   128;
          worksetSizeMap[standardKey][4] =    16;
          worksetSizeMap[standardKey][5] =     4;
          worksetSizeMap[standardKey][6] =     1;
          worksetSizeMap[standardKey][7] =     1;
          worksetSizeMap[standardKey][8] =     1;
          
          worksetSizeMap[nonAffineTensorKey][1] = 32768;
          worksetSizeMap[nonAffineTensorKey][2] = 32768;
          worksetSizeMap[nonAffineTensorKey][3] = 16384;
          worksetSizeMap[nonAffineTensorKey][4] =  8192;
          worksetSizeMap[nonAffineTensorKey][5] =  4096;
          worksetSizeMap[nonAffineTensorKey][6] =  2048;
          worksetSizeMap[nonAffineTensorKey][7] =   256;
          worksetSizeMap[nonAffineTensorKey][8] =   128;
          
          worksetSizeMap[affineTensorKey][1] = 8192;
          worksetSizeMap[affineTensorKey][2] = 8192;
          worksetSizeMap[affineTensorKey][3] = 8192;
          worksetSizeMap[affineTensorKey][4] = 8192;
          worksetSizeMap[affineTensorKey][5] = 4096;
          worksetSizeMap[affineTensorKey][6] = 2048;
          worksetSizeMap[affineTensorKey][7] =  256;
          worksetSizeMap[affineTensorKey][8] =  128;
        } // L^2 formulation
        {
          // VectorWeightedPoisson
          // TODO: set this with some actual calibration result values.  For now, we just borrow from Poisson
          
          FormulationChoice formulation = VectorWeightedPoisson;
          tuple<Mode,FormulationChoice,AlgorithmChoice> standardKey {mode,formulation,Standard};
          tuple<Mode,FormulationChoice,AlgorithmChoice> nonAffineTensorKey {mode,formulation,NonAffineTensor};
          tuple<Mode,FormulationChoice,AlgorithmChoice> affineTensorKey {mode,formulation,AffineTensor};
          
          tuple<Mode,FormulationChoice,AlgorithmChoice> standardKey_Poisson {mode,Poisson,Standard};
          tuple<Mode,FormulationChoice,AlgorithmChoice> nonAffineTensorKey_Poisson {mode,Poisson,NonAffineTensor};
          tuple<Mode,FormulationChoice,AlgorithmChoice> affineTensorKey_Poisson {mode,Poisson,AffineTensor};
          
          worksetSizeMap[standardKey]        = worksetSizeMap[standardKey_Poisson];
          worksetSizeMap[nonAffineTensorKey] = worksetSizeMap[nonAffineTensorKey_Poisson];
          worksetSizeMap[affineTensorKey]    = worksetSizeMap[affineTensorKey_Poisson];
        }
    } // BestCuda case
        break;
      case Precalibrated:
      {
        if (precalibrationFile != "")
        {
          std::ifstream calibrationFileStream (precalibrationFile, std::ios::in);
          
          while (calibrationFileStream.good())
          {
            std::string formulationString, algorithmString;
            int polyOrder, worksetSize;
            
            std::string line;
            std::getline(calibrationFileStream, line, '\n');
            if (line == "") continue;
            std::istringstream linestream(line);
            linestream >> formulationString;
            linestream >> algorithmString;
            linestream >> polyOrder;
            linestream >> worksetSize;
            
            AlgorithmChoice algorithm     = algorithm_from_string(algorithmString);
            FormulationChoice formulation = formulation_from_string(formulationString);
            
            tuple<Mode,FormulationChoice,AlgorithmChoice> key { mode, formulation, algorithm };
            worksetSizeMap[key][polyOrder] = worksetSize;
          }
        }
      }
        break;
      default:
        std::cout << "WARNING: Unhandled mode.\n";
    } // mode switch
  } // mode for loop
  return worksetSizeMap;
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
    // Here, we support various 3D "formulations": Poisson, H^1, H(div), H(curl), and L^2 norms.  (Poisson is the default.)
    // The geometry we use is an axis-aligned, uniform hypercube mesh, but in the non-affine tensor case, the CellGeometry object does
    // not express the axis-aligned or uniform structure, and the computations are as they would be in a more general hypercube mesh.
    // Similarly, the Affine case only assumes an affine mesh; it need not be uniform or axis-aligned.
   
    string timingsFilePath;
    
    const int spaceDim = 3;
    
    Mode mode;
    
    vector<AlgorithmChoice> allAlgorithmChoices {Standard, NonAffineTensor, AffineTensor, Uniform};
    vector<FormulationChoice> allFormulationChoices {Poisson, Hgrad, Hdiv, Hcurl, L2};
    
    Teuchos::CommandLineProcessor cmdp(false,true); // false: don't throw exceptions; true: do return errors for unrecognized options
    
    string algorithmChoiceString = "All"; // alternatives: Standard, NonAffineTensor, AffineTensor, Uniform
    string formulationChoiceString = "Poisson";
    string basisFamilyChoiceString = "Nodal";
    
    int polyOrderFixed = -1;
    int polyOrderMin = 1;
    int polyOrderMax = 8;
    
    string modeChoiceString = "Test"; // alternatives: Calibration, BestSerial, BestCuda
    
    bool saveTimingsToFile = false;
    string outputDir = ".";
    string precalibrationFile = "";
    
    cmdp.setOption("algorithm", &algorithmChoiceString, "Options: All, Standard, NonAffineTensor, AffineTensor, Uniform");
    cmdp.setOption("formulation", &formulationChoiceString, "Options: Poisson, Hgrad, Hdiv, Hcurl, L2");
    cmdp.setOption("polyOrder", &polyOrderFixed, "Single polynomial degree to run at");
    cmdp.setOption("minPolyOrder", &polyOrderMin, "Starting polynomial degree to run at");
    cmdp.setOption("maxPolyOrder", &polyOrderMax, "Maximum polynomial degree to run at");
    cmdp.setOption("mode", &modeChoiceString);
    cmdp.setOption("basisFamily", &basisFamilyChoiceString, "Options: Nodal, Hierarchical, Serendipity");
    cmdp.setOption("saveTimings", "dontSaveTimings", &saveTimingsToFile, "Save timings to a file in outputDir.");
    cmdp.setOption("outputDir", &outputDir, "Directory for saving timings file");
    cmdp.setOption("precalibrationFile", &precalibrationFile, "File into which calibration output has been written");
    
    Teuchos::RCP<std::ofstream> timingsFileStream;
    
    bool success = true;
    
    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
    {
  #ifdef HAVE_MPI
      MPI_Finalize();
  #endif
      return -1;
    }
    
    if (saveTimingsToFile)
    {
      std::ostringstream fileNameStream;
      fileNameStream << outputDir << "/";
      fileNameStream << "timings_" << algorithmChoiceString << "_";
      fileNameStream << formulationChoiceString << "_";
      if (polyOrderFixed != -1)
      {
        fileNameStream << "p" << polyOrderFixed;
      }
      else
      {
        fileNameStream << "p" << polyOrderMin << "_to_p" << polyOrderMax << "_";
      }
      fileNameStream << modeChoiceString << "_";
      fileNameStream << basisFamilyChoiceString;
      fileNameStream << ".dat";
      
      timingsFilePath = fileNameStream.str();
      
      timingsFileStream = Teuchos::rcp( new std::ofstream(timingsFilePath, std::ios::out) );
      
      *timingsFileStream << "Algorithm\t";
      *timingsFileStream << "p\t";
      *timingsFileStream << "Element Count\t";
      *timingsFileStream << "Workset Size\t";
      *timingsFileStream << "mode\t";
      *timingsFileStream << "Basis Family\t";
      *timingsFileStream << "Core Integration Timing\t";
      *timingsFileStream << "Core Integration Flops\t";
      *timingsFileStream << "Core Integration Throughput\t";
      *timingsFileStream << "Jac. Timing\t";
      *timingsFileStream << "Jac. Flops\t";
      *timingsFileStream << "Jac. Throughput\t";
      *timingsFileStream << "Initialization Timing\t";
      *timingsFileStream << "Other Timing\t";
      *timingsFileStream << "Total Time\t";
      *timingsFileStream << "Total Flops\t";
      *timingsFileStream << "Total Throughput";
      *timingsFileStream << std::endl;
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
    
    Teuchos::RCP<Kokkos::Array<double,spaceDim>> vectorWeight1, vectorWeight2; // used for VectorWeightedPoisson
    vector<FormulationChoice> formulationChoices;
    if (formulationChoiceString == "All")
    {
      formulationChoices = allFormulationChoices;
    }
    else if (formulationChoiceString == "Poisson")
    {
      formulationChoices = vector<FormulationChoice>{Poisson};
    }
    else if (formulationChoiceString == "Hgrad")
    {
      formulationChoices = vector<FormulationChoice>{Hgrad};
    }
    else if (formulationChoiceString == "Hdiv")
    {
      formulationChoices = vector<FormulationChoice>{Hdiv};
    }
    else if (formulationChoiceString == "Hcurl")
    {
      formulationChoices = vector<FormulationChoice>{Hcurl};
    }
    else if (formulationChoiceString == "L2")
    {
      formulationChoices = vector<FormulationChoice>{L2};
    }
    else if (formulationChoiceString == "VectorWeightedPoisson")
    {
      formulationChoices = vector<FormulationChoice>{VectorWeightedPoisson};
      vectorWeight1 = Teuchos::rcp( new Kokkos::Array<double, spaceDim>() );
      vectorWeight2 = Teuchos::rcp( new Kokkos::Array<double, spaceDim>() );
      for (int d=0; d<spaceDim; d++)
      {
        (*vectorWeight1)[d] = 1.0;
        (*vectorWeight2)[d] = 1.0;
      }
    }
    else
    {
      cout << "Unrecognized formulation choice: " << formulationChoiceString << endl;
#ifdef HAVE_MPI
      MPI_Finalize();
#endif
      return -1;
    }
    
    vector<BasisFamilyChoice> basisFamilyChoices;
    if (basisFamilyChoiceString == "Nodal")
    {
      basisFamilyChoices = vector<BasisFamilyChoice>{Nodal};
    }
    else if (basisFamilyChoiceString == "Hierarchical")
    {
      basisFamilyChoices = vector<BasisFamilyChoice>{Hierarchical};
    }
    else if (basisFamilyChoiceString == "Serendipity")
    {
      basisFamilyChoices = vector<BasisFamilyChoice>{Serendipity};
    }
    else
    {
      cout << "Unrecognized basis family choice: " << basisFamilyChoiceString << endl;
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
    else if (modeChoiceString == "BestOpenMP_16")
    {
      mode = BestOpenMP_16;
    }
    else if (modeChoiceString == "Test")
    {
      mode = Test;
    }
    else if (modeChoiceString == "Precalibrated")
    {
      mode = Precalibrated;
    }
    else
    {
      cout << "Unrecognized mode choice: " << modeChoiceString << endl;
#ifdef HAVE_MPI
      MPI_Finalize();
#endif
      return -1;
    }
    
    if (mode == Test)
    {
      // then overwrite poly orders to match the test cases we define in getWorksetSizeMap
      polyOrderFixed = -1;
      polyOrderMin = 1;
      polyOrderMax = 3;
    }
    
    using Scalar = double;
    using ExecutionSpace = Kokkos::DefaultExecutionSpace;
    using DeviceType = Kokkos::DefaultExecutionSpace::device_type;
    
    using std::vector;
    using std::map;
    using std::pair;
    using std::make_pair;
    using std::tuple;
    using std::cout;
    using std::endl;
    using std::setw;
    
    using WorksetForAlgorithmChoice = map<AlgorithmChoice, int>;
    
    vector< tuple<int,FormulationChoice,Kokkos::Array<int,spaceDim>,WorksetForAlgorithmChoice> > polyOrderGridDimsWorksetTestCases;
    
    const int maxStiffnessGB = 2;
    const int maxEntryCount = maxStiffnessGB * 1024 * (1024 * 1024 / sizeof(Scalar));
    
    const int maxCellCount = 32768;
    
    map<int,int> cellCountForPolyOrder;
    
    map<int, Kokkos::Array<int,spaceDim> > gridCellCountsForPolyOrder;
    for (int p=polyOrderMin; p<=polyOrderMax; p++)
    {
      if (mode != Test)
      {
        int maxBasisCardinality = 1;
        for (auto formulationChoice : formulationChoices)
        {
          for (auto basisFamilyChoice : basisFamilyChoices)
          {
            auto basis = getHypercubeBasisForFormulation<Scalar, Scalar, DeviceType, spaceDim>(formulationChoice, basisFamilyChoice, p);
            maxBasisCardinality = std::max(basis->getCardinality(), maxBasisCardinality);
          }
        }
        gridCellCountsForPolyOrder[p] = getMeshWidths<spaceDim>(maxBasisCardinality, maxEntryCount, maxCellCount);
      }
      else // mode == Test; set the gridCellCounts according to mesh widths indicated in getWorksetSizes above.
      {
        Kokkos::Array<int,spaceDim> meshWidths_8, meshWidths_4;
        for (int d=0; d<spaceDim; d++)
        {
          meshWidths_4[d] = 4;
          meshWidths_8[d] = 8;
        }
        gridCellCountsForPolyOrder[1] = meshWidths_8;
        gridCellCountsForPolyOrder[2] = meshWidths_8;
        gridCellCountsForPolyOrder[3] = meshWidths_4;
      }
      
      int cellCount = 1;
      for (int d=0; d<spaceDim; d++)
      {
        cellCount *= gridCellCountsForPolyOrder[p][d];
      }
      
      cellCountForPolyOrder[p] = cellCount;
    }
    
    auto worksetSizeMap = getWorksetSizeMap(precalibrationFile, cellCountForPolyOrder, polyOrderMin, polyOrderMax);
    
    vector<int> worksetSizes {1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768};
    
    map<int,int> minWorksetSizeForPolyOrder;
    minWorksetSizeForPolyOrder[1] = 256;
    minWorksetSizeForPolyOrder[2] = 64;
    minWorksetSizeForPolyOrder[3] = 16;
    minWorksetSizeForPolyOrder[4] = 4;
    minWorksetSizeForPolyOrder[5] = 1;
    minWorksetSizeForPolyOrder[6] = 1;
    minWorksetSizeForPolyOrder[7] = 1;
    minWorksetSizeForPolyOrder[8] = 1;
    
    if (mode == Calibration)
    {
      for (auto formulationChoice : formulationChoices)
      {
        for (int polyOrder=polyOrderMin; polyOrder<=polyOrderMax; polyOrder++)
        {
          for (int worksetSize : worksetSizes)
          {
            if (worksetSize > cellCountForPolyOrder[polyOrder])
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
            auto gridDims = gridCellCountsForPolyOrder[polyOrder];
            polyOrderGridDimsWorksetTestCases.push_back(tuple<int,FormulationChoice,Kokkos::Array<int,spaceDim>,WorksetForAlgorithmChoice>{polyOrder,formulationChoice,gridDims,worksetForAlgorithmChoice} );
          }
        }
      }
    }
    else
    {
     std::vector<AlgorithmChoice> algorithmChoices {Standard,NonAffineTensor,AffineTensor,Uniform};
      WorksetForAlgorithmChoice worksetForAlgorithmChoice;
      for (auto formulationChoice : formulationChoices)
      {
        for (int polyOrder=polyOrderMin; polyOrder<=polyOrderMax; polyOrder++)
        {
          for (auto algorithmChoice : algorithmChoices)
          {
            worksetForAlgorithmChoice[algorithmChoice] = worksetSizeMap[{mode,formulationChoice,algorithmChoice}][polyOrder];
          }
          const auto & gridDims = gridCellCountsForPolyOrder[polyOrder];
          
          polyOrderGridDimsWorksetTestCases.push_back(tuple<int,FormulationChoice,Kokkos::Array<int,spaceDim>,WorksetForAlgorithmChoice>{polyOrder,formulationChoice,gridDims,worksetForAlgorithmChoice} );
        }
      }
    }
    
    cout << std::setprecision(2) << std::scientific;
    
    map< pair<FormulationChoice,AlgorithmChoice>, map<int, pair<double,int> > > maxAlgorithmThroughputForPolyOrderCore;  // values are (throughput in GFlops/sec, worksetSize)
    map< pair<FormulationChoice,AlgorithmChoice>, map<int, pair<double,int> > > maxAlgorithmThroughputForPolyOrderTotal; // values are (throughput in GFlops/sec, worksetSize)
    
    const int charWidth = 15;
    
    FormulationChoice previousFormulation = UnknownFormulation;
    for (auto basisFamilyChoice : basisFamilyChoices)
    {
      for (auto & testCase : polyOrderGridDimsWorksetTestCases)
      {
          int polyOrder       = std::get<0>(testCase);
          auto formulation    = std::get<1>(testCase);
          auto gridDims       = std::get<2>(testCase);
          auto worksetSizeMap = std::get<3>(testCase);
          if (formulation != previousFormulation)
          {
            std::cout << "\n\n***** Formulation: " << to_string(formulation) << " *******\n";
            previousFormulation = formulation;
          }
          std::cout << "\n\n";
          std::cout << "Running with polyOrder = " << polyOrder << ", mesh dims = ";
          for (int d=0; d<spaceDim; d++)
          {
            std::cout << gridDims[d];
            if (d < spaceDim - 1) std::cout << " x ";
            else                  std::cout << std::endl;
          }
          
          std::map<AlgorithmChoice, Intrepid2::ScalarView<Scalar,DeviceType> > assembledMatrices;
          for (auto algorithmChoice : algorithmChoices)
          {
            int worksetSize = 1;
            if (worksetSizeMap.find(algorithmChoice) != worksetSizeMap.end())
              worksetSize = worksetSizeMap[algorithmChoice];
            if (mode == Calibration)
            {
              // if this workset size is bigger than the optimal for p-1, skip it -- it's highly
              // unlikely that for a larger p, the optimal workset size will be *larger*.
              const auto & bestThroughputs = maxAlgorithmThroughputForPolyOrderCore[{formulation,algorithmChoice}];
              if (bestThroughputs.find(polyOrder-1) != bestThroughputs.end() )
              {
                int bestWorksetSize = bestThroughputs.find(polyOrder-1)->second.second;
                if (bestWorksetSize < worksetSize)
                {
                  continue;
                }
              }
            }
            auto geometry = getMesh<Scalar, spaceDim, ExecutionSpace>(algorithmChoice, gridDims);
            
            // timers recorded in performStructuredQuadratureGRADGRAD, performStandardQuadratureGRADGRAD
            auto jacobianAndCellMeasureTimer = Teuchos::TimeMonitor::getNewTimer("Jacobians");
            auto fstIntegrateCall = Teuchos::TimeMonitor::getNewTimer("transform + integrate()");
            auto initialSetupTimer = Teuchos::TimeMonitor::getNewTimer("Initial Setup");

            jacobianAndCellMeasureTimer->reset();
            fstIntegrateCall->reset();
            initialSetupTimer->reset();
            
            double elapsedTimeSeconds = 0;
            double jacobianCellMeasureFlopCount = 0;
            double transformIntegrateFlopCount = 0;
            
            Intrepid2::ScalarView<Scalar,DeviceType> assembledMatrix;
            if (algorithmChoice == Standard)
            {
              // each cell needs on the order of polyOrder^N quadrature points, each of which has a Jacobian of size N * N.
              auto timer = Teuchos::TimeMonitor::getNewTimer("Standard Integration");
              timer->start();
              switch (basisFamilyChoice)
              {
                case Nodal:
                {
                  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
                  assembledMatrix = performStandardQuadrature<Scalar,BasisFamily>(formulation, geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount, vectorWeight1, vectorWeight2);
                }
                  break;
                case Hierarchical:
                {
                  using BasisFamily = HierarchicalBasisFamily<DeviceType>;
                  assembledMatrix = performStandardQuadrature<Scalar,BasisFamily>(formulation, geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount, vectorWeight1, vectorWeight2);
                }
                  break;
                case Serendipity:
                  INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "basis family choice not yet implemented");
              }
              timer->stop();
              elapsedTimeSeconds = timer->totalElapsedTime();
              
              cout << "Standard, workset size:          " << setw(charWidth) << worksetSize << endl;
              
              timer->reset();
            }
            else if (algorithmChoice == AffineTensor)
            {
              auto timer = Teuchos::TimeMonitor::getNewTimer("Affine tensor Integration");
              timer->start();
              switch (basisFamilyChoice)
              {
                case Nodal:
                {
                  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
                  assembledMatrix = performStructuredQuadrature<Scalar,BasisFamily>(formulation, geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount, vectorWeight1, vectorWeight2);
                }
                  break;
                case Hierarchical:
                {
                  using BasisFamily = HierarchicalBasisFamily<DeviceType>;
                  assembledMatrix = performStructuredQuadrature<Scalar,BasisFamily>(formulation, geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount, vectorWeight1, vectorWeight2);
                }
                  break;
                case Serendipity:
                  INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "basis family choice not yet implemented");
              }
              timer->stop();
              
              elapsedTimeSeconds = timer->totalElapsedTime();
              
              cout << "Affine Tensor, workset size:     " << setw(charWidth) << worksetSize << endl;
                        
              timer->reset();
            }
            else if (algorithmChoice == NonAffineTensor)
            {
              auto timer = Teuchos::TimeMonitor::getNewTimer("Non-affine tensor Integration");
              timer->start();
              switch (basisFamilyChoice)
              {
                case Nodal:
                {
                  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
                  assembledMatrix = performStructuredQuadrature<Scalar,BasisFamily>(formulation, geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount, vectorWeight1, vectorWeight2);
                }
                  break;
                case Hierarchical:
                {
                  using BasisFamily = HierarchicalBasisFamily<DeviceType>;
                  assembledMatrix = performStructuredQuadrature<Scalar,BasisFamily>(formulation, geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount, vectorWeight1, vectorWeight2);
                }
                  break;
                case Serendipity:
                  INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "basis family choice not yet implemented");
              }
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
                numCells *= gridDims[d];
              }
              auto timer = Teuchos::TimeMonitor::getNewTimer("Uniform Integration");
              timer->start();
              switch (basisFamilyChoice)
              {
                case Nodal:
                {
                  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
                  assembledMatrix = performStructuredQuadrature<Scalar,BasisFamily>(formulation, geometry, polyOrder, numCells, transformIntegrateFlopCount, jacobianCellMeasureFlopCount, vectorWeight1, vectorWeight2);
                }
                  break;
                case Hierarchical:
                {
                  using BasisFamily = HierarchicalBasisFamily<DeviceType>;
                  assembledMatrix = performStructuredQuadrature<Scalar,BasisFamily>(formulation, geometry, polyOrder, numCells, transformIntegrateFlopCount, jacobianCellMeasureFlopCount, vectorWeight1, vectorWeight2);
                }
                  break;
                case Serendipity:
                  INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "basis family choice not yet implemented");
              }
              timer->stop();
              
              elapsedTimeSeconds = timer->totalElapsedTime();
              
              cout << "Uniform, workset size:           " << setw(charWidth) << worksetSize << endl;
              
              timer->reset();
            }
            
            assembledMatrices[algorithmChoice] = assembledMatrix;
            
            const double approximateFlopCountTotal = transformIntegrateFlopCount + jacobianCellMeasureFlopCount;
            const double overallThroughputInGFlops = approximateFlopCountTotal / elapsedTimeSeconds / 1.0e9;
            
            const double previousMaxThroughput = maxAlgorithmThroughputForPolyOrderTotal[{formulation,algorithmChoice}][polyOrder].first;
            if (overallThroughputInGFlops > previousMaxThroughput)
            {
              maxAlgorithmThroughputForPolyOrderTotal[{formulation,algorithmChoice}][polyOrder] = make_pair(overallThroughputInGFlops,worksetSize);
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
            
            const double previousMaxThroughputCore = maxAlgorithmThroughputForPolyOrderCore[{formulation,algorithmChoice}][polyOrder].first;
            if (transformIntegrateThroughputInGFlops > previousMaxThroughputCore)
            {
              maxAlgorithmThroughputForPolyOrderCore[{formulation,algorithmChoice}][polyOrder] = make_pair(transformIntegrateThroughputInGFlops,worksetSize);
            }
            
            cout << "Time (core integration)      " << setw(charWidth) << std::scientific << integrateCallTime << " seconds (" << std::fixed << integrateCallPercentage << "%)." << endl;
            cout << "flop estimate (core):        " << setw(charWidth) << std::scientific << transformIntegrateFlopCount << endl;
            cout << "estimated throughput (core): " << setw(charWidth) << std::scientific << transformIntegrateThroughputInGFlops << " GFlops" << endl;
            cout << "************************************************" << endl;
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
            
            if (saveTimingsToFile)
            {
              *timingsFileStream << std::scientific;
              
              *timingsFileStream << to_string(algorithmChoice) << "\t";
              *timingsFileStream << polyOrder << "\t";
              *timingsFileStream << geometry.numCells() << "\t";
              *timingsFileStream << worksetSize << "\t";
              *timingsFileStream << modeChoiceString << "\t";
              *timingsFileStream << basisFamilyChoiceString << "\t";
              *timingsFileStream << integrateCallTime << "\t";
              *timingsFileStream << transformIntegrateFlopCount << "\t";
              *timingsFileStream << transformIntegrateThroughputInGFlops << "\t";
              *timingsFileStream << jacobianTime << "\t";
              *timingsFileStream << jacobianCellMeasureFlopCount << "\t";
              *timingsFileStream << jacobiansThroughputInGFlops << "\t";
              *timingsFileStream << initialSetupTime << "\t";
              *timingsFileStream << remainingTime << "\t";
              *timingsFileStream << elapsedTimeSeconds << "\t";
              *timingsFileStream << approximateFlopCountTotal << "\t";
              *timingsFileStream << overallThroughputInGFlops;
              *timingsFileStream << std::endl;
            }
          }
          
          if (assembledMatrices.size() > 1)
          {
            // if we have multiple, then let's compare values to make sure they agree.
            Teuchos::basic_FancyOStream<char> out(Teuchos::rcp(&std::cout,false));
            const double relTol = 1e-10; // pretty loose tolerances are required, especially for higher-order hierarchical comparisons to Standard
            const double absTol = 1e-8;
            auto firstAlgorithm = assembledMatrices.begin()->first;
            auto firstMatrix = assembledMatrices.begin()->second;
            std::string algorithmName1 = to_string(firstAlgorithm);
            
            for (const auto &entry : assembledMatrices)
            {
              auto secondAlgorithm = entry.first;
              auto secondMatrix    = entry.second;
              std::string algorithmName2 = to_string(secondAlgorithm);
              testViewFloatingEquality(firstMatrix, secondMatrix, relTol, absTol, out, success, algorithmName1, algorithmName2);
//              printFunctor3(firstMatrix, std::cout, algorithmName1);
//              printFunctor3(secondMatrix, std::cout, algorithmName2);
            }
          }
        } // testCase for loop
        if (mode == Calibration)
        {
          std::ostringstream fileNameStream;
          fileNameStream << outputDir << "/";
          fileNameStream << "bestWorkset_";
          if (polyOrderFixed != -1)
          {
            fileNameStream << "p" << polyOrderFixed;
          }
          else
          {
            fileNameStream << "p" << polyOrderMin << "_to_p" << polyOrderMax << "_";
          }
          fileNameStream << basisFamilyChoiceString;
          fileNameStream << ".dat";

          auto calibrationFilePath = fileNameStream.str();
          auto calibrationFileStream = Teuchos::rcp( new std::ofstream(calibrationFilePath, std::ios::out) );

          cout << "Best workset sizes (as determined by 'core integration' throughput, which includes basis transforms, but not setup and/or Jacobian computations):\n";
          for (auto & formulation : formulationChoices)
          {
            for (auto & algorithmChoice : algorithmChoices)
            {
              if (algorithmChoice == Uniform) continue; // workset size is not meaningful for uniform (workset is always effectively one cell, or all cells, depending on how you choose to frame it).
              for (auto & maxThroughputEntry : maxAlgorithmThroughputForPolyOrderCore[{formulation,algorithmChoice}])
              {
                int polyOrder   = maxThroughputEntry.first;
                int worksetSize = maxThroughputEntry.second.second;
                *calibrationFileStream << to_string(formulation) << "\t" << to_string(algorithmChoice) << "\t" << polyOrder << "\t" << worksetSize << std::endl;
                
                double throughput = maxThroughputEntry.second.first;
                cout << to_string(formulation) << "\t" << to_string(algorithmChoice) << "\t" << polyOrder << "\t" << worksetSize << " (" << throughput << " GFlops/sec)\n";
              }
            }
          } // formulation loop
          
          calibrationFileStream->close();
        } // if (Calibration)
    } // basisFamilyChoices
    
    if (success)
    {
      return 0;
    }
    else
    {
      std::cout << "ERROR: Assembled matrices did *NOT* match across algorithms.\n";
      return -1;
    }
  }
}
