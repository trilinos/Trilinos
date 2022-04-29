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

enum FormulationChoice
{
  Poisson, // (grad, grad)
  Hgrad,   // (grad, grad) + (value, value)
  Hdiv,    // (div, div)   + (value, value)
  Hcurl,   // (curl, curl) + (value, value)
  L2       // (value, value)
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
    case Poisson: return "Poisson";
    case Hgrad:   return "Hgrad";
    case Hdiv:    return "Hdiv";
    case Hcurl:   return "Hcurl";
    case L2:      return "L2";
    
    default:      return "Unknown FormulationChoice";
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
      return Intrepid2::ScalarView<Scalar,DeviceType>();
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
      return Intrepid2::ScalarView<Scalar,DeviceType>();
  }
}

template<class Scalar, class BasisFamily, class PointScalar, int spaceDim>
typename BasisFamily::BasisPtr getBasisForFormulation(FormulationChoice formulation, shards::CellTopology &cellTopo, const int polyOrder)
{
  Intrepid2::EFunctionSpace fs;
  switch (formulation)
  {
    case Poisson: fs = FUNCTION_SPACE_HGRAD; break;
    case Hgrad:   fs = FUNCTION_SPACE_HGRAD; break;
    case Hdiv:    fs = FUNCTION_SPACE_HDIV;  break;
    case Hcurl:   fs = FUNCTION_SPACE_HCURL; break;
    case L2:      fs = FUNCTION_SPACE_HVOL;  break;
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
    
    enum Mode
    {
      Calibration,
      Test,
      BestSerial,
      BestOpenMP_16,
      BestCuda
    };
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
    
    cmdp.setOption("algorithm", &algorithmChoiceString, "Options: All, Standard, NonAffineTensor, AffineTensor, Uniform");
    cmdp.setOption("formulation", &formulationChoiceString, "Options: Poisson, Hgrad, Hdiv, Hcurl, L2");
    cmdp.setOption("polyOrder", &polyOrderFixed, "Single polynomial degree to run at");
    cmdp.setOption("minPolyOrder", &polyOrderMin, "Starting polynomial degree to run at");
    cmdp.setOption("maxPolyOrder", &polyOrderMax, "Maximum polynomial degree to run at");
    cmdp.setOption("mode", &modeChoiceString);
    cmdp.setOption("basisFamily", &basisFamilyChoiceString, "Options: Nodal, Hierarchical, Serendipity");
    cmdp.setOption("saveTimings", "dontSaveTimings", &saveTimingsToFile, "Save timings to a file in outputDir.");
    cmdp.setOption("outputDir", &outputDir, "Directory for saving timings file");
    
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
    
    vector< tuple<int,Kokkos::Array<int,spaceDim>,WorksetForAlgorithmChoice> > polyOrderGridDimsWorksetTestCases;
    
    const int maxStiffnessGB = 2;
    const int maxEntryCount = maxStiffnessGB * 1024 * (1024 * 1024 / sizeof(Scalar));
    
    const int maxCellCount = 32768;
    
    map<int,int> cellCountForPolyOrder;
    
    map<int, Kokkos::Array<int,spaceDim> > gridCellCountsForPolyOrder;
    for (int p=polyOrderMin; p<=polyOrderMax; p++)
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
      
      int cellCount = 1;
      for (int d=0; d<spaceDim; d++)
      {
        cellCount *= gridCellCountsForPolyOrder[p][d];
      }
      
      cellCountForPolyOrder[p] = cellCount;
    }
    
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
    
    switch (mode) {
      case Calibration:
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
            polyOrderGridDimsWorksetTestCases.push_back(tuple<int,Kokkos::Array<int,spaceDim>,WorksetForAlgorithmChoice>{polyOrder,gridDims,worksetForAlgorithmChoice} );
          }
        }
      }
      break;
      case Test:
      {
        // DEBUGGING -- for ease of debugging, a single test case.
//        vector< tuple<int,int,int> > testCases { tuple<int,int,int> {1,1,1} };
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
          Kokkos::Array<int,spaceDim> gridDims;
          for (int d=0; d<spaceDim; d++)
          {
            gridDims[d] = meshWidth;
          }
          polyOrderGridDimsWorksetTestCases.push_back(tuple<int,Kokkos::Array<int,spaceDim>,WorksetForAlgorithmChoice>{polyOrder,gridDims,worksetForAlgorithmChoice} );
        }
      }
      break;
      case BestSerial:
      {
        if (formulationChoices.size() != 1)
        {
          std::cout << "BestSerial mode is not supported when running multiple formulations.\n";
          exit(-1);
        }
        
        auto formulationChoice = formulationChoices[0];
        
        // manually calibrated workset sizes on Mac Pro (2.5 GHz Xeon W, 28-core, running in serial)
        
        map<int,int> standardWorksetForPolyOrder;
        map<int,int> nonAffineTensorWorksetForPolyOrder;
        map<int,int> affineTensorWorksetForPolyOrder;
        
        switch(formulationChoice)
        {
          case Poisson:
          {
            // best for Poisson - these are for meshes that range from 32768 for p=1 to 256 for p=8
            standardWorksetForPolyOrder[1] = 8192;
            standardWorksetForPolyOrder[2] = 4096;
            standardWorksetForPolyOrder[3] =   64;
            standardWorksetForPolyOrder[4] =   16;
            standardWorksetForPolyOrder[5] =   16;
            standardWorksetForPolyOrder[6] =    1;
            standardWorksetForPolyOrder[7] =    1;
            standardWorksetForPolyOrder[8] =    1;
            
            nonAffineTensorWorksetForPolyOrder[1] = 2048;
            nonAffineTensorWorksetForPolyOrder[2] =  256;
            nonAffineTensorWorksetForPolyOrder[3] =  128;
            nonAffineTensorWorksetForPolyOrder[4] =   16;
            nonAffineTensorWorksetForPolyOrder[5] =    2;
            nonAffineTensorWorksetForPolyOrder[6] =    1;
            nonAffineTensorWorksetForPolyOrder[7] =    1;
            nonAffineTensorWorksetForPolyOrder[8] =    1;
            
            affineTensorWorksetForPolyOrder[1] = 4096;
            affineTensorWorksetForPolyOrder[2] =   64;
            affineTensorWorksetForPolyOrder[3] =   32;
            affineTensorWorksetForPolyOrder[4] =    4;
            affineTensorWorksetForPolyOrder[5] =    2;
            affineTensorWorksetForPolyOrder[6] =    1;
            affineTensorWorksetForPolyOrder[7] =    1;
            affineTensorWorksetForPolyOrder[8] =    1;
          }
            break;
          case Hgrad:
          {
            // best for Hgrad - these are for meshes that range from 32768 for p=1 to 256 for p=8
            standardWorksetForPolyOrder[1] = 32768;
            standardWorksetForPolyOrder[2] = 16384;
            standardWorksetForPolyOrder[3] =   512;
            standardWorksetForPolyOrder[4] =   512;
            standardWorksetForPolyOrder[5] =   512;
            standardWorksetForPolyOrder[6] =     2;
            standardWorksetForPolyOrder[7] =     1;
            standardWorksetForPolyOrder[8] =     1;
            
            nonAffineTensorWorksetForPolyOrder[1] = 4096;
            nonAffineTensorWorksetForPolyOrder[2] =  512;
            nonAffineTensorWorksetForPolyOrder[3] =  128;
            nonAffineTensorWorksetForPolyOrder[4] =   32;
            nonAffineTensorWorksetForPolyOrder[5] =   16;
            nonAffineTensorWorksetForPolyOrder[6] =    1;
            nonAffineTensorWorksetForPolyOrder[7] =    1;
            nonAffineTensorWorksetForPolyOrder[8] =    1;
            
            affineTensorWorksetForPolyOrder[1] = 8192;
            affineTensorWorksetForPolyOrder[2] =  512;
            affineTensorWorksetForPolyOrder[3] =  128;
            affineTensorWorksetForPolyOrder[4] =   64;
            affineTensorWorksetForPolyOrder[5] =   16;
            affineTensorWorksetForPolyOrder[6] =    1;
            affineTensorWorksetForPolyOrder[7] =    1;
            affineTensorWorksetForPolyOrder[8] =    1;
          }
            break;
          case Hdiv:
          {
            // best for Hdiv - these are for meshes that range from 32768 for p=1 to 64 for p=8
            standardWorksetForPolyOrder[1] = 256;
            standardWorksetForPolyOrder[2] =  64;
            standardWorksetForPolyOrder[3] =  64;
            standardWorksetForPolyOrder[4] =  16;
            standardWorksetForPolyOrder[5] =   4;
            standardWorksetForPolyOrder[6] =   1;
            standardWorksetForPolyOrder[7] =   1;
            standardWorksetForPolyOrder[8] =   1;
            
            nonAffineTensorWorksetForPolyOrder[1] = 4096;
            nonAffineTensorWorksetForPolyOrder[2] =  256;
            nonAffineTensorWorksetForPolyOrder[3] =   64;
            nonAffineTensorWorksetForPolyOrder[4] =   16;
            nonAffineTensorWorksetForPolyOrder[5] =    4;
            nonAffineTensorWorksetForPolyOrder[6] =    1;
            nonAffineTensorWorksetForPolyOrder[7] =    1;
            nonAffineTensorWorksetForPolyOrder[8] =    1;
            
            affineTensorWorksetForPolyOrder[1] = 8192;
            affineTensorWorksetForPolyOrder[2] =  512;
            affineTensorWorksetForPolyOrder[3] =   64;
            affineTensorWorksetForPolyOrder[4] =   16;
            affineTensorWorksetForPolyOrder[5] =    8;
            affineTensorWorksetForPolyOrder[6] =    1;
            affineTensorWorksetForPolyOrder[7] =    1;
            affineTensorWorksetForPolyOrder[8] =    1;
          }
            break;
          case Hcurl:
          {
            // best for Hcurl - these are for meshes that range from 32768 for p=1 to 64 for p=8
            standardWorksetForPolyOrder[1] = 1024;
            standardWorksetForPolyOrder[2] =  512;
            standardWorksetForPolyOrder[3] =  256;
            standardWorksetForPolyOrder[4] =    4;
            standardWorksetForPolyOrder[5] =    1;
            standardWorksetForPolyOrder[6] =    1;
            standardWorksetForPolyOrder[7] =    1;
            standardWorksetForPolyOrder[8] =    1;
            
            nonAffineTensorWorksetForPolyOrder[1] = 512;
            nonAffineTensorWorksetForPolyOrder[2] =  64;
            nonAffineTensorWorksetForPolyOrder[3] =  16;
            nonAffineTensorWorksetForPolyOrder[4] =   4;
            nonAffineTensorWorksetForPolyOrder[5] =   1;
            nonAffineTensorWorksetForPolyOrder[6] =   1;
            nonAffineTensorWorksetForPolyOrder[7] =   1;
            nonAffineTensorWorksetForPolyOrder[8] =   1;
            
            affineTensorWorksetForPolyOrder[1] = 1024;
            affineTensorWorksetForPolyOrder[2] =  128;
            affineTensorWorksetForPolyOrder[3] =   16;
            affineTensorWorksetForPolyOrder[4] =    4;
            affineTensorWorksetForPolyOrder[5] =    1;
            affineTensorWorksetForPolyOrder[6] =    1;
            affineTensorWorksetForPolyOrder[7] =    1;
            affineTensorWorksetForPolyOrder[8] =    1;
          }
            break;
          case L2:
          {
            // best for L^2 - these are for meshes that range from 32768 for p=1 to 256 for p=8
            standardWorksetForPolyOrder[1] = 1024;
            standardWorksetForPolyOrder[2] =  256;
            standardWorksetForPolyOrder[3] =   64;
            standardWorksetForPolyOrder[4] =   16;
            standardWorksetForPolyOrder[5] =   16;
            standardWorksetForPolyOrder[6] =   16;
            standardWorksetForPolyOrder[7] =    1;
            standardWorksetForPolyOrder[8] =    1;
            
            nonAffineTensorWorksetForPolyOrder[1] = 16384;
            nonAffineTensorWorksetForPolyOrder[2] =   512;
            nonAffineTensorWorksetForPolyOrder[3] =   256;
            nonAffineTensorWorksetForPolyOrder[4] =    64;
            nonAffineTensorWorksetForPolyOrder[5] =    16;
            nonAffineTensorWorksetForPolyOrder[6] =     8;
            nonAffineTensorWorksetForPolyOrder[7] =     2;
            nonAffineTensorWorksetForPolyOrder[8] =     1;
            
            affineTensorWorksetForPolyOrder[1] = 32768;
            affineTensorWorksetForPolyOrder[2] =  1024;
            affineTensorWorksetForPolyOrder[3] =   256;
            affineTensorWorksetForPolyOrder[4] =   128;
            affineTensorWorksetForPolyOrder[5] =    16;
            affineTensorWorksetForPolyOrder[6] =     8;
            affineTensorWorksetForPolyOrder[7] =     1;
            affineTensorWorksetForPolyOrder[8] =     1;
          }
            break;
        }
        
        // for the cases that we have not tried yet (polyOrder > 8), we try to choose sensible guesses for workset size:
        // 1 is best, we think, for polyOrder 8, so it'll be the best for the rest.
        int worksetSize = 1;
        for (int polyOrder=9; polyOrder <= polyOrderMax; polyOrder++)
        {
          nonAffineTensorWorksetForPolyOrder[polyOrder] = worksetSize;
          affineTensorWorksetForPolyOrder[polyOrder]    = worksetSize;
          standardWorksetForPolyOrder[polyOrder]        = worksetSize;
        }
        
        for (int polyOrder=polyOrderMin; polyOrder<=polyOrderMax; polyOrder++)
        {
          WorksetForAlgorithmChoice worksetForAlgorithmChoice;
          worksetForAlgorithmChoice[Standard]        = standardWorksetForPolyOrder       [polyOrder];
          worksetForAlgorithmChoice[NonAffineTensor] = nonAffineTensorWorksetForPolyOrder[polyOrder];
          worksetForAlgorithmChoice[AffineTensor]    = affineTensorWorksetForPolyOrder[polyOrder];
          worksetForAlgorithmChoice[Uniform]         = cellCountForPolyOrder[polyOrder];
          
          const auto & gridDims = gridCellCountsForPolyOrder[polyOrder];
          
          polyOrderGridDimsWorksetTestCases.push_back(tuple<int,Kokkos::Array<int,spaceDim>,WorksetForAlgorithmChoice>{polyOrder,gridDims,worksetForAlgorithmChoice} );
        }
      }
        break;
      
      case BestOpenMP_16:
      {
        if (formulationChoices.size() != 1)
        {
          std::cout << "BestOpenMP_16 mode is not supported when running multiple formulations.\n";
          exit(-1);
        }
        
        auto formulationChoice = formulationChoices[0];
        
        // manually calibrated workset sizes on Mac Pro (2.5 GHz Xeon W, 28-core, running with OpenMP, OMP_NUM_THREADS=16)
        // Calibration for sum factorization cases was run while usePointCacheForRank3Tensor = true.
        
        map<int,int> standardWorksetForPolyOrder;
        map<int,int> nonAffineTensorWorksetForPolyOrder;
        map<int,int> affineTensorWorksetForPolyOrder;
        
        switch(formulationChoice)
        {
          case Poisson:
          {
            // best for Poisson - these are for meshes that range from 32768 for p=1 to 256 for p=8
            standardWorksetForPolyOrder[1] = 4096;
            standardWorksetForPolyOrder[2] = 2048;
            standardWorksetForPolyOrder[3] = 2048;
            standardWorksetForPolyOrder[4] = 2048;
            standardWorksetForPolyOrder[5] = 2048;
            standardWorksetForPolyOrder[6] = 2048;
            standardWorksetForPolyOrder[7] =    4;
            standardWorksetForPolyOrder[8] =    2;
            
            nonAffineTensorWorksetForPolyOrder[1] = 2048;
            nonAffineTensorWorksetForPolyOrder[2] =  512;
            nonAffineTensorWorksetForPolyOrder[3] =  256;
            nonAffineTensorWorksetForPolyOrder[4] =  128;
            nonAffineTensorWorksetForPolyOrder[5] =   64;
            nonAffineTensorWorksetForPolyOrder[6] =   32;
            nonAffineTensorWorksetForPolyOrder[7] =   16;
            nonAffineTensorWorksetForPolyOrder[8] =   16;
            
            affineTensorWorksetForPolyOrder[1] = 8192;
            affineTensorWorksetForPolyOrder[2] = 4096;
            affineTensorWorksetForPolyOrder[3] = 1024;
            affineTensorWorksetForPolyOrder[4] =  256;
            affineTensorWorksetForPolyOrder[5] =   64;
            affineTensorWorksetForPolyOrder[6] =   32;
            affineTensorWorksetForPolyOrder[7] =   16;
            affineTensorWorksetForPolyOrder[8] =   16;
          }
            break;
          case Hgrad:
          {
            // best for Hgrad - these are for meshes that range from 32768 for p=1 to 256 for p=8
            standardWorksetForPolyOrder[1] = 16384;
            standardWorksetForPolyOrder[2] =  8192;
            standardWorksetForPolyOrder[3] =  8192;
            standardWorksetForPolyOrder[4] =  2048;
            standardWorksetForPolyOrder[5] =   512;
            standardWorksetForPolyOrder[6] =   512;
            standardWorksetForPolyOrder[7] =   512;
            standardWorksetForPolyOrder[8] =     1;
            
            nonAffineTensorWorksetForPolyOrder[1] = 16384;
            nonAffineTensorWorksetForPolyOrder[2] =  8192;
            nonAffineTensorWorksetForPolyOrder[3] =   256;
            nonAffineTensorWorksetForPolyOrder[4] =   256;
            nonAffineTensorWorksetForPolyOrder[5] =    64;
            nonAffineTensorWorksetForPolyOrder[6] =    32;
            nonAffineTensorWorksetForPolyOrder[7] =    16;
            nonAffineTensorWorksetForPolyOrder[8] =    16;
            
            affineTensorWorksetForPolyOrder[1] =  8192;
            affineTensorWorksetForPolyOrder[2] =  4096;
            affineTensorWorksetForPolyOrder[3] =  1024;
            affineTensorWorksetForPolyOrder[4] =   256;
            affineTensorWorksetForPolyOrder[5] =    64;
            affineTensorWorksetForPolyOrder[6] =    32;
            affineTensorWorksetForPolyOrder[7] =    16;
            affineTensorWorksetForPolyOrder[8] =    16;
          }
            break;
          case Hdiv:
          {
            // best for Hdiv - these are for meshes that range from 32768 for p=1 to 64 for p=8
            standardWorksetForPolyOrder[1] = 32768;
            standardWorksetForPolyOrder[2] = 32768;
            standardWorksetForPolyOrder[3] =   512;
            standardWorksetForPolyOrder[4] =   256;
            standardWorksetForPolyOrder[5] =    64;
            standardWorksetForPolyOrder[6] =     2;
            standardWorksetForPolyOrder[7] =     2;
            standardWorksetForPolyOrder[8] =     1;
            
            nonAffineTensorWorksetForPolyOrder[1] = 32768;
            nonAffineTensorWorksetForPolyOrder[2] = 16384;
            nonAffineTensorWorksetForPolyOrder[3] =  8192;
            nonAffineTensorWorksetForPolyOrder[4] =    64;
            nonAffineTensorWorksetForPolyOrder[5] =    16;
            nonAffineTensorWorksetForPolyOrder[6] =    16;
            nonAffineTensorWorksetForPolyOrder[7] =    16;
            nonAffineTensorWorksetForPolyOrder[8] =    16;
            
            affineTensorWorksetForPolyOrder[1] = 16384;
            affineTensorWorksetForPolyOrder[2] =  4096;
            affineTensorWorksetForPolyOrder[3] =   256;
            affineTensorWorksetForPolyOrder[4] =   128;
            affineTensorWorksetForPolyOrder[5] =    64;
            affineTensorWorksetForPolyOrder[6] =    16;
            affineTensorWorksetForPolyOrder[7] =    16;
            affineTensorWorksetForPolyOrder[8] =    16;
          }
            break;
          case Hcurl:
          {
            // best for Hcurl - these are for meshes that range from 32768 for p=1 to 64 for p=8
            standardWorksetForPolyOrder[1] = 4096;
            standardWorksetForPolyOrder[2] =  128;
            standardWorksetForPolyOrder[3] =  128;
            standardWorksetForPolyOrder[4] =   32;
            standardWorksetForPolyOrder[5] =    4;
            standardWorksetForPolyOrder[6] =    1;
            standardWorksetForPolyOrder[7] =    1;
            standardWorksetForPolyOrder[8] =    1;
            
            nonAffineTensorWorksetForPolyOrder[1] = 16384;
            nonAffineTensorWorksetForPolyOrder[2] =   512;
            nonAffineTensorWorksetForPolyOrder[3] =   128;
            nonAffineTensorWorksetForPolyOrder[4] =    64;
            nonAffineTensorWorksetForPolyOrder[5] =    32;
            nonAffineTensorWorksetForPolyOrder[6] =    16;
            nonAffineTensorWorksetForPolyOrder[7] =    16;
            nonAffineTensorWorksetForPolyOrder[8] =    16;
            
            affineTensorWorksetForPolyOrder[1] = 32768;
            affineTensorWorksetForPolyOrder[2] =  4096;
            affineTensorWorksetForPolyOrder[3] =   128;
            affineTensorWorksetForPolyOrder[4] =    64;
            affineTensorWorksetForPolyOrder[5] =    16;
            affineTensorWorksetForPolyOrder[6] =    16;
            affineTensorWorksetForPolyOrder[7] =    16;
            affineTensorWorksetForPolyOrder[8] =    16;
          }
            break;
          case L2:
          {
            // best for L^2 - these are for meshes that range from 32768 for p=1 to 256 for p=8
            standardWorksetForPolyOrder[1] = 8192;
            standardWorksetForPolyOrder[2] =  512;
            standardWorksetForPolyOrder[3] =   32;
            standardWorksetForPolyOrder[4] =   32;
            standardWorksetForPolyOrder[5] =   32;
            standardWorksetForPolyOrder[6] =    1;
            standardWorksetForPolyOrder[7] =    1;
            standardWorksetForPolyOrder[8] =    1;
            
            nonAffineTensorWorksetForPolyOrder[1] = 16384;
            nonAffineTensorWorksetForPolyOrder[2] =  4096;
            nonAffineTensorWorksetForPolyOrder[3] =  1024;
            nonAffineTensorWorksetForPolyOrder[4] =   256;
            nonAffineTensorWorksetForPolyOrder[5] =    64;
            nonAffineTensorWorksetForPolyOrder[6] =    32;
            nonAffineTensorWorksetForPolyOrder[7] =    16;
            nonAffineTensorWorksetForPolyOrder[8] =    16;
            
            affineTensorWorksetForPolyOrder[1] = 32768;
            affineTensorWorksetForPolyOrder[2] =  4096;
            affineTensorWorksetForPolyOrder[3] =  1024;
            affineTensorWorksetForPolyOrder[4] =   256;
            affineTensorWorksetForPolyOrder[5] =   128;
            affineTensorWorksetForPolyOrder[6] =    32;
            affineTensorWorksetForPolyOrder[7] =    16;
            affineTensorWorksetForPolyOrder[8] =    16;
          }
            break;
        }
        
        // for the cases that we have not tried yet (polyOrder > 8), we try to choose sensible guesses for workset size:
        // Standard: 1 is best for polyOrder 8, so it'll be the best for the rest.
        // NonAffineTensor, AffineTensor: we seem to bottom out at the number of OpenMP threads: here, 16.
        int standardWorksetSize = 1;
        int tensorWorksetSize = 16;
        for (int polyOrder=9; polyOrder <= polyOrderMax; polyOrder++)
        {
          nonAffineTensorWorksetForPolyOrder[polyOrder] = tensorWorksetSize;
          affineTensorWorksetForPolyOrder[polyOrder]    = tensorWorksetSize;
          standardWorksetForPolyOrder[polyOrder]        = standardWorksetSize;
        }
        
        for (int polyOrder=polyOrderMin; polyOrder<=polyOrderMax; polyOrder++)
        {
          WorksetForAlgorithmChoice worksetForAlgorithmChoice;
          worksetForAlgorithmChoice[Standard]        = standardWorksetForPolyOrder       [polyOrder];
          worksetForAlgorithmChoice[NonAffineTensor] = nonAffineTensorWorksetForPolyOrder[polyOrder];
          worksetForAlgorithmChoice[AffineTensor]    = affineTensorWorksetForPolyOrder[polyOrder];
          worksetForAlgorithmChoice[Uniform]         = cellCountForPolyOrder[polyOrder];
          const auto & gridDims = gridCellCountsForPolyOrder[polyOrder];
          
          polyOrderGridDimsWorksetTestCases.push_back(tuple<int,Kokkos::Array<int,spaceDim>,WorksetForAlgorithmChoice>{polyOrder,gridDims,worksetForAlgorithmChoice} );
        }
      }
        break;
      case BestCuda:
      {
        {
          if (formulationChoices.size() != 1)
          {
            std::cout << "BestCuda mode is not supported when running multiple formulations.\n";
            exit(-1);
          }
          
          auto formulationChoice = formulationChoices[0];
          
          // STANDARD
          // manually calibrated workset size on P100 (weaver)
          map<int,int> standardWorksetForPolyOrder;
          map<int,int> nonAffineTensorWorksetForPolyOrder;
          map<int,int> affineTensorWorksetForPolyOrder;
          
          switch(formulationChoice)
          {
            case Poisson:
            {
              // best for Poisson - these are for meshes that range from 32768 for p=1 to 256 for p=8
              standardWorksetForPolyOrder[1] = 16384;
              standardWorksetForPolyOrder[2] =   512;
              standardWorksetForPolyOrder[3] =   128;
              standardWorksetForPolyOrder[4] =     8;
              standardWorksetForPolyOrder[5] =     4;
              standardWorksetForPolyOrder[6] =     1;
              standardWorksetForPolyOrder[7] =     1;
              standardWorksetForPolyOrder[8] =     1;
              
              nonAffineTensorWorksetForPolyOrder[1] = 32768;
              nonAffineTensorWorksetForPolyOrder[2] = 32768;
              nonAffineTensorWorksetForPolyOrder[3] = 16384;
              nonAffineTensorWorksetForPolyOrder[4] =  8192;
              nonAffineTensorWorksetForPolyOrder[5] =  4096;
              nonAffineTensorWorksetForPolyOrder[6] =  2048;
              nonAffineTensorWorksetForPolyOrder[7] =   256;
              nonAffineTensorWorksetForPolyOrder[8] =   256;
              
              affineTensorWorksetForPolyOrder[1] = 8192;
              affineTensorWorksetForPolyOrder[2] = 8192;
              affineTensorWorksetForPolyOrder[3] = 8192;
              affineTensorWorksetForPolyOrder[4] = 8192;
              affineTensorWorksetForPolyOrder[5] = 4096;
              affineTensorWorksetForPolyOrder[6] = 2048;
              affineTensorWorksetForPolyOrder[7] =  256;
              affineTensorWorksetForPolyOrder[8] =  128;
            }
              break;
            case Hgrad:
            {
              // best for Hgrad - these are for meshes that range from 32768 for p=1 to 256 for p=8
              standardWorksetForPolyOrder[1] = 32768;
              standardWorksetForPolyOrder[2] =   512;
              standardWorksetForPolyOrder[3] =   128;
              standardWorksetForPolyOrder[4] =    16;
              standardWorksetForPolyOrder[5] =     4;
              standardWorksetForPolyOrder[6] =     1;
              standardWorksetForPolyOrder[7] =     1;
              standardWorksetForPolyOrder[8] =     1;
              
              nonAffineTensorWorksetForPolyOrder[1] = 32768;
              nonAffineTensorWorksetForPolyOrder[2] = 32768;
              nonAffineTensorWorksetForPolyOrder[3] = 16384;
              nonAffineTensorWorksetForPolyOrder[4] =  8192;
              nonAffineTensorWorksetForPolyOrder[5] =  4096;
              nonAffineTensorWorksetForPolyOrder[6] =  2048;
              nonAffineTensorWorksetForPolyOrder[7] =   256;
              nonAffineTensorWorksetForPolyOrder[8] =   256;
              
              affineTensorWorksetForPolyOrder[1] = 32768;
              affineTensorWorksetForPolyOrder[2] = 32768;
              affineTensorWorksetForPolyOrder[3] =  8192;
              affineTensorWorksetForPolyOrder[4] =  8192;
              affineTensorWorksetForPolyOrder[5] =  4096;
              affineTensorWorksetForPolyOrder[6] =  2048;
              affineTensorWorksetForPolyOrder[7] =   256;
              affineTensorWorksetForPolyOrder[8] =   256;
            }
              break;
            case Hdiv:
            {
              // best for Hdiv - these are for meshes that range from 32768 for p=1 to 64 for p=8
              standardWorksetForPolyOrder[1] = 32768;
              standardWorksetForPolyOrder[2] =   512;
              standardWorksetForPolyOrder[3] =    32;
              standardWorksetForPolyOrder[4] =     4;
              standardWorksetForPolyOrder[5] =     1;
              standardWorksetForPolyOrder[6] =     1;
              standardWorksetForPolyOrder[7] =     1;
              standardWorksetForPolyOrder[8] =     1;
              
              nonAffineTensorWorksetForPolyOrder[1] = 32768;
              nonAffineTensorWorksetForPolyOrder[2] = 32768;
              nonAffineTensorWorksetForPolyOrder[3] = 16384;
              nonAffineTensorWorksetForPolyOrder[4] =  4096;
              nonAffineTensorWorksetForPolyOrder[5] =  1024;
              nonAffineTensorWorksetForPolyOrder[6] =   256;
              nonAffineTensorWorksetForPolyOrder[7] =   128;
              nonAffineTensorWorksetForPolyOrder[8] =    64;
              
              affineTensorWorksetForPolyOrder[1] = 32768;
              affineTensorWorksetForPolyOrder[2] = 32768;
              affineTensorWorksetForPolyOrder[3] = 16384;
              affineTensorWorksetForPolyOrder[4] =  4096;
              affineTensorWorksetForPolyOrder[5] =  1024;
              affineTensorWorksetForPolyOrder[6] =   256;
              affineTensorWorksetForPolyOrder[7] =   128;
              affineTensorWorksetForPolyOrder[8] =    64;
            }
              break;
            case Hcurl:
            {
              standardWorksetForPolyOrder[1] = 1024;
              standardWorksetForPolyOrder[2] =  128;
              standardWorksetForPolyOrder[3] =   16;
              standardWorksetForPolyOrder[4] =    4;
              standardWorksetForPolyOrder[5] =    1;
              standardWorksetForPolyOrder[6] =    1;
              standardWorksetForPolyOrder[7] =    1;
              standardWorksetForPolyOrder[8] =    1;
              
              nonAffineTensorWorksetForPolyOrder[1] = 32768;
              nonAffineTensorWorksetForPolyOrder[2] = 32768;
              nonAffineTensorWorksetForPolyOrder[3] =  8192;
              nonAffineTensorWorksetForPolyOrder[4] =  2048;
              nonAffineTensorWorksetForPolyOrder[5] =   512;
              nonAffineTensorWorksetForPolyOrder[6] =   256;
              nonAffineTensorWorksetForPolyOrder[7] =   128;
              nonAffineTensorWorksetForPolyOrder[8] =    64;
              
              affineTensorWorksetForPolyOrder[1] = 32768;
              affineTensorWorksetForPolyOrder[2] = 32768;
              affineTensorWorksetForPolyOrder[3] =  8192;
              affineTensorWorksetForPolyOrder[4] =  2048;
              affineTensorWorksetForPolyOrder[5] =   512;
              affineTensorWorksetForPolyOrder[6] =   256;
              affineTensorWorksetForPolyOrder[7] =   128;
              affineTensorWorksetForPolyOrder[8] =    64;
            }
              break;
            case L2:
            {
              standardWorksetForPolyOrder[1] = 32768;
              standardWorksetForPolyOrder[2] =  1024;
              standardWorksetForPolyOrder[3] =   128;
              standardWorksetForPolyOrder[4] =    16;
              standardWorksetForPolyOrder[5] =     4;
              standardWorksetForPolyOrder[6] =     1;
              standardWorksetForPolyOrder[7] =     1;
              standardWorksetForPolyOrder[8] =     1;
              
              nonAffineTensorWorksetForPolyOrder[1] = 32768;
              nonAffineTensorWorksetForPolyOrder[2] = 32768;
              nonAffineTensorWorksetForPolyOrder[3] = 16384;
              nonAffineTensorWorksetForPolyOrder[4] =  8192;
              nonAffineTensorWorksetForPolyOrder[5] =  4096;
              nonAffineTensorWorksetForPolyOrder[6] =  2048;
              nonAffineTensorWorksetForPolyOrder[7] =   256;
              nonAffineTensorWorksetForPolyOrder[8] =   128;
              
              affineTensorWorksetForPolyOrder[1] = 8192;
              affineTensorWorksetForPolyOrder[2] = 8192;
              affineTensorWorksetForPolyOrder[3] = 8192;
              affineTensorWorksetForPolyOrder[4] = 8192;
              affineTensorWorksetForPolyOrder[5] = 4096;
              affineTensorWorksetForPolyOrder[6] = 2048;
              affineTensorWorksetForPolyOrder[7] =  256;
              affineTensorWorksetForPolyOrder[8] =  128;
            }
              break;
          }
          
          // for the cases that we have not tried yet (polyOrder > 8), we try to choose sensible guesses for workset size:
          int standardWorksetSize  = 1;  // 1 is best for polyOrder 8, so it'll be the best for the rest.
          // for the rest under CUDA, we observe that in most cases, the optimal workset size for non-affine tensor is the cell count.  For affine, it's lower by a factor of 2 or 4 in most cases.
          for (int polyOrder=9; polyOrder <= polyOrderMax; polyOrder++)
          {
            nonAffineTensorWorksetForPolyOrder[polyOrder] = cellCountForPolyOrder[polyOrder];
            affineTensorWorksetForPolyOrder[polyOrder]    = cellCountForPolyOrder[polyOrder] / 2;
            standardWorksetForPolyOrder[polyOrder] = standardWorksetSize;
          }
          
          for (int polyOrder=polyOrderMin; polyOrder<=polyOrderMax; polyOrder++)
          {
            WorksetForAlgorithmChoice worksetForAlgorithmChoice;
            worksetForAlgorithmChoice[Standard]        = standardWorksetForPolyOrder       [polyOrder];
            worksetForAlgorithmChoice[NonAffineTensor] = nonAffineTensorWorksetForPolyOrder[polyOrder];
            worksetForAlgorithmChoice[AffineTensor]    = affineTensorWorksetForPolyOrder[polyOrder];
            worksetForAlgorithmChoice[Uniform]         = cellCountForPolyOrder[polyOrder];
            
            const auto & gridDims = gridCellCountsForPolyOrder[polyOrder];
            
            polyOrderGridDimsWorksetTestCases.push_back(tuple<int,Kokkos::Array<int,spaceDim>,WorksetForAlgorithmChoice>{polyOrder,gridDims,worksetForAlgorithmChoice} );
          }
        }
        break;
        
      default:
        break;
    }
    }
    
    cout << std::setprecision(2) << std::scientific;
    
    map< AlgorithmChoice, map<int, pair<double,int> > > maxAlgorithmThroughputForPolyOrderCore;  // values are (throughput in GFlops/sec, worksetSize)
    map< AlgorithmChoice, map<int, pair<double,int> > > maxAlgorithmThroughputForPolyOrderTotal; // values are (throughput in GFlops/sec, worksetSize)
    
    const int charWidth = 15;
    
    for (auto basisFamilyChoice : basisFamilyChoices)
    {
      for (auto formulation : formulationChoices)
      {
        std::cout << "\n\n***** Formulation: " << to_string(formulation) << " *******\n";
        for (auto & testCase : polyOrderGridDimsWorksetTestCases)
        {
          int polyOrder       = std::get<0>(testCase);
          auto gridDims       = std::get<1>(testCase);
          auto worksetSizeMap = std::get<2>(testCase);
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
            int worksetSize = worksetSizeMap[algorithmChoice];
            if (mode == Calibration)
            {
              // if this workset size is bigger than the optimal for p-1, skip it -- it's highly
              // unlikely that for a larger p, the optimal workset size will be *larger*.
              const auto & bestThroughputs = maxAlgorithmThroughputForPolyOrderCore[algorithmChoice];
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
                  assembledMatrix = performStandardQuadrature<Scalar,BasisFamily>(formulation, geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
                }
                  break;
                case Hierarchical:
                {
                  using BasisFamily = HierarchicalBasisFamily<DeviceType>;
                  assembledMatrix = performStandardQuadrature<Scalar,BasisFamily>(formulation, geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
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
                  assembledMatrix = performStructuredQuadrature<Scalar,BasisFamily>(formulation, geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
                }
                  break;
                case Hierarchical:
                {
                  using BasisFamily = HierarchicalBasisFamily<DeviceType>;
                  assembledMatrix = performStructuredQuadrature<Scalar,BasisFamily>(formulation, geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
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
                  assembledMatrix = performStructuredQuadrature<Scalar,BasisFamily>(formulation, geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
                }
                  break;
                case Hierarchical:
                {
                  using BasisFamily = HierarchicalBasisFamily<DeviceType>;
                  assembledMatrix = performStructuredQuadrature<Scalar,BasisFamily>(formulation, geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
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
                  assembledMatrix = performStructuredQuadrature<Scalar,BasisFamily>(formulation, geometry, polyOrder, numCells, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
                }
                  break;
                case Hierarchical:
                {
                  using BasisFamily = HierarchicalBasisFamily<DeviceType>;
                  assembledMatrix = performStructuredQuadrature<Scalar,BasisFamily>(formulation, geometry, polyOrder, numCells, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
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
            
            const double previousMaxThroughput = maxAlgorithmThroughputForPolyOrderTotal[algorithmChoice][polyOrder].first;
            if (overallThroughputInGFlops > previousMaxThroughput)
            {
              maxAlgorithmThroughputForPolyOrderTotal[algorithmChoice][polyOrder] = make_pair(overallThroughputInGFlops,worksetSize);
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
            
            const double previousMaxThroughputCore = maxAlgorithmThroughputForPolyOrderCore[algorithmChoice][polyOrder].first;
            if (transformIntegrateThroughputInGFlops > previousMaxThroughputCore)
            {
              maxAlgorithmThroughputForPolyOrderCore[algorithmChoice][polyOrder] = make_pair(transformIntegrateThroughputInGFlops,worksetSize);
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
        }
      }
    } // basisFamilyChoices
    
    if (mode == Calibration)
    {
      cout << "Best workset sizes (as determined by 'core integration' throughput, which includes basis transforms, but not setup and/or Jacobian computations):\n";
      for (auto & algorithmChoice : algorithmChoices)
      {
        if (algorithmChoice == Uniform) continue; // workset size is not meaningful for uniform (workset is always effectively one cell, or all cells, depending on how you choose to frame it).
        
        cout << "Best workset sizes for " << to_string(algorithmChoice) << ":" << endl;
        
        for (auto & maxThroughputEntry : maxAlgorithmThroughputForPolyOrderCore[algorithmChoice])
        {
          int polyOrder   = maxThroughputEntry.first;
          int worksetSize = maxThroughputEntry.second.second;
          double throughput = maxThroughputEntry.second.first;
          cout << "p = " << polyOrder << ":" << setw(5) << worksetSize << " (" << throughput << " GFlops/sec)\n";
        }
      }
    }
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
