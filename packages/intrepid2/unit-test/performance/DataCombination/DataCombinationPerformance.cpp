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

/** \file   DataCombinationPerformance.cpp
    \brief  Main for performance tests comparing performance when combining Intrepid2 Data objects (as sums and products) with the performance of (expanded) Kokkos View objects.
 
 Specifically, we consider a few use cases, each with nominal shape (C,P):
 1. Constant data.  This case favors Data objects most heavily, since redundancy in the Views will be maximal.
 2. "Affine" data.  This has shape (C,P), but only varies in the cell dimension.
 3. General data.  There is no redundancy in the data.  This case favors the View objects most heavily, and will maximally expose overhead from the Data implementation.
 
 We fix the cell count at 16,000, and allow the point count to vary.  For the first two cases, we expect there to be some threshold point count at which the data redundancy is sufficient to
 overcome the overhead of Data; for all higher point counts, we expect performance savings for Data.  For the third case, we expect (with a "smart" Data combination implementation) that the
 relative overhead of Data will diminish as the point count goes up.
 
 In addition to combinations of "like" Data (e.g., constant plus constant), we can also test combinations of "unlike" Data (e.g., constant plus affine).  We expect these tests to have performance
 characteristics somewhere between the corresponding "like" tests.
 
 After measuring timings, we also confirm that the two algorithms agree on the results.
 */

#include "Teuchos_GlobalMPISession.hpp"

#include "Teuchos_StackedTimer.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_DefaultComm.hpp"

#include "Kokkos_Core.hpp"

#include "Intrepid2_Data.hpp"
#include "Intrepid2_TestUtils.hpp"
#include "Intrepid2_Types.hpp"

enum CaseChoice
{
  Constant,
  Affine,
  General
};

std::string to_string(CaseChoice choice)
{
  switch (choice) {
    case Constant: return "Constant";
    case Affine:   return "Affine";
    case General:  return "General";
    
    default:       return "Unknown CaseChoice";
  }
}

using namespace Intrepid2;

static const int NUM_CELLS = 16000;

template< typename Scalar, typename DeviceType >
inline
Data<Scalar, DeviceType> getData(CaseChoice caseChoice, const int numPoints, const double baseValue)
{
  using ExecutionSpace = typename DeviceType::execution_space;
  const int numCells = NUM_CELLS;
  Kokkos::Array<ordinal_type,2> extents {numCells, numPoints};
  Kokkos::Array<DataVariationType,2> variationTypes {GENERAL,GENERAL};
  
  switch (caseChoice) {
    case Constant:
      return Data<Scalar, DeviceType>(baseValue,extents);
    case Affine:
    {
      // (C,P); varies in C dimension
      variationTypes[1] = CONSTANT;
      Kokkos::View<Scalar*,DeviceType> cellView("affine case - underlying view",numCells);
      Kokkos::RangePolicy<ExecutionSpace> policy(ExecutionSpace(), 0, numCells);
      Kokkos::parallel_for("initialize underlying view data", policy,
      KOKKOS_LAMBDA (const int &i0) {
        cellView(i0) = i0 * baseValue;
      });
      return Data<Scalar, DeviceType>(cellView,extents,variationTypes);
    }
    case General:
    {
      // (C,P); varies in C and P dimensions
      variationTypes[1] = GENERAL;
      Kokkos::View<Scalar**,DeviceType> cellView("affine case - underlying view",numCells,numPoints);
      Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>> policy({0,0},{numCells,numPoints});
      Kokkos::parallel_for("initialize underlying view data", policy,
      KOKKOS_LAMBDA (const int &i0, const int &i1) {
        cellView(i0,i1) = i0 * baseValue + i1;
      });
      return Data<Scalar, DeviceType>(cellView,extents,variationTypes);
    }
    default:
      break;
  }
}

double theoreticalSpeedup(CaseChoice caseChoice, const int numPoints)
{
  switch (caseChoice) {
    case Constant:
      return NUM_CELLS * numPoints;
    case Affine:
      return numPoints;
    case General:
      return 1.0;
    default:
      break;
  }
}

template< typename Scalar, typename DeviceType >
Kokkos::View<Scalar**, DeviceType> allocateView(const int numPoints)
{
  Kokkos::View<Scalar**,DeviceType> view("DataCombinationPerformance - View", NUM_CELLS, numPoints);
  return view;
}

template< typename Scalar, typename DeviceType >
inline
void fillView(CaseChoice caseChoice, Kokkos::View<Scalar**,DeviceType> view, const double baseValue)
{
  using ExecutionSpace = typename DeviceType::execution_space;
  const int numCells  = view.extent_int(0);
  const int numPoints = view.extent_int(1);
  
  switch (caseChoice) {
    case Constant:
      Kokkos::deep_copy(view, baseValue);
    case Affine:
    {
      Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>> policy({0,0},{numCells,numPoints});
      // (C,P); varies in C dimension
      Kokkos::parallel_for("initialize underlying view data", policy,
      KOKKOS_LAMBDA (const int &i0, const int &i1) {
        view(i0,i1) = i0 * baseValue;
      });
    }
    case General:
    {
      Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>> policy({0,0},{numCells,numPoints});
      // (C,P); varies in C and P dimensions
      Kokkos::parallel_for("initialize underlying view data", policy,
      KOKKOS_LAMBDA (const int &i0, const int &i1) {
        view(i0,i1) = i0 * baseValue + i1;
      });
    }
    default:
      break;
  }
  ExecutionSpace().fence();
}

template< typename Scalar, typename DeviceType >
void sumViews(Kokkos::View<Scalar**,DeviceType> resultView,
              Kokkos::View<Scalar**,DeviceType> view1, Kokkos::View<Scalar**,DeviceType> view2)
{
  using ExecutionSpace = typename DeviceType::execution_space;
  const int numCells  = resultView.extent_int(0);
  const int numPoints = resultView.extent_int(1);
  Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>> policy({0,0},{numCells,numPoints});
  
  Kokkos::parallel_for("initialize underlying view data", policy,
  KOKKOS_LAMBDA (const int &i0, const int &i1) {
    resultView(i0,i1) = view1(i0,i1) + view2(i0,i1);
  });
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
    vector<CaseChoice> allCaseChoices {Constant, Affine, General};
    
    Teuchos::CommandLineProcessor cmdp(false,true); // false: don't throw exceptions; true: do return errors for unrecognized options
    
    string caseChoiceString = "All"; // alternatives: Standard, NonAffineTensor, AffineTensor, Uniform
    
    int pointCountFixed = -1;
    int pointCountMin = 64;
    int pointCountMax = 4096;
    
    cmdp.setOption("case", &caseChoiceString, "Options: All, Constant, Affine, General");
    cmdp.setOption("pointCount", &pointCountFixed, "Single point count to run with");
    cmdp.setOption("minPointCount", &pointCountMin, "Starting point count (will double until max count is reached)");
    cmdp.setOption("maxPointCount", &pointCountMax, "Maximum point count");
    
    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
    {
  #ifdef HAVE_MPI
      MPI_Finalize();
  #endif
      return -1;
    }

    vector<CaseChoice> caseChoices;
    if (caseChoiceString == "All")
    {
      caseChoices = allCaseChoices;
    }
    else if (caseChoiceString == "Constant")
    {
      caseChoices = vector<CaseChoice>{Constant};
    }
    else if (caseChoiceString == "Affine")
    {
      caseChoices = vector<CaseChoice>{Affine};
    }
    else if (caseChoiceString == "General")
    {
      caseChoices = vector<CaseChoice>{General};
    }
    else
    {
      cout << "Unrecognized case choice: " << caseChoiceString << endl;
#ifdef HAVE_MPI
      MPI_Finalize();
#endif
      return -1;
    }
    
    if (pointCountFixed > 0)
    {
      pointCountMin = pointCountFixed;
      pointCountMax = pointCountFixed;
    }
    
    using Scalar = double;
    using DeviceType = Kokkos::DefaultExecutionSpace::device_type;
    
    using DataType = Data<Scalar, DeviceType>;
    
    const int charWidth = 15;
    using std::vector;
    using std::map;
    using std::pair;
    using std::make_pair;
    using std::tuple;
    using std::cout;
    using std::endl;
    using std::setw;
    using std::scientific;
    using std::fixed;
    
    for (CaseChoice caseChoice : caseChoices)
    {
      cout << "\n\n***************************************************\n";
      cout <<     "******   " << setw(20) << to_string(caseChoice) << setw(23) << "   ******\n";
      cout << "***************************************************\n";
      for (int pointCount=pointCountMin; pointCount<=pointCountMax; pointCount *= 2)
      {
        double baseValue1 = M_PI;
        auto data1 = getData<Scalar, DeviceType>(caseChoice, pointCount, baseValue1);
        
        double baseValue2 = 1.0;
        auto data2 = getData<Scalar, DeviceType>(caseChoice, pointCount, baseValue2);
        
        auto result = DataType::allocateInPlaceCombinationResult(data1, data2);
        
        DeviceType::execution_space().fence();
        auto dataTimer = Teuchos::TimeMonitor::getNewTimer("Data sum");
        dataTimer->start();
        result.storeInPlaceSum(data1, data2);
        DeviceType::execution_space().fence();
        dataTimer->stop();
        double dataElapsedTimeSeconds = dataTimer->totalElapsedTime();
        
        cout << "Point count:          " << setw(charWidth) << pointCount << endl;
        cout << "Time (sum - data):    " << setw(charWidth) << std::setprecision(2) << scientific << dataElapsedTimeSeconds << endl;
        
        dataTimer->reset();
        
        auto viewTimer = Teuchos::TimeMonitor::getNewTimer("View sum");
        auto view1 = allocateView<Scalar, DeviceType>(pointCount);
        auto view2 = allocateView<Scalar, DeviceType>(pointCount);
        auto resultView = allocateView<Scalar, DeviceType>(pointCount);
        
        fillView(caseChoice, view1, baseValue1);
        fillView(caseChoice, view2, baseValue2);
        
        DeviceType::execution_space().fence();
        viewTimer->start();
        sumViews(resultView, view1, view2);
        DeviceType::execution_space().fence();
        viewTimer->stop();
        double viewElapsedTimeSeconds = viewTimer->totalElapsedTime();
        cout << "Time (sum - view):    " << setw(charWidth) << std::setprecision(2) << scientific << viewElapsedTimeSeconds << endl;
        
        viewTimer->reset();
        
        const double maxSpeedup = theoreticalSpeedup(caseChoice, pointCount);
        const double actualSpeedup = viewElapsedTimeSeconds / dataElapsedTimeSeconds;
        const double percentage = actualSpeedup / maxSpeedup * 100.0;
        cout << "Ideal speedup:        " << setw(charWidth) << std::setprecision(2) << scientific << maxSpeedup << endl;
        cout << "Actual speedup:       " << setw(charWidth) << std::setprecision(2) << scientific << actualSpeedup << endl;
        cout << "Percentage of ideal:  " << setw(charWidth) << std::setprecision(2) << fixed << percentage << "%" << endl;
        cout << endl;
      }
    }
    
//
//    using std::vector;
//    using std::map;
//    using std::pair;
//    using std::make_pair;
//    using std::tuple;
//    using std::cout;
//    using std::endl;
//    using std::setw;
//
//    using WorksetForAlgorithmChoice = map<AlgorithmChoice, int>;
//
//    vector< tuple<int,int,WorksetForAlgorithmChoice> > polyOrderMeshWidthWorksetTestCases;
//
//    const int meshWidth = 16;
//    vector<int> worksetSizes {1,2,4,8,16,32,64,128,256,512,1024,2048,4096};
//
//    // due to memory constraints, restrict the workset size for higher orders
//    map<int,int> maxWorksetSizeForPolyOrder;
//    maxWorksetSizeForPolyOrder[1] = 4096;
//    maxWorksetSizeForPolyOrder[2] = 4096;
//    maxWorksetSizeForPolyOrder[3] = 4096;
//    maxWorksetSizeForPolyOrder[4] = 4096;
//    maxWorksetSizeForPolyOrder[5] = 4096;
//    maxWorksetSizeForPolyOrder[6] = 2048;
//    maxWorksetSizeForPolyOrder[7] = 1024;
//    maxWorksetSizeForPolyOrder[8] = 512;
//
//    map<int,int> minWorksetSizeForPolyOrder;
//    minWorksetSizeForPolyOrder[1] = 1;
//    minWorksetSizeForPolyOrder[2] = 1;
//    minWorksetSizeForPolyOrder[3] = 1;
//    minWorksetSizeForPolyOrder[4] = 1;
//    minWorksetSizeForPolyOrder[5] = 1;
//    minWorksetSizeForPolyOrder[6] = 1;
//    minWorksetSizeForPolyOrder[7] = 1;
//    minWorksetSizeForPolyOrder[8] = 1;
//
//    switch (mode) {
//      case Calibration:
//      {
//        for (int polyOrder=polyOrderMin; polyOrder<=polyOrderMax; polyOrder++)
//        {
//          for (int worksetSize : worksetSizes)
//          {
//            if (worksetSize > maxWorksetSizeForPolyOrder[polyOrder])
//            {
//              continue;
//            }
//            if (worksetSize < minWorksetSizeForPolyOrder[polyOrder])
//            {
//              continue;
//            }
//            // use the same worksetSize for all AlgorithmChoice's
//            WorksetForAlgorithmChoice worksetForAlgorithmChoice;
//            for (auto algorithmChoice : algorithmChoices)
//            {
//              worksetForAlgorithmChoice[algorithmChoice] = worksetSize;
//            }
//            polyOrderMeshWidthWorksetTestCases.push_back(tuple<int,int,WorksetForAlgorithmChoice>{polyOrder,meshWidth,worksetForAlgorithmChoice} );
//          }
//        }
//      }
//      break;
//      case Test:
//      {
//        // for test run, use the same modestly-sized tuples for each AlgorithmChoice
//        // (note that meshWidth varies here)
//        vector< tuple<int,int,int> > testCases { tuple<int,int,int> {1,8,512},
//                                                 tuple<int,int,int> {2,8,256},
//                                                 tuple<int,int,int> {3,4,64},
//                                                 tuple<int,int,int> {3,4,9} // test case whose workset size does not evenly divide the cell count
//        };
//
//        for (auto testCase : testCases )
//        {
//          int polyOrder   = std::get<0>(testCase);
//          int meshWidth   = std::get<1>(testCase);
//
//          int numCells = 1;
//          for (int d=0; d<spaceDim; d++)
//          {
//            numCells *= meshWidth;
//          }
//
//          WorksetForAlgorithmChoice worksetForAlgorithmChoice;
//          for (auto algorithmChoice : algorithmChoices)
//          {
//            worksetForAlgorithmChoice[algorithmChoice] = std::get<2>(testCase);
//          }
//          worksetForAlgorithmChoice[Uniform] = numCells;
//          polyOrderMeshWidthWorksetTestCases.push_back(tuple<int,int,WorksetForAlgorithmChoice>{polyOrder,meshWidth,worksetForAlgorithmChoice} );
//        }
//      }
//      break;
//      case BestSerial:
//      {
//        // manually calibrated workset sizes on iMac Pro (2.3 GHz Xeon W, 18-core, running in serial)
//        // (these were calibrated without much tuning for the affine tensor case; if/when that happens, will want to recalibrate.)
//
//        map<int,int> standardWorksetForPolyOrder;
//        standardWorksetForPolyOrder[1] = 64;
//        standardWorksetForPolyOrder[2] = 64;
//        standardWorksetForPolyOrder[3] = 128;
//        standardWorksetForPolyOrder[4] = 64;
//        standardWorksetForPolyOrder[5] = 16;
//        standardWorksetForPolyOrder[6] = 2;
//        standardWorksetForPolyOrder[7] = 2;
//        standardWorksetForPolyOrder[8] = 1;
//
//        // Non-Affine Tensor
//        // it turns out 4096 is the best choice for the PointValueCache algorithm for any polyOrder from 1 to 5
//        // this likely means we're not exposing enough parallelism within the cell.
//        map<int,int> nonAffineTensorWorksetForPolyOrder;
//        nonAffineTensorWorksetForPolyOrder[1] = 256;
//        nonAffineTensorWorksetForPolyOrder[2] = 256;
//        nonAffineTensorWorksetForPolyOrder[3] = 128;
//        nonAffineTensorWorksetForPolyOrder[4] = 64;
//        nonAffineTensorWorksetForPolyOrder[5] = 16;
//        nonAffineTensorWorksetForPolyOrder[6] = 8;
//        nonAffineTensorWorksetForPolyOrder[7] = 2;
//        nonAffineTensorWorksetForPolyOrder[8] = 2;
//
//        map<int,int> affineTensorWorksetForPolyOrder;
//        affineTensorWorksetForPolyOrder[1] = 256;
//        affineTensorWorksetForPolyOrder[2] = 128;
//        affineTensorWorksetForPolyOrder[3] = 16;
//        affineTensorWorksetForPolyOrder[4] = 4;
//        affineTensorWorksetForPolyOrder[5] = 2;
//        affineTensorWorksetForPolyOrder[6] = 1;
//        affineTensorWorksetForPolyOrder[7] = 1;
//        affineTensorWorksetForPolyOrder[8] = 1;
//
//        // for the cases that we have not tried yet, we try to choose sensible guesses for workset size:
//        // 1 is best, we think, for polyOrder 8, so it'll be the best for the rest.
//        int worksetSize = 1;
//        for (int polyOrder=9; polyOrder <= polyOrderMax; polyOrder++)
//        {
//          nonAffineTensorWorksetForPolyOrder[polyOrder] = worksetSize;
//          affineTensorWorksetForPolyOrder[polyOrder]    = worksetSize;
//          standardWorksetForPolyOrder[polyOrder]        = worksetSize;
//        }
//
//        int numCells = 1;
//        for (int d=0; d<spaceDim; d++)
//        {
//          numCells *= meshWidth;
//        }
//
//        for (int polyOrder=polyOrderMin; polyOrder<=polyOrderMax; polyOrder++)
//        {
//          WorksetForAlgorithmChoice worksetForAlgorithmChoice;
//          worksetForAlgorithmChoice[Standard]        = standardWorksetForPolyOrder       [polyOrder];
//          worksetForAlgorithmChoice[NonAffineTensor] = nonAffineTensorWorksetForPolyOrder[polyOrder];
//          worksetForAlgorithmChoice[AffineTensor]    = nonAffineTensorWorksetForPolyOrder[polyOrder];
//          worksetForAlgorithmChoice[Uniform]         = numCells;
//
//          polyOrderMeshWidthWorksetTestCases.push_back(tuple<int,int,WorksetForAlgorithmChoice>{polyOrder,meshWidth,worksetForAlgorithmChoice} );
//        }
//      }
//        break;
//      case BestCuda:
//      {
//        {
//          // STANDARD
//          // manually calibrated workset size on P100 (ride) - best for Standard
//          // these are for 4096-element meshes
//          map<int,int> standardWorksetForPolyOrder;
//          standardWorksetForPolyOrder[1] = 4096;
//          standardWorksetForPolyOrder[2] = 1024;
//          standardWorksetForPolyOrder[3] = 128;
//          standardWorksetForPolyOrder[4] = 64;
//          standardWorksetForPolyOrder[5] = 8;
//          standardWorksetForPolyOrder[6] = 4;
//          standardWorksetForPolyOrder[7] = 2;
//          standardWorksetForPolyOrder[8] = 1;
//
//          // Non-Affine Tensor
//          // it turns out 4096 is the best choice for the PointValueCache algorithm for any polyOrder from 1 to 5
//          // this likely means we're not exposing enough parallelism within the cell.
//          map<int,int> nonAffineTensorWorksetForPolyOrder;
//          nonAffineTensorWorksetForPolyOrder[1] = 4096;
//          nonAffineTensorWorksetForPolyOrder[2] = 4096;
//          nonAffineTensorWorksetForPolyOrder[3] = 4096;
//          nonAffineTensorWorksetForPolyOrder[4] = 4096;
//          nonAffineTensorWorksetForPolyOrder[5] = 4096;
//          nonAffineTensorWorksetForPolyOrder[6] = 2048;
//          nonAffineTensorWorksetForPolyOrder[7] = 512;
//          nonAffineTensorWorksetForPolyOrder[8] = 512;
//
//          // for the cases that we have not tried yet, we try to choose sensible guesses for workset size:
//          int nonAffineWorksetSize = 256; // divide by 2 for each polyOrder beyond 8
//          int standardWorksetSize  = 1;  // 1 is best, we think, for polyOrder 8, so it'll be the best for the rest.
//          for (int polyOrder=9; polyOrder <= polyOrderMax; polyOrder++)
//          {
//            nonAffineTensorWorksetForPolyOrder[polyOrder] = nonAffineWorksetSize;
//            nonAffineWorksetSize = (nonAffineWorksetSize > 1) ? nonAffineWorksetSize / 2 : 1;
//            standardWorksetForPolyOrder[polyOrder] = standardWorksetSize;
//          }
//
//          int numCells = 1;
//          for (int d=0; d<spaceDim; d++)
//          {
//            numCells *= meshWidth;
//          }
//
//          for (int polyOrder=polyOrderMin; polyOrder<=polyOrderMax; polyOrder++)
//          {
//            WorksetForAlgorithmChoice worksetForAlgorithmChoice;
//            worksetForAlgorithmChoice[Standard]        = standardWorksetForPolyOrder       [polyOrder];
//            worksetForAlgorithmChoice[NonAffineTensor] = nonAffineTensorWorksetForPolyOrder[polyOrder];
//            worksetForAlgorithmChoice[AffineTensor]    = nonAffineTensorWorksetForPolyOrder[polyOrder];
//            worksetForAlgorithmChoice[Uniform]         = numCells;
//
//            polyOrderMeshWidthWorksetTestCases.push_back(tuple<int,int,WorksetForAlgorithmChoice>{polyOrder,meshWidth,worksetForAlgorithmChoice} );
//          }
//        }
//        break;
//
//      default:
//        break;
//    }
//    }
//
//    cout << std::setprecision(2) << std::scientific;
//
//    map< AlgorithmChoice, map<int, pair<double,int> > > maxAlgorithmThroughputForPolyOrder; // values are (throughput in GFlops/sec, worksetSize)
//
//    const int charWidth = 15;
//
//    for (auto & testCase : polyOrderMeshWidthWorksetTestCases)
//    {
//      int polyOrder       = std::get<0>(testCase);
//      int meshWidth       = std::get<1>(testCase);
//      auto worksetSizeMap = std::get<2>(testCase);
//      std::cout << "\n\n";
//      std::cout << "Running with polyOrder = " << polyOrder << ", meshWidth = " << meshWidth << std::endl;
//      for (auto algorithmChoice : algorithmChoices)
//      {
//        int worksetSize = worksetSizeMap[algorithmChoice];
//        auto geometry = getMesh<Scalar, spaceDim, ExecutionSpace>(algorithmChoice, meshWidth);
//
//        // timers recorded in performStructuredQuadratureHypercubeGRADGRAD, performStandardQuadratureHypercubeGRADGRAD
//        auto jacobianAndCellMeasureTimer = Teuchos::TimeMonitor::getNewTimer("Jacobians");
//        auto fstIntegrateCall = Teuchos::TimeMonitor::getNewTimer("transform + integrate()");
//        auto initialSetupTimer = Teuchos::TimeMonitor::getNewTimer("Initial Setup");
//
//        jacobianAndCellMeasureTimer->reset();
//        fstIntegrateCall->reset();
//        initialSetupTimer->reset();
//
//        double elapsedTimeSeconds = 0;
//        double jacobianCellMeasureFlopCount = 0;
//        double transformIntegrateFlopCount = 0;
//
//        if (algorithmChoice == Standard)
//        {
//          // each cell needs on the order of polyOrder^N quadrature points, each of which has a Jacobian of size N * N.
//          auto timer = Teuchos::TimeMonitor::getNewTimer("Standard Integration");
//          timer->start();
//          performStandardQuadratureHypercubeGRADGRAD<Scalar,Scalar,spaceDim,ExecutionSpace>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
//          timer->stop();
//          elapsedTimeSeconds = timer->totalElapsedTime();
//
//          cout << "Standard, workset size:          " << setw(charWidth) << worksetSize << endl;
//
//          timer->reset();
//        }
//        else if (algorithmChoice == AffineTensor)
//        {
//          auto timer = Teuchos::TimeMonitor::getNewTimer("Affine tensor Integration");
//          timer->start();
//          performStructuredQuadratureHypercubeGRADGRAD<Scalar>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
//          timer->stop();
//
//          elapsedTimeSeconds = timer->totalElapsedTime();
//
//          cout << "Affine Tensor, workset size:     " << setw(charWidth) << worksetSize << endl;
//
//          timer->reset();
//        }
//        else if (algorithmChoice == NonAffineTensor)
//        {
//          auto timer = Teuchos::TimeMonitor::getNewTimer("Non-affine tensor Integration");
//          timer->start();
//          performStructuredQuadratureHypercubeGRADGRAD<Scalar>(geometry, polyOrder, worksetSize, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
//          timer->stop();
//
//          elapsedTimeSeconds = timer->totalElapsedTime();
//
//          cout << "Non-Affine Tensor, workset size: " << setw(charWidth) << worksetSize << endl;
//
//          timer->reset();
//        }
//        else if (algorithmChoice == Uniform)
//        {
//          // for uniform, override worksetSize: no loss in taking maximal worksetSize
//          int numCells = 1;
//          for (int d=0; d<spaceDim; d++)
//          {
//            numCells *= meshWidth;
//          }
//          auto timer = Teuchos::TimeMonitor::getNewTimer("Uniform Integration");
//          timer->start();
//          performStructuredQuadratureHypercubeGRADGRAD<Scalar>(geometry, polyOrder, numCells, transformIntegrateFlopCount, jacobianCellMeasureFlopCount);
//          timer->stop();
//
//          elapsedTimeSeconds = timer->totalElapsedTime();
//
//          cout << "Uniform, workset size:           " << setw(charWidth) << worksetSize << endl;
//
//          timer->reset();
//        }
//
//        const double approximateFlopCountTotal = transformIntegrateFlopCount + jacobianCellMeasureFlopCount;
//        const double overallThroughputInGFlops = approximateFlopCountTotal / elapsedTimeSeconds / 1.0e9;
//
//        const double previousMaxThroughput = maxAlgorithmThroughputForPolyOrder[algorithmChoice][polyOrder].first;
//        if (overallThroughputInGFlops > previousMaxThroughput)
//        {
//          maxAlgorithmThroughputForPolyOrder[algorithmChoice][polyOrder] = make_pair(overallThroughputInGFlops,worksetSize);
//        }
//
//        // timing details
//        double integrateCallTime       = fstIntegrateCall->totalElapsedTime();
//        double integrateCallPercentage = integrateCallTime / elapsedTimeSeconds * 100.0;
//        double jacobianTime            = jacobianAndCellMeasureTimer->totalElapsedTime();
//        double jacobianPercentage      = jacobianTime / elapsedTimeSeconds * 100.0;
//        double initialSetupTime        = initialSetupTimer->totalElapsedTime();
//        double initialSetupPercentage  = initialSetupTime / elapsedTimeSeconds * 100.0;
//        double remainingTime           = elapsedTimeSeconds - (integrateCallTime + jacobianTime + initialSetupTime);
//        double remainingPercentage     = remainingTime / elapsedTimeSeconds * 100.0;
//
//        const double transformIntegrateThroughputInGFlops = transformIntegrateFlopCount  / integrateCallTime / 1.0e9;
//        const double jacobiansThroughputInGFlops          = jacobianCellMeasureFlopCount / jacobianTime      / 1.0e9;
//        cout << "Time (core integration)          " << setw(charWidth) << std::scientific << integrateCallTime << " seconds (" << std::fixed << integrateCallPercentage << "%)." << endl;
//        cout << "flop estimate (core):            " << setw(charWidth) << std::scientific << transformIntegrateFlopCount << endl;
//        cout << "estimated throughput (core):     " << setw(charWidth) << std::scientific << transformIntegrateThroughputInGFlops << " GFlops" << endl;
//        cout << std::fixed;
//        cout << "Time (Jacobians)                 " << setw(charWidth) << std::scientific << jacobianTime      << " seconds (" << std::fixed << jacobianPercentage      << "%)." << endl;
//        cout << "flop estimate (Jacobians):       " << setw(charWidth) << std::scientific << jacobianCellMeasureFlopCount << endl;
//        cout << "estimated throughput (Jac.):     " << setw(charWidth) << std::scientific << jacobiansThroughputInGFlops << " GFlops" << endl;
//        cout << "Time (initial setup)             " << setw(charWidth) << std::scientific << initialSetupTime  << " seconds (" << std::fixed << initialSetupPercentage  << "%)." << endl;
//        cout << "Time (other)                     " << setw(charWidth) << std::scientific << remainingTime     << " seconds (" << std::fixed << remainingPercentage     << "%)." << endl;
//        cout << "Time (total):                    " << setw(charWidth) << std::scientific << elapsedTimeSeconds   << " seconds.\n";
//        cout << "flop estimate (total):           " << setw(charWidth) << std::scientific << approximateFlopCountTotal << endl;
//        cout << "estimated throughput (total):    " << setw(charWidth) << std::scientific << overallThroughputInGFlops << " GFlops" << endl;
//
//        cout << endl;
//      }
//    }
//
//    if (mode == Calibration)
//    {
//      for (auto & algorithmChoice : algorithmChoices)
//      {
//        cout << "Best workset sizes for " << to_string(algorithmChoice) << ":" << endl;
//
//        for (auto & maxThroughputEntry : maxAlgorithmThroughputForPolyOrder[algorithmChoice])
//        {
//          int polyOrder   = maxThroughputEntry.first;
//          int worksetSize = maxThroughputEntry.second.second;
//          double throughput = maxThroughputEntry.second.first;
//          cout << "p = " << polyOrder << ":" << setw(5) << worksetSize << " (" << throughput << " GFlops/sec)\n";
//        }
//      }
//    }
  }
  
  return 0;
}
