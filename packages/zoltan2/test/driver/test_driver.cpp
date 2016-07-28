// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/* \file test_driver.cpp
 * \brief Test driver for Zoltan2. Facilitates generation of test problem via
 * a simple .xml input interface
 */

// taking headers from existing driver template
// will keep or remove as needed
#include <UserInputForTests.hpp>
#include <Zoltan2_Typedefs.hpp>
#include <AdapterForTests.hpp>
#include <Zoltan2_ComparisonHelper.hpp>
#include <Zoltan2_MetricAnalyzer.hpp>

#include <Zoltan2_ProblemFactory.hpp>
#include <Zoltan2_EvaluatePartitionFactory.hpp>
#include <Zoltan2_BasicIdentifierAdapter.hpp>
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>

#include <Zoltan2_Parameters.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Teuchos_FileInputSource.hpp>

#include <sstream>
#include <string>
#include <map>
#include <iostream>
#include <queue>

//#include <BDD_PamgenUtils.hpp>

using Teuchos::ParameterList;
using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::XMLObject;
using namespace Zoltan2_TestingFramework;

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::map;
using std::pair;
using std::exception;
using std::ostringstream;
using std::queue;

#define ERRMSG(msg) if (rank == 0){ cerr << "FAIL: " << msg << endl; }
#define EXC_ERRMSG(msg, e) \
if (rank==0){ cerr << "FAIL: " << msg << endl << e.what() << endl;}

void xmlToModelPList(const Teuchos::XMLObject &xml, Teuchos::ParameterList & plist)
{
  // This method composes a plist for the problem
  Teuchos::XMLParameterListReader reader;
  plist = reader.toParameterList(xml);
  
  //  Get list of valid Zoltan2 Parameters
  // Zoltan 2 parameters appear in the input file
  // Right now we have default values stored in
  // the parameter list, we would like to apply
  // the options specified by the user in their
  // input file
  Teuchos::ParameterList zoltan2Parameters;
  Zoltan2::createAllParameters(zoltan2Parameters);
  
  if (plist.isSublist("Zoltan2Parameters")) {
    // Apply user specified zoltan2Parameters
    ParameterList &sub = plist.sublist("Zoltan2Parameters");
    zoltan2Parameters.setParameters(sub);
  }
}

bool getParameterLists(const string &inputFileName,
                       queue<ParameterList> &problems,
                       queue<ParameterList> &comparisons,
                       const RCP<const Teuchos::Comm<int> > & comm)
{
  int rank = comm->getRank();
  // return a parameter list of problem definitions
  // and a parameter list for solution comparisons
  Teuchos::FileInputSource inputSource(inputFileName);
  if(rank == 0) {
    cout << "Input file source: " << inputFileName << endl;
  }
  XMLObject xmlInput;
  
  // Try to get xmlObject from inputfile
  try
  {
    xmlInput = inputSource.getObject();
  }
  catch(exception &e)
  {
    EXC_ERRMSG("Test Driver error: reading", e); // error reading input
    return false;
  }

  // get the parameter lists for each model
  for(int i = 0; i < xmlInput.numChildren(); i++)
  {
    ParameterList plist;
    xmlToModelPList(xmlInput.getChild(i), plist);

    if(plist.name() == "Comparison") {
      comparisons.emplace(plist);
    }
    else {
      problems.emplace(plist);
    }
  }

  return true;
}

bool run(const UserInputForTests &uinput,
        const ParameterList &problem_parameters,
        bool bHasComparisons,
        RCP<ComparisonHelper> & comparison_helper,
        const RCP<const Teuchos::Comm<int> > & comm)
{
  // Major steps in running a problem in zoltan 2
  // 1. get an input adapter
  // 2. construct the problem
  // 3. solve the problem
  // 4. analyze metrics
  // 5. clean up

  int rank = comm->getRank();
  if(!problem_parameters.isParameter("kind"))
  {
    if(rank == 0) {
      std::cout << "Problem kind not provided" << std::endl;
    }
    return false;
  }
  if(!problem_parameters.isParameter("InputAdapterParameters"))
  {
    if(rank == 0) {
      std::cout << "Input adapter parameters not provided" << std::endl;
    }
    return false;
  }
  if(!problem_parameters.isParameter("Zoltan2Parameters"))
  {
    if(rank == 0) {
      std::cout << "Zoltan2 problem parameters not provided" << std::endl;
    }
    return false;
  }

  if(rank == 0) {
    cout << "\n\nRunning test: " << problem_parameters.name() << endl;
  }
  
  ////////////////////////////////////////////////////////////
  // 0. add comparison source
  ////////////////////////////////////////////////////////////
  ComparisonSource * comparison_source = new ComparisonSource;
  comparison_helper->AddSource(problem_parameters.name(), comparison_source);
  comparison_source->addTimer("adapter construction time");
  comparison_source->addTimer("problem construction time");
  comparison_source->addTimer("solve time");

  ////////////////////////////////////////////////////////////
  // 1. get basic input adapter
  ////////////////////////////////////////////////////////////
  const ParameterList &adapterPlist =
                      problem_parameters.sublist("InputAdapterParameters");
  comparison_source->timers["adapter construction time"]->start();

  // a pointer to a basic type
  AdapterWithOptionalCoordinateAdapter adapters = AdapterForTests::getAdapterForInput(
                                        const_cast<UserInputForTests*>(&uinput),
                                        adapterPlist,comm); 
  comparison_source->timers["adapter construction time"]->stop();

  if(adapters.mainAdapter == nullptr)
  {
    if(rank == 0) {
      cout << "Get adapter for input failed" << endl;
    }
    return false;
  }
  RCP<basic_id_t> iaRCP = rcp(reinterpret_cast<basic_id_t *>
    (adapters.mainAdapter), true);

  RCP<Zoltan2::VectorAdapter<tMVector_t>> coordinateAdapterRCP = 
    rcp(adapters.coordinateAdapter, true);

  ////////////////////////////////////////////////////////////
  // 2. construct a Zoltan2 problem
  ////////////////////////////////////////////////////////////
  // If we are here we have an input adapter, no need to check for one.
  string adapter_name = adapterPlist.get<string>("input adapter"); 
  // get Zoltan2 partition parameters
  ParameterList zoltan2_parameters = 
   const_cast<ParameterList &>(problem_parameters.sublist("Zoltan2Parameters"));
  
  if(rank == 0) {
    cout << endl;
  }

  comparison_source->timers["problem construction time"]->start();
  std::string problem_kind = problem_parameters.get<std::string>("kind"); 
  if (rank == 0) {
    std::cout << "Creating a new " << problem_kind << " problem." << std::endl;
  }
#ifdef HAVE_ZOLTAN2_MPI
  base_problem_t * problem = 
    Zoltan2_TestingFramework::ProblemFactory::newProblem(problem_kind,
                                                         adapter_name,
                                                         adapters.mainAdapter,
                                                         &zoltan2_parameters,
                                                         MPI_COMM_WORLD);
#else
  base_problem_t * problem = 
    Zoltan2_TestingFramework::ProblemFactory::newProblem(problem_kind,
                                                         adapter_name,
                                                         adapters.mainAdapter,
                                                         &zoltan2_parameters);
#endif

  if (problem == nullptr) {
    if (rank == 0) {
      std::cerr << "Input adapter type: " << adapter_name 
                << ", is unavailable, or misspelled." << std::endl;
    }
    return false;
  }
  else if(rank == 0) {
    std::cout << "Using input adapter type: " + adapter_name << std::endl;
  }
  RCP<base_problem_t> problemRCP = rcp(problem, true);

  ////////////////////////////////////////////////////////////
  // 3. Solve the problem
  ////////////////////////////////////////////////////////////
  comparison_source->timers["solve time"]->start();
  if (problem_kind == "partitioning") {
    reinterpret_cast<partitioning_problem_t *>(problem)->solve();
  } else if (problem_kind == "ordering") {
    reinterpret_cast<ordering_problem_t *>(problem)->solve();
  } else if (problem_kind == "coloring") {
    reinterpret_cast<coloring_problem_t *>(problem)->solve();
  }

  comparison_source->timers["solve time"]->stop();
  if (rank == 0) {
    cout << problem_kind + " problem solved." << endl;
  }
 
#define KDDKDD
#ifdef KDDKDD
  {
  const base_adapter_t::gno_t *kddIDs = NULL;
  adapters.mainAdapter->getIDsView(kddIDs);
    for (size_t i = 0; i < adapters.mainAdapter->getLocalNumIDs(); i++) {
      std::cout << rank << " LID " << i
                << " GID " << kddIDs[i]
                << " PART " 
                << reinterpret_cast<partitioning_problem_t *>
                               (problem)->getSolution().getPartListView()[i]
                << std::endl;
    }
  }
  if (adapter_name == "XpetraCrsGraph") {
    typedef xcrsGraph_adapter::lno_t lno_t;
    typedef xcrsGraph_adapter::gno_t gno_t;
    typedef xcrsGraph_adapter::scalar_t scalar_t;
    int ewgtDim = 
        reinterpret_cast<const xcrsGraph_adapter *>(adapters.mainAdapter)->
          getNumWeightsPerEdge();
    lno_t localNumObj = 
        reinterpret_cast<const xcrsGraph_adapter *>(adapters.mainAdapter)->
          getLocalNumVertices();
    const gno_t *vertexIds;
    reinterpret_cast<const xcrsGraph_adapter *>(adapters.mainAdapter)->
      getVertexIDsView(vertexIds);
    const lno_t *offsets;
    const gno_t *adjIds;
    reinterpret_cast<const xcrsGraph_adapter *>(adapters.mainAdapter)->
      getEdgesView(offsets, adjIds);
    for (int edim = 0; edim < ewgtDim; edim++) {
      const scalar_t *weights;
      int stride=0;
      reinterpret_cast<xcrsGraph_adapter *>(adapters.mainAdapter)->
        getEdgeWeightsView(weights, stride, edim);
      for (lno_t i=0; i < localNumObj; i++)
        for (lno_t j=offsets[i]; j < offsets[i+1]; j++)
          std::cout << edim << " " << vertexIds[i] << " " 
                    << adjIds[j] << " " << weights[stride*j] << std::endl;
    }
  }
#endif

  ////////////////////////////////////////////////////////////
  // 4. Print problem metrics
  ////////////////////////////////////////////////////////////
  // An environment.  This is usually created by the problem.
  // BDD unused, only applicable to partitioning problems
  // RCP<const Zoltan2::Environment> env =
  //   reinterpret_cast<partitioning_problem_t *>(problem)->getEnvironment();

  // get metric object
  // this is not the most beautiful thing, but comparison parameters is checked 
  // as well because it's possible we are checking comparisons of metrics but 
  // not individual metrics
  // we want to only load the EvaluatePartition when Metrics is requested, or 
  // some comparison is requested

  bool bSuccess = true;

  if(problem_parameters.isSublist("Metrics") || bHasComparisons) { 
    // the specification is that we don't create anything unless 
    // the Metrics list exists
    RCP<EvaluatePartition<basic_id_t> > metricObject = rcp(
       Zoltan2_TestingFramework::EvaluatePartitionFactory::newEvaluatePartition(
               reinterpret_cast<partitioning_problem_t*> (problem), 
               adapter_name, adapters.mainAdapter, &zoltan2_parameters));

    std::ostringstream msgSummary;
    metricObject->printMetrics(msgSummary, true); //
    if(rank == 0) {
      cout << msgSummary.str();
    }

    std::ostringstream msgResults;
    if (!MetricAnalyzer::analyzeMetrics(metricObject, 
                                        problem_parameters.sublist("Metrics"), 
                                        msgResults)) 
    { 
     // Note the MetricAnalyzer only cares about the data found in the 
     // "Metrics" sublist
      bSuccess = false;
      if (rank == 0) {
	std::cout << "MetricAnalyzer::analyzeMetrics() "
                  << "returned false and the test is FAILED." << std::endl;
      }
    }
    if(rank == 0) {
      cout << msgResults.str();
    }

//#define BDD
#ifdef BDD 
    if (problem_kind == "ordering") {
      std::cout << "\nLet's examine the solution..." << std::endl;
      auto solution = reinterpret_cast<ordering_problem_t *>
                                       (problem)->getSolution();
      if (solution->haveSeparators() ) {
      
        std::ostringstream sol;
        sol << "Number of column blocks: " << solution->getNumSeparatorBlocks() 
            << std::endl;
        if (solution->getPermutationSize() < 100) {
          if (solution->havePerm()) {
            sol << "permutation: {";
            for (auto &x : solution->getPermutationRCPConst(false)) 
              sol << " " << x;
            sol << "}" << std::endl;
          }
       
         if (solution->haveInverse()) { 
            sol << "inverse permutation: {";
            for (auto &x : solution->getPermutationRCPConst(true)) 
              sol << " " << x;
            sol << "}" << std::endl;
         }
        
         if (solution->haveSeparatorRange()) {
            sol << "separator range: {";
            for (auto &x : solution->getSeparatorRangeRCPConst()) 
              sol << " " << x;
            sol << "}" << std::endl;
         }
         
          if (solution->haveSeparatorTree()) { 
            sol << "separator tree: {";
            for (auto &x : solution->getSeparatorTreeRCPConst()) 
              sol << " " << x;
            sol << "}" << std::endl;
          }
        }

        std::cout << sol.str() << std::endl;
      }
    }
#endif
    // 4b. timers
    //  if(zoltan2_parameters.isParameter("timer_output_stream"))
    //    reinterpret_cast<partitioning_problem_t *>(problem)->printTimers();

    ////////////////////////////////////////////////////////////
    // 5. Add solution to map for possible comparison testing
    ////////////////////////////////////////////////////////////

    comparison_source->adapter = iaRCP;
    comparison_source->coordinateAdapterRCP = coordinateAdapterRCP;
    comparison_source->problem = problemRCP;
    comparison_source->metricObject = metricObject;
    comparison_source->problem_kind = (problem_parameters.isParameter("kind") ? 
                                       problem_parameters.get<string>("kind") :
                                       "?");
    comparison_source->adapter_kind = adapter_name;

  
    // write mesh solution
    //  auto sol = reinterpret_cast<partitioning_problem_t *>(problem)->getSolution();
    //  MyUtils::writePartionSolution(sol.getPartListView(), ia->getLocalNumIDs(), comm);

    ////////////////////////////////////////////////////////////
    // 6. Clean up
    ////////////////////////////////////////////////////////////
  }

  return bSuccess;
}

bool mainExecute(int argc, char *argv[], RCP<const Comm<int> > &comm) 
{
  ////////////////////////////////////////////////////////////
  // (0) Set up MPI environment and timer
  ////////////////////////////////////////////////////////////
  int rank = comm->getRank(); // get rank

  ////////////////////////////////////////////////////////////
  // (1) Get and read the input file
  // the input file defines tests to be run
  ////////////////////////////////////////////////////////////
  string inputFileName(""); 
  if(argc > 1)
    inputFileName = argv[1]; // user has provided an input file
  else{
    if(rank == 0){
      std::cout << "\nFAILED to specify xml input file!" << std::endl;
      ostringstream msg;
      msg << "\nStandard use of test_driver.cpp:\n";
      msg << "mpiexec -n <procs> ./Zoltan2_test_driver.exe <input_file.xml>\n";
      std::cout << msg.str() << std::endl;
    }
    return false;
  }

  ////////////////////////////////////////////////////////////
  // (2) Get All Input Parameter Lists
  ////////////////////////////////////////////////////////////
  queue<ParameterList> problems, comparisons;
  if( !getParameterLists(inputFileName,problems, comparisons, comm) ) {
    return false;
  }
  
  ////////////////////////////////////////////////////////////
  // (3) Get Input Data Parameters
  ////////////////////////////////////////////////////////////

  // assumes that first block will always be the input block
  const ParameterList inputParameters = problems.front();
  if(inputParameters.name() != "InputParameters")
  {
    if(rank == 0)
      cout << "InputParameters not defined. Testing FAILED." << endl;
    return false;
  }
  
  // get the user input for all tests
  UserInputForTests uinput(inputParameters,comm);

  problems.pop();
  comm->barrier();

  bool bPass = true;
  if(uinput.hasInput())
  {
    ////////////////////////////////////////////////////////////
    // (4) Perform all tests
    ////////////////////////////////////////////////////////////
//     pamgen debugging
//    MyUtils::writeMesh(uinput,comm);
//    MyUtils::getConnectivityGraph(uinput, comm);
        
    RCP<ComparisonHelper> comparison_manager = rcp(new ComparisonHelper);
    while (!problems.empty()) {
      if (!run(uinput, problems.front(), !comparisons.empty(),
               comparison_manager, comm)) {
        std::cout << "Problem run returned false" << std::endl;
        bPass = false;
      }
      problems.pop();
    }

    ////////////////////////////////////////////////////////////
    // (5) Compare solutions
    ////////////////////////////////////////////////////////////

    while (!comparisons.empty()) {
      if (!comparison_manager->Compare(comparisons.front(),comm)) {
        if (rank == 0) {
          std::cout << "Comparison manager returned false so the "
                    << "test should fail." << std::endl;
        }
        bPass = false;
      }
      comparisons.pop();
    }
  }
  else {
    if(rank == 0) {
      cout << "\nFAILED to load input data source. Skipping all tests." << endl;
      return false;
    }
  }

  return bPass;
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv); 
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  int result = 0;
  int rank = comm->getRank();
  try {
    result = mainExecute(argc, argv, comm) ? 0 : 1; // code 0 is ok,
                                                    // 1 is a failed test
  }
  catch(std::logic_error &e) { 
    // logic_error exceptions can be thrown by EvaluatePartition or 
    // MetricAnalyzer if any problem is detected in the formatting of the 
    // input xml
    if (rank == 0) {
      std::cout << "Test driver for rank " << rank 
                << " caught the following exception: " << e.what() << std::endl;
    }
    result = 1; // fail for any exception
  }
  catch(std::runtime_error &e) { 
    std::cout << "Test driver for rank " << rank 
              << " caught the following exception: " << e.what() << std::endl;
    result = 1; // fail for any exception
  }
  catch(std::exception &e) { 
    std::cout << "Test driver for rank " << rank 
              << " caught the following exception: " << e.what() << std::endl;
    result = 1; // fail for any exception
  }

  // clean up - reduce the result codes
  comm->barrier();
  int resultReduced;

  // for a passed test all of these values should return 0 - 
  // if any result is 1 this will reduce to 1 and the test will fail
  reduceAll<int,int>(*comm, Teuchos::EReductionType::REDUCE_MAX, 1, 
                     &result, &resultReduced); 

  // provide a final message which guarantees that the test will fail 
  // if any of the processes failed
  if (rank == 0) {
    std::cout << "Test Driver with " << comm->getSize() 
              << " processes has reduced result code " << resultReduced
              << " and is exiting in the " 
              << ((resultReduced == 0 ) ? "PASSED" : "FAILED") << " state." 
              << std::endl;
  }

  return result;
}
