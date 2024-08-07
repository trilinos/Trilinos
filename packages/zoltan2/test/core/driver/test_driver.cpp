// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include <Zoltan2_EvaluateFactory.hpp>

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

using Teuchos::ParameterList;
using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::XMLObject;
using namespace Zoltan2_TestingFramework;

using std::string;
using std::map;
using std::pair;
using std::exception;
using std::ostringstream;
using std::queue;

#define ERRMSG(msg) if (rank == 0){ std::cerr << "FAIL: " << msg << std::endl; }
#define EXC_ERRMSG(msg, e) \
if (rank==0){ std::cerr << "FAIL: " << msg << std::endl << e.what() << std::endl;}

void xmlToModelPList(const Teuchos::XMLObject &xml,
  Teuchos::ParameterList & plist)
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
    std::cout << "Input file source: " << inputFileName << std::endl;
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

// Utility function for safe type conversion of adapter
bool analyzeMetrics(RCP<EvaluateFactory> evaluateFactory,
                    std::ostringstream & msg,
                    const ParameterList &problem_parameters) {
  #define ANALYZE_METRICS(adapterClass, metricAnalyzerClass)               \
    RCP<EvaluateBaseClass<adapterClass>> pCast =                           \
      rcp_dynamic_cast<EvaluateBaseClass<adapterClass>>(                   \
         evaluateFactory->getEvaluateClass());                             \
    if(pCast == Teuchos::null) throw std::logic_error(                     \
      "Bad evaluate class cast in analyzeMetrics!"  );                     \
    metricAnalyzerClass analyzer(pCast);                                   \
    return analyzer.analyzeMetrics(                                        \
      problem_parameters.sublist("Metrics"), msg);

  #define ANALYZE_METRICS_PARTITIONING(adapterClass)                       \
    ANALYZE_METRICS(adapterClass,                                          \
      MetricAnalyzerEvaluatePartition<adapterClass>)

  #define ANALYZE_METRICS_ORDERING(adapterClass)                           \
    ANALYZE_METRICS(adapterClass,                                          \
      MetricAnalyzerEvaluateOrdering<adapterClass>)

  if(evaluateFactory->getProblemName() == "partitioning") {
    Z2_TEST_UPCAST(evaluateFactory->getAdapterType(), ANALYZE_METRICS_PARTITIONING)
  }
  else if(evaluateFactory->getProblemName() == "ordering") {
    Z2_TEST_UPCAST(evaluateFactory->getAdapterType(), ANALYZE_METRICS_ORDERING)
  }
  else {
    throw std::logic_error(
      "analyzeMetrics not implemented for this problem type!"  );
  }
}

// Utility function for safe type conversion of adapter
LocalOrderingSolution<zlno_t> * getLocalOrderingSolution(
  RCP<ProblemFactory> problemFactory) {
  #define GET_LOCAL_ORDERING(adapterClass)                                 \
      return (rcp_dynamic_cast<OrderingProblem<adapterClass>>(             \
        problemFactory->getProblem()))->getLocalOrderingSolution();
  Z2_TEST_UPCAST(problemFactory->getAdapterType(), GET_LOCAL_ORDERING)
}

// Utility function for safe type conversion of adapter
const zpart_t * getPartListView(RCP<ProblemFactory> problemFactory) {
  #define GET_PROBLEM_PARTS(adapterClass)                                  \
      return (rcp_dynamic_cast<PartitioningProblem<adapterClass>>(         \
        problemFactory->getProblem()))->getSolution().getPartListView();
  Z2_TEST_UPCAST(problemFactory->getAdapterType(), GET_PROBLEM_PARTS)
}

// Utility function for safe type conversion of adapter
void getIDsView(RCP<AdapterFactory> adapterFactory, const zgno_t *&Ids) {
    #define GET_IDS_VIEW(adapterClass)                                       \
        return dynamic_cast<adapterClass*>(                                  \
          adapterFactory->getMainAdapter())->getIDsView(Ids);
    Z2_TEST_UPCAST(adapterFactory->getMainAdapterType(), GET_IDS_VIEW);
    throw std::logic_error( "getIDsView() failed to match adapter name" );
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
    std::cout << "\n\nRunning test: " << problem_parameters.name() << std::endl;
  }

  ////////////////////////////////////////////////////////////
  // 0. add comparison source
  ////////////////////////////////////////////////////////////
  RCP<ComparisonSource> comparison_source = rcp(new ComparisonSource);

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
  RCP<AdapterFactory> adapterFactory = rcp(new AdapterFactory(
    const_cast<UserInputForTests*>(&uinput), adapterPlist, comm));

  comparison_source->timers["adapter construction time"]->stop();

  comparison_source->adapterFactory = adapterFactory; // saves until done

  ////////////////////////////////////////////////////////////
  // 2. construct a Zoltan2 problem
  ////////////////////////////////////////////////////////////
  // If we are here we have an input adapter, no need to check for one.
  string adapter_name = adapterPlist.get<string>("input adapter");
  // get Zoltan2 Parameters
  ParameterList zoltan2_parameters =
   const_cast<ParameterList &>(problem_parameters.sublist("Zoltan2Parameters"));
  if(rank == 0) {
    std::cout << std::endl;
  }

  comparison_source->timers["problem construction time"]->start();
  std::string problem_kind = problem_parameters.get<std::string>("kind");
  if (rank == 0) {
    std::cout << "Creating a new " << problem_kind << " problem." << std::endl;
  }

  RCP<ProblemFactory> problemFactory = rcp(new ProblemFactory(
                                      problem_kind,
                                      adapterFactory,
                                      &zoltan2_parameters
                                    #ifdef HAVE_ZOLTAN2_MPI
                                      ,MPI_COMM_WORLD
                                    #endif
                                      ));

  if(rank == 0) {
    std::cout << "Using input adapter type: " + adapter_name << std::endl;
  }

  comparison_source->problemFactory = problemFactory; // saves until done

  ////////////////////////////////////////////////////////////
  // 3. Solve the problem
  ////////////////////////////////////////////////////////////
  comparison_source->timers["solve time"]->start();

  problemFactory->getProblem()->solve();

  comparison_source->timers["solve time"]->stop();
  if (rank == 0) {
    std::cout << problem_kind + " problem solved." << std::endl;
  }

#undef KDDKDD
#ifdef KDDKDD
  if(problem_kind == "partitioning") {
    const base_adapter_t::gno_t *kddIDs = NULL;
    getIDsView(adapterFactory, kddIDs);
    for (size_t i = 0;
      i < adapterFactory->getMainAdapter()->getLocalNumIDs(); i++) {
      std::cout << rank << " LID " << i
                << " GID " << kddIDs[i]
                << " PART "
                << getPartListView(problemFactory)[i]
                << std::endl;
    }
  }
  if (adapter_name == "XpetraCrsGraph") {
    typedef xCG_xCG_t::lno_t lno_t;
    typedef xCG_xCG_t::gno_t gno_t;
    typedef xCG_xCG_t::scalar_t scalar_t;
    const xCG_xCG_t * xscrsGraphAdapter =
      dynamic_cast<const xCG_xCG_t *>(adapterFactory->getMainAdapter());

    int ewgtDim = xscrsGraphAdapter->getNumWeightsPerEdge();
    lno_t localNumObj = xscrsGraphAdapter->getLocalNumVertices();
    const gno_t *vertexIds;
    xscrsGraphAdapter->getVertexIDsView(vertexIds);
    const offset_t *offsets;
    const gno_t *adjIds;
    xscrsGraphAdapter->getEdgesView(offsets, adjIds);
    for (int edim = 0; edim < ewgtDim; edim++) {
      const scalar_t *weights;
      int stride=0;
      xscrsGraphAdapter->getEdgeWeightsView(weights, stride, edim);
      for (lno_t i=0; i < localNumObj; i++)
        for (offset_t j=offsets[i]; j < offsets[i+1]; j++)
          std::cout << edim << " " << vertexIds[i] << " "
                    << adjIds[j] << " " << weights[stride*j] << std::endl;
    }
  }
#endif

  ////////////////////////////////////////////////////////////
  // 4. Print problem metrics
  ////////////////////////////////////////////////////////////
  bool bSuccess = true;
  if(problem_parameters.isSublist("Metrics") || bHasComparisons) {
    RCP<EvaluateFactory> evaluateFactory = rcp(new EvaluateFactory(
                                        problem_kind,
                                        adapterFactory,
                                        &zoltan2_parameters,
                                        problemFactory));

    if(rank == 0) {
      std::cout << "Create evaluate class for: " + problem_kind << std::endl;
    }

    // must add for proper deletion
    comparison_source->evaluateFactory = evaluateFactory; // saves until done

    std::ostringstream msgSummary;

    evaluateFactory->getEvaluateClass()->printMetrics(msgSummary);
    if(rank == 0) {
      std::cout << msgSummary.str();
    }

    std::ostringstream msgResults;
    if (!analyzeMetrics(evaluateFactory, msgResults, problem_parameters)) {
      bSuccess = false;
      std::cout << "MetricAnalyzer::analyzeMetrics() "
                << "returned false and the test is FAILED." << std::endl;
    }
    if(rank == 0) {
      std::cout << msgResults.str();
    }

//#define BDD
#ifdef BDD
    if (problem_kind == "ordering") {
      std::cout << "\nLet's examine the solution..." << std::endl;
      LocalOrderingSolution<zlno_t> * localOrderingSolution =
         getLocalOrderingSolution(problemFactory);
      if (localOrderingSolution->haveSeparators() ) {
        std::cout << "Number of column blocks: "
          << localOrderingSolution->getNumSeparatorBlocks() << std::endl;
        {
          if (localOrderingSolution->havePerm()) {
            std::cout << "permutation: {";
            for (auto &x : localOrderingSolution->getPermutationRCPConst(false))
              std::cout << " " << x;
            std::cout << "}" << std::endl;
          }

          if (localOrderingSolution->haveInverse()) {
            std::cout << "inverse permutation: {";
            for (auto &x : localOrderingSolution->getPermutationRCPConst(true))
              std::cout << " " << x;
            std::cout << "}" << std::endl;
          }

          if (localOrderingSolution->haveSeparatorRange()) {
            std::cout << "separator range: {";
            for (auto &x : localOrderingSolution->getSeparatorRangeRCPConst())
              std::cout << " " << x;
            std::cout << "}" << std::endl;
          }

          if (localOrderingSolution->haveSeparatorTree()) {
            std::cout << "separator tree: {";
            for (auto &x : localOrderingSolution->getSeparatorTreeRCPConst())
              std::cout << " " << x;
            std::cout << "}" << std::endl;
          }
        }
      }
    }
#endif

    comparison_source->printTimers();

    // write mesh solution
    // if(problem_kind == "partitioning") {
    //  auto sol = reinterpret_cast<partitioning_problem_t *>(problem)->getSolution();
    //  MyUtils::writePartionSolution(sol.getPartListView(), ia->getLocalNumIDs(), comm);
    // }
  }

  return bSuccess;
}

bool mainExecute(int narg, char *arg[], RCP<const Comm<int> > &comm)
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
  if(narg > 1)
    inputFileName = arg[1]; // user has provided an input file
  else{
    if(rank == 0){
      std::cout << "\nFAILED to specify xml input file!" << std::endl;
      std::ostringstream msg;
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
  if( !getParameterLists(inputFileName, problems, comparisons, comm) ) {
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
      std::cout << "InputParameters not defined. Testing FAILED." << std::endl;
    return false;
  }

  // get the user input for all tests
  // UserInputForTests uinput(inputParameters,comm);

  problems.pop();
  comm->barrier();

  bool bPass = true;
  if(true)
  {
    ////////////////////////////////////////////////////////////
    // (4) Perform all tests
    ////////////////////////////////////////////////////////////
//     pamgen debugging
//    MyUtils::writeMesh(uinput,comm);
//    MyUtils::getConnectivityGraph(uinput, comm);

    RCP<ComparisonHelper> comparison_manager = rcp(new ComparisonHelper);
    while (!problems.empty()) {
      UserInputForTests uinput(inputParameters,comm);

      if(uinput.hasInput()) {
      if (!run(uinput, problems.front(), !comparisons.empty(),
               comparison_manager, comm)) {
        std::cout << "Problem run returned false" << std::endl;
        bPass = false;
      }
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
      std::cout << "\nFAILED to load input data source. Skipping "
        "all tests." << std::endl;
      return false;
    }
  }

  return bPass;
}

int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int result = 0;
  int rank = comm->getRank();
  try {
    result = mainExecute(narg, arg, comm) ? 0 : 1; // code 0 is ok,
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
