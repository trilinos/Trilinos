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

// taking headers from existing driver template
// will keep or remove as needed
#include <UserInputForTests.hpp>
#include <AdapterForTests.hpp>

#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_PartitioningSolutionQuality.hpp>
#include <Zoltan2_BasicIdentifierAdapter.hpp>
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>

#include <Zoltan2_Parameters.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Teuchos_FileInputSource.hpp>

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <sstream>
#include <string>
#include <iostream>
#include <queue>

using Teuchos::ParameterList;
using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::XMLObject;

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::exception;
using std::ostringstream;
using std::queue;

#define ERRMSG(msg) if (rank == 0){ cerr << "FAIL: " << msg << endl; }
#define EXC_ERRMSG(msg, e) \
if (rank==0){ cerr << "FAIL: " << msg << endl << e.what() << endl;}

// temporary methods for debugging and leanring
void readXML(const XMLObject &xml, const string &title);
void readPList(const ParameterList &plist,
               const string &title,
               bool doc = false,
               bool unused = false);

typedef Zoltan2::MetricValues<zscalar_t> metric_t; // typedef metric_type

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

  zoltan2Parameters.set("compute_metrics", "true");
  // update zoltan 2 parameters
  //        plist.remove("Zoltan2Parameters");
  //        plist.set("Zoltan2Parameters", zoltan2Parameters);
}

// this method takes the user input file and extracts all
// pLists definitions
queue<ParameterList> getParameterLists(const string &inputFileName, int rank)
{
  Teuchos::FileInputSource inputSource(inputFileName);
  cout << "input file source: " << inputFileName << endl;
  XMLObject xmlInput;
  queue<ParameterList> pLists;
  // Try to get xmlObject from inputfile
  try{
    xmlInput = inputSource.getObject();
  }
  catch(exception &e)
  {
    EXC_ERRMSG("Test Driver error: reading", e); // error reading input
  }

  // get the parameter lists for each model
  for(int i = 0; i < xmlInput.numChildren(); i++)
  {
    ParameterList plist;
    xmlToModelPList(xmlInput.getChild(i), plist);
    pLists.emplace(plist);
  }

  return pLists;
}

bool minMaxTest(const metric_t & metric,
                const Teuchos::ParameterList & metricPlist,
                ostringstream &msg)
{
  // run a comparison of min and max agains a given metric
  // return an error message on failure
  bool pass = true;
  string test_name = metric.getName() + " test";
  if (metricPlist.isParameter("lower"))
  {
    double min = metricPlist.get<double>("lower");

    if(metric.getMinImbalance() < min)
    {
      msg << test_name << " FAILED: Minimum imbalance per part, "
      << metric.getMinImbalance() <<
      ", less than specified allowable minimum, " << min;
      pass = false;
    }
  }

  if(metricPlist.isParameter("upper" ) && pass != false) {
    double max = metricPlist.get<double>("upper");
    if (metric.getMaxImbalance() > max)
    {
      msg << test_name << " FAILED: Maximum imbalance per part, "
      << metric.getMaxImbalance() <<
      ", greater than specified allowable maximum, " << max;
      pass = false;
    }

  }

  if(pass){
    msg << test_name << " PASSED.";
    pass = true;
  }

  return pass;
}

void run(const UserInputForTests &uinput, const ParameterList &problem_parameters,const RCP<const Teuchos::Comm<int> > & comm)
{
  // Major steps in running a problem in zoltan 2
  // 1. get an input adapter
  // 2. construct the problem
  // 3. solve the problem
  // 4. analyze metrics
  // 5. clean up

  typedef AdapterForTests::base_adapter_t base_t;
  typedef AdapterForTests::basic_id_t basic_id_t; // basic_identifier_type
  typedef AdapterForTests::xpetra_mv_adapter xpetra_mv_t; // xpetra_mv_type
  typedef AdapterForTests::xcrsGraph_adapter xcrsGraph_t;
  typedef AdapterForTests::xcrsMatrix_adapter xcrsMatrix_t;
  typedef AdapterForTests::basic_vector_adapter basic_vector_t;

  typedef Zoltan2::Problem<base_t> problem_t;
  typedef Zoltan2::PartitioningProblem<base_t> partioning_problem_t; // base abstract type
  typedef Zoltan2::PartitioningProblem<basic_id_t> basic_problem_t; // basic id problem type
  typedef Zoltan2::PartitioningProblem<xpetra_mv_t> xpetra_mv_problem_t; // xpetra_mb problem type
  typedef Zoltan2::PartitioningProblem<xcrsGraph_t> xcrsGraph_problem_t; // xpetra_mb problem type
  typedef Zoltan2::PartitioningProblem<xcrsMatrix_t> xcrsMatrix_problem_t; // xpetra_mb problem type
  typedef Zoltan2::PartitioningProblem<basic_vector_t> basicVector_problem_t; // xpetra_mb problem type


  int rank = comm->getRank();
  if(rank == 0)
    cout << "\nPeforming test: " << problem_parameters.get<string>("Name") << endl;


  ////////////////////////////////////////////////////////////
  // 1. get basic input adapter
  ////////////////////////////////////////////////////////////
  if(!problem_parameters.isParameter("InputAdapterParameters"))
    throw std::runtime_error("Input adapter parameters not provided");
  if(!problem_parameters.isParameter("Zoltan2Parameters"))
    throw std::runtime_error("Zoltan2 probnlem parameters not provided");



  const ParameterList &adapterPlist = problem_parameters.sublist("InputAdapterParameters");
  base_t * ia = AdapterForTests::getAdapterForInput(const_cast<UserInputForTests *>(&uinput), adapterPlist); // a pointer to a basic type
  if(ia == nullptr)
  {
    if(rank == 0)
      cout << "Get adapter for input failed" << endl;

    return;
  }

  ////////////////////////////////////////////////////////////
  // 2. construct partitioning problem
  ////////////////////////////////////////////////////////////
  problem_t * problem;
  string adapter_name = adapterPlist.get<string>("inputAdapter"); // If we are here we have an input adapter, no need to check for one.
  // get Zoltan2 partion parameters
  ParameterList zoltan2_parameters = const_cast<ParameterList &>(problem_parameters.sublist("Zoltan2Parameters"));
  zoltan2_parameters.set("num_global_parts", comm->getSize());

  if(rank == 0){
    readPList(zoltan2_parameters, "Zoltan 2 Params:\n");
    cout <<"\n\n"<<endl;}

#ifdef HAVE_ZOLTAN2_MPI

  if(adapter_name == "BasicIdentifier"){
    problem = reinterpret_cast<problem_t * >(new basic_problem_t(reinterpret_cast<basic_id_t *>(ia),
                                                                 &zoltan2_parameters,
                                                                 MPI_COMM_WORLD));
  }else if(adapter_name == "XpetraMultiVector")
  {
    problem = reinterpret_cast<problem_t * >(new xpetra_mv_problem_t(reinterpret_cast<xpetra_mv_t *>(ia),
                                                                     &zoltan2_parameters,
                                                                     MPI_COMM_WORLD));
  }else if(adapter_name == "XpetraCrsGraph"){
    problem = reinterpret_cast<problem_t * >(new xcrsGraph_problem_t(reinterpret_cast<xcrsGraph_t *>(ia),
                                                                     &zoltan2_parameters,
                                                                     MPI_COMM_WORLD));
  }
  else if(adapter_name == "XpetraCrsMatrix")
  {
    problem = reinterpret_cast<problem_t * >(new xcrsMatrix_problem_t(reinterpret_cast<xcrsMatrix_t *>(ia),
                                                                      &zoltan2_parameters,
                                                                      MPI_COMM_WORLD));
  }  else if(adapter_name == "BasicVector")
  {
    problem = reinterpret_cast<problem_t * >(new basicVector_problem_t(reinterpret_cast<basic_vector_t *>(ia),
                                                                      &zoltan2_parameters,
                                                                      MPI_COMM_WORLD));
  }
  else
    throw std::runtime_error("Input adapter type not avaible, or misspelled.");


#else
  if(adapter_name == "BasicIdentifier"){
    problem = reinterpret_cast<problem_t * >(new basic_problem_t(reinterpret_cast<basic_id_t *>(ia),
                                                                 &zoltan2_parameters));
  }else if(adapter_name == "XpetraMultiVector")
  {
    problem = reinterpret_cast<problem_t * >(new xpetra_mv_problem_t(reinterpret_cast<xpetra_mv_t *>(ia),
                                                                     &zoltan2_parameters));
  }else if(adapter_name == "XpetraCrsGraph"){
    problem = reinterpret_cast<problem_t * >(new xcrsGraph_problem_t(reinterpret_cast<xcrsGraph_t *>(ia),
                                                                     &zoltan2_parameters));
  }
  else if(adapter_name == "XpetraCrsMatrix")
  {
    problem = reinterpret_cast<problem_t * >(new xcrsMatrix_problem_t(reinterpret_cast<xcrsMatrix_t *>(ia),
                                                                      &zoltan2_parameters));
  } else if(adapter_name == "BasicVector")
  {
    problem = reinterpret_cast<problem_t * >(new basicVector_problem_t(reinterpret_cast<basic_vector_t *>(ia),
                                                                       &zoltan2_parameters));
  }
  else
    throw std::runtime_error("Input adapter type not avaible, or misspelled.");
#endif

  ////////////////////////////////////////////////////////////
  // 3. Solve the problem
  ////////////////////////////////////////////////////////////
  if(rank == 0) cout << "...Solving problem..." << endl;
  reinterpret_cast<basic_problem_t *>(problem)->solve();
  if (rank == 0)
    cout << "Problem solved" << endl;

  ////////////////////////////////////////////////////////////
  // 4. Print problem metrics
  ////////////////////////////////////////////////////////////
  if (comm->getRank() == 0)
  {
    // calculate pass fail based on imbalance
    reinterpret_cast<basic_problem_t *>(problem)->printMetrics(cout);

    if(problem_parameters.isParameter("Metrics"))
    {

      ArrayRCP<const metric_t> metrics
      = reinterpret_cast<basic_problem_t *>(problem)->getMetrics();

      // get metric plist
      const ParameterList &metricsPlist = problem_parameters.sublist("Metrics");

      string test_name;
      bool all_tests_pass = true;
      for(int i = 0; i < metrics.size(); i++)
      {
        // print their names...
        ostringstream msg;
        test_name = metrics[i].getName();
        if(metricsPlist.isSublist(test_name))
        {
          if(!minMaxTest(metrics[i], metricsPlist.sublist(test_name), msg))
            all_tests_pass = false;
          cout << msg.str() << endl;

        }
      }

      if(all_tests_pass) cout << "All tests PASSED." << endl;
      else cout << "Testing FAILED." << endl;

    }else{
      cout << "No test metrics provided." << endl;
      reinterpret_cast<basic_problem_t *>(problem)->printMetrics(cout);
    }
  }
  // 4a. timers
  if(zoltan2_parameters.isParameter("timer_output_stream"))
    reinterpret_cast<basic_problem_t *>(problem)->printTimers();

  ////////////////////////////////////////////////////////////
  // 5. Clean up
  ////////////////////////////////////////////////////////////

  if(adapter_name == "XpetraCrsGraph")
    delete reinterpret_cast<xcrsGraph_t *>(ia)->getCoordinateInput();
  if(adapter_name == "XpetraCrsMatrix")
    delete reinterpret_cast<xcrsMatrix_t *>(ia)->getCoordinateInput();

  delete ia;
  delete reinterpret_cast<basic_problem_t *>(problem);
}

int main(int argc, char *argv[])
{

  ////////////////////////////////////////////////////////////
  // (0) Set up MPI environment
  ////////////////////////////////////////////////////////////
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank(); // get rank
  //    if(rank == 0) // exit for rank 0
  //    {
  //        cout << "PASS" << endl;
  //        cout << "FINISHING TEST DRIVER (RANK 0 EXIT)...." << endl;
  //        return 0;
  //    }

  ////////////////////////////////////////////////////////////
  // (1) Get and read the input file
  // the input file defines tests to be run
  ////////////////////////////////////////////////////////////
  string inputFileName("driver.xml"); // assumes a default input file exists
  if(argc > 1)
    inputFileName = argv[1]; // user has provided an input file

  ////////////////////////////////////////////////////////////
  // (2) Get All Input Parameter Lists
  ////////////////////////////////////////////////////////////
  queue<ParameterList>pLists = getParameterLists(inputFileName,rank);

  ////////////////////////////////////////////////////////////
  // (3) Get Input Data Parameters
  ////////////////////////////////////////////////////////////
  const ParameterList inputParameters = pLists.front();
  if(inputParameters.name() != "InputParameters")
  {
    if(rank == 0)
      cout << "InputParameters not defined" << endl;

    return 1;
  }

  // get the user input for all tests
  UserInputForTests uinput(inputParameters, comm,true,true);
  pLists.pop();

  ////////////////////////////////////////////////////////////
  // (4) Perform all tests
  ////////////////////////////////////////////////////////////
  while (!pLists.empty()) {
    run(uinput, pLists.front(), comm);
    pLists.pop();
  }

  return 0;
}

// helper functions

void readXML(const XMLObject &xml, const string &title)
{
  cout << "\nReading XML object " << title << " ...." << endl;
  xml.print(cout , 5);
}

void readPList(const ParameterList &plist, const string &title, bool doc, bool unused)
{
  cout << "\nReading parameter list: " << title << " ...." << endl;
  plist.print(cout, ParameterList::PrintOptions().showDoc(doc).indent(3).showTypes(true));

  if(unused)
  {
    cout << "\nUnused fields: " << title << " ...." << endl;
    plist.unused(cout);
  }
}
