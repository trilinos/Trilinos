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
#include <vector>

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
using std::vector;

#define ERRMSG(msg) if (rank == 0){ cerr << "FAIL: " << msg << endl; }
#define EXC_ERRMSG(msg, e) \
if (rank==0){ cerr << "FAIL: " << msg << endl << e.what() << endl;}

// temporary methods for debugging and leanring
void readXML(const XMLObject &xml, const string &title);
void readPList(const ParameterList &plist,
               const string &title,
               bool doc = false,
               bool unused = false);



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
// problem definitions
vector<ParameterList> getProblems(const string &inputFileName, int rank)
{
    Teuchos::FileInputSource inputSource(inputFileName);
    XMLObject xmlInput;
    vector<ParameterList> problems;
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
        problems.push_back(plist);
    }
    
    return problems;
}

void problemWithBasicInputAdapter(const ParameterList &problem_parameters,
                                const RCP<const Teuchos::Comm<int> > & comm)
{
    // Major steps in running a problem in zoltan 2
    // 1. get a user input object
    // 2. get an input adapter
    // 3. construct the problem
    // 4. solve the problem
    // 5. analyze metrics
    int rank = comm->getRank();
    if(rank == 0)
        cout << "\nPerforming test: " << problem_parameters.get<string>("Name") << endl;
    
    // 1. get appropriate problem
    //    if(!problem.isParameter("Name"))
    //        std::runtime_error("no test name provided");
    
    // 1. get uinput for this problem
    if(!problem_parameters.isParameter("TestParameters"))
        std::runtime_error("Test parameters not provided");
    
    const ParameterList &input = problem_parameters.sublist("TestParameters");
    UserInputForTests uinput(input,comm,true, true);
    
    // 2. get Zoltan2 partion parameters
    // get Zoltan2 partion parameters
    const ParameterList &zoltan2params = problem_parameters.sublist("Zoltan2Parameters");
    
    // 3. get basic input adapter
    auto ia = AdapterForTests::getBasicIdentiferAdapterForInput(&uinput, input);
    
    //4. construct partitioning problem
#ifdef HAVE_ZOLTAN2_MPI
    Zoltan2::PartitioningProblem<decltype(ia)> problem(&ia,
                                                          const_cast<ParameterList *>(&zoltan2params),
                                                          MPI_COMM_WORLD);
#else
    Zoltan2::PartitioningProblem<adapter_t> problem(ia,
                                                          const_cast<ParameterList *>(&zoltan2params));
#endif
    if (rank == 0)
        cout << "Problem constructed" << endl;
    
    // 5. Solve the problem
    problem.solve();
    if (rank == 0)
        cout << "Problem solved" << endl;
    
    // Print problem metrics
    if (comm->getRank() == 0)
        problem.printMetrics(cout);
    
    if (rank == 0)
        cout << "\n" << endl;
}

void problemWithXpetraMVAdapter(const ParameterList &problem_parameters,
                                const RCP<const Teuchos::Comm<int> > & comm)
{

    // Major steps in running a problem in zoltan 2
    // 1. get a user input object
    // 2. get an input adapter
    // 3. construct the problem
    // 4. solve the problem
    // 5. analyze metrics
    int rank = comm->getRank();
    if(rank == 0)
        cout << "\nPeforming test: " << problem_parameters.get<string>("Name") << endl;
    
    // 1. get appropriate problem
    //    if(!problem.isParameter("Name"))
    //        std::runtime_error("no test name provided");
    
    // 1. get uinput for this problem
    if(!problem_parameters.isParameter("TestParameters"))
        std::runtime_error("Test parameters not provided");
    
    const ParameterList &input = problem_parameters.sublist("TestParameters");
    UserInputForTests uinput(input,comm,true, true);
    
    // 2. get Zoltan2 partion parameters
    // get Zoltan2 partion parameters
    const ParameterList &zoltan2params = problem_parameters.sublist("Zoltan2Parameters");

    // 3. get basic input adapter
    auto ia = AdapterForTests::getXpetraMVAdapterForInput(&uinput, input);
    
    //4. construct partitioning problem
#ifdef HAVE_ZOLTAN2_MPI
    Zoltan2::PartitioningProblem<decltype(ia)> problem(&ia,
                                                       const_cast<ParameterList *>(&zoltan2params),
                                                       MPI_COMM_WORLD);
#else
    Zoltan2::PartitioningProblem<adapter_t> problem(ia,
                                                          const_cast<ParameterList *>(&zoltan2params));
#endif
    if (rank == 0)
        cout << "Problem constructed" << endl;
    
    // 5. Solve the problem
    problem.solve();
    if (rank == 0)
        cout << "Problem solved" << endl;
    
    // Print problem metrics
    if (comm->getRank() == 0)
        problem.printMetrics(cout);
    
    if (rank == 0)
        cout << "\n" << endl;
}

//template<typename base_adapter_t,typename data_t>
//void problemWithXpetraCRSGraphAdapter(base_adapter_t * adapter,
//                                      const ParameterList &zoltan2_params,
//                                      const RCP<const Teuchos::Comm<int> > & comm)
//{
//    int rank = comm->getRank();
//    typedef Zoltan2::XpetraCrsGraphAdapter<data_t> exact_adapter_t;
//    auto tmp = dynamic_cast<exact_adapter_t *>(adapter);
//    Zoltan2::PartitioningProblem<exact_adapter_t> problem(tmp,
//                                                          const_cast<ParameterList *>(&zoltan2_params));
//    if (rank == 0)
//        cout << "Problem constructed" << endl;
//    
//    // 5. Solve the problem
//    problem.solve();
//    if (rank == 0)
//        cout << "Problem solved" << endl;
//    
//    // Print problem metrics
//    if (comm->getRank() == 0)
//        problem.printMetrics(cout);
//}
//
//template<typename base_adapter_t,typename data_t>
//void problemWithXpetraCRSMatrixAdapter(base_adapter_t * adapter,
//                                       const ParameterList &zoltan2_params,
//                                       const RCP<const Teuchos::Comm<int> > & comm)
//{
//    int rank = comm->getRank();
//    typedef Zoltan2::XpetraCrsMatrixAdapter<data_t> exact_adapter_t;
//    auto tmp = dynamic_cast<exact_adapter_t *>(adapter);
//    Zoltan2::PartitioningProblem<exact_adapter_t> problem(tmp,
//                                                          const_cast<ParameterList *>(&zoltan2_params));
//    
//    
//    if (rank == 0)
//        cout << "Problem constructed" << endl;
//    
//    // 5. Solve the problem
//    problem.solve();
//    if (rank == 0)
//        cout << "Problem solved" << endl;
//    
//    // Print problem metrics
//    if (comm->getRank() == 0)
//        problem.printMetrics(cout);
//}

void run(const ParameterList &problem_parameters,const RCP<const Teuchos::Comm<int> > & comm)
{
    
   // Run Test for specified input adapter
    if(!problem_parameters.isParameter("TestParameters"))
        throw std::runtime_error("Test parameters not provided");
    
    const ParameterList &input = problem_parameters.sublist("TestParameters");
    if(!input.isParameter("inputAdapter"))
        throw std::runtime_error("Input adapter not specified");
    
    // pick method for chosen adapter
    string adapter_name = input.get<string>("inputAdapter");
    if(adapter_name == "BasicIdentifier")
        problemWithBasicInputAdapter(problem_parameters, comm);
    else if(adapter_name == "XpetraMultiVector")
        problemWithXpetraMVAdapter(problem_parameters,comm);
//    else if(adapter_name == "XpetraCRSGraph")
//        problemWithXpetraCRSGraphAdapter<inputAdapter_t, data_t>(ia,zoltan2params,comm);
//    else if(adapter_name == "XpetraCRSMatrix")
//        problemWithXpetraCRSMatrixAdapter<inputAdapter_t, data_t>(ia,zoltan2params,comm);
    else
        throw std::runtime_error("Input adapter type not avaible, or misspelled.");
    
    
    //    Zoltan2Test * test = getZoltan2Test(problem.get<string>("Name"));
    
    //    if(test == nullptr)
    //            std::runtime_error("no such test");
    
    // 2. solve
    //    test->Run(problem, comm);
    
    // 3. pass
    //    if (rank == 0)
    //        std::cout << (test->didPass()? "PASS\n" : "FAIL\n") << std::endl;
}


int main(int argc, char *argv[])
{
    //------------------------------------------------>>
    // (0) Set up MPI environment
    //------------------------------------------------>>
    Teuchos::GlobalMPISession session(&argc, &argv);
    RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    
    int rank = comm->getRank(); // get rank
    //    if(rank == 0) // exit for rank 0
    //    {
    //        cout << "PASS" << endl;
    //        cout << "FINISHING TEST DRIVER (RANK 0 EXIT)...." << endl;
    //        return 0;
    //    }
    
    //------------------------------------------------>>
    // (1) Get and read the input file
    // the input file defines tests to be run
    //------------------------------------------------>>
    string inputFileName("driver.xml"); // assumes a default input file exists
    if(argc > 1)
        inputFileName = argv[1]; // user has provided an input file
    
    //------------------------------------------------>>
    // (2) Get Model Parameter Lists
    //------------------------------------------------>>
    vector<ParameterList>tests = getProblems(inputFileName,rank);
    
    //------------------------------------------------>>
    // (3) Loop over all tests and execute them
    //------------------------------------------------>>
    for(auto i = tests.begin(); i != tests.end(); ++i) run(*i, comm);
    
    
    if(rank == 0)
        cout << "....finished tests\n" << endl;
    
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
