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


void run(const ParameterList &problem_parameters,const RCP<const Teuchos::Comm<int> > & comm)
{
    // Major steps in running a problem in zoltan 2
    // 1. get a user input object
    // 2. get an input adapter
    // 3. construct the problem
    // 4. solve the problem
    // 5. analyze metrics
    // 6. clean up
    
    typedef AdapterForTests::base_adapter_t base_t;
    typedef AdapterForTests::basic_id_t basic_id_t; // basic_identifier_type
    typedef AdapterForTests::xpetra_mv_adapter xpetra_mv_t; // xpetra_mv_type
    
    typedef Zoltan2::PartitioningProblem<base_t> problem_t; // base abstract type
    typedef Zoltan2::PartitioningProblem<basic_id_t> basic_problem_t; // basic id problem type
    typedef Zoltan2::PartitioningProblem<xpetra_mv_t> xpetra_mv_problem_t; // xpetra_mb problem type
    
    int rank = comm->getRank();
    if(rank == 0)
        cout << "\nPeforming test: " << problem_parameters.get<string>("Name") << endl;
    
    ////////////////////////////////////////////////////////////
    // 1. get uinput for this problem
    ////////////////////////////////////////////////////////////
    if(!problem_parameters.isParameter("TestParameters"))
        std::runtime_error("Test parameters not provided");
    
    const ParameterList &input = problem_parameters.sublist("TestParameters");
    UserInputForTests uinput(input,comm,true, true);
    
    ////////////////////////////////////////////////////////////
    // 2. get basic input adapter
    ////////////////////////////////////////////////////////////
    base_t * ia = AdapterForTests::getAdapterForInput(&uinput, input); // a pointer to a basic type
    
    ////////////////////////////////////////////////////////////
    // 3. construct partitioning problem
    ////////////////////////////////////////////////////////////
    problem_t * problem;
    string adapter_name = input.get<string>("inputAdapter"); // If we are here we have an input adapter, no need to check for one.
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
            
        }
        else if(adapter_name == "XpetraCrsMatrix")
        {
            
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
        
    }
    else if(adapter_name == "XpetraCrsMatrix")
    {
        
    }
    else
        throw std::runtime_error("Input adapter type not avaible, or misspelled.");
#endif

    ////////////////////////////////////////////////////////////
    // 4. Solve the problem
    ////////////////////////////////////////////////////////////
    reinterpret_cast<basic_problem_t *>(problem)->solve();
    if (rank == 0)
        cout << "Problem solved" << endl;
    
    ////////////////////////////////////////////////////////////
    // 5. Print problem metrics
    ////////////////////////////////////////////////////////////
    if (comm->getRank() == 0)
        reinterpret_cast<basic_problem_t *>(problem)->printMetrics(cout);
    
    if (rank == 0)
        cout << "\n" << endl;
    
    ////////////////////////////////////////////////////////////
    // 6. Clean up
    ////////////////////////////////////////////////////////////
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
    // (2) Get Model Parameter Lists
    ////////////////////////////////////////////////////////////
    vector<ParameterList>tests = getProblems(inputFileName,rank);
    
    ////////////////////////////////////////////////////////////
    // (3) Loop over all tests and execute them
    ////////////////////////////////////////////////////////////
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
