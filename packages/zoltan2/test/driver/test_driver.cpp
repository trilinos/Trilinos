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
#include "Zoltan2_Tests.hpp"
#include <Zoltan2_Parameters.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_PartitioningSolutionQuality.hpp>
#include <Zoltan2_BasicIdentifierAdapter.hpp>
#include <Zoltan2_BasicVectorAdapter.hpp>
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Teuchos_FileInputSource.hpp>

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <stack>

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
using std::stack;

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
    
//    // update zoltan 2 parameters
//    plist.remove("Zoltan2Parameters");
//    plist.set("Zoltan2Parameters", zoltan2Parameters);
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

void run(const ParameterList &problem,const RCP<const Teuchos::Comm<int> > & comm)
{
    // Major steps in running a problem in zoltan 2
    // 1. Get a test object
    // 2. Solve
    // 3. Check for pass/fail if test metrics are provided
    int rank = comm->getRank();
    if(rank == 0)
        cout << "\nPERFORMING TEST: " << problem.get<string>("Name") << endl;
//    

    // 1. get appropriate problem
    if(!problem.isParameter("Name"))
        std::runtime_error("no test name provided");
    
    Zoltan2Test * test = getZoltan2Test(problem.get<string>("Name"));
    
    if(test == nullptr)
            std::runtime_error("no such test");
    
    // 2. solve
    test->Run(problem, comm);
    
    // 3. pass
    if (rank == 0)
        std::cout << (test->didPass()? "PASS\n" : "FAIL\n") << std::endl;
    
    delete test;
}


int main(int argc, char *argv[])
{
    //------------------------------------------------>>
    // (0) Set up MPI environment
    //------------------------------------------------>>
    Teuchos::GlobalMPISession session(&argc, &argv);
    RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    
    int rank = comm->getRank(); // get rank
    if(rank ==0)
        cout << "\nBEGINNING TEST DRIVER...." << endl;
    
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
        cout << "....FINISHING TEST DRIVER\n" << endl;
    
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
