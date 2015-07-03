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
//#include <Zoltan2_TestHelpers.hpp>
//#include <Zoltan2_Parameters.hpp>
//#include <Zoltan2_PartitioningProblem.hpp>
//#include <Zoltan2_PartitioningSolutionQuality.hpp>
//#include <Zoltan2_BasicIdentifierAdapter.hpp>
//#include <Zoltan2_BasicVectorAdapter.hpp>
//#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
//#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
//#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
//
//#include <Teuchos_DefaultComm.hpp>
//#include <Teuchos_XMLObject.hpp>
//#include <Teuchos_FileInputSource.hpp>
//
//#include <Tpetra_MultiVector.hpp>
//#include <Tpetra_CrsMatrix.hpp>

#include <sstream>
#include <string>
#include <iostream>
#include <vector>

//using Teuchos::ParameterList;
//using Teuchos::Comm;
//using Teuchos::RCP;
//using Teuchos::ArrayRCP;
//using Teuchos::XMLObject;

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

int main(int argc, char *argv[])
{
//------------------------------------------------>>
// (0) Set up MPI environment
//------------------------------------------------>>
//    Teuchos::GlobalMPISession session(&argc, &argv);
//    RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
//    
//    int rank = comm->getRank(); // get rank
//    if(rank == 0){
//        cout << "PASS" << endl;
//        return 0;
//    }
//------------------------------------------------>>
// (1) Get and read the input file
// the input file defines tests to be run
//------------------------------------------------>>
//string inputFileName("driver.xml"); // assumes a default input file exists
//if(argc > 1)
//    inputFileName = argv[1]; // user has provided an input file
//
//Teuchos::FileInputSource inputFile(inputFileName);
//XMLObject xmlInput;
//
//// Try to get xmlObject from inputfile
//try{
//    xmlInput = inputFile.getObject();
//}
//catch(exception &e)
//{
////        EXC_ERRMSG("Test Driver error: reading", e); // error reading input
//}


//------------------------------------------------>>
// (2) Get list of valid Zoltan2 Parameters
//------------------------------------------------>>
//    Teuchos::ParameterList zoltan2Parameters;
//    Zoltan2::createAllParameters(zoltan2Parameters);

//------------------------------------------------>>
// (3) Loop over all tests and execute them
//------------------------------------------------>>
//

return 0;
}
