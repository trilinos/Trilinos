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
//

/*! \file VerifyVectorGenerateFiles.hpp
 *  \brief Shared code that verifies the GenerateFiles method.
 */

#ifndef VERIFYVECTORGENERATEFILES_HPP_
#define VERIFYVECTORGENERATEFILES_HPP_

#include <string>

#include <Zoltan2_TestHelpers.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>


template <typename User>
int verifyGenerateFiles(
  Zoltan2::VectorAdapter<User> &ia, 
  const char *fileprefixInput,
  const Teuchos::Comm<int> &comm
)
{
  int fail = 0, gfail=0;

  const char *fileprefixGen = "unitTestOutput";
  ia.generateFiles(fileprefixGen, comm);

  // Only rank zero needs to check the resulting files
  if (comm.getRank() == 0) {

    size_t nIDsGen, nIDsInput;
    size_t nEdgeGen, nEdgeInput;
    char codeGen[4], codeInput[4];

    std::ifstream fpGen, fpInput;
    std::string graphFilenameGen = fileprefixGen;
    graphFilenameGen = graphFilenameGen + ".graph";
    std::string graphFilenameInput = fileprefixInput;
    graphFilenameInput = graphFilenameInput + ".graph";

    // Read header info from generated file
    fpGen.open(graphFilenameGen.c_str(), std::ios::in);
    std::string lineGen;
    std::getline(fpGen, lineGen);
    std::istringstream issGen(lineGen);
    issGen >> nIDsGen >> nEdgeGen >> codeGen;

    // Read header info from input file
    fpInput.open(graphFilenameInput.c_str(), std::ios::in);
    std::string lineInput;
    std::getline(fpInput, lineInput);
    while (lineInput[0]=='#') std::getline(fpInput, lineInput); // skip comments
    std::istringstream issInput(lineGen);
    issInput >> nIDsInput >> nEdgeInput >> codeInput;

    // input file and generated file should have same number of IDs
    if (nIDsGen != nIDsInput) {
      std::cout << "GenerateFiles:  nIDsGen " << nIDsGen
                << " != nIDsInput " << nIDsInput << std::endl;
      fail = 222;
    }

    // Vector adapters don't have edges
    if (!fail && nEdgeGen != 0) {
      std::cout << "GenerateFiles:  nEdgeGen " << nEdgeGen << " != 0" << std::endl;
      fail = 222;
    }

    // Check the weights, if any
    if (!fail && !strcmp(codeGen, "010")) {
      // TODO
      // If input file has weights, compare weights
      // Otherwise, just make sure there are enough weights in file
    }

    fpGen.close();
    fpInput.close();
    
    // check coordinate files
    if (!fail) {
      std::string coordsFilenameGen = fileprefixGen;
      coordsFilenameGen = coordsFilenameGen + ".coords";
      std::string coordsFilenameInput = fileprefixInput;
      coordsFilenameInput = coordsFilenameInput + ".coords";

      fpGen.open(coordsFilenameGen.c_str(), std::ios::in);
      fpInput.open(coordsFilenameInput.c_str(), std::ios::in);

      size_t cnt = 0;
      for (; std::getline(fpGen,lineGen) && std::getline(fpInput,lineInput);) { 

        cnt++;

        // Check each token
        issGen.str(lineGen);
        issInput.str(lineInput);

        while (issGen && issInput) {
          double xGen, xInput;
          issGen >> xGen;
          issInput >> xInput;

          if (xGen != xInput) {
            std::cout << "Coordinates " << xGen << " != " << xInput 
                      << std::endl;
            fail = 333;
          }
        }

        // Check same number of tokens:  are there any left in either line?
        if (issGen || issInput) {
          std::cout << "Dimension of generated file != dimension of input file"
                    << std::endl;
          fail = 334;
        }
      }

      // Did we have the correct number of coordinates?
      if (!fail && cnt != nIDsGen) {
        std::cout << "Number of coordinates read " << cnt 
                  << " != number of IDs " << nIDsGen << std::endl;
        fail = 444;
      }

      fpGen.close();
      fpInput.close();
    }
  }

  gfail = globalFail(comm, fail);
  return gfail;
}

#endif
