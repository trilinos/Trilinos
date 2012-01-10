/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/

#include <fei_iostream.hpp>
#include <fei_sstream.hpp>
#include <fei_fstream.hpp>

#include <test_utils/fei_test_utils.hpp>

#include <test_utils/SolnCheck.hpp>

//==============================================================================
int SolnCheck::readSoln(const char* baseName, int np, fei::FillableMat& solution)
{
  for(int i=0; i<np; i++) {
    FEI_OSTRINGSTREAM osstr;
    osstr << baseName << "." << np << "." << i;
    FEI_IFSTREAM infile(osstr.str().c_str());
    if (!infile || infile.bad()) return(-1);

    int node, numDOF;
    double tmpValue;
    infile >> node;
    while(!infile.eof()) {
      infile >> numDOF;

      for(int j=0; j<numDOF; j++) {
        infile >> tmpValue;
        solution.putCoef(node,j,tmpValue);
      }
      infile >> node;
    }
  }

  return(0);
}

//==============================================================================
int SolnCheck::compareSoln(fei::FillableMat& solution1, fei::FillableMat& solution2,
			   double tol)
{
  return(fei_test_utils::compareMatrices(solution1, solution2, tol) );
}

//==============================================================================
int SolnCheck::readMatrix(const char* baseName, int np, fei::FillableMat& matrix)
{
  return( fei_test_utils::readMatrix(baseName, np, matrix) );
}

//==============================================================================
int SolnCheck::compareMatrices(fei::FillableMat& mat1, fei::FillableMat& mat2)
{
  return( fei_test_utils::compareMatrices(mat1, mat2) );
}

//----------------------------------------------------------------------------
int SolnCheck::checkSolution(int localProc, int numProcs,
			     const char* solnFileName,
			     const char* checkFileName,
			     const char* extension,
			     int solveCounter)
{
  if (localProc == 0) {
    fei::FillableMat soln, correctSoln;
    FEI_OSTRINGSTREAM fullSolnFileName;
    FEI_OSTRINGSTREAM fullCheckFileName;

    fullSolnFileName << solnFileName<<"."<<extension<<"."<<solveCounter;
    fullCheckFileName<< checkFileName<<"."<<extension<<".correct."<<solveCounter;

    std::string fullCheck_str = fullCheckFileName.str();
    const char* check_c_str = fullCheck_str.c_str();
    int err = SolnCheck::readSoln(check_c_str, 1, correctSoln);
    if (err != 0) {
      //If we failed to read the data for the "correct" solution, assume that
      //this is simply a portion of the solution (e.g., lagrange multipliers)
      //that this test isn't supposed to compare.
      //FEI_COUT << "FEI_tester: checkSolution: no check-file for '"<<extension
      //    << "' portion of solution, skipping..." << FEI_ENDL;
      return(0);
    }

    std::string fullSoln_str = fullSolnFileName.str();
    const char* soln_c_str = fullSoln_str.c_str();
    err = SolnCheck::readSoln(soln_c_str, numProcs, soln);
    if (err != 0) return(err);

    FEI_COUT << "FEI_tester:checkSolution: checking '"<<extension<<"' solution...";
    int solnCheckCode = SolnCheck::compareSoln(soln, correctSoln);

    if (solnCheckCode != 0) {
      FEI_COUT << "soln file-name: " << soln_c_str << FEI_ENDL;
      FEI_COUT << "soln-check failed, checkFileName="<<checkFileName<<FEI_ENDL;
      FEI_COUT << "soln: " << FEI_ENDL;
      fei::print(FEI_COUT, soln);
      FEI_COUT << "correctSoln file-name: " << check_c_str << FEI_ENDL;
      FEI_COUT << "correctSoln: " << FEI_ENDL;
      fei::print(FEI_COUT, correctSoln);
      return(-1);
    }
    FEI_COUT << " ok"<<FEI_ENDL;
  }
  return(0);
}

