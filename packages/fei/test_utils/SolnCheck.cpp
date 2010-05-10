/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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

