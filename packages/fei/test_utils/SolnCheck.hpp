#ifndef _SolnCheck_hpp_
#define _SolnCheck_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_FillableMat.hpp"

namespace SolnCheck {
  int readSoln(const char* baseName, int np, fei::FillableMat& solution);

  int compareSoln(fei::FillableMat& solution1, fei::FillableMat& solution2,
			 double tol=1.e-3);

  int readMatrix(const char* baseName, int np, fei::FillableMat& matrix);

  int compareMatrices(fei::FillableMat& mat1, fei::FillableMat& mat2);

  int checkSolution(int localProc, int numProcs,
		    const char* solnFileName, const char* checkFileName,
		    const char* extension,
		    int solveCounter);
}

#endif // _SolnCheck_hpp_

