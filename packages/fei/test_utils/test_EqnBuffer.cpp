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

#include <fei_macros.hpp>

#include <test_utils/test_EqnBuffer.hpp>
#include <cmath>
#include <fei_CSVec.hpp>
#include <fei_EqnBuffer.hpp>

#undef fei_file
#define fei_file "test_EqnBuffer.cpp"
#include <fei_ErrMacros.hpp>

test_EqnBuffer::test_EqnBuffer(MPI_Comm comm)
 : tester(comm)
{
}

test_EqnBuffer::~test_EqnBuffer()
{
}

int test_EqnBuffer::runtests()
{
  CHK_ERR( test1() );
  CHK_ERR( test2() );
  CHK_ERR( test3() );
  CHK_ERR( test4() );
  return(0);
}

int test_EqnBuffer::test1()
{
  FEI_COUT << "testing EqnBuffer...";

  EqnBuffer eqns;

  std::vector<double> eqnCoefs(3);
  std::vector<int> colIndices(3);

  eqnCoefs[0] = 0.5;  colIndices[0] = 2;
  eqnCoefs[1] = 0.5;  colIndices[1] = 4;
  eqnCoefs[2] = 0.5;   colIndices[2] = 6;

  CHK_ERR( eqns.addIndices(1, &colIndices[0], eqnCoefs.size()) );

  CHK_ERR( eqns.addEqn(1, &eqnCoefs[0], &colIndices[0],
		       eqnCoefs.size(), false) );

  eqnCoefs[0] = 0.5;  colIndices[0] = 1;
  eqnCoefs[1] = 0.5;  colIndices[1] = 3;
  eqnCoefs[2] = 0.5;   colIndices[2] = 5;

  CHK_ERR( eqns.addEqn(7, &eqnCoefs[0], &colIndices[0],
		       eqnCoefs.size(), true) );

  eqnCoefs[0] = 0.25;  colIndices[0] = 2;
  eqnCoefs[1] = 0.25;  colIndices[1] = 3;
  eqnCoefs[2] = 0.5;   colIndices[2] = 6;

  CHK_ERR( eqns.addEqn(8, &eqnCoefs[0], &colIndices[0],
		       eqnCoefs.size(), false) );

  eqns.setNumRHSs(1);

  double coef = 0.0;
  CHK_ERR( eqns.getCoef(7, 3, coef) );
  if (std::abs(coef - 0.5) > 1.e-49) ERReturn(-1);

  CHK_ERR( eqns.removeIndex(7,3) );
  int err = eqns.getCoef(7, 3, coef);
  if (err != -1) ERReturn(-1);

  CHK_ERR( eqns.getCoefAndRemoveIndex(1, 6, coef) );
  if (std::abs(coef - 0.5) > 1.e-49) ERReturn(-1);

  err = eqns.getCoef(1, 6, coef);
  if (err != -1) ERReturn(-1);


  std::vector<int>& eqnNumbers = eqns.eqnNumbers();
  std::vector<fei::CSVec*>& rows = eqns.eqns();

  EqnBuffer* eqnsCopy = eqns.deepCopy();

  std::vector<int>& eqnNumbersCopy = eqnsCopy->eqnNumbers();

  if (eqnNumbersCopy != eqnNumbers) {
    ERReturn(-1);
  }

  CHK_ERR( eqnsCopy->addRHS(1, 0, 1.0, true) );
  CHK_ERR( eqnsCopy->addEqns(eqns, true) );

  std::vector<double> tempCoefs;
  std::vector<int> tempIndices;

  int levelsOfCoupling = 0;
  bool finished = false;
  while(!finished) {
    bool foundCoupling = false;
    for(size_t i=0; i<eqnNumbers.size(); i++) {
      int rowIndex = eqns.isInIndices(eqnNumbers[i]);

      while(rowIndex >= 0) {
	foundCoupling = true;
	coef = 0.0;
	CHK_ERR( eqns.getCoefAndRemoveIndex( eqnNumbers[rowIndex], eqnNumbers[i],
					     coef) );

	std::vector<int>& indicesRef = rows[i]->indices();
	std::vector<double>& coefsRef = rows[i]->coefs();

	int len = indicesRef.size();
	tempCoefs.resize(len);
	tempIndices.resize(len);

	double* tempCoefsPtr = &tempCoefs[0];
	int* tempIndicesPtr = &tempIndices[0];
	double* coefsPtr = &coefsRef[0];
	int* indicesPtr = &indicesRef[0];

	for(int j=0; j<len; ++j) {
	  tempIndicesPtr[j] = indicesPtr[j];
	  tempCoefsPtr[j] = coef*coefsPtr[j];
	}

	CHK_ERR( eqns.addEqn(eqnNumbers[rowIndex], tempCoefsPtr, tempIndicesPtr,
			     len, true) );

	rowIndex = eqns.isInIndices(eqnNumbers[i]);
      }
    }
    if (foundCoupling) ++levelsOfCoupling;
    else finished = true;
  }

  if (levelsOfCoupling != 1) {
    return(-1);
  }

  eqns.resetCoefs();

  delete eqnsCopy;

  FEI_COUT << "ok"<<FEI_ENDL;

  return(0);
}

int test_EqnBuffer::test2()
{

  return(0);
}

int test_EqnBuffer::test3()
{
  return(0);
}

int test_EqnBuffer::test4()
{
  return(0);
}
