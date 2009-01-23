/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <test_utils/test_EqnBuffer.hpp>
#include <cmath>
#include <feiArray.hpp>
#include <fei_SSVec.hpp>
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

  feiArray<double> eqnCoefs(3);
  feiArray<int> colIndices(3);

  eqnCoefs[0] = 0.5;  colIndices[0] = 2;
  eqnCoefs[1] = 0.5;  colIndices[1] = 4;
  eqnCoefs[2] = 0.5;   colIndices[2] = 6;

  CHK_ERR( eqns.addIndices(1, colIndices.dataPtr(), eqnCoefs.length()) );

  CHK_ERR( eqns.addEqn(1, eqnCoefs.dataPtr(), colIndices.dataPtr(),
		       eqnCoefs.length(), false) );

  eqnCoefs[0] = 0.5;  colIndices[0] = 1;
  eqnCoefs[1] = 0.5;  colIndices[1] = 3;
  eqnCoefs[2] = 0.5;   colIndices[2] = 5;

  CHK_ERR( eqns.addEqn(7, eqnCoefs.dataPtr(), colIndices.dataPtr(),
		       eqnCoefs.length(), true) );

  eqnCoefs[0] = 0.25;  colIndices[0] = 2;
  eqnCoefs[1] = 0.25;  colIndices[1] = 3;
  eqnCoefs[2] = 0.5;   colIndices[2] = 6;

  CHK_ERR( eqns.addEqn(8, eqnCoefs.dataPtr(), colIndices.dataPtr(),
		       eqnCoefs.length(), false) );

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


  feiArray<int>& eqnNumbers = eqns.eqnNumbersPtr();
  feiArray<SSVec*>& rows = eqns.eqns();

  EqnBuffer* eqnsCopy = eqns.deepCopy();

  feiArray<int>& eqnNumbersCopy = eqnsCopy->eqnNumbersPtr();

  if (eqnNumbersCopy != eqnNumbers) {
    ERReturn(-1);
  }

  CHK_ERR( eqnsCopy->addRHS(1, 0, 1.0, true) );
  CHK_ERR( eqnsCopy->addEqns(eqns, true) );

  feiArray<double> tempCoefs;
  feiArray<int> tempIndices;

  int levelsOfCoupling = 0;
  bool finished = false;
  while(!finished) {
    bool foundCoupling = false;
    for(int i=0; i<eqnNumbers.length(); i++) {
      int rowIndex = eqns.isInIndices(eqnNumbers[i]);

      while(rowIndex >= 0) {
	foundCoupling = true;
	coef = 0.0;
	CHK_ERR( eqns.getCoefAndRemoveIndex( eqnNumbers[rowIndex], eqnNumbers[i],
					     coef) );

	feiArray<int>& indicesRef = rows[i]->indices();
	feiArray<double>& coefsRef = rows[i]->coefs();

	int len = indicesRef.length();
	tempCoefs.resize(len);
	tempIndices.resize(len);

	double* tempCoefsPtr = tempCoefs.dataPtr();
	int* tempIndicesPtr = tempIndices.dataPtr();
	double* coefsPtr = coefsRef.dataPtr();
	int* indicesPtr = indicesRef.dataPtr();

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
