/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _FEI_tester_h_
#define _FEI_tester_h_

#include <fei_macros.hpp>
#include <fei_SharedPtr.hpp>
#include <fei_defs.h>

class LibraryWrapper;
class LinearSystemCore;
class FiniteElementData;

#include <test_utils/feitester.hpp>
#include <test_utils/DataReader.hpp>

#include <FEI_Implementation.hpp>

class FEI_tester : public feitester {
 public:
  FEI_tester(fei::SharedPtr<DataReader> data_reader,
	     MPI_Comm comm, int localProc, int numProcs, bool useNewFEI=false);
  ~FEI_tester();

  const char* getName()
    {
      static const char name[] = "FEI_tester";
      return((const char*)name);
    }

  int testInitialization();

  int testLoading();

  int testSolve();

  int testCheckResult();

  void dumpMatrixFiles();

  void setParameter(const char* param);

 private:
  int createFEIinstance(const char* solverName);
  int setIDlists();
  int initializationPhase();
  int normalLoadPhase();
  int aggregateLoadPhase();
  int exerciseResidualNorm();
  int exercisePutFunctions();

  int save_block_node_soln(DataReader& data, FEI& fei,
			   const char* solnFileName, int numProcs,
			   int localProc, int solveCounter);

  int save_block_elem_soln(DataReader& data, FEI& fei,
			   const char* solnFileName,
			   int numProcs, int localProc, int solveCounter);

  int save_multiplier_soln(DataReader& data, FEI& fei,
			   const char* solnFileName,
			   int numProcs, int localProc, int solveCounter);

  int checkSolution(int localProc, int numProcs,
		  const char* solnFileName, const char* checkFileName,
		  const char* extension, int solveCounter);

  int lsc_matrix_check();

  MPI_Comm comm_;

  fei::SharedPtr<FEI> fei_;

  fei::SharedPtr<LibraryWrapper> wrapper_;

  fei::SharedPtr<DataReader> data_;

  int localProc_, numProcs_;

  int numMatrices;
  int* matrixIDs;
  int numRHSs;
  int* rhsIDs;
  bool useNewFEI_;
};

#endif // _FEI_tester_h_
