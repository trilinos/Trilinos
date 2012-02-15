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
