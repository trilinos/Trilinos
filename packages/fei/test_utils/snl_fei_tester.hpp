/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _snl_fei_tester_h_
#define _snl_fei_tester_h_

#include <fei_macros.hpp>
#include <fei_mpi.h>
#include <fei_SharedPtr.hpp>

#include <test_utils/feitester.hpp>
#include <test_utils/DataReader.hpp>

#include <fei_fwd.hpp>

class snl_fei_tester : public feitester {
 public:
  snl_fei_tester(fei::SharedPtr<DataReader> data_reader,
		 MPI_Comm comm, int localProc, int numProcs);
  ~snl_fei_tester();

  const char* getName()
    {
      static const char name[] = "snl_fei_tester";
      return((const char*)name);
    }

  int testInitialization();

  int testLoading();

  int testSolve();

  int testCheckResult();

  void dumpMatrixFiles();

  void setParameter(const char* param);

 private:
  void defineFieldsAndIDTypes();
  int initElemBlocks();
  int loadElemBlocks();
  int initConstraints();
  int loadConstraints();
  void definePattern(ElemBlock& eb, int& patternID);
  int createLibraryInstance(const char* solverName);

  int save_block_node_soln(DataReader& data, fei::Vector* vec,
			   const char* solnFileName, int numProcs,
			   int localProc, int solveCounter);

  int save_block_elem_soln(DataReader& data, fei::Vector* vec,
			   const char* solnFileName,
			   int numProcs, int localProc, int solveCounter);

  int save_multiplier_soln(DataReader& data, fei::Vector* vec,
			   const char* solnFileName,
			   int numProcs, int localProc, int solveCounter);

  int checkSolution(int localProc, int numProcs,
		  const char* solnFileName, const char* checkFileName,
		  const char* extension, int solveCounter);

  MPI_Comm comm_;

  fei::SharedPtr<fei::Factory> factory_;

  fei::SharedPtr<fei::VectorSpace> vecSpace_;
  fei::SharedPtr<fei::MatrixGraph> matrixGraph_;

  fei::SharedPtr<fei::Matrix> A_;
  fei::SharedPtr<fei::Vector> x_;
  fei::SharedPtr<fei::Vector> b_;

  fei::SharedPtr<fei::LinearSystem> linSys_;

  LinearSystemCore* linSysCore_;
  FiniteElementData* feData_;

  fei::SharedPtr<DataReader> data_;

  std::vector<int> idTypes_;
  int numPatterns_;
  int nodeTypeOffset_, elemTypeOffset_, constraintTypeOffset_;

  int localProc_, numProcs_;
};

#endif // _snl_fei_tester_h_
