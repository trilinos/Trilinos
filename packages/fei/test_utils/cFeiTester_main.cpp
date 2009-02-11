/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_fstream.hpp>
#include <fei_utils.hpp>
#include <stdio.h>

#include <feiArray.hpp>
#include <snl_fei_ArrayUtils.hpp>

//A user can enable all solver-support with the macro FEI_ALL_SOLVERS.
//Otherwise, to enable support for any particular solver, build with
//-DHAVE_FEI_<solver> where <solver> is one of the libraries listed below.

#if defined(FEI_ALL_SOLVERS)
#define HAVE_FEI_FETI

#define HAVE_FEI_PETSC

#define HAVE_FEI_AZTECOO
#endif

//ok, so we need to include a solver-support header. We need the prototype
//for the <solver>_LinSysCore_create(...) function.
//

#if defined(HAVE_FEI_AZTECOO)
#include <cfei_aztec.h>
#endif

#if defined(HAVE_FEI_PETSC)
#include <cfei_petsc.h>
#endif

//And we also need the FEI prototypes. e.g., FEI_create(...), etc.
//

#include <cfei.h>

//The following headers are only used in this driver program -- a C
//application would not use/need these.

#include <BCNodeSet.hpp>
#include <CRSet.hpp>
#include <CommNodeSet.hpp>
#include <ElemBlock.hpp>
#include <DataReader.hpp>
#include <fei_test_utils.hpp>
#include <snl_fei_Utils.hpp>
#include <SolnCheck.hpp>

#define fei_file "cFeiTester_main"
#include <fei_ErrMacros.hpp>

//==============================================================================
int cfei_handle_args(int argc, char** argv, int localProc,
		     char*& solverName,
		     char*& paramFileName,
		     char*& inputFileName,
		     char*& solnFileName,
		     char*& checkFileName,
		     int& numSolves);

void cfei_print_data(DataReader& data);

int cfei_normalLoadPhase(DataReader& data, CFEI* cfei);

int cfei_aggregateInitSolveStep(DataReader& data, CFEI* cfei,
				int*& matrixIDs, int& numMatrices,
				int*& rhsIDs, int& numRHSs);

int cfei_aggregateLoadPhase(DataReader& data, CFEI* cfei,
			    int*& matrixIDs, int& numMatrices,
			    int*& rhsIDs, int& numRHSs);

int cfei_exercisePutFunctions(DataReader& data, CFEI* cfei);

void cfei_deleteAggregateStuff(int*& matrixIDs, int& numMatrices,
			       int*& rhsIDs, int& numRHSs);

int cfei_save_block_node_soln(DataReader& data, CFEI* cfei,
			      char* solnFileName,
			      int numProcs, int localProc, int solveCounter);

int cfei_save_block_elem_soln(DataReader& data, CFEI* cfei,
			      char* solnFileName,
			      int numProcs, int localProc, int solveCounter);

int cfei_save_multiplier_soln(DataReader& data, CFEI* cfei,
			      char* solnFileName,
			      int numProcs, int localProc, int solveCounter);

int cfei_checkSolution(int localProc, int numProcs,
		       const char* solnFileName, const char* checkFileName,
		       const char* extension, int solveCounter);

bool cfei_namesMatch(const char* name1, const char* name2);

//==============================================================================
int cFeiTester_main(int argc, char** argv,
		    MPI_Comm comm, int numProcs, int localProc)
{
  int i, masterProc = 0, numSolves = 0, err = 0;
  int numMatrices = 1;
  int* matrixIDs = NULL;
  int numRHSs = 1;
  int* rhsIDs = NULL;

#ifdef SIERRA_BUILD_DATE
  double cpu_time = fei::utils::cpu_time();
#endif

  char* solverName = NULL, *inputFileName = NULL;
  char* paramFileName = NULL;
  char* solnFileName = NULL, *checkFileName = NULL;

  CHK_ERR( cfei_handle_args(argc, argv, localProc, solverName,
			    paramFileName, inputFileName,
			    solnFileName, checkFileName, numSolves) );

  DataReader data;

  char* fullName = new char[strlen(inputFileName)+32];
  sprintf(fullName, "%s.%d.%d", inputFileName, numProcs, localProc);

  //let's read all the data out of the input file.
  CHK_ERR(data.readData(fullName));

  //now we'll get the parameters out of another file.
  sprintf(fullName, "%s.parameters", inputFileName);
  CHK_ERR( data.readData(fullName));

  delete [] fullName;

  //the first line of the input file is a banner, let's print that now.
  if (localProc==masterProc) {
    char* name = new char[strlen(inputFileName)+32];
    sprintf(name, "%s.%d.%d", inputFileName, numProcs, localProc);
    FILE* file = fopen(name, "r");
    delete [] name;

    char line[256];
    fgets(line, 256, file);
    FEI_COUT << line << FEI_ENDL;
    fclose(file);
  }

  if (localProc==masterProc) {
    //      sprintf(fei_version_string, "FEI %s Built: %s", fei_version_number,
    //                                                      fei_date);
    //      FEI_COUT << fei_version_string << FEI_ENDL;
    cfei_print_data(data);
  }

  //ok, all the data is in the 'data' object, so we're ready to start
  //handing it all over to an instantiation of the FEI.

  //first, we have to create an instance of a LinSysCore object...

  LinSysCore* linSys = NULL;
  CFEI* cfei = NULL;


   if (cfei_namesMatch(solverName, "Aztec")) {
#ifdef HAVE_FEI_AZTECOO
     CHK_ERR( Aztec_LinSysCore_create(&linSys, MPI_COMM_WORLD));
#endif
   }
   if (cfei_namesMatch(solverName, "PETSc")) {
#ifdef HAVE_FEI_PETSC
     CHK_ERR( PETSc_LinSysCore_create(&linSys, MPI_COMM_WORLD));
#endif
   }


   //now we're ready to create an fei instance...

   CHK_ERR( FEI_create(&cfei, linSys, MPI_COMM_WORLD, 0));

   CHK_ERR( FEI_parameters(cfei, data.numParams_, data.paramStrings_));

   if (data.solveType_ == 0) {
     CHK_ERR( FEI_setSolveType(cfei, data.solveType_));
   }
   else if (data.solveType_ == 2) {
      CHK_ERR( cfei_aggregateInitSolveStep(data, cfei,
					   matrixIDs, numMatrices, rhsIDs, numRHSs));
   }
   else {
      FEI_COUT << "cFeiTester: bad solveType: " << data.solveType_ << FEI_ENDL;
      return(1);
   }

   CHK_ERR( FEI_initFields(cfei, data.numFields_, data.fieldSizes_,
                           data.fieldIDs_));

   for(i=0; i<data.numElemBlocks_; i++) {
      ElemBlock& block = data.elemBlocks_[i];

      CHK_ERR( FEI_initElemBlock(cfei, block.blockID_,
                              block.numElements_,
                              block.numNodesPerElement_,
                              block.numFieldsPerNode_,
                              block.nodalFieldIDs_,
                              block.numElemDOF_,
                              block.elemDOFFieldIDs_,
				 block.interleaveStrategy_) );

      for(int el=0; el<block.numElements_; el++) {
         CHK_ERR(FEI_initElem(cfei,
                              block.blockID_,
                              block.elemIDs_[el],
                              block.elemConn_[el]));
      }
   }

   for(i=0; i<data.numSharedNodeSets_; i++) {
      CommNodeSet& shNodeSet = data.sharedNodeSets_[i];

      CHK_ERR(FEI_initSharedNodes(cfei, shNodeSet.numNodes_, shNodeSet.nodeIDs_,
                                  shNodeSet.procsPerNode_, shNodeSet.procs_));
   }

   for(i=0; i<data.numCRMultSets_; i++) {
      CRSet& crSet = data.crMultSets_[i];

      for(int j=0; j<1; j++) {
         CHK_ERR( FEI_initCRMult(cfei, crSet.numNodes_,
                               crSet.nodeIDs_[j],
                               crSet.fieldIDs_,
				 &(crSet.crID_)));
      }
   }

   for(i=0; i<data.numCRPenSets_; i++) {
      CRSet& crSet = data.crPenSets_[i];

      for(int j=0; j<1; j++) {
         CHK_ERR( FEI_initCRPen(cfei, crSet.numNodes_,
                              crSet.nodeIDs_[j],
                              crSet.fieldIDs_,
				&(crSet.crID_)));
      }
   }

   CHK_ERR( FEI_initComplete(cfei));

   //************** Initialization Phase is now complete *****************

   for(int solveCounter=1; solveCounter<=numSolves; solveCounter++) {

     CHK_ERR( FEI_resetSystem(cfei, 0.0));
     CHK_ERR( FEI_resetMatrix(cfei, 0.0)); //try out all reset functions,
     CHK_ERR( FEI_resetRHSVector(cfei, 0.0));//just for the sake of coverage...

     // ************ Load Phase (delegated to a function) ***************

     if (data.solveType_ == 0) {
       CHK_ERR( cfei_normalLoadPhase(data, cfei));
     }
     if (data.solveType_ == 2) {
       CHK_ERR( cfei_aggregateLoadPhase(data, cfei,
					matrixIDs, numMatrices, rhsIDs, numRHSs));
     }
     
     // **** let's try out the residualNorm function.
     FEI_COUT << "numFields: " << data.numFields_ << FEI_ENDL;
     double* norms = new double[data.numFields_];
     int n, *fields = new int[data.numFields_];
     CHK_ERR( FEI_residualNorm(cfei, 1, data.numFields_, fields, norms) );
     for(n=0; n<data.numFields_; n++) {
       FEI_COUT << " field["<<n<<"]: " << fields[n]
		<< ", 1-norm: " << norms[n] << FEI_ENDL;
     }

     // **** let's try out the 'put' functions. ******************
     CHK_ERR( cfei_exercisePutFunctions(data, cfei) );


     // ************ Finally do the Solve ************************

     int status;
     err = FEI_solve(cfei, &status);

     FEI_COUT << "FEI_solve, err = " << err << ", status = " 
	      << status << FEI_ENDL;

     if (err) return(err);

     //NOTE: normally a non-zero err return from iterateToSolve is not a
     //fatal error, and probably just indicates non-convergence. But since
     //this program is a unit-testing program, non-convergence probably means
     //that something's wrong with the test anyway.

     double iTime, lTime, sTime, rTime;
     CHK_ERR( FEI_cumulative_cpu_times(cfei, &iTime, &lTime, &sTime, &rTime) );

     if (localProc==masterProc) {
       FEI_COUT << "Cumulative cpu-times (seconds):" << FEI_ENDL
		<< "    init time:        " << iTime << FEI_ENDL
		<< "    load time:        " << lTime << FEI_ENDL
		<< "    solve time:       " << sTime << FEI_ENDL
		<< "    soln return time: " << rTime << FEI_ENDL;

       int itersTaken;
       CHK_ERR( FEI_iterations(cfei, &itersTaken));
       FEI_COUT << "iterations: " << itersTaken << FEI_ENDL;
     }

     CHK_ERR( cfei_save_block_node_soln(data, cfei, solnFileName,
					numProcs, localProc, solveCounter));

      CHK_ERR( cfei_save_block_elem_soln(data, cfei, solnFileName,
					 numProcs, localProc, solveCounter));

      CHK_ERR( cfei_save_multiplier_soln(data, cfei, solnFileName,
					 numProcs, localProc, solveCounter));

      CHK_ERR( FEI_residualNorm(cfei, 2, data.numFields_, fields, norms) );;
      for(n=0; n<data.numFields_; n++) FEI_COUT << " field["<<n<<"]: " << fields[n]
<< ", 2-norm: " << norms[n] << FEI_ENDL;

      delete [] fields;
      delete [] norms;

   } //end of 'for solveCounter<=numSolves loop'

   if (data.solveType_ == FEI_AGGREGATE_SUM) {
      cfei_deleteAggregateStuff(matrixIDs, numMatrices, rhsIDs, numRHSs);
   }

#ifndef FEI_SER
   MPI_Barrier(MPI_COMM_WORLD);
#endif

   if (localProc == 0) {
     CHK_ERR( SolnCheck::checkSolution(localProc, numProcs,
				       solnFileName, checkFileName,
				       "node", 1));

     FEI_COUT << "cFeiTester: TEST PASSED" << FEI_ENDL;

     //This is something the SIERRA runtest tool looks for in test output...
     FEI_COUT << "SIERRA execution successful" << FEI_ENDL;
#ifdef SIERRA_BUILD_DATE
     double elapsedTime = snl_fei::utils::cpu_time() - cpu_time;
     FEI_COUT.setf(IOS_FIXED, IOS_FLOATFIELD);
     FEI_COUT << "Maximum CPU  time: " << elapsedTime << " seconds." << FEI_ENDL;
#endif
   }

   CHK_ERR( FEI_destroy(&cfei))

   LinSysCore_destroy(&linSys);

   delete [] paramFileName;
   delete [] checkFileName;
   delete [] inputFileName;
   delete [] solnFileName;
   delete [] solverName;

   return(err);
}

//==============================================================================
int cfei_handle_args(int argc, char** argv, int localProc,
		char*& solverName,
		char*& paramFileName,
		char*& inputFileName,
		char*& solnFileName,
		char*& checkFileName,
		int& numSolves)
{
  if (localProc == 0) {
    int dashDarg = fei_test_utils::whichArg(argc, argv, "-d");

    if (dashDarg < 0) {
      FEI_CERR << "cFeiTester: argument '-d' not found." << FEI_ENDL;
      return(-1);
    }

    int dashIarg = fei_test_utils::whichArg(argc, argv, "-i");

    if (dashIarg < 0) {
      FEI_CERR << "cFeiTester: argument '-i' not found." << FEI_ENDL;
      return(-1);
    }

    char* path = argv[dashDarg+1];
    int pathLength = strlen(path)+1;
    char* newpath = new char[pathLength+1];
    if (path[pathLength-1] == '/') {
      sprintf(newpath, "%s", path);
    }
    else {
      sprintf(newpath, "%s/", path);
    }
    char* filename = new char[strlen(argv[dashIarg+1])+2+pathLength];
    sprintf(filename, "%s/%s", newpath,argv[dashIarg+1]);

    //we're going to assume that argv[2] is the name of the file in which we can
    //find the following lines (in this order):
    //  SOLVER_LIBRARY <library-name>
    //  PARAM_FILE <param-file-name>
    //  INPUT_FILE <input-file-base-name>
    //  SOLN_FILE <output-file-base-name>
    //  CHECK_FILE <check-file-base-name>
    //  NUM_SOLVES <number-of-solves-to-perform>
    //
    if (argv[dashIarg+1] == NULL) {
      FEI_CERR << "cFeiTester ERROR, argv["<<(dashIarg+1)<<"] == NULL. " << FEI_ENDL;
      return(-1);
    }

    FILE* infile = fopen(filename, "r");
    if (infile == NULL) {
      FEI_CERR << "cFeiTester ERROR, couldn't open file '"<<filename<<"'."<<FEI_ENDL;
      return(-1);
    }

    delete [] filename;

    char* keyString = new char[512];
    char* line = new char[512];
    char* temp = new char[512];

    solverName = new char[512];
    paramFileName = new char[512];
    inputFileName = new char[512];
    solnFileName = new char[512];
    checkFileName = new char[512];

    fgets(line, 512, infile);

    sscanf(line, "%s %s", keyString, solverName);

    if (strcmp("SOLVER_LIBRARY", keyString) != 0) {
      FEI_CERR << "cFeiTester ERROR parsing input file, found key-string "
	   << "'" << keyString << "', expected 'SOLVER_LIBRARY'." << FEI_ENDL;
      return(-1);
    }

    fgets(line, 512, infile);
    sscanf(line, "%s %s", keyString, temp);

    if (strcmp("PARAM_FILE", keyString) != 0) {
      FEI_CERR << "cFeiTester ERROR parsing input file, found key-string "
	   << "'" << keyString << "', expected 'PARAM_FILE'." << FEI_ENDL;
      return(-1);
    }

    sprintf(paramFileName, "%s%s", newpath, temp);

    fgets(line, 512, infile);
    sscanf(line, "%s %s", keyString, temp);

    if (strcmp("INPUT_FILE", keyString) != 0) {
      FEI_CERR << "cFeiTester ERROR parsing input file, found key-string "
	   << "'" << keyString << "', expected 'INPUT_FILE'." << FEI_ENDL;
      return(-1);
    }

    sprintf(inputFileName, "%s%s", newpath, temp);

    fgets(line, 512, infile);
    sscanf(line, "%s %s", keyString, temp);

    if (strcmp("SOLN_FILE", keyString) != 0) {
      FEI_CERR << "cFeiTester ERROR parsing input file, found key-string "
	   << "'" << keyString << "', expected 'SOLN_FILE'." << FEI_ENDL;
      return(-1);
    }

    sprintf(solnFileName, "%s%s", newpath, temp);

    fgets(line, 512, infile);
    sscanf(line, "%s %s", keyString, temp);

    if (strcmp("CHECK_FILE", keyString) != 0) {
      FEI_CERR << "cFeiTester ERROR parsing input file, found key-string "
	   << "'" << keyString << "', expected 'CHECK_FILE'." << FEI_ENDL;
      return(-1);
    }

    sprintf(checkFileName, "%s%s", newpath, temp);

    fgets(line, 512, infile);
    numSolves = -1;
    sscanf(line, "%s %d", keyString, &numSolves);

    if (numSolves < 0) {
      FEI_CERR << "cFeiTester ERROR parsing num-solves." << FEI_ENDL;
      return(-1);
    }

    delete [] newpath;
    delete [] keyString;
    delete [] line;
    delete [] temp;

#ifndef FEI_SER
    int len = strlen(solverName)+1;
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(solverName, len, MPI_CHAR, 0, MPI_COMM_WORLD);

    len = strlen(paramFileName)+1;
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(paramFileName, len, MPI_CHAR, 0, MPI_COMM_WORLD);

    len = strlen(inputFileName)+1;
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(inputFileName, len, MPI_CHAR, 0, MPI_COMM_WORLD);

    len = strlen(solnFileName)+1;
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(solnFileName, len, MPI_CHAR, 0, MPI_COMM_WORLD);

    len = strlen(checkFileName)+1;
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(checkFileName, len, MPI_CHAR, 0, MPI_COMM_WORLD);

    MPI_Bcast(&numSolves, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  }
#ifndef FEI_SER
  else {
    int len = -1;
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    solverName = new char[len];
    MPI_Bcast(solverName, len, MPI_CHAR, 0, MPI_COMM_WORLD);

    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    paramFileName = new char[len];
    MPI_Bcast(paramFileName, len, MPI_CHAR, 0, MPI_COMM_WORLD);

    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    inputFileName = new char[len];
    MPI_Bcast(inputFileName, len, MPI_CHAR, 0, MPI_COMM_WORLD);

    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    solnFileName = new char[len];
    MPI_Bcast(solnFileName, len, MPI_CHAR, 0, MPI_COMM_WORLD);

    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    checkFileName = new char[len];
    MPI_Bcast(checkFileName, len, MPI_CHAR, 0, MPI_COMM_WORLD);

    MPI_Bcast(&numSolves, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
#endif

  FEI_CERR << "cFeiTester: solverName: " << solverName << FEI_ENDL;
  FEI_CERR << "cFeiTester: paramFileName: " << paramFileName << FEI_ENDL;
  FEI_CERR << "cFeiTester: inputFileName: " << inputFileName << FEI_ENDL;
  FEI_CERR << "cFeiTester: solnFileName: " << solnFileName << FEI_ENDL;
  FEI_CERR << "cFeiTester: checkFileName: " << checkFileName << FEI_ENDL;
  FEI_CERR << "cFeiTester: numSolves: " << numSolves << FEI_ENDL;

  return(0);
}

//==============================================================================
void cfei_print_data(DataReader& data) {
   int i;
   FEI_COUT << "data.solveType_: " << data.solveType_ << FEI_ENDL << FEI_ENDL;

   FEI_COUT << "data.numParams_: " << data.numParams_ << ", parameters:" << FEI_ENDL;
   for(i=0; i<data.numParams_; i++) {
      FEI_COUT << "    <" << data.paramStrings_[i] << ">" << FEI_ENDL;
   }

   FEI_COUT << FEI_ENDL << "data.numFields_: " << data.numFields_ << FEI_ENDL;

   for(i=0; i<data.numFields_; i++) { 
      FEI_COUT << "   field " << i << ", fieldID: " << data.fieldIDs_[i]
          << ", fieldSize: " << data.fieldSizes_[i] << FEI_ENDL;
   }

   int neb = data.numElemBlocks_;
   FEI_COUT << "data.numElemBlocks_: " << neb << FEI_ENDL;

   for(i=0; i<neb; i++) {
      ElemBlock& eb = data.elemBlocks_[i];
      FEI_COUT << "elemBlock " << i << ", blockID: " << (int)eb.blockID_ << FEI_ENDL;
      FEI_COUT << "         numElements: " << eb.numElements_ << FEI_ENDL;
      FEI_COUT << "         numNumNodesPerElement: " << eb.numNodesPerElement_
           << FEI_ENDL;

      int j;
      for(j=0; j<eb.numNodesPerElement_; j++) {
         FEI_COUT << "    numFieldsPerNode["<<j<<"]: " << eb.numFieldsPerNode_[j]
             << " fieldIDs: ";
         for(int k=0; k<eb.numFieldsPerNode_[j]; k++) {
            FEI_COUT << eb.nodalFieldIDs_[j][k] << " ";
         }
         FEI_COUT << FEI_ENDL;
      }

      for(j=0; j<eb.numElements_; j++) {
         FEI_COUT << "element " << j << ", elemID " << eb.elemIDs_[j] << FEI_ENDL;
         FEI_COUT << "connectivity: ";
         for(int k=0; k<eb.numNodesPerElement_; k++) {
            FEI_COUT << eb.elemConn_[j][k] << " ";
         }
         FEI_COUT << FEI_ENDL << "stiffness: " << FEI_ENDL;

         int row;
         for(row=0; row<eb.numStiffRows_; row++) {
            FEI_COUT << "   ";
            for(int col=0; col<eb.numStiffRows_; col++) {
               FEI_COUT << eb.elemStiff_[j][row][col] << " ";
            }
            FEI_COUT << FEI_ENDL;
         }

         FEI_COUT << "load: " << FEI_ENDL << "   ";
         for(row=0; row<eb.numStiffRows_; row++) {
            FEI_COUT << eb.elemLoad_[j][row] << " ";
         }
         FEI_COUT << FEI_ENDL << FEI_ENDL;
      }
   }

   FEI_COUT << "numSharedNodeSets_: " << data.numSharedNodeSets_ << FEI_ENDL;
   for(i=0; i<data.numSharedNodeSets_; i++) {
      CommNodeSet& shNodeSet = data.sharedNodeSets_[i];

      FEI_COUT << "set " << i << ", numNodes: " << shNodeSet.numNodes_ << FEI_ENDL;
      for(int j=0; j<shNodeSet.numNodes_; j++) {
         FEI_COUT << "node " << (int)shNodeSet.nodeIDs_[j] << ", procs: ";
         for(int k=0; k<shNodeSet.procsPerNode_[j]; k++) {
            FEI_COUT << shNodeSet.procs_[j][k] << " ";
         }
         FEI_COUT << FEI_ENDL;
      }
   }

   FEI_COUT << "numBCNodeSets_: " << data.numBCNodeSets_ << FEI_ENDL;
   for(i=0; i<data.numBCNodeSets_; i++) {
      BCNodeSet& bcSet = data.bcNodeSets_[i];

      FEI_COUT << "set " << i << ", numNodes: " << bcSet.numNodes_ << FEI_ENDL;
      for(int j=0; j<bcSet.numNodes_; j++) {
         FEI_COUT << " nodeID " << (int)bcSet.nodeIDs_[j] << ", fieldID "
              << bcSet.fieldID_
              << ", offsetIntoField " << bcSet.offsetsIntoField_[j]
              << ", value " << bcSet.prescribed_values_[j] << FEI_ENDL;
      }
   }

   FEI_COUT << "numCRMultSets_: " << data.numCRMultSets_ << FEI_ENDL;
   for(i=0; i<data.numCRMultSets_; i++) {
      CRSet& crSet = data.crMultSets_[i];

      int j;
      FEI_COUT << " set " << i << ", numCRs: 1, numNodes: "
           << crSet.numNodes_ << FEI_ENDL;
      FEI_COUT << " nodes: " << FEI_ENDL;
      for(j=0; j<1; j++) {
         for(int k=0; k<crSet.numNodes_; k++) {
            FEI_COUT << (int)crSet.nodeIDs_[j][k] << " ";
         }
         FEI_COUT << FEI_ENDL;
      }

      FEI_COUT << " weights: " << FEI_ENDL;
      int offset = 0;
      for(j=0; j<crSet.numNodes_; j++) {
         int size = data.getFieldSize(crSet.fieldIDs_[j]);
         for(int k=0; k<size; k++) {
            FEI_COUT << crSet.weights_[offset++] << " ";
         }
         FEI_COUT << FEI_ENDL;
      }
      FEI_COUT << " values: " << FEI_ENDL;
      for(j=0; j<1; j++) {
         FEI_COUT << crSet.values_[j] << " ";
      }
      FEI_COUT << FEI_ENDL;
   }

   FEI_COUT << "numCRPenSets_: " << data.numCRPenSets_ << FEI_ENDL;
   for(i=0; i<data.numCRPenSets_; i++) {
      CRSet& crSet = data.crPenSets_[i];

      int j;
      FEI_COUT << " set " << i << ", numCRs: 1, numNodes: "
           << crSet.numNodes_ << FEI_ENDL;
      FEI_COUT << " nodes: " << FEI_ENDL;
      for(j=0; j<1; j++) {
         for(int k=0; k<crSet.numNodes_; k++) {
            FEI_COUT << (int)crSet.nodeIDs_[j][k] << " ";
         }
         FEI_COUT << FEI_ENDL;
      }

      FEI_COUT << " weights: " << FEI_ENDL;
      int offset = 0;
      for(j=0; j<crSet.numNodes_; j++) {
         int size = data.getFieldSize(crSet.fieldIDs_[j]);
         for(int k=0; k<size; k++) {
            FEI_COUT << crSet.weights_[offset++] << " ";
         }
         FEI_COUT << FEI_ENDL;
      }

      FEI_COUT << " values: " << FEI_ENDL;
      for(j=0; j<1; j++) {
         FEI_COUT << crSet.values_[j] << " ";
      }
      FEI_COUT << FEI_ENDL;

      FEI_COUT << " penValues: " << FEI_ENDL;
      for(j=0; j<1; j++) {
         FEI_COUT << crSet.penValues_[j] << " ";
      }
      FEI_COUT << FEI_ENDL;
   }
}

//------------------------------------------------------------------------------
int cfei_normalLoadPhase(DataReader& data, CFEI* cfei) {
   int i;

   for(i=0; i<data.numBCNodeSets_; i++) {
      BCNodeSet& bcSet = data.bcNodeSets_[i];

      CHK_ERR(FEI_loadNodeBCs(cfei, bcSet.numNodes_,
                     bcSet.nodeIDs_,
                     bcSet.fieldID_,
                     bcSet.offsetsIntoField_,
                     bcSet.prescribed_values_))
   }

   for(i=0; i<data.numElemBlocks_; i++) {
      ElemBlock& block = data.elemBlocks_[i];

      for(int el=0; el<block.numElements_; el++) {

         CHK_ERR(FEI_sumInElem(cfei, block.blockID_,
                          block.elemIDs_[el],
                          block.elemConn_[el],
                          block.elemStiff_[el],
                          block.elemLoad_[el],
                          block.elemFormat_))
      }
   }

   //******** Load Constraint Relation Equations ***********************

   for(i=0; i<data.numCRMultSets_; i++) {
      CRSet& crSet = data.crMultSets_[i];

      for(int j=0; j<1; j++) {
         CHK_ERR(FEI_loadCRMult(cfei, crSet.crID_,
                         crSet.numNodes_,
                         crSet.nodeIDs_[j],
                         crSet.fieldIDs_,
                         crSet.weights_,
                         crSet.values_[j]))
      }
   }

   for(i=0; i<data.numCRPenSets_; i++) {
      CRSet& crSet = data.crPenSets_[i];

      for(int j=0; j<1; j++) {
         CHK_ERR(FEI_loadCRPen(cfei, crSet.crID_,
                        crSet.numNodes_,
                        crSet.nodeIDs_[j],
                        crSet.fieldIDs_,
                        crSet.weights_,
                        crSet.values_[j],
                        crSet.penValues_[j]))
      }
   }

   return(0);
}

//------------------------------------------------------------------------------
int cfei_aggregateInitSolveStep(DataReader& data, CFEI* cfei,
                           int*& matrixIDs, int& numMatrices,
                           int*& rhsIDs, int& numRHSs) {
   int i;

   const char* param = snl_fei::getParam("numMatrices", data.numParams_,
						data.paramStrings_);

   if (param != NULL) {
      sscanf(param+12, "%d", &numMatrices);
   }

   matrixIDs = new int[numMatrices];
   numRHSs = 1;
   rhsIDs = new int[1];
   rhsIDs[0] = 0;

   if (!matrixIDs || !numRHSs || !rhsIDs) return(1);

   for(i=0; i<numMatrices; i++) {
      matrixIDs[i] = i;
   }

   CHK_ERR(FEI_setIDLists(cfei, numMatrices, matrixIDs, numRHSs, rhsIDs))
   return(FEI_setSolveType(cfei, data.solveType_));
}

//------------------------------------------------------------------------------
int cfei_aggregateLoadPhase(DataReader& data, CFEI* cfei,
                       int*& matrixIDs, int& numMatrices,
                       int*& rhsIDs, int& numRHSs) {
   int i;

   for(i=0; i<numMatrices; i++) {
      CHK_ERR(FEI_setCurrentMatrix(cfei, matrixIDs[i]))

      for(int j=0; j<data.numElemBlocks_; j++) {
         ElemBlock& block = data.elemBlocks_[j];

         for(int el=0; el<block.numElements_; el++) {

            CHK_ERR(FEI_sumInElemMatrix(cfei, block.blockID_,
                                   block.elemIDs_[el],
                                   block.elemConn_[el],
                                   block.elemStiff_[el],
                                   block.elemFormat_))
         }
      }
   }

   for(i=0; i<numRHSs; i++) {
     CHK_ERR(FEI_setCurrentRHS(cfei, rhsIDs[i]))

     for(int j=0; j<data.numElemBlocks_; j++) {
        ElemBlock& block = data.elemBlocks_[j];

        for(int el=0; el<block.numElements_; el++) {
	  CHK_ERR(FEI_sumInElemRHS(cfei, block.blockID_,
                                   block.elemIDs_[el],
                                   block.elemConn_[el],
                                   block.elemLoad_[el]))
         }
      }
   }

   //*************** Boundary Condition Node Sets *************************

   for(i=0; i<data.numBCNodeSets_; i++) {
      BCNodeSet& bcSet = data.bcNodeSets_[i];

      CHK_ERR(FEI_loadNodeBCs(cfei, bcSet.numNodes_,
                     bcSet.nodeIDs_,
                     bcSet.fieldID_,
                     bcSet.offsetsIntoField_,
                     bcSet.prescribed_values_))
   }

   double* matScalars = new double[numMatrices];
   for(i=0; i<numMatrices; i++) {
      matScalars[i] = 1.0;
   }

   int rhsScaleID = rhsIDs[0];
   double rhsScalar = 1.0;

   CHK_ERR(FEI_setMatScalars(cfei, numMatrices, matrixIDs, matScalars))
   CHK_ERR(FEI_setRHSScalars(cfei, 1, &rhsScaleID, &rhsScalar))

   delete [] matScalars;

   return(0);
}

//------------------------------------------------------------------------------
int cfei_exercisePutFunctions(DataReader& data, CFEI* cfei) {
   int i, returnValue = 0;

   //first let's exercise the putBlockNodeSolution function.

   for(i=0; i<data.numElemBlocks_; i++) {
      if (returnValue != 0) break;

      ElemBlock& eb = data.elemBlocks_[i];

      GlobalID blockID = eb.blockID_;
      int numActNodes;
      CHK_ERR( FEI_getNumBlockActNodes(cfei, blockID, &numActNodes))
      int numActEqns;
      CHK_ERR( FEI_getNumBlockActEqns(cfei, blockID, &numActEqns))

      GlobalID* nodeList = new GlobalID[numActNodes];

      returnValue = FEI_getBlockNodeIDList(cfei, blockID, numActNodes,
                                           nodeList);

      double* solnValues = new double[numActEqns];
      int* offsets = new int[numActNodes+1];

      if (returnValue == 0) {
         int counter = 0;
         for(int j=0; j<numActNodes; j++) {
            int numValues;
            CHK_ERR( FEI_getNumSolnParams(cfei, nodeList[j], &numValues))

            offsets[j] = counter;
            for(int k=0; k<numValues; k++) {
               solnValues[counter++] = 0.0001;
            }
         }
         offsets[numActNodes] = counter;

         returnValue = FEI_putBlockNodeSolution(cfei, blockID, numActNodes,
                                                nodeList,
                                                offsets, solnValues);
      }

      delete [] nodeList;
      delete [] solnValues;
      delete [] offsets;
   }

   //and now let's exercise the putBlockFieldNodeSolution function.
   int fieldID = data.fieldIDs_[0];
   int fieldSize = data.fieldSizes_[0];

   for(i=0; i<data.numElemBlocks_; i++) {
      if (returnValue != 0) break;

      ElemBlock& eb = data.elemBlocks_[i];

      GlobalID blockID = eb.blockID_;
      int numActNodes;
      CHK_ERR( FEI_getNumBlockActNodes(cfei, blockID, &numActNodes))

      GlobalID* nodeList = new GlobalID[numActNodes];
      double* ddata = new double[numActNodes*fieldSize];

      returnValue = FEI_getBlockNodeIDList(cfei, blockID, numActNodes,
                                           nodeList);

      if (returnValue == 0) {
         int counter = 0;

         for(int j=0; j<numActNodes; j++) {

            for(int k=0; k<fieldSize; k++) {
               ddata[counter++] = 0.0001;
            }
         }

         returnValue = FEI_putBlockFieldNodeSolution(cfei, blockID, fieldID,
                                                      numActNodes, nodeList,
                                                      ddata);
      }

      delete [] nodeList;
      delete [] ddata;
   }

   return(returnValue);
}

//------------------------------------------------------------------------------
void cfei_deleteAggregateStuff(int*& matrixIDs, int& numMatrices,
                       int*& rhsIDs, int& numRHSs)
{
   if (numRHSs > 0) delete [] rhsIDs;
   if (numMatrices > 0) delete [] matrixIDs;

   numMatrices = 0;
}

//------------------------------------------------------------------------------
int cfei_save_block_node_soln(DataReader& data, CFEI* cfei, char* solnFileName,
                          int numProcs, int localProc, int solveCounter)
{
  (void)solveCounter;
   int returnValue = 0;

   std::vector<GlobalID> allNodes;
   feiArray<feiArray<double> > allSolns;

   for(int i=0; i<data.numElemBlocks_; i++) {
      if (returnValue != 0) break;

      ElemBlock& eb = data.elemBlocks_[i];

      GlobalID blockID = eb.blockID_;
      int numActNodes;
      CHK_ERR( FEI_getNumBlockActNodes(cfei, blockID, &numActNodes));

      GlobalID* nodeList = new GlobalID[numActNodes];

      int err = FEI_getBlockNodeIDList(cfei, blockID, numActNodes, nodeList);
      if (err) returnValue = 1;

      for(int jj=0; jj<numActNodes; jj++) {
	snl_fei::sortedListInsert(nodeList[jj], allNodes);
      }

      delete [] nodeList;
   }

   allSolns.resize(allNodes.size());

   for(int i=0; i<data.numElemBlocks_; i++) {
      if (returnValue != 0) break;

      ElemBlock& eb = data.elemBlocks_[i];

      GlobalID blockID = eb.blockID_;
      int numActNodes;
      CHK_ERR( FEI_getNumBlockActNodes(cfei, blockID, &numActNodes));
      int numActEqns;
      CHK_ERR( FEI_getNumBlockActEqns(cfei, blockID, &numActEqns));

      double* solnValues = new double[numActEqns];

      GlobalID* nodeList = new GlobalID[numActNodes];
      int* offsets = new int[numActNodes+1];

      int err = FEI_getBlockNodeIDList(cfei, blockID, numActNodes, nodeList);
      if (err) returnValue = 1;
   
      err = FEI_getBlockNodeSolution(cfei, blockID, numActNodes, nodeList,
                                          offsets, solnValues);
      if (err) returnValue = 1;
      
      if (!err) {
         for(int j=0; j<numActNodes; j++) {

            int ndof = offsets[j+1]-offsets[j];
            int index = snl_fei::binarySearch(nodeList[j], allNodes);

            for(int k=0; k<ndof; k++) {
               allSolns[index].append(solnValues[offsets[j]+k]);
            }
         }
      }

      delete [] nodeList;
      delete [] solnValues;
      delete [] offsets;
   }

   char* fileName = new char[strlen(solnFileName)+32];
   sprintf(fileName, "%s.node.%d.%d.%d", solnFileName, solveCounter,
           numProcs, localProc);
   FEI_OFSTREAM outfile(fileName);

   if (!outfile || outfile.bad()) {
      FEI_COUT << "ERROR opening solution output file " << fileName << FEI_ENDL;
      return(1);
   }

   delete [] fileName;

   for(size_t i=0; i<allNodes.size(); i++) {
     outfile << allNodes[i] << " " << allSolns[i].length() << FEI_ENDL;
     for(int j=0; j<allSolns[i].length(); j++) {
       outfile << allSolns[i][j] << " ";
     }
     outfile << FEI_ENDL;
   }

   return(0);
}

//------------------------------------------------------------------------------
int cfei_save_block_elem_soln(DataReader& data, CFEI* cfei, char* solnFileName,
                          int numProcs, int localProc, int solveCounter)
{
   int returnValue = 0;
   char* fileName = new char[strlen(solnFileName)+32];
   sprintf(fileName, "%s.elem.elem.%d.%d.%d", solnFileName, solveCounter,
           numProcs, localProc);
   FEI_OFSTREAM outfile(fileName);

   if (!outfile || outfile.bad()) {
      FEI_COUT << "ERROR opening elem-solution output file " << fileName << FEI_ENDL;
      return(1);
   }

   delete [] fileName;

   for(int i=0; i<data.numElemBlocks_; i++) {
      if (returnValue != 0) break;

      ElemBlock& eb = data.elemBlocks_[i];

      GlobalID blockID = eb.blockID_;
      int numElems;
      CHK_ERR( FEI_getNumBlockElements(cfei, blockID, &numElems))

      GlobalID* elemIDs = new GlobalID[numElems];
      if (elemIDs==NULL) return(-1);

      int err = FEI_getBlockElemIDList(cfei, blockID, numElems, elemIDs);
      if (err) returnValue = 1;

      int dofPerElem;
      CHK_ERR( FEI_getNumBlockElemDOF(cfei, blockID, &dofPerElem))
      int totalNumElemDOF = numElems*dofPerElem;

      int* offsets = new int[numElems+1];
      if (offsets == NULL) return(-1);

      if (totalNumElemDOF > 0) {
         double* solnValues = new double[totalNumElemDOF];
         if (solnValues == NULL) return(-1);

         err = FEI_getBlockElemSolution(cfei, blockID, numElems, elemIDs,
                                         &dofPerElem, solnValues);
         if (err) returnValue = 1;

         if (!err) {
            for(int j=0; j<numElems; j++) {

               outfile << (int)elemIDs[j] << " " << dofPerElem << FEI_ENDL << "  ";
               for(int k=0; k<dofPerElem; k++) {
                  outfile << solnValues[j*dofPerElem + k] << " ";
               }
               outfile << FEI_ENDL;
            }
         }

         delete [] solnValues;
      }

      delete [] elemIDs;
      delete [] offsets;
   }

   return(0);
}

//------------------------------------------------------------------------------
int cfei_save_multiplier_soln(DataReader& data, CFEI* cfei, char* solnFileName,
                          int numProcs, int localProc, int solveCounter)
{
   int numCRs = 0;

   CHK_ERR( FEI_getNumCRMultipliers(cfei, &numCRs) )

   if (numCRs == 0) return(0);

   char* fileName = new char[strlen(solnFileName)+32];
   sprintf(fileName, "%s.mult.mult.%d.%d.%d", solnFileName, solveCounter,
           numProcs, localProc);
   FEI_OFSTREAM outfile(fileName);

   if (!outfile || outfile.bad()) {
      FEI_COUT << "ERROR opening mult-solution output file " << fileName << FEI_ENDL;
      return(1);
   }

   delete [] fileName;

   int* CRIDs = new int[numCRs];
   double* results = new double[numCRs];

   if (CRIDs==NULL || results==NULL) return(-1);

   CHK_ERR( FEI_getCRMultIDList(cfei, numCRs, CRIDs) )

   CHK_ERR( FEI_getCRMultipliers(cfei, numCRs, CRIDs, results))

   for(int i=0; i<numCRs; i++) {
      outfile << CRIDs[i] << " " << 1 << FEI_ENDL;

      outfile << "   " << results[i] << FEI_ENDL;
   }

   delete [] CRIDs;
   delete [] results;

   return(0);
}

//------------------------------------------------------------------------------
bool cfei_namesMatch(const char* name1, const char* name2)
{
  int len1 = strlen(name1);
  int len2 = strlen(name2);

  if (len1 != len2) return(false);

  if (strncmp(name1, name2, len1)) return(false);

  return(true);
}
