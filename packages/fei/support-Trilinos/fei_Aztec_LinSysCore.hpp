#ifndef _fei_Aztec_LinSysCore_hpp_
#define _fei_Aztec_LinSysCore_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_defs.h>
#include <fei_Data.hpp>
#include <fei_LinearSystemCore.hpp>
#include <fei_SharedPtr.hpp>

#include <string>
#include <map>

#include <az_aztec.h>

//
//This is the Aztec implementation of LinSysCore.
//
namespace fei_trilinos {

class Aztec_Map;
class Aztec_BlockMap;
class Aztec_LSVector;
class AztecDMSR_Matrix;
class AztecDVBR_Matrix;

class Aztec_LinSysCore: public LinearSystemCore {
 public:
   Aztec_LinSysCore(MPI_Comm comm);
   virtual ~Aztec_LinSysCore();

   //for creating another instance of LinearSystemCore without knowing
   //the run-time type of 'this' one.
   LinearSystemCore* clone();

   //int parameters:
   //for setting generic argc/argv style parameters.

   int parameters(int numParams, const char*const * params);

   int setLookup(Lookup& lookup);

   int setGlobalOffsets(int len, int* nodeOffsets, int* eqnOffsets,
                                 int* blkEqnOffsets);

   int setConnectivities(GlobalID elemBlock,
                          int numElements,
                          int numNodesPerElem,
                          const GlobalID* elemIDs,
                          const int* const* connNodes) ;

   int setStiffnessMatrices(GlobalID,
			    int,
			    const GlobalID*,
			    const double *const *const *,
			    int,
			    const int *const *)
     { return(0); }

   int setLoadVectors(GlobalID,
		      int,
		      const GlobalID*,
		      const double *const *,
		      int,
		      const int *const *)
     { return(0); }

   int setMatrixStructure(int** ptColIndices,
                           int* ptRowLengths,
                           int** blkColIndices,
                           int* blkRowLengths,
                           int* ptRowsPerBlkRow) ;

   int setMultCREqns(int,
		     int, int,
		     int**, int**,
		     int*,
		     int*)
     { return(0); }

   int setPenCREqns(int,
		     int, int,
		     int**, int**,
		     int*)
     { return(0); }

   //int resetMatrixAndVector:
   //don't destroy the structure of the matrix, but set the value 's'
   //throughout the matrix and vectors.

   int resetMatrixAndVector(double s);
   int resetMatrix(double s);
   int resetRHSVector(double s);

   //int sumIntoSystemMatrix:
   //this is the primary assembly function. The coefficients 'values'
   //are to be accumumlated into (added to any values already in place)
   //global (0-based) equations in 'ptRows' of the matrix.

   int sumIntoSystemMatrix(int numPtRows, const int* ptRows,
                            int numPtCols, const int* ptCols,
                            int numBlkRows, const int* blkRows,
                            int numBlkCols, const int* blkCols,
                            const double* const* values);

   int sumIntoSystemMatrix(int numPtRows, const int* ptRows,
                            int numPtCols, const int* ptCols,
                            const double* const* values);
   int putIntoSystemMatrix(int numPtRows, const int* ptRows,
                            int numPtCols, const int* ptCols,
                            const double* const* values);

   int getMatrixRowLength(int row, int& length);

   int getMatrixRow(int row, double* coefs, int* indices,
                            int len, int& rowLength);

   //int sumIntoRHSVector:
   //this is the rhs vector equivalent to sumIntoSystemMatrix above.

   int sumIntoRHSVector(int num,
                          const double* values,
                          const int* indices);
   int putIntoRHSVector(int num,
                          const double* values,
                          const int* indices);
   int getFromRHSVector(int num,
                          double* values,
                          const int* indices);

   //int matrixLoadComplete:
   //do any internal synchronization/communication.

   int matrixLoadComplete();
   
   //functions for enforcing boundary conditions.
   int enforceEssentialBC(int* globalEqn,
                           double* alpha,
                           double* gamma, int len);

   int enforceBlkEssentialBC(int* blkEqn, int* blkOffset,
                              double* alpha, double* gamma,
                              int len);

   int enforceRemoteEssBCs(int numEqns, int* globalEqns,
                            int** colIndices, int* colIndLen,
                            double** coefs);

   int enforceBlkRemoteEssBCs(int numEqns, int* blkEqns,
                               int** blkColInds, int** blkColOffsets,
                               int* blkColLens, double** remEssBCCoefs);

   //functions for getting/setting matrix or vector pointers.

   //getMatrixPtr:
   //obtain a pointer to the 'A' matrix. This should be considered a
   //constant pointer -- i.e., this class remains responsible for the
   //matrix (e.g., de-allocation upon destruction). 
   int getMatrixPtr(Data& data);

   //copyInMatrix:
   //replaces the internal matrix with a copy of the input argument, scaled
   //by the coefficient 'scalar'.

   int copyInMatrix(double scalar, const Data& data);

   //copyOutMatrix:
   //passes out a copy of the internal matrix, scaled by the coefficient
   //'scalar'.

   int copyOutMatrix(double scalar, Data& data);

   //sumInMatrix:
   //accumulate (sum) a copy of the input argument into the internal
   //matrix, scaling the input by the coefficient 'scalar'.

   int sumInMatrix(double scalar, const Data& data);

   //get/setRHSVectorPtr:
   //the same semantics apply here as for the matrixPtr functions above.

   int getRHSVectorPtr(Data& data);

   //copyInRHSVector/copyOutRHSVector/sumInRHSVector:
   //the same semantics apply here as for the matrix functions above.

   int copyInRHSVector(double scalar, const Data& data);
   int copyOutRHSVector(double scalar, Data& data);
   int sumInRHSVector(double scalar, const Data& data);

   //destroyMatrixData/destroyVectorData:
   //Utility function for destroying the matrix (or vector) in Data

   int destroyMatrixData(Data& data);
   int destroyVectorData(Data& data);

   //functions for managing multiple rhs vectors
   int setNumRHSVectors(int numRHSs, const int* rhsIDs);

   //int setRHSID:
   //set the 'current' rhs context, assuming multiple rhs vectors.
   int setRHSID(int rhsID);

   //int putInitialGuess:
   //function for setting (a subset of) the initial-guess
   //solution values (i.e., in the 'x' vector).

   int putInitialGuess(const int* eqnNumbers, const double* values,
                        int len);

   //function for getting all of the answers ('x' vector).
   int getSolution(double* answers, int len);

   //function for getting the (single) entry at equation
   //number 'eqnNumber'.
   int getSolnEntry(int eqnNumber, double& answer);

   //function for obtaining the residual vector (using current solution or
   //initial guess)
   int formResidual(double* values, int len);

   //function for launching the linear solver
   int launchSolver(int& solveStatus, int& iterations);

   int putNodalFieldData(int, int, int*, int, const double*)
     { return(0); }

   int writeSystem(const char* name);

 private:        //functions

   int createMiscStuff();

   int allocateMatrix(int** ptColIndices, int* ptRowLengths,
                       int** blkColIndices, int* blkRowLengths,
                       int* ptRowsPerBlkRow);

   int VBRmatPlusScaledMat(AztecDVBR_Matrix* A, double scalar,
                            AztecDVBR_Matrix* source);

   int MSRmatPlusScaledMat(AztecDMSR_Matrix* A, double scalar,
                            AztecDMSR_Matrix* source);

   int createBlockMatrix(int** blkColIndices,
                          int* blkRowLengths,
                          int* ptRowsPerBlkRow);

   int sumIntoBlockRow(int numBlkRows, const int* blkRows,
                        int numBlkCols, const int* blkCols,
                        const double* const* values,
                        int numPtCols,
		       bool overwriteInsteadOfAccumulate);

   int copyBlockRow(int i, const int* blkRows,
                     int numBlkCols, const int* blkCols,
                     const double* const* values,
                     double* coefs);

   int modifyRHSforBCs();

   int explicitlySetDirichletBCs();

   int blockRowToPointRow(int blkRow);

   int getBlockRow(int blkRow, double*& val, int& valLen,
                   int*& blkColInds, int& blkColIndLen,
                   int& numNzBlks, int& numNNZ);

   int getBlkEqnsAndOffsets(int* ptEqns, int* blkEqns, int* blkOffsets,
                             int numEqns);

   int getBlockSize(int blkInd);

   int sumIntoPointRow(int numPtRows, const int* ptRows,
		       int numPtCols, const int* ptColIndices,
		       const double* const* values,
		       bool overwriteInsteadOfAccumulate);

   int sumPointIntoBlockRow(int blkRow, int rowOffset,
			    int blkCol, int colOffset, double value);

   int setMatrixType(const char* name);
   int selectSolver(const char* name);
   int selectPreconditioner(const char* name);
   void setSubdomainSolve(const char* name);
   void setScalingOption(const char* param);
   void setConvTest(const char* param);
   void setPreCalc(const char* param);
   void setTypeOverlap(const char* param);
   void setOverlap(const char* param);
   void setOrthog(const char* param);
   void setAuxVec(const char* param);
   void setAZ_output(const char* param);

   void recordUserParams();

   void checkForParam(const char* paramName, int numParams_,
                      char** paramStrings,
                      double& param);

   void recordUserOptions();

   void checkForOption(const char* paramName, int numParams_,
                       char** paramStrings,
                       int& param);

   int blkRowEssBCMod(int blkEqn, int blkOffset, double* val, int* blkCols,
                        int numCols, int numPtNNZ, double alpha, double gamma);

   int blkColEssBCMod(int blkRow, int blkEqn, int blkOffset, double* val,
                      int* blkCols,
                      int numCols, int numPtNNZ, double alpha, double gamma);

   void setDebugOutput(const char* path, const char* name);

   void debugOutput(const char* msg) const;

   int writeA(const char* name);
   int writeVec(Aztec_LSVector* v, const char* name);

   int messageAbort(const char* msg) const;

 private:            //variables

   MPI_Comm comm_;

   Lookup* lookup_;
   bool haveLookup_;

   int numProcs_;
   int thisProc_;
   int masterProc_;

   int* update_;
   fei::SharedPtr<Aztec_Map> map_;
   AztecDMSR_Matrix *A_;
   AztecDMSR_Matrix *A_ptr_;
   Aztec_LSVector *x_, **b_, *bc_;
   int* essBCindices_;
   int numEssBCs_;
   bool bcsLoaded_;
   bool explicitDirichletBCs_;
   bool BCenforcement_no_column_mod_;
   Aztec_LSVector *b_ptr_;
   bool matrixAllocated_;
   bool vectorsAllocated_;
   bool blkMatrixAllocated_;
   bool matrixLoaded_;
   bool rhsLoaded_;
   bool needNewPreconditioner_;

   bool tooLateToChooseBlock_;
   bool blockMatrix_;
   fei::SharedPtr<Aztec_BlockMap> blkMap_;
   AztecDVBR_Matrix *blkA_;
   AztecDVBR_Matrix *blkA_ptr_;
   int* blkUpdate_;

   AZ_MATRIX *azA_;
   AZ_PRECOND *azP_;
   bool precondCreated_;
   AZ_SCALING *azS_;
   bool scalingCreated_;

   int *aztec_options_;
   double *aztec_params_;
   double *aztec_status_;

   double* tmp_x_;
   bool tmp_x_touched_;
   double** tmp_b_;
   double* tmp_bc_;
   bool tmp_b_allocated_;

   bool ML_Vanek_;
   int numLevels_;

   int* rhsIDs_;
   int numRHSs_;

   int currentRHS_;

   int numGlobalEqns_;
   int localOffset_;
   int numLocalEqns_;

   int numGlobalEqnBlks_;
   int localBlkOffset_;
   int numLocalEqnBlks_;
   int* localBlockSizes_;

   int numNonzeroBlocks_;

   int outputLevel_;
   int numParams_;
   char** paramStrings_;

   std::string name_;
   int debugOutput_;
   int debugFileCounter_;
   char* debugPath_;
   char* debugFileName_;
   FILE* debugFile_;

   std::map<std::string,unsigned>& named_solve_counter_;
};

}//namespace fei_trilinos

#endif

