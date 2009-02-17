#ifndef _fei_PETSc_LinSysCore_hpp_
#define _fei_PETSc_LinSysCore_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

//
//This is the PETSc implementation of LinearSystemCore.
//


#include <mpi.h>
class PromCRVec;
class PromCRMat;
class PETSc_EssBCData;
class PETSc_OthBCData;
class PETSc_ZeroEquations;
class PETSc_EssBCData;

#include <fei_defs.h>
#include <fei_Data.hpp>
#include <fei_LinearSystemCore.hpp>

class PETSc_LinSysCore : public LinearSystemCore 
{
public:
  PETSc_LinSysCore(MPI_Comm comm);
  virtual ~PETSc_LinSysCore();
  
  //for creating another one, without knowing the run-time type
  //of 'this' one.
  LinearSystemCore* clone();
  
  //int parameters:
  //for setting generic argc/argv style parameters.
  
  int parameters(int numParams, const char*const* params);
  
  int setLookup(Lookup& lookup);
  
  int setGlobalOffsets(int len, int* nodeOffsets, int* eqnOffsets,
			int* blkEqnOffsets);
  
  int setConnectivities(GlobalID elemBlock,
			 int numElements,
			 int numNodesPerElem,
			 const GlobalID* elemIDs,
			 const int* const* connNodes) ;
  
  int setStiffnessMatrices(GlobalID elemBlock,
			    int numElems,
			    const GlobalID* elemIDs,
			    const double *const *const *stiff,
			    int numEqnsPerElem,
			    const int *const * eqnIndices);
  
  int setLoadVectors(GlobalID elemBlock,
		      int numElems,
		      const GlobalID* elemIDs,
		      const double *const * load,
		      int numEqnsPerElem,
		      const int *const * eqnIndices);
  
  int setMatrixStructure(int** ptColIndices, int* ptRowLengths,
			  int** blkColIndices, int* blkRowLengths,
			  int* ptRowsPerBlkRow) ;
  
  int setMultCREqns(int multCRSetID,
		     int numCRs, int numNodesPerCR,
		     int** nodeNumbers, int** eqnNumbers,
		     int* fieldIDs,
		     int* multiplierEqnNumbers);
  
  int setPenCREqns(int penCRSetID,
		    int numCRs, int numNodesPerCR,
		    int** nodeNumbers, int** eqnNumbers,
		    int* fieldIDs);
  
  //int resetMatrixAndVector:
  //don't destroy the structure of the matrix, but set the value 's'
  //throughout the matrix and vectors.
  
  int resetMatrixAndVector(double s);
  int resetMatrix(double s);
  int resetRHSVector(double s);
  
  //int sumIntoSystemMatrix:
  //this is the primary assembly function. The coefficients 'values'
  //are to be accumumlated into (added to any values already in place)
  //global (0-based) equation 'row' of the matrix.
  
  int sumIntoSystemMatrix(int numPtRows, const int* ptRows,
			   int numPtCols, const int* ptCols,
			   int numBlkRows, const int* blkRows,
			   int numBlkCols, const int* blkCols,
			   const double* const* values);
  
  int sumIntoSystemMatrix_private(int numPtRows, const int* ptRows,
			   int numPtCols, const int* ptCols,
			   const double* const* values, int add_type);

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
  
  int sumIntoRHSVector_private(int num,
			const double* values,
			const int* indices, int add_type); 
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
  
  int putNodalFieldData(int fieldID, int fieldSize, int* nodeNumbers,
			 int numNodes, const double* data);
  
  //functions for enforcing boundary conditions.
  int enforceEssentialBC(int* globalEqn,
			  double* alpha,
			  double* gamma, int len);
  
  int enforceRemoteEssBCs(int numEqns, int* globalEqns,
			   int** colIndices, int* colIndLen,
			   double** coefs);
  
  int enforceOtherBC(int* globalEqn, double* alpha,
		      double* beta, double* gamma,
		      int len);
  
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
  
  int writeSystem(const char* name);
 private:
  // mapping for Langrange Multiplier problems
  int eq_map_private( const int geq, int &myeq, int &prim, int &proc ) const;
  int getProc_private( const int geq, int &proc ) const;
  int applyBCs_private();
  int enforceEssentialBC_private( int* globalEqn, double* alpha,
				  double* gamma, int len);
  int enforceOtherBC_private( int* globalEqn, double* alpha,
			      double* beta, double* gamma, int len);

  // data      
  int numProcs_;
  int thisProc_;
  int masterProc_;
  int ndf_;
  int dirty_guess_;   // keep track of initial solution guess in "work_"
  int dirty_system_; // used to rebuild coarse grids. etc. w/ new matrix
  static int nActiveLSC_; 
  static int petsc_init_called_; 
  
  PromCRMat *K_;
  PromCRVec *x_;
  PromCRVec **b_;
  PromCRVec *init_guess_;
  void      *Pen_stiff_void_;
  
  int* rhsIDs_;
  int numRHSs_;
  
  int currentRHS_;
  
  int *proc_gnode_;

  int *proc_globalEq_;
  int localStartRow_;
  int numLocalRows_;
  int localEndRow_;
  int numGlobalRows_;
  
  int *proc_lmEq_;
  int localStartRowLm_;
  int numLocalRowsLm_;
  int localEndRowLm_;
  int numGlobalRowsLm_;
  int maxNodesPerCR_;
  
  int *proc_primEq_;
  int localStartRowPrim_;
  int numLocalRowsPrim_;
  int localEndRowPrim_;
  int numGlobalRowsPrim_;

  MPI_Comm FEI_PETSc_Comm_;
  PETSc_EssBCData *essbcdata_;
  PETSc_OthBCData *othbcdata_;
  PETSc_ZeroEquations *zeroEqs_;

  int verbose_;
};

#endif

