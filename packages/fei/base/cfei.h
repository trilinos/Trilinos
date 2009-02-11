#ifndef _cfei_h_
#define _cfei_h_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
   This is the header for the prototypes of the C interface
   of the Finite Element Interface to Linear Solvers (FEI).

   For explanations of parameters and semantics, see the C++ FEI.h header, or
   doxygen output created there-from.
   With the exception of added create/destroy functions, all FEI functions
   in this header mirror those in the C++ header except as follows:
   - 'FEI_' is pre-pended to function names
   - 'CFEI* cfei' is the first argument for every function
   -  arguments which are references in C++ (e.g., an output int) are
      pointers in this C interface.

   NOTE: ALL functions return an error code which is 0 if successful,
         non-zero if un-successful.

   Noteworthy special case: the solve function may return non-zero
   if the solver failed to converge. This is, of course, a non-fatal 
   situation, and the caller should then check the 'status' argument for
   possible further information (solver-specific/solver-dependent).
------------------------------------------------------------------------------*/

#include "fei_LinSysCore_struct.h"


/*------------------------------------------------------------------------------
   Next, define an opaque CFEI object which will be an FEI context, and will
   be the first argument to all of the C FEI functions which follow in this
   header.
------------------------------------------------------------------------------*/

struct CFEI_struct {
   void* cfei_;
};
typedef struct CFEI_struct CFEI;


/*------------------------------------------------------------------------------
   And now, the function prototypes...
------------------------------------------------------------------------------*/

/* include fei_defs.h for the #defines of parameters such as FEI_LOCAL_TIMES,
  FEI_NODE_MAJOR, etc. */
#include "fei_defs.h"
#include "fei_mpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
   Initialization function. Creates an FEI instance, wrapped in a CFEI pointer.
*/
int FEI_create(CFEI** cfei,
               LinSysCore* lsc,
               MPI_Comm FEI_COMM_WORLD,
               int masterRank);

/* A function to destroy allocated memory. */
int FEI_destroy(CFEI** cfei);

/* A function to destroy LinSysCore objects. (Note that the create function
 is specific to the implementation and so is not provided here.) */
int LinSysCore_destroy(LinSysCore** lsc);

  /* A function to query a named property. */
int LinSysCore_getProperty_double(LinSysCore* lsc,
				  const char* name, double* value);

/*                                     */
/* And now all of the FEI functions... */
/*                                     */

int FEI_parameters(CFEI* cfei, 
                   int numParams, 
                   char **paramStrings);

int FEI_setIDLists(CFEI* cfei,
                   int numMatrices,
                   const int* matrixIDs,
                   int numRHSs,
                   const int* rhsIDs);

int FEI_setSolveType(CFEI* cfei, 
                     int solveType);

int FEI_initFields(CFEI* cfei, 
                   int numFields, 
                   int *fieldSizes, 
                   int *fieldIDs); 

int FEI_initElemBlock(CFEI* cfei, 
                      GlobalID elemBlockID, 
                      int numElements, 
                      int numNodesPerElement, 
                      int* numFieldsPerNode,
                      int** nodalFieldIDs,
                      int numElemDofFieldsPerElement,
                      int* elemDOFFieldIDs,
                      int interleaveStrategy); 

int FEI_initElem(CFEI* cfei, 
                 GlobalID elemBlockID, 
                 GlobalID elemID, 
                 GlobalID *elemConn);

int FEI_initSharedNodes(CFEI* cfei,
                        int numSharedNodes, 
                        GlobalID *sharedNodeIDs,
                        int* numProcsPerNode,
                        int** sharingProcIDs);

int FEI_initCRMult(CFEI* cfei,
                   int numCRNodes,
                   GlobalID* CRNodes,
                   int *CRFields,
                   int* CRID);

int FEI_initCRPen(CFEI* cfei, 
                  int numCRNodes,
                  GlobalID* CRNodes, 
                  int *CRFields,
                  int* CRID); 

int FEI_initSlaveVariable(CFEI* cfei,
			  GlobalID slaveNodeID,
			  int slaveFieldID,
			  int offsetIntoSlaveField,
			  int numMasterNodes,
			  const GlobalID* masterNodeIDs,
			  const int* masterFieldIDs,
			  const double* weights,
			  double rhsValue);

int FEI_initComplete(CFEI* cfei);

int FEI_resetSystem(CFEI* cfei, double s);
int FEI_resetMatrix(CFEI* cfei, double s);
int FEI_resetRHSVector(CFEI* cfei, double s);
int FEI_resetInitialGuess(CFEI* cfei, double s);

int FEI_deleteMultCRs(CFEI* cfei);

int FEI_setCurrentMatrix(CFEI* cfei, int matID);
int FEI_setCurrentRHS(CFEI* cfei, int rhsID);

int FEI_loadNodeBCs(CFEI* cfei,
                    int numNodes,
                    GlobalID *BCNodes,
                    int fieldID,
                    int* offsetsIntoField,
                    double* prescribed_values);

int FEI_loadElemBCs( CFEI* cfei,
                     int numElems,
                     GlobalID *elemIDs,
                     int fieldID,
                     double **alpha,  
                     double **beta,  
                     double **gamma );

int FEI_sumInElem(CFEI* cfei, 
                  GlobalID elemBlockID, 
                  GlobalID elemID,
                  GlobalID* elemConn,
                  double **elemStiffness,
                  double *elemLoad,
                  int elemFormat);

int FEI_sumInElemMatrix(CFEI* cfei, 
                        GlobalID elemBlockID, 
                        GlobalID elemID,
                        GlobalID* elemConn,
                        double **elemStiffness,
                        int elemFormat);

int FEI_sumInElemRHS(CFEI* cfei, 
                     GlobalID elemBlockID, 
                     GlobalID elemID,
                     GlobalID* elemConn,
                     double *elemLoad);

int FEI_loadCRMult(CFEI* cfei, 
                   int CRID, 
                   int numCRNodes,
                   GlobalID *CRNodes,  
                   int *CRFields,
                   double *CRWeights,
                   double CRValue);

int FEI_loadCRPen(CFEI* cfei, 
                  int CRID,
                  int numCRNodes, 
                  GlobalID *CRNodes,
                  int *CRFields,
                  double *CRWeights,  
                  double CRValue,
                  double penValue);

int FEI_setMatScalars(CFEI* cfei,
                      int numScalars,
                      int* IDs,
                      double* scalars);

int FEI_setRHSScalars(CFEI* cfei,
                      int numScalars,
                      int* IDs,
                      double* scalars);

int FEI_loadComplete(CFEI* cfei);

int FEI_residualNorm(CFEI* cfei,
                      int whichNorm,
                     int numFields,
                     int* fieldIDs,
                     double* norms);

int FEI_solve(CFEI* cfei, int* status);

int FEI_iterations(CFEI* cfei, int* itersTaken);

int FEI_getFieldSize(CFEI* cfei, int fieldID, int* numScalars);

int FEI_getEqnNumbers(CFEI* cfei,
		      GlobalID ID,
		      int idType, 
		      int fieldID,
		      int* numEqns,
		      int* eqnNumbers);

int FEI_getNodalFieldSolution(CFEI* cfei,
			      int fieldID,
			      int numNodes,
			      GlobalID* nodeIDs,
			      double* results);

int FEI_getNumLocalNodes(CFEI* cfei, int* numNodes);

int FEI_getLocalNodeIDList(CFEI* cfei,
			   int* numNodes,
			   GlobalID* nodeIDs,
			   int lenNodeIDs);

int FEI_version(CFEI* cfei, const char** versionStringPtr);

int FEI_cumulative_cpu_times(CFEI* cfei,
                              double* initTime,
                              double* loadTime,
                              double* solveTime,
                              double* solnReturnTime);

int FEI_allocatedSize(CFEI* cfei,
                      int* bytes);

int FEI_getBlockNodeSolution(CFEI* cfei, 
                             GlobalID elemBlockID,
                             int numNodes, 
                             GlobalID* nodeIDs, 
                             int *offsets,
                             double *results);

int FEI_getBlockFieldNodeSolution(CFEI* cfei, 
                                  GlobalID elemBlockID,
                                  int fieldID,
                                  int numNodes, 
                                  GlobalID* nodeIDs, 
                                  double *results);

int FEI_getBlockElemSolution(CFEI* cfei, 
                             GlobalID elemBlockID,
                             int numElems, 
                             GlobalID *elemIDs, 
                             int* numElemDOFPerElement,
                             double *results);

int FEI_getNumCRMultipliers(CFEI* cfei, 
                            int* numMultCRs);

int FEI_getCRMultIDList(CFEI* cfei,
                        int numMultCRs,
                        int* multIDs);

int FEI_getCRMultipliers(CFEI* cfei,
                         int numMultCRs,
                         int* CRIDs,
                         double* multipliers);

int FEI_putBlockNodeSolution(CFEI* cfei, 
                             GlobalID elemBlockID,  
                             int numNodes, 
                             GlobalID *nodeIDs, 
                             int *offsets,  
                             double *estimates);

int FEI_putNodalFieldData(CFEI* cfei,
			  int fieldID, 
			  int numNodes, 
			  GlobalID* nodeIDs,
			  double* nodeData);

int FEI_putBlockFieldNodeSolution(CFEI* cfei, 
                                  GlobalID elemBlockID,  
                                  int fieldID, 
                                  int numNodes, 
                                  GlobalID *nodeIDs, 
                                  double *estimates);
         
int FEI_putBlockElemSolution(CFEI* cfei, 
                             GlobalID elemBlockID,  
                             int numElems, 
                             GlobalID *elemIDs, 
                             int dofPerElem,
                             double *estimates);
 
int FEI_putCRMultipliers(CFEI* cfei, 
                         int numMultCRs, 
                         int* CRIDs,
                         double *multEstimates);
 
int FEI_getBlockNodeIDList(CFEI* cfei, 
                           GlobalID elemBlockID,
                           int numNodes, 
                           GlobalID* nodeIDs);

int FEI_getBlockElemIDList(CFEI* cfei, 
                           GlobalID elemBlockID,
                           int numElems,
                           GlobalID* elemIDs);

int FEI_getNumSolnParams(CFEI* cfei,
                         GlobalID nodeID,
                         int* numSolnParams);

int FEI_getNumElemBlocks(CFEI* cfei, int* numElemBlocks);

int FEI_getNumBlockActNodes(CFEI* cfei,
                            GlobalID blockID,
                            int* numNodes);

int FEI_getNumBlockActEqns(CFEI* cfei,
                           GlobalID blockID,
                           int* numEqns);

int FEI_getNumNodesPerElement(CFEI* cfei,
                              GlobalID blockID,
                              int* nodesPerElem);

int FEI_getNumEqnsPerElement(CFEI* cfei,
                             GlobalID blockID,
                             int* numEqns);
 
int FEI_getNumBlockElements(CFEI* cfei,
                            GlobalID blockID,
                            int* numElems);

int FEI_getNumBlockElemDOF(CFEI* cfei,
                           GlobalID blockID,
                           int* DOFPerElem);

#ifdef __cplusplus
}
#endif

#endif
