/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_mpi.h>

#include <fei_LinearSystemCore.hpp>
#include <fei_FiniteElementData.hpp>
#include <fei_LibraryWrapper.hpp>
#include <FEI_Implementation.hpp>

#include <cfei.h>

#ifdef CHK_CFEI_NULL
#undef CHK_CFEI_NULL
#endif

/* Define a macro to use on entry to function wrappers, to check whether
   the cfei struct holds a valid pointer to a non-NULL FEI instance. */

#define CHK_CFEI_NULL(cfei) \
 { if (cfei == NULL) return(FEI_FATAL_ERROR); \
   if (cfei->cfei_ == NULL) return(FEI_FATAL_ERROR); }

/*
   To see descriptions of parameters and semantics for these functions,
   consult the C++ header fei.h, and/or doxygen output thereof.
*/

/*============================================================================*/
/* Create function. This is not in the C++ interface. Must be called before
   any of the other functions, to instantiate the CFEI* object.
*/
extern "C" int FEI_create(CFEI** cfei,
                          LinSysCore* lsc, 
                          MPI_Comm comm, 
                          int masterRank)
{
   //
   // First let's try to get a LinearSystemCore out of 'lsm'.
   //

   if (lsc == NULL) return(FEI_FATAL_ERROR);

   fei::SharedPtr<LinearSystemCore>* linSys =
     (fei::SharedPtr<LinearSystemCore>*)lsc->lsc_;

   if (linSys->get() == NULL) return(FEI_FATAL_ERROR);

   fei::SharedPtr<LibraryWrapper> libWrap(new LibraryWrapper(*linSys));

   //
   // Now create an FEI instance.
   //

   FEI_Implementation* fei = new FEI_Implementation(libWrap, comm, masterRank);

   //
   // now create the cfei holder.
   //

   *cfei = new CFEI;

   (*cfei)->cfei_ = (void*)fei;

   return(0);
}

/*============================================================================*/
/* Destroy function, to get rid of allocated memory.
*/
extern "C" int FEI_destroy(CFEI** cfei)
{
   if (*cfei == NULL) return(0);

   FEI_Implementation* fei = (FEI_Implementation*)((*cfei)->cfei_);

   delete fei;

   delete *cfei;
   *cfei = NULL;

   return(0);
}

/*============================================================================*/
/* Destroy function for the LinSysCore objects.
*/
extern "C" int LinSysCore_destroy(LinSysCore** lsc) {

   if (*lsc == NULL) return(0);

   fei::SharedPtr<LinearSystemCore>* linSys =
     (fei::SharedPtr<LinearSystemCore>*)((*lsc)->lsc_);

   if (linSys->get() == NULL) return(0);

   delete linSys;
   delete *lsc;
   *lsc = NULL;

   return(0);
}

/*============================================================================*/
/* A function to query a named property.
*/
extern "C" int LinSysCore_getProperty_double(LinSysCore* lsc,
					     const char* name,
					     double* value)
{
  if (lsc == NULL) return(-1);

  fei::SharedPtr<LinearSystemCore>* linSys =
    (fei::SharedPtr<LinearSystemCore>*)(lsc->lsc_);

  if (linSys->get() == NULL) return(-1);

  (void)name;
  (void)value;
  //   return(linSys->getProperty(name, *value));
  return(-1);
}

/*============================================================================*/
extern "C" int FEI_setIDLists(CFEI* cfei,
                              int numMatrices,
                              const int* matrixIDs,
                              int numRHSs,
                              const int* rhsIDs)
{
  CHK_CFEI_NULL(cfei);

  return(
	 ((FEI_Implementation*)(cfei->cfei_))->setIDLists(numMatrices, matrixIDs,
							  numRHSs, rhsIDs)
	 );
}

/*============================================================================*/
extern "C" int FEI_setSolveType(CFEI* cfei, 
                                int solveType)
{
  CHK_CFEI_NULL(cfei);

  return( ((FEI_Implementation*)(cfei->cfei_))->setSolveType(solveType) );
}

/*============================================================================*/
extern "C" int FEI_initFields(CFEI* cfei, 
                              int numFields, 
                              int *fieldSizes, 
                              int *fieldIDs)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->initFields(numFields, 
							  fieldSizes,
							  fieldIDs));
}

/*============================================================================*/
extern "C" int FEI_initElemBlock(CFEI* cfei, 
                                 GlobalID elemBlockID,
                                 int numElements, 
                                 int numNodesPerElement, 
                                 int* numFieldsPerNode, 
                                 int** nodalFieldIDs,
                                 int numElemDOFPerElement, 
                                 int* elemDOFFieldIDs,
                                 int interleaveStrategy)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->initElemBlock(elemBlockID,
		       				     numElements, 
		       				     numNodesPerElement, 
		       				     numFieldsPerNode, 
		       				     nodalFieldIDs,
		       				     numElemDOFPerElement, 
		       				     elemDOFFieldIDs, 
		       				     interleaveStrategy));
}

/*============================================================================*/
extern "C" int FEI_initElem(CFEI* cfei, 
                            GlobalID elemBlockID, 
                            GlobalID elemID,
                            GlobalID *elemConn)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->initElem(elemBlockID, 
							elemID,
							elemConn));
}
 
/*============================================================================*/
extern "C" int FEI_initSharedNodes(CFEI* cfei, 
                                   int numSharedNodes, 
                                   GlobalID *sharedNodeIDs,
                                   int* numProcsPerNode,
                                   int** sharingProcIDs)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->initSharedNodes(numSharedNodes,
							       sharedNodeIDs,
							       numProcsPerNode, 
							       sharingProcIDs));
}
 
/*============================================================================*/
extern "C" int FEI_initCRMult(CFEI* cfei, 
                              int numCRNodes,  
                              GlobalID *CRNodes,  
                              int *CRFields,
                              int* CRID)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->initCRMult(numCRNodes,
							  CRNodes,
							  CRFields, 
							  *CRID));
}

/*============================================================================*/
extern "C" int FEI_initCRPen(CFEI* cfei, 
                             int numCRNodes, 
                             GlobalID *CRNodes, 
                             int *CRFields,
                             int* CRID)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->initCRPen(numCRNodes,
							 CRNodes,
							 CRFields, 
							 *CRID));
}
 
/*============================================================================*/
extern "C" int FEI_initSlaveVariable(CFEI* cfei,
			  GlobalID slaveNodeID,
			  int slaveFieldID,
			  int offsetIntoSlaveField,
			  int numMasterNodes,
			  const GlobalID* masterNodeIDs,
			  const int* masterFieldIDs,
			  const double* weights,
			  double rhsValue)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->initSlaveVariable(slaveNodeID,
								 slaveFieldID,
							    offsetIntoSlaveField,
								 numMasterNodes,
								 masterNodeIDs,
								 masterFieldIDs,
								 weights,
								 rhsValue));
}

/*============================================================================*/
extern "C" int FEI_initCoefAccessPattern(CFEI* cfei,
                               int patternID,
                               int numRowIDs,
                               int* numFieldsPerRow,
                               int** rowFieldIDs,
                               int numColIDsPerRow,
                               int* numFieldsPerCol,
                               int** colFieldIDs,
                               int interleaveStrategy )
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->
            initCoefAccessPattern(patternID, numRowIDs, 
                                  numFieldsPerRow, rowFieldIDs,
                                  numColIDsPerRow,
                                  numFieldsPerCol, colFieldIDs,
                                  interleaveStrategy));
}

/*============================================================================*/
extern "C" int FEI_initCoefAccess( CFEI* cfei,
				   int patternID, int* rowIDTypes,
				   GlobalID* rowIDs,
				   int* colIDTypes,
				   GlobalID* colIDs )
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->
             initCoefAccess(patternID, rowIDTypes, rowIDs, colIDTypes, colIDs));
}

/*============================================================================*/
extern "C" int FEI_initComplete(CFEI* cfei)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->initComplete());
}

/*============================================================================*/
extern "C" int FEI_resetSystem(CFEI* cfei, double s)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->resetSystem(s));
}

/*============================================================================*/
extern "C" int FEI_resetMatrix(CFEI* cfei, double s)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->resetMatrix(s));
}

/*============================================================================*/
extern "C" int FEI_resetRHSVector(CFEI* cfei, double s)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->resetRHSVector(s));
}

/*============================================================================*/
extern "C" int FEI_resetInitialGuess(CFEI* cfei, double s)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->resetInitialGuess(s));
}

/*============================================================================*/
extern "C" int FEI_deleteMultCRs(CFEI* cfei)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->deleteMultCRs());
}

/*============================================================================*/
extern "C" int FEI_setCurrentMatrix(CFEI* cfei, int matID)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->setCurrentMatrix(matID));
}

/*============================================================================*/
extern "C" int FEI_setCurrentRHS(CFEI* cfei, int rhsID)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->setCurrentRHS(rhsID));
}


/*============================================================================*/
extern "C" int FEI_loadNodeBCs(CFEI* cfei, 
                               int numNodes,
                               GlobalID *BCNodes,
                               int fieldID,
                               int* offsetsIntoField,
                               double* prescribed_values)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->loadNodeBCs(numNodes,
                                             BCNodes, fieldID,
                                        offsetsIntoField, prescribed_values));
}

/*============================================================================*/
extern "C" int FEI_loadElemBCs( CFEI* cfei,
                     int numElems,
                     GlobalID *elemIDs,
                     int fieldID,
                     double **alpha,  
                     double **beta,  
                     double **gamma )
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->loadElemBCs(numElems,
                elemIDs, fieldID, alpha, beta, gamma));
}
 
/*============================================================================*/
extern "C" int FEI_sumInElem(CFEI* cfei, 
                             GlobalID elemBlockID, 
                             GlobalID elemID,  
                             GlobalID *elemConn,
                             double **elemStiffness,
                             double *elemLoad,
                             int elemFormat)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->sumInElem(elemBlockID, 
                                           elemID,  
                                           elemConn,
                                           elemStiffness,
                                           elemLoad,
                                           elemFormat));
}

/*============================================================================*/
extern "C" int FEI_sumInElemMatrix(CFEI* cfei,
                                   GlobalID elemBlockID,
                                   GlobalID elemID,
                                   GlobalID *elemConn,
                                   double **elemStiffness,
                                   int elemFormat)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->sumInElemMatrix(elemBlockID,
                                                 elemID,
                                                 elemConn,
                                                 elemStiffness,
                                                 elemFormat));
}

/*============================================================================*/
extern "C" int FEI_sumInElemRHS(CFEI* cfei,
                                GlobalID elemBlockID,
                                GlobalID elemID,
                                GlobalID *elemConn,
                                double *elemLoad)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->sumInElemRHS(elemBlockID,
                                              elemID,
                                              elemConn,
                                              elemLoad));
}

/*============================================================================*/
extern "C" int FEI_loadCRMult(CFEI* cfei, 
                              int CRID, 
                              int numCRNodes,
                              GlobalID *CRNodes,  
                              int *CRFields,
                              double *CRWeights,
                              double CRValue)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->loadCRMult(CRID, 
                                          numCRNodes,
                                          CRNodes,  
                                          CRFields,
                                          CRWeights, 
                                          CRValue));
}

/*============================================================================*/
extern "C" int FEI_loadCRPen(CFEI* cfei, 
                             int CRID,
                             int numCRNodes, 
                             GlobalID *CRNodes,
                             int *CRFields,
                             double *CRWeights,  
                             double CRValue,
                             double penValue)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->loadCRPen(CRID, 
                                         numCRNodes,
                                         CRNodes,  
                                         CRFields,
                                         CRWeights, 
                                         CRValue, 
                                         penValue));
}

/*============================================================================*/
extern "C" int FEI_sumIntoMatrix(CFEI* cfei,
				 int patternID,
				 int* rowIDTypes,
				 GlobalID* rowIDs,
				 int* colIDTypes,
				 GlobalID* colIDs,
				 double** matrixEntries)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->sumIntoMatrix(patternID,
                   rowIDTypes, rowIDs, colIDTypes, colIDs, matrixEntries));
}

/*============================================================================*/
extern "C" int FEI_getFromMatrix(CFEI* cfei,
				 int patternID,
				 int* rowIDTypes,
				 GlobalID* rowIDs,
				 int* colIDTypes,
				 GlobalID* colIDs,
				 double** matrixEntries)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->getFromMatrix(patternID,
           rowIDTypes, rowIDs, colIDTypes, colIDs, matrixEntries));
}

/*============================================================================*/
extern "C" int FEI_putIntoMatrix(CFEI* cfei, int patternID,
				 int* rowIDTypes,
				 GlobalID* rowIDs,
				 int* colIDTypes,
				 GlobalID* colIDs,
				 double* * matrixEntries)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->putIntoMatrix(patternID,
                   rowIDTypes, rowIDs, colIDTypes, colIDs, matrixEntries));
}

/*============================================================================*/
extern "C" int FEI_sumIntoRHS(CFEI* cfei, int patternID,
			      int* IDTypes,
                               GlobalID* IDs,
                               double* vectorEntries)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->sumIntoRHS(patternID,
                IDTypes, IDs, vectorEntries));
}

/*============================================================================*/
extern "C" int FEI_getFromRHS(CFEI* cfei, int patternID,
			      int* IDTypes,
                               GlobalID* IDs,
                              double* vectorEntries)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->getFromRHS(patternID,
                 IDTypes, IDs, vectorEntries));
}

/*============================================================================*/
extern "C" int FEI_putIntoRHS(CFEI* cfei, int patternID,
			      int* IDTypes,
                           GlobalID* IDs,
                           double* vectorEntries)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->putIntoRHS(patternID,
					       IDTypes, IDs, vectorEntries));
}

/*============================================================================*/
extern "C" int FEI_parameters(CFEI* cfei, 
                              int numParams, 
                              char **paramStrings)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->parameters(numParams,
							  paramStrings));
}
 
/*============================================================================*/
extern "C" int FEI_setMatScalars(CFEI* cfei,
                                 int numScalars,
                                 int* IDs,
                                 double* scalars)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->setMatScalars(numScalars,
                                               IDs,
                                               scalars));
}

/*============================================================================*/
extern "C" int FEI_setRHSScalars(CFEI* cfei,
                                 int numScalars,
                                 int* IDs,
                                 double* scalars)
{
  CHK_CFEI_NULL(cfei);

   return(((FEI_Implementation*)(cfei->cfei_))->setRHSScalars(numScalars,
                                               IDs,
                                               scalars));
}

/*============================================================================*/
extern "C" int FEI_residualNorm(CFEI* cfei,
                                int whichNorm,
                                int numFields,
                                int* fieldIDs,
                                double* norms)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->residualNorm(whichNorm,
                                              numFields,
                                              fieldIDs,
                                              norms));
}

/*============================================================================*/
extern "C" int FEI_loadComplete(CFEI* cfei)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->loadComplete());
}

/*============================================================================*/
extern "C" int FEI_solve(CFEI* cfei, int* status)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->solve(*status));
}

/*============================================================================*/
extern "C" int FEI_iterations(CFEI* cfei, int* itersTaken)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->iterations(*itersTaken));
}

/*============================================================================*/
extern "C" int FEI_getFieldSize(CFEI* cfei, int fieldID, int* numScalars)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->getFieldSize(fieldID, 
							    *numScalars));
}

/*============================================================================*/
extern "C" int FEI_getEqnNumbers(CFEI* cfei,
				 GlobalID ID,
				 int idType,
				 int fieldID,
				 int* numEqns,
				 int* eqnNumbers)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->getEqnNumbers(ID,
							     idType,
							     fieldID,
							     *numEqns,
							     eqnNumbers));
}

/*============================================================================*/
extern "C" int FEI_getNodalFieldSolution(CFEI* cfei,
					 int fieldID,
					 int numNodes,
					 GlobalID* nodeIDs,
					 double* results)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->getNodalFieldSolution(
						     fieldID,
						     numNodes,
						     nodeIDs,
						     results));
}

/*============================================================================*/
extern "C" int FEI_getNumLocalNodes(CFEI* cfei, int* numNodes)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->getNumLocalNodes(*numNodes));
}

/*============================================================================*/
extern "C" int FEI_getLocalNodeIDList(CFEI* cfei,
				      int* numNodes,
				      GlobalID* nodeIDs,
				      int lenNodeIDs)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->getLocalNodeIDList(*numNodes,
								  nodeIDs,
								  lenNodeIDs));
}

/*============================================================================*/
extern "C" int FEI_version(CFEI* cfei, const char** versionStringPtr)
{
  return(((FEI_Implementation*)(cfei->cfei_))->version(*versionStringPtr));
}

/*============================================================================*/
extern "C" int FEI_cumulative_cpu_times(CFEI* cfei,
					double* initTime,
					double* loadTime,
					double* solveTime,
					double* solnReturnTime)
{
   CHK_CFEI_NULL(cfei);

   return(((FEI_Implementation*)(cfei->cfei_))->cumulative_cpu_times(*initTime,
								     *loadTime,
								     *solveTime,
								     *solnReturnTime));
}

/*============================================================================*/
extern "C" int FEI_getBlockNodeSolution(CFEI* cfei, 
                                        GlobalID elemBlockID,
                                        int numNodes, 
                                        GlobalID *nodeIDs, 
                                        int *offsets,
                                        double *results)
{
  CHK_CFEI_NULL(cfei);

   return(((FEI_Implementation*)(cfei->cfei_))->getBlockNodeSolution(elemBlockID,
                                                      numNodes, 
                                                      nodeIDs, 
                                                      offsets, 
                                                      results));
}

/*============================================================================*/
extern "C" int FEI_getBlockFieldNodeSolution(CFEI* cfei, 
                                         GlobalID elemBlockID,
                                         int fieldID,
                                         int numNodes, 
                                         GlobalID *nodeIDs, 
                                         double *results)
{
  CHK_CFEI_NULL(cfei);

   return(((FEI_Implementation*)(cfei->cfei_))->getBlockFieldNodeSolution(elemBlockID,
                                                         fieldID, 
                                                         numNodes, 
                                                         nodeIDs, 
                                                         results));
}
 
/*============================================================================*/
extern "C" int FEI_getBlockElemSolution(CFEI* cfei, 
                                        GlobalID elemBlockID,
                                        int numElems, 
                                        GlobalID *elemIDs, 
                                        int* numElemDOFPerElement,
                                        double *results)
{
  CHK_CFEI_NULL(cfei);

   return(((FEI_Implementation*)(cfei->cfei_))->getBlockElemSolution(elemBlockID,
                                                      numElems, 
                                                      elemIDs, 
                                                     *numElemDOFPerElement,
                                                      results));
}

/*============================================================================*/
extern "C" int FEI_getNumCRMultipliers(CFEI* cfei, 
                                       int* numMultCRs) 
{
  CHK_CFEI_NULL(cfei);

   return(((FEI_Implementation*)(cfei->cfei_))->getNumCRMultipliers(*numMultCRs));
}

/*============================================================================*/
extern "C" int FEI_getCRMultIDList(CFEI* cfei, 
                                   int numMultCRs, 
                                   int* multIDs)
{
  CHK_CFEI_NULL(cfei);

   return(((FEI_Implementation*)(cfei->cfei_))->getCRMultIDList(numMultCRs,
                                                 multIDs));
}
 
/*============================================================================*/
extern "C" int FEI_getCRMultipliers(CFEI* cfei, 
                                    int numMultCRs, 
                                    int* CRIDs, 
                                    double* multipliers)
{
  CHK_CFEI_NULL(cfei);

   return(((FEI_Implementation*)(cfei->cfei_))->getCRMultipliers(numMultCRs, 
                                                                 CRIDs,
                                                                 multipliers));
}

/*============================================================================*/
extern "C" int FEI_putBlockNodeSolution(CFEI* cfei,
                                        GlobalID elemBlockID,
                                        int numNodes,
                                        GlobalID* nodeIDs,
                                        int* offsets,
                                        double* estimates)
{
  CHK_CFEI_NULL(cfei);

   return(((FEI_Implementation*)(cfei->cfei_))->putBlockNodeSolution(elemBlockID,
                                                      numNodes, 
                                                      nodeIDs,
                                                      offsets,
                                                      estimates));
}


/*============================================================================*/
extern "C" int FEI_putNodalFieldData(CFEI* cfei,
				     int fieldID, 
				     int numNodes, 
				     GlobalID* nodeIDs,
				     double* nodeData)
{
  CHK_CFEI_NULL(cfei);

   return(((FEI_Implementation*)(cfei->cfei_))->putNodalFieldData(fieldID, 
								  numNodes, 
								  nodeIDs,
								  nodeData));
}

/*============================================================================*/
extern "C" int FEI_putBlockFieldNodeSolution(CFEI* cfei,
					     GlobalID elemBlockID,
				     int fieldID, 
				     int numNodes, 
				     GlobalID* nodeIDs,
				     double* estimates)
{
  CHK_CFEI_NULL(cfei);

   return(((FEI_Implementation*)(cfei->cfei_))->putBlockFieldNodeSolution(elemBlockID,
								  fieldID, 
								  numNodes, 
								  nodeIDs,
								  estimates));
}

/*============================================================================*/
extern "C" int FEI_putBlockElemSolution(CFEI* cfei,
                                        GlobalID elemBlockID,
                                        int numElems, 
                                        GlobalID* elemIDs,
                                        int dofPerElem,
                                        double* estimates)
{
  CHK_CFEI_NULL(cfei);

   return(((FEI_Implementation*)(cfei->cfei_))->putBlockElemSolution(elemBlockID,
                                                      numElems,
                                                      elemIDs, 
                                                      dofPerElem,
                                                      estimates));
}

/*============================================================================*/
extern "C" int FEI_putCRMultipliers(CFEI* cfei,
                                    int numMultCRs,
                                    int* CRIDs, 
                                    double* multEstimates)
{
  CHK_CFEI_NULL(cfei);

   return(((FEI_Implementation*)(cfei->cfei_))->putCRMultipliers(numMultCRs,
                                                  CRIDs, 
                                                  multEstimates));
}

/*============================================================================*/
extern "C" int FEI_getBlockNodeIDList(CFEI* cfei,
                                      GlobalID elemBlockID,
                                      int numNodes,
                                      GlobalID *nodeIDs)
{
  CHK_CFEI_NULL(cfei);

   return(((FEI_Implementation*)(cfei->cfei_))->getBlockNodeIDList(elemBlockID,
                                                    numNodes,
                                                    nodeIDs));
}

/*============================================================================*/
extern "C" int FEI_getBlockElemIDList(CFEI* cfei, 
                                      GlobalID elemBlockID, 
                                      int numElems, 
                                      GlobalID* elemIDs)
{
  CHK_CFEI_NULL(cfei);

   return(((FEI_Implementation*)(cfei->cfei_))->getBlockElemIDList(elemBlockID,
                                                    numElems, 
                                                    elemIDs));
}

/*============================================================================*/
extern "C" int FEI_getNumSolnParams(CFEI* cfei, 
                                    GlobalID nodeID,
                                    int* numSolnParams)
{
   CHK_CFEI_NULL(cfei);

   return(((FEI_Implementation*)(cfei->cfei_))->getNumSolnParams(nodeID, *numSolnParams));
}

/*============================================================================*/
extern "C" int FEI_getNumElemBlocks(CFEI* cfei,
                                    int* numElemBlocks)
{
   CHK_CFEI_NULL(cfei);

   return(((FEI_Implementation*)(cfei->cfei_))->getNumElemBlocks(*numElemBlocks));
}

/*============================================================================*/
extern "C" int FEI_getNumBlockActNodes(CFEI* cfei, 
                                       GlobalID blockID,
                                       int* numNodes)
{
  CHK_CFEI_NULL(cfei);

   return(((FEI_Implementation*)(cfei->cfei_))->getNumBlockActNodes(blockID, *numNodes));
}

/*============================================================================*/
extern "C" int FEI_getNumBlockActEqns(CFEI* cfei, 
                                      GlobalID blockID,
                                      int* numEqns)
{
  CHK_CFEI_NULL(cfei);

   return(((FEI_Implementation*)(cfei->cfei_))->getNumBlockActEqns(blockID, *numEqns));
}
  
/*============================================================================*/
extern "C" int FEI_getNumNodesPerElement(CFEI* cfei, 
                                         GlobalID blockID,
                                         int* nodesPerElem)
{
  CHK_CFEI_NULL(cfei);

   return(((FEI_Implementation*)(cfei->cfei_))->getNumNodesPerElement(blockID, *nodesPerElem));
}

/*============================================================================*/
extern "C" int FEI_getNumEqnsPerElement(CFEI* cfei, 
                                        GlobalID blockID,
                                        int* numEqns)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->getNumEqnsPerElement(blockID,
								    *numEqns));
}
 
/*============================================================================*/
extern "C" int FEI_getNumBlockElements(CFEI* cfei,
                                       GlobalID blockID,
                                       int* numElems)
{
  CHK_CFEI_NULL(cfei);

  return(((FEI_Implementation*)(cfei->cfei_))->getNumBlockElements(blockID,
								   *numElems));
}

/*============================================================================*/
extern "C" int FEI_getNumBlockElemDOF(CFEI* cfei,
                                      GlobalID blockID,
                                      int* DOFPerElem)
{
  CHK_CFEI_NULL(cfei);

   return(((FEI_Implementation*)(cfei->cfei_))->getNumBlockElemDOF(blockID, *DOFPerElem));
}

