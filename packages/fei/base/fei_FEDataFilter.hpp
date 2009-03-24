#ifndef _fei_FEDataFilter_hpp_
#define _fei_FEDataFilter_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_fwd.hpp"
#include "fei_defs.h"
#include "fei_Filter.hpp"

/**
FEI_Implementation manages one or several instances of the Filter class in the
process of assembling and solving a linear-system. Many of the public FEI
function calls are simply forwarded from FEI_Implementation to this class.
This class then performs the "filtering" process of turning nodeIDs and solution
fields into equation numbers and then dropping the data through to the
underlying solver library by way of the FiniteElementData interface that is
implemented for the solver library in question.
 */
 
class FEDataFilter : public Filter {

 public:
   // Constructor.
   FEDataFilter(FEI_Implementation* owner, MPI_Comm comm,
                SNL_FEI_Structure* probStruct,
                LibraryWrapper* wrapper,
                int masterRank=0);

   //Destructor
   virtual ~FEDataFilter();


   // set a value (usually zeros) throughout the linear system
   int resetSystem(double s);
   int resetMatrix(double s);
   int resetRHSVector(double s);
   int resetInitialGuess(double s);

   int deleteMultCRs();

   int loadNodeBCs(int numNodes,
                   const GlobalID *nodeIDs,
                   int fieldID,
                   const int* offsetsIntoField,
                   const double* prescribedValues);

   int loadElemBCs(int numElems,
                   const GlobalID *elemIDs,
                   int fieldID,
                   const double *const *alpha,  
                   const double *const *beta,  
                   const double *const *gamma);

   int sumInElem(GlobalID elemBlockID,
                 GlobalID elemID,
                 const GlobalID* elemConn,
                 const double* const* elemStiffness,
                 const double* elemLoad,
                 int elemFormat);

   int sumInElemMatrix(GlobalID elemBlockID,
                       GlobalID elemID,
                       const GlobalID* elemConn,
                       const double* const* elemStiffness,
                       int elemFormat);

   int sumInElemRHS(GlobalID elemBlockID,
                    GlobalID elemID,
                    const GlobalID* elemConn,
                    const double* elemLoad);

    int loadCRMult(int CRMultID, 
                   int numCRNodes,
                   const GlobalID* CRNodes, 
                   const int* CRFields,
                   const double* CRWeights,
                   double CRValue);

    int loadCRPen(int CRPenID, 
                  int numCRNodes, 
                  const GlobalID* CRNodes,
                  const int *CRFields,
                  const double* CRWeights,
                  double CRValue,
                  double penValue);

   int putIntoRHS(int IDType,
                  int fieldID,
                          int numIDs,
                  const GlobalID* IDs,
                  const double* rhsEntries);

   int sumIntoRHS(int IDType,
                  int fieldID,
                          int numIDs,
                  const GlobalID* IDs,
                  const double* rhsEntries);

   int sumIntoMatrixDiagonal(int  IDType,
                             int  fieldID,
                             int  numIDs,
                             const GlobalID*  IDs,
                             const double*  coefficients);

   int loadComplete();

    // set parameters associated with solver choice, etc.
    int parameters(int numParams, const char *const* paramStrings);

    //get residual norms
    int residualNorm(int whichNorm, int numFields,
                     int* fieldIDs, double* norms, double& residTime);

    // start iterative solution
    int solve(int& status, double& sTime);

    // query function iterations performed.
    int iterations() const {return(iterations_);};

// Solution return services.......................................
 
    // return all nodal solution params on a block-by-block basis 
    int getBlockNodeSolution(GlobalID elemBlockID,  
                             int numNodes, 
                             const GlobalID *nodeIDs, 
                             int *offsets,
                             double *results);
 
    int getNodalSolution(int numNodes, 
                         const GlobalID *nodeIDs, 
                         int *offsets,
                         double *results);

    // return nodal solution for one field on a block-by-block basis 
    int getBlockFieldNodeSolution(GlobalID elemBlockID,
                                  int fieldID,
                                  int numNodes, 
                                  const GlobalID *nodeIDs, 
                                  double *results);
         
    // return element solution params on a block-by-block basis 
    int getBlockElemSolution(GlobalID elemBlockID,  
                             int numElems, 
                             const GlobalID *elemIDs,
                             int& numElemDOFPerElement,
                             double *results);

   int getCRMultipliers(int numCRs, const int* CRIDs, double* multipliers);

// associated "puts" paralleling the solution return services.
// 
// the int sizing parameters are passed for error-checking purposes, so
// that the interface implementation can tell if the passed estimate
// vectors make sense -before- an attempt is made to utilize them as
// initial guesses by unpacking them into the solver's native solution
// vector format (these parameters include lenNodeIDList, lenElemIDList,
// numElemDOF, and numMultCRs -- all other passed params are either 
// vectors or block/constraint-set IDs)

   // put nodal-based solution guess on a block-by-block basis 
   int putBlockNodeSolution(GlobalID elemBlockID,
                            int numNodes,
                            const GlobalID *nodeIDs, 
                            const int *offsets,
                            const double *estimates);

    // put nodal-based guess for one field on a block-by-block basis 
    int putBlockFieldNodeSolution(GlobalID elemBlockID, 
                                  int fieldID, 
                                  int numNodes, 
                                  const GlobalID *nodeIDs, 
                                  const double *estimates);
         
    // put element-based solution guess on a block-by-block basis  
    int putBlockElemSolution(GlobalID elemBlockID,  
                             int numElems, 
                             const GlobalID *elemIDs, 
                             int dofPerElem,
                             const double *estimates);
  
    int putCRMultipliers(int numMultCRs, 
                         const int* CRIDs,
                         const double *multEstimates);

//===== a couple of public non-FEI functions... ================================
//These are intended to be used by an 'outer-layer' class like 
//FEI_Implementation.
//
  public:
    int getNodalFieldSolution(int fieldID,
                              int numNodes,
                              const GlobalID* nodeIDs,
                              double* results);

    int putNodalFieldData(int fieldID,
                          int numNodes,
                          const GlobalID* nodeIDs,
                          const double* nodeData);

    int putNodalFieldSolution(int fieldID,
                              int numNodes,
                              const GlobalID* nodeIDs,
                              const double* nodeData);

    int unpackSolution();

    void setEqnCommMgr(EqnCommMgr* eqnCommMgr);

   EqnCommMgr* getEqnCommMgr() {return(eqnCommMgr_);};

   int setNumRHSVectors(int numRHSs, int* rhsIDs);
   int setCurrentRHS(int rhsID);

   int enforceEssentialBCs(const int* eqns, const double* alpha,
                                  const double* gamma, int numEqns);

   int initialize();

//==============================================================================
//private functions for internal implementation of FEDataFilter.
//==============================================================================
  private:
   FEDataFilter(const FEDataFilter& src);
   FEDataFilter& operator=(const FEDataFilter& src);

   int initLinSysCore();

   int loadFEDataMultCR(int CRID,
                        int numCRNodes,
                        const GlobalID* CRNodes, 
                        const int* CRFields,
                        const double* CRWeights,
                        double CRValue);

   int loadFEDataPenCR(int CRID,
                       int numCRNodes,
                       const GlobalID* CRNodes, 
                       const int* CRFields,
                       const double* CRWeights,
                       double CRValue,
                       double penValue);

   int generalElemInput(GlobalID elemBlockID,
                        GlobalID elemID,
                        const double* const* elemStiffness,
                        const double* elemLoad,
                        int elemFormat);

   int generalElemInput(GlobalID elemBlockID,
                        GlobalID elemID,
                        const GlobalID* elemConn,
                        const double* const* elemStiffness,
                        const double* elemLoad,
                        int elemFormat);

   void allocElemStuff();

   int giveToMatrix(int numPtRows, const int* ptRows,
                    int numPtCols, const int* ptCols,
                    const double* const* values,
                    int mode);

   int giveToLocalReducedMatrix(int numPtRows, const int* ptRows,
                                int numPtCols, const int* ptCols,
                                const double* const* values,
                                int mode);

   int getFromMatrix(int numPtRows, const int* ptRows,
                     const int* rowColOffsets, const int* ptCols,
                     int numColsPerRow, double** values);

   int getEqnsFromMatrix(ProcEqns& procEqns, EqnBuffer& eqnData);

   int getEqnsFromRHS(ProcEqns& procEqns, EqnBuffer& eqnData);

   int giveToRHS(int num, const double* values,
                 const int* indices, int mode);

   int giveToLocalReducedRHS(int num, const double* values,
                             const int* indices, int mode);

   int getFromRHS(int num, double* values, const int* indices);

   int getEqnSolnEntry(int eqnNumber, double& solnValue);

   int getSharedRemoteSolnEntry(int eqnNumber, double& solnValue);

   int getReducedSolnEntry(int eqnNumber, double& solnValue);

   int formResidual(double* residValues, int numLocalEqns);

   int getRemoteSharedEqns(int numPtRows, const int* ptRows,
                           ProcEqns& remoteProcEqns);

   int resetTheMatrix(double s);
   int resetTheRHSVector(double s);

   int assembleEqns(int numPtRows, 
                    int numPtCols,
                    const int* rowNumbers,
                    const int* colIndices,
                    const double* const* coefs,
                    bool structurallySymmetric,
                    int mode);

   int assembleRHS(int numValues, const int* indices, const double* coefs, int mode);

   void debugOutput(const char* mesg);

   int createEqnCommMgr_put();

//==============================================================================
//private FEDataFilter variables
//==============================================================================
  private:

    LibraryWrapper* wrapper_;
    fei::SharedPtr<FiniteElementData> feData_;
    bool useLookup_;

    int internalFei_;

    bool newData_;

    int localStartRow_, localEndRow_, numGlobalEqns_;
    int reducedStartRow_, reducedEndRow_, numReducedRows_;

    int iterations_;
    int numRHSs_;
    int currentRHS_;
    std::vector<int> rhsIDs_;

    int outputLevel_;

    MPI_Comm comm_;
    int masterRank_;

    SNL_FEI_Structure* problemStructure_;

    std::vector<GlobalID> penCRIDs_;

    std::vector<int> rowIndices_;
    std::vector<int> rowColOffsets_, colIndices_;

    EqnCommMgr* eqnCommMgr_; //equation communication manager
    EqnCommMgr* eqnCommMgr_put_; //only created if users call
                                 // the 'put' functions

    int maxElemRows_;

    double** eStiff_;
    double* eStiff1D_;
    double* eLoad_;

    int numRegularElems_;
    std::vector<int> constraintBlocks_;
    std::vector<int> constraintNodeOffsets_;
    std::vector<int> packedFieldSizes_;
};

#endif

