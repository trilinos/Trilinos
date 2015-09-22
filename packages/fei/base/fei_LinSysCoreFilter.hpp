#ifndef _LinSysCoreFilter_hpp_
#define _LinSysCoreFilter_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_macros.hpp"
#include "fei_defs.h"
#include "fei_fwd.hpp"
#include "fei_Filter.hpp"
#include <fei_CSRMat.hpp>
#include <fei_CSVec.hpp>

namespace fei {
  class DirichletBCManager;
}

/**
FEI_Implementation manages one or several instances of this class in the process
of assembling and solving a linear-system. Many of the public FEI function calls
are simply forwarded from FEI_Implementation to this class. This class then
performs the "filtering" process of turning nodeIDs and solution fields into
equation numbers and then dropping the data on through to the underlying solver
library by way of the LinearSystemCore
interface that is implemented for the solver library in question.
 */

class LinSysCoreFilter : public Filter {

 public:
   // Constructor.
   LinSysCoreFilter(FEI_Implementation* owner, MPI_Comm comm,
		    SNL_FEI_Structure* probStruct,
		    LinearSystemCore* lsc,
		    int masterRank=0);

   //Destructor
   virtual ~LinSysCoreFilter();


   // set a value (usually zeros) throughout the linear system
   virtual int resetSystem(double s);
   virtual int resetMatrix(double s);
   virtual int resetRHSVector(double s);
   virtual int resetInitialGuess(double s);

   virtual int deleteMultCRs();

   virtual int loadNodeBCs(int numNodes,
                   const GlobalID *nodeIDs,
                   int fieldID,
                   const int* offsetsIntoField,
                   const double* prescribedValues);

   virtual int loadElemBCs(int numElems,
                   const GlobalID *elemIDs,
                   int fieldID,
                   const double *const *alpha,  
                   const double *const *beta,
                   const double *const *gamma);

   virtual int sumInElem(GlobalID elemBlockID,
                 GlobalID elemID,
                 const GlobalID* elemConn,
                 const double* const* elemStiffness,
                 const double* elemLoad,
                 int elemFormat);

   virtual int sumInElemMatrix(GlobalID elemBlockID,
                       GlobalID elemID,
                       const GlobalID* elemConn,
                       const double* const* elemStiffness,
                       int elemFormat);

   virtual int sumInElemRHS(GlobalID elemBlockID,
                    GlobalID elemID,
                    const GlobalID* elemConn,
                    const double* elemLoad);

    virtual int loadCRMult(int CRMultID, 
                   int numCRNodes,
                   const GlobalID* CRNodes, 
                   const int* CRFields,
                   const double* CRWeights,
                   double CRValue);

    virtual int loadCRPen(int CRPenID, 
                  int numCRNodes, 
                  const GlobalID* CRNodes,
                  const int *CRFields,
                  const double* CRWeights,
                  double CRValue,
                  double penValue);

   virtual int putIntoRHS(int IDType,
		  int fieldID,
			  int numIDs,
		  const GlobalID* IDs,
		  const double* rhsEntries);

   virtual int sumIntoRHS(int IDType,
		  int fieldID,
			  int numIDs,
		  const GlobalID* IDs,
		  const double* rhsEntries);

   virtual int loadComplete();

    // set parameters associated with solver choice, etc.
    virtual int parameters(int numParams, const char *const* paramStrings);

    //get residual norms
    virtual int residualNorm(int whichNorm, int numFields,
                     int* fieldIDs, double* norms, double& residTime);

    // start iterative solution
    virtual int solve(int& status, double& sTime);

    // query function iterations performed.
    virtual int iterations() const {return(iterations_);};

// Solution return services.......................................
 
    // return all nodal solution params on a block-by-block basis 
    virtual int getBlockNodeSolution(GlobalID elemBlockID,  
                             int numNodes, 
                             const GlobalID *nodeIDs, 
                             int *offsets,
                             double *results);
 
    virtual int getNodalSolution(int numNodes, 
				 const GlobalID *nodeIDs, 
				 int *offsets,
				 double *results);

    // return nodal solution for one field on a block-by-block basis 
    virtual int getBlockFieldNodeSolution(GlobalID elemBlockID,
                                  int fieldID,
                                  int numNodes, 
                                  const GlobalID *nodeIDs, 
                                  double *results);
         
    // return element solution params on a block-by-block basis 
    virtual int getBlockElemSolution(GlobalID elemBlockID,  
                             int numElems, 
                             const GlobalID *elemIDs,
                             int& numElemDOFPerElement,
                             double *results);

   virtual int getCRMultipliers(int numCRs, const int* CRIDs, double* multipliers);

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
   virtual int putBlockNodeSolution(GlobalID elemBlockID,
                            int numNodes,
                            const GlobalID *nodeIDs, 
                            const int *offsets,
                            const double *estimates);

    // put nodal-based guess for one field on a block-by-block basis 
    virtual int putBlockFieldNodeSolution(GlobalID elemBlockID, 
                                  int fieldID, 
                                  int numNodes, 
                                  const GlobalID *nodeIDs, 
                                  const double *estimates);
         
    // put element-based solution guess on a block-by-block basis
    virtual int putBlockElemSolution(GlobalID elemBlockID,  
                             int numElems, 
                             const GlobalID *elemIDs, 
                             int dofPerElem,
                             const double *estimates);
  
    virtual int putCRMultipliers(int numMultCRs, 
                         const int* CRIDs,
                         const double *multEstimates);

//===== a couple of public non-FEI functions... ================================
//These are intended to be used by an 'outer-layer' class like 
//FEI_Implementation.
//
  public:
    virtual int getNodalFieldSolution(int fieldID,
			      int numNodes,
			      const GlobalID* nodeIDs,
			      double* results);

    virtual int putNodalFieldData(int fieldID,
			  int numNodes,
			  const GlobalID* nodeIDs,
			  const double* nodeData);

    virtual int putNodalFieldSolution(int fieldID,
			      int numNodes,
			      const GlobalID* nodeIDs,
			      const double* nodeData);

    virtual int unpackSolution();

    void setEqnCommMgr(EqnCommMgr* eqnCommMgr);

   EqnCommMgr* getEqnCommMgr() {return(eqnCommMgr_);};

   virtual int setNumRHSVectors(int numRHSs, int* rhsIDs);
   virtual int setCurrentRHS(int rhsID);

   virtual int exchangeRemoteEquations();
   virtual int exchangeRemoteBCs(std::vector<int>& essEqns,
                                 std::vector<double>& essAlpha,
                                 std::vector<double>& essGamma);

   virtual int implementAllBCs();

   virtual int enforceEssentialBCs(const int* eqns, const double* alpha,
                                  const double* gamma, int numEqns);

   virtual int enforceRemoteEssBCs(int numEqns, const int* eqns, 
			   const int* const* colIndices, const int* colIndLens,
			   const double* const* coefs);

   virtual int initialize();

//==============================================================================
//private functions for internal implementation of LinSysCoreFilter.
//==============================================================================
  private:
   int initLinSysCore();
   void setLinSysCoreCREqns();

   int unpackRemoteContributions(EqnCommMgr& eqnCommMgr,
				 int assemblyMode);

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

   int storeNodalColumnCoefs(int eqn, const NodeDescriptor& node,
			     int fieldID, int fieldSize,
			     double* coefs);

   int storeNodalRowCoefs(const NodeDescriptor& node,
			  int fieldID, int fieldSize,
			  double* coefs, int eqn);

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

   void storeNodalSendIndex(const NodeDescriptor& node, int fieldID, int col);
   void storeNodalSendEqn(const NodeDescriptor& node, int fieldID, int col,
			  double* coefs);
   void storeNodalSendIndices(const NodeDescriptor& iNode, int iField,
			      const NodeDescriptor& jNode, int jField);

   void storePenNodeSendData(const NodeDescriptor& iNode,
			     int iField, int iFieldSize,
			     double* iCoefs,
			     const NodeDescriptor& jNode,
			     int jField, int jFieldSize,
			     double* jCoefs,
			     double penValue, double CRValue);

   int storePenNodeData(const NodeDescriptor& iNode,
			int iField, int iFieldSize,
			double* iCoefs,
			const NodeDescriptor& jNode,
			int jField, int jFieldSize,
			double* jCoefs,
			double penValue, double CRValue);

   void allocElemStuff();

   int resolveConflictingCRs(EqnBuffer& bcEqns);

   int giveToMatrix_symm_noSlaves(int numPtRows,
				  const int* ptRowNumbers,
				  const double* const* coefs,
				  int mode);

   int giveToBlkMatrix_symm_noSlaves(int numPtRows, const int* ptRows,
				     int numBlkRows, const int* blkRowNumbers,
				     const int* blkRowSizes,
				     const double* const* coefs,
				     int mode);

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

   int sumIntoMatrix(fei::CSRMat& mat);

   int getEqnsFromMatrix(ProcEqns& procEqns, EqnBuffer& eqnData);

   int getEqnsFromRHS(ProcEqns& procEqns, EqnBuffer& eqnData);

   int giveToRHS(int num, const double* values,
		 const int* indices, int mode);

   int giveToLocalReducedRHS(int num, const double* values,
			     const int* indices, int mode);

   int getFromRHS(int num, double* values, const int* indices);

   int sumIntoRHS(fei::CSVec& vec);

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
		    int numBlkEqns, int* blkEqns, int* blkSizes,
		    bool useBlkEqns, int mode);

   int assembleReducedEqns();

   int assembleRHS(int numValues, const int* indices, const double* coefs, int mode);

   int assembleReducedRHS();
 
   void debugOutput(const char* mesg);

   int createEqnCommMgr_put();

//==============================================================================
//private LinSysCoreFilter variables
//==============================================================================
  private:

   int timesInitializeCalled_;

    LinearSystemCore* lsc_;
    bool useLookup_;

    int internalFei_;

    bool newMatrixData_;
    bool newVectorData_;
    bool newConstraintData_;
    bool newBCData_;
    bool connectivitiesInitialized_;
    bool firstRemEqnExchange_;
    bool needToCallMatrixLoadComplete_;

    bool resolveConflictRequested_;

    int localStartRow_, localEndRow_, numLocalEqns_, numGlobalEqns_;
    int reducedStartRow_, reducedEndRow_, numReducedRows_;
    int numLocallyOwnedNodes_, numGlobalNodes_, firstLocalNodeNumber_;

    bool blockMatrix_;
    bool tooLateToChooseBlock_;

    int numLocalEqnBlks_;
    int localReducedBlkOffset_, numLocalReducedEqnBlks_;

    int iterations_;
    int numRHSs_;
    int currentRHS_;
    std::vector<int> rhsIDs_;

    int outputLevel_;

    MPI_Comm comm_;
    int masterRank_;

    SNL_FEI_Structure* problemStructure_;
    bool matrixAllocated_;

    std::vector<int> rowIndices_;
    std::vector<int> rowColOffsets_, colIndices_;

    fei::FillableMat *Kid_, *Kdi_, *Kdd_;
    fei::CSRMat csrD, csrKid, csrKdi, csrKdd, tmpMat1_, tmpMat2_;
    fei::CSVec fd_, tmpVec1_;
    int reducedEqnCounter_, reducedRHSCounter_;
    std::vector<int> rSlave_, cSlave_;

    int nodeIDType_;
    fei::DirichletBCManager* bcManager_; //Boundary condition manager

    EqnCommMgr* eqnCommMgr_; //equation communication manager
    EqnCommMgr* eqnCommMgr_put_; //only created if users call the 'put'
                                 // functions

    int maxElemRows_;
    std::vector<int> scatterIndices_;
    std::vector<int> blkScatterIndices_;
    std::vector<int> iworkSpace_, iworkSpace2_;
    std::vector<double> dworkSpace_;
    std::vector<const double*> dworkSpace2_;

    double** eStiff_;
    double* eStiff1D_;
    double* eLoad_;
};

#endif

