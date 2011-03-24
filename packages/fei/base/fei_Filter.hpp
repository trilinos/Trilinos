#ifndef _Filter_hpp_
#define _Filter_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_fwd.hpp>
#include <fei_defs.h>
#include <fei_macros.hpp>
#include <fei_iostream.hpp>

/**
FEI_Implementation manages one or several instances of this class in the process
of assembling and solving a linear-system. Many of the public FEI function calls
are simply forwarded from FEI_Implementation to this class. This class then
performs the "filtering" process of turning nodeIDs and solution fields into
equation numbers and then dropping the data on through to the underlying solver
library by way of a LinearSystemCore or FiniteElementData
interface that is implemented for the solver library in question.
 */
 
class Filter {

 public:
  /** Constructor */
  Filter(SNL_FEI_Structure* probStruct);

  /** Destructor */
  virtual ~Filter();

  virtual int initialize() = 0;

   // set a value (usually zeros) throughout the linear system
   virtual int resetSystem(double s) = 0;
   virtual int resetMatrix(double s) = 0;
   virtual int resetRHSVector(double s) = 0;
   virtual int resetInitialGuess(double s) = 0;

   virtual int deleteMultCRs() = 0;

   virtual int loadNodeBCs(int /*numNodes*/,
                   const GlobalID* /*nodeIDs*/,
                   int /*fieldID*/,
                   const int* /*offsetsIntoField*/,
                   const double* /*prescribedValues*/)
   {
      fei::console_out() << "fei ERROR, Filter::loadNodeBCs not overridden."<<FEI_ENDL;
      return -1;
   }

   virtual int loadElemBCs(int numElems,
                   const GlobalID *elemIDs,
                   int fieldID,
                   const double *const *alpha,  
                   const double *const *beta,  
                   const double *const *gamma) = 0;

   virtual int sumInElem(GlobalID /*elemBlockID*/,
                         GlobalID /*elemID*/,
                         const GlobalID* /*elemConn*/,
                         const double* const* /*elemStiffness*/,
                         const double* /*elemLoad*/,
                           int /*elemFormat*/) { return(0); }

   virtual int sumInElemMatrix(GlobalID /*elemBlockID*/,
                               GlobalID /*elemID*/,
                               const GlobalID* /*elemConn*/,
                               const double* const* /*elemStiffness*/,
                               int /*elemFormat*/) { return(0); }

   virtual int sumInElemRHS(GlobalID /*elemBlockID*/,
                            GlobalID /*elemID*/,
                            const GlobalID* /*elemConn*/,
                            const double* /*elemLoad*/) { return(0); }

    virtual int loadCRMult(int CRMultID, 
                   int numCRNodes,
                   const GlobalID* CRNodes, 
                   const int* CRFields,
                   const double* CRWeights,
                   double CRValue) = 0;

    virtual int loadCRPen(int CRPenID, 
                  int numCRNodes, 
                  const GlobalID* CRNodes,
                  const int *CRFields,
                  const double* CRWeights,
                  double CRValue,
                  double penValue) = 0;

   virtual int putIntoRHS(int IDType,
                  int fieldID,
                          int numIDs,
                  const GlobalID* IDs,
                  const double* rhsEntries) = 0;

   virtual int sumIntoRHS(int IDType,
                  int fieldID,
                          int numIDs,
                  const GlobalID* IDs,
                  const double* rhsEntries) = 0;

   virtual int sumIntoMatrixDiagonal(int /* IDType*/,
                                     int /* fieldID*/,
                                     int /* numIDs*/,
                                     const GlobalID* /* IDs*/,
                                     const double* /* coefficients*/)
   { return -1; }

   virtual int loadComplete() = 0;

    // set parameters associated with solver choice, etc.
    virtual int parameters(int numParams, const char *const* paramStrings);

    //get residual norms
    virtual int residualNorm(int whichNorm, int numFields,
                     int* fieldIDs, double* norms, double& residTime) = 0;

    // start iterative solution
    virtual int solve(int& status, double& sTime) = 0;

    // query function iterations performed.
    virtual int iterations() const = 0;

// Solution return services.......................................
 
    // return all nodal solution params on a block-by-block basis 
    virtual int getBlockNodeSolution(GlobalID elemBlockID,  
                             int numNodes, 
                             const GlobalID *nodeIDs, 
                             int *offsets,
                             double *results) = 0;
 
    virtual int getNodalSolution(int numNodes, 
                                 const GlobalID *nodeIDs, 
                                 int *offsets,
                                 double *results) = 0;

    // return nodal solution for one field on a block-by-block basis 
    virtual int getBlockFieldNodeSolution(GlobalID elemBlockID,
                                  int fieldID,
                                  int numNodes, 
                                  const GlobalID *nodeIDs, 
                                  double *results) = 0;
         
    // return element solution params on a block-by-block basis 
    virtual int getBlockElemSolution(GlobalID elemBlockID,  
                             int numElems, 
                             const GlobalID *elemIDs,
                             int& numElemDOFPerElement,
                             double *results) = 0;

   virtual int getCRMultipliers(int numCRs, const int* CRIDs,
                                double* multipliers) = 0;

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
                            const double *estimates) = 0;

    // put nodal-based guess for one field on a block-by-block basis 
    virtual int putBlockFieldNodeSolution(GlobalID elemBlockID, 
                                  int fieldID, 
                                  int numNodes, 
                                  const GlobalID *nodeIDs, 
                                  const double *estimates) = 0;
  
    virtual int putBlockElemSolution(GlobalID elemBlockID,  
                                     int numElems, 
                                     const GlobalID *elemIDs, 
                                     int dofPerElem,
                                     const double *estimates) = 0;
  
    virtual int putCRMultipliers(int numMultCRs, 
                         const int* CRIDs,
                         const double *multEstimates) = 0;

//===== a couple of public non-FEI functions... ================================
//These are intended to be used by an 'outer-layer' class like 
//FEI_Implementation.
//
  public:

    virtual int getNodalFieldSolution(int fieldID,
                              int numNodes,
                              const GlobalID* nodeIDs,
                              double* results) = 0;

    virtual int putNodalFieldData(int fieldID,
                          int numNodes,
                          const GlobalID* nodeIDs,
                          const double* nodeData) = 0;

    virtual int putNodalFieldSolution(int fieldID,
                              int numNodes,
                              const GlobalID* nodeIDs,
                              const double* nodeData) = 0;

    virtual int unpackSolution() = 0;

    virtual void setEqnCommMgr(EqnCommMgr* eqnCommMgr) = 0;

   virtual EqnCommMgr* getEqnCommMgr() = 0;

   virtual int setNumRHSVectors(int numRHSs, int* rhsIDs) = 0;
   virtual int setCurrentRHS(int rhsID) = 0;

   virtual int exchangeRemoteEquations() { return 0; }

   virtual int enforceEssentialBCs(const int* eqns, const double* alpha,
                                  const double* gamma, int numEqns) = 0;

   static void copyStiffness(const double* const* elemStiff, int numRows,
                             int elemFormat, double** copy);

   void setLogStream(std::ostream* logstrm);
   std::ostream* logStream();

 protected:
   virtual int generalElemInput(GlobalID /*elemBlockID*/,
                                GlobalID /*elemID*/,
                                const GlobalID* /*elemConn*/,
                                const double* const* /*elemStiffness*/,
                                const double* /*elemLoad*/,
                                  int /*elemFormat*/) {return(-1);}

   int generalCoefInput(int /*patternID*/,
                        const int* /*rowIDTypes*/,
                        const GlobalID* /*rowIDs*/,
                        const int* /*colIDTypes*/,
                        const GlobalID* /*colIDs*/,
                        const double* const* /*matrixEntries*/,
                        const double* /*vectorEntries*/,
                          int /*numRows*/, int /*numCols*/) {return(-1);}

   int calculateResidualNorms(int whichNorm, int numFields,
                              int* fieldIDs, double* norms,
                              std::vector<double>& residValues);

   const NodeDescriptor* findNode(GlobalID nodeID) const;
   const NodeDescriptor& findNodeDescriptor(GlobalID nodeID) const;

   SNL_FEI_Structure* problemStructure_;

   bool logInput_;
   std::ostream* logInputStream_;

   int outputLevel_;

   int numProcs_;
   int localRank_;

 private:
   Filter(const Filter& /*src*/)
     : problemStructure_(NULL), logInput_(false), logInputStream_(NULL),
     outputLevel_(0), numProcs_(0), localRank_(0)
     {}

   Filter& operator=(const Filter& /*src*/)
     {
       return(*this);
     }
};

#endif

