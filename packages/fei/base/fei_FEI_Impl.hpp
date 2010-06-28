/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_FEI_Impl_hpp_
#define _fei_FEI_Impl_hpp_

#include <fei_macros.hpp>
#include <fei_mpi.h>
#include <FEI.hpp>
#include <fei_Factory.hpp>

#include <set>
#include <vector>
#include <map>

namespace fei {

  /** Implementation of the original abstract FEI interface. This FEI_Impl
      class plays the same role as the FEI_Implementation class, except that
      this class internally is based on the newer more modular fei:: classes.
      The result should be the same functionality as FEI_Implementation but
      with improved performance.

      The ultimate goal is for this class to replace the FEI_Implementation
      class, thereby greatly reducing the total body of code that constitutes
      the fei code-distribution.
  */
  class FEI_Impl : public FEI, private fei::Logger {
  public:
    /** Constructor
     */
    FEI_Impl(fei::SharedPtr<LibraryWrapper> wrapper,
          MPI_Comm comm,
          int masterRank=0);

    /** Constructor
     */
    FEI_Impl(const fei::Factory* factory,
          MPI_Comm comm,
          int masterRank=0);

    /** Destructor
     */
    virtual ~FEI_Impl();

   /** This function isn't in the original FEI 2.1 interface. But it is useful
      to users who may be mixing old and new fei objects and need access
      to underlying stuff.
   */
   fei::SharedPtr<fei::LinearSystem> getLinearSystem();

    int parameters(int numParams, 
                   const char *const* paramStrings);

   int setIDLists(int numMatrices,
                  const int* matrixIDs,
                  int numRHSs,
                  const int* rhsIDs);

   int setSolveType(int solveType);

   int initFields(int numFields, 
                  const int *fieldSizes, 
                  const int *fieldIDs,
                  const int *fieldTypes=NULL);

   int initElemBlock(GlobalID elemBlockID,
                     int numElements,
                     int numNodesPerElement,
                     const int* numFieldsPerNode,
                     const int* const* nodalFieldIDs,
                     int numElemDofFieldsPerElement,
                     const int* elemDOFFieldIDs,
                     int interleaveStrategy);

   int initElem(GlobalID elemBlockID,
                GlobalID elemID,
                const GlobalID* elemConn);

   int initSlaveVariable(GlobalID slaveNodeID, 
                         int slaveFieldID,
                         int offsetIntoSlaveField,
                         int numMasterNodes,
                         const GlobalID* masterNodeIDs,
                         const int* masterFieldIDs,
                         const double* weights,
                         double rhsValue);

   int deleteMultCRs();

   int initSharedNodes(int numSharedNodes,
                       const GlobalID *sharedNodeIDs,  
                       const int* numProcsPerNode, 
                       const int *const *sharingProcIDs);

   int initCRMult(int numCRNodes,
                  const GlobalID* CRNodeIDs,
                  const int *CRFieldIDs,
                  int& CRID); 

   int initCRPen(int numCRNodes,
                 const GlobalID* CRNodeIDs, 
                 const int *CRFieldIDs,
                 int& CRID); 

   /** indicate that overall initialization sequence is complete */
   int initComplete();

   int setCurrentMatrix(int matrixID);

   int setCurrentRHS(int rhsID);

   int resetSystem(double s=0.0);

   int resetMatrix(double s=0.0);

   int resetRHSVector(double s=0.0);

  int resetInitialGuess(double s=0.0);

    int loadNodeBCs(int numNodes,
                    const GlobalID *nodeIDs,
                    int fieldID,
                    const int* offsetsIntoField,
                    const double* prescribedValues);

    int loadElemBCs(int numElems,
                    const GlobalID* elemIDs,  
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
                  const GlobalID* CRNodeIDs,
                  const int* CRFieldIDs,
                  const double* CRWeights,
                  double CRValue);

   int loadCRPen(int CRPenID,
                 int numCRNodes,
                 const GlobalID* CRNodeIDs,
                 const int* CRFieldIDs,
                 const double* CRWeights,
                 double CRValue,
                 double penValue);

   /**
   */
   int putIntoRHS(int IDType,
                  int fieldID,
                  int numIDs,
                  const GlobalID* IDs,
                  const double* coefficients);

   /**
   */
   int sumIntoRHS(int IDType,
                  int fieldID,
                  int numIDs,
                  const GlobalID* IDs,
                  const double* coefficients);

   int setMatScalars(int numScalars,
                     const int* IDs, 
                     const double* scalars);

   int setRHSScalars(int numScalars,
                     const int* IDs,
                     const double* scalars);

   int loadComplete(bool applyBCs=true,
                    bool globalAssemble=true);

   /** get residual norms */
   int residualNorm(int whichNorm,
                    int numFields,
                    int* fieldIDs,
                    double* norms);

   /** launch underlying solver */
   int solve(int& status);

   /** Query number of iterations taken for last solve.
      @param itersTaken Iterations performed during any previous solve.
   */
   int iterations(int& itersTaken) const;

   int version(const char*& versionString);

   /** query for some accumulated timing information. Collective function. */
   int cumulative_cpu_times(double& initTime,
                            double& loadTime,
                            double& solveTime,
                            double& solnReturnTime);

   /** return all nodal solution params on a block-by-block basis */
    int getBlockNodeSolution(GlobalID elemBlockID,  
                             int numNodes, 
                             const GlobalID *nodeIDs, 
                             int *offsets,
                             double *results);

    /** return all nodal solution params for an arbitrary list of nodes */
    int getNodalSolution(int numNodes,
                         const GlobalID* nodeIDs,
                         int* offsets,
                         double* results);

    /** return nodal solution for one field on a block-by-block basis */
    int getBlockFieldNodeSolution(GlobalID elemBlockID,
                                  int fieldID,
                                  int numNodes, 
                                  const GlobalID *nodeIDs, 
                                  double *results);
         
    /** return element solution params on a block-by-block basis */
    int getBlockElemSolution(GlobalID elemBlockID,  
                             int numElems, 
                             const GlobalID *elemIDs,
                             int& numElemDOFPerElement,
                             double *results);

    /** Query number of lagrange-multiplier constraints on local processor */
   int getNumCRMultipliers(int& numMultCRs);

   /** Obtain list of lagrange-multiplier IDs */
   int getCRMultIDList(int numMultCRs, int* multIDs);

   /** get Lagrange Multipliers */
   int getCRMultipliers(int numCRs,
                        const int* CRIDs,
                        double *multipliers);

// associated "puts" paralleling the solution return services.
// 
// the int sizing parameters are passed for error-checking purposes, so
// that the interface implementation can tell if the passed estimate
// vectors make sense -before- an attempt is made to utilize them as
// initial guesses by unpacking them into the solver's native solution
// vector format (these parameters include lenNodeIDList, lenElemIDList,
// numElemDOF, and numMultCRs -- all other passed params are either 
// vectors or block/constraint-set IDs)

   /** put nodal-based solution guess on a block-by-block basis */
    int putBlockNodeSolution(GlobalID elemBlockID, 
                             int numNodes, 
                             const GlobalID *nodeIDs, 
                             const int *offsets,
                             const double *estimates);

    /** put nodal-based guess for one field on a block-by-block basis */
    int putBlockFieldNodeSolution(GlobalID elemBlockID, 
                                  int fieldID, 
                                  int numNodes, 
                                  const GlobalID *nodeIDs, 
                                  const double *estimates);
         
    /** put element-based solution guess on a block-by-block basis */
    int putBlockElemSolution(GlobalID elemBlockID,  
                             int numElems, 
                             const GlobalID *elemIDs, 
                             int dofPerElem,
                             const double *estimates);

    /** put Lagrange solution to FE analysis on a constraint-set basis */
    int putCRMultipliers(int numMultCRs, 
                         const int* CRIDs,
                         const double* multEstimates);

// utility functions that aid in integrating the FEI calls..............

// support methods for the "gets" and "puts" of the soln services.


    /** return info associated with blocked nodal solution */
    int getBlockNodeIDList(GlobalID elemBlockID,
                           int numNodes,
                           GlobalID *nodeIDs);

    /** return info associated with blocked element solution */
   int getBlockElemIDList(GlobalID elemBlockID, 
                          int numElems, 
                          GlobalID* elemIDs);
 
// miscellaneous self-explanatory "read-only" query functions............ 
 
   /** Query number of degrees-of-freedom for specified node */
    int getNumSolnParams(GlobalID nodeID, int& numSolnParams) const;

    /** Query number of element-blocks on local processor */
    int getNumElemBlocks(int& numElemBlocks) const;

    /**  return the number of active nodes in a given element block */
    int getNumBlockActNodes(GlobalID blockID, int& numNodes) const;

    /**  return the number of active equations in a given element block */
    int getNumBlockActEqns(GlobalID blockID, int& numEqns) const;

    /**  return the number of nodes associated with elements of a
         given block ID */
    int getNumNodesPerElement(GlobalID blockID, int& nodesPerElem) const;
    
    /**  return the number of equations (including element eqns)
         associated with elements of a given block ID */
    int getNumEqnsPerElement(GlobalID blockID, int& numEqns) const;

    /**  return the number of elements associated with this blockID */
    int getNumBlockElements(GlobalID blockID, int& numElems) const;

    /**  return the number of elements eqns for elems w/ this blockID */
    int getNumBlockElemDOF(GlobalID blockID, int& DOFPerElem) const;


    /** return the parameters that have been set so far. The caller should
     NOT delete the paramStrings pointer.
    */
    int getParameters(int& numParams, char**& paramStrings);

    //And now a couple of non-FEI query functions that Sandia applications
    //need to augment the matrix-access functions. I (Alan Williams) will
    //argue to have these included in the FEI 2.1 specification update.

    /** Query the size of a field. This info is supplied to the FEI (initFields)
        by the application, but may not be easily obtainable on the app side at
        all times. Thus, it would be nice if the FEI could answer this query.
    */
    int getFieldSize(int fieldID, int& numScalars);

    /**Since the ultimate intent for matrix-access is to bypass the FEI and go
     straight to data objects, we need a translation function to map between
     the IDs that the FEI deals in, and equation numbers that linear algebra
     objects deal in.
     @param ID Identifier of either a node or an element.
     @param idType Can take either of the values FEI_NODE or FEI_ELEMENT.
     @param fieldID Identifies a particular field at this [node||element].
     @param numEqns Output. Number of equations associated with this
     node/field (or element/field) pair.
     @param eqnNumbers Caller-allocated array. On exit, this is filled with the
     above-described equation-numbers. They are global 0-based numbers.
    */
    int getEqnNumbers(GlobalID ID,
                      int idType, 
                      int fieldID,
                      int& numEqns,
                      int* eqnNumbers);

    /**Get the solution data for a particular field, on an arbitrary set of
       nodes.
       @param fieldID Input. field identifier for which solution data is being
       requested.
       @param numNodes Input. Length of the nodeIDs list.
       @param nodeIDs Input. List specifying the nodes on which solution
       data is being requested.
       @param results Allocated by caller, but contents are output.
       Solution data for the i-th node/element starts in position i*fieldSize,
       where fieldSize is the number of scalar components that make up 
       'fieldID'.
       @return error-code 0 if successful
    */
    int getNodalFieldSolution(int fieldID,
                              int numNodes,
                              const GlobalID* nodeIDs,
                              double* results);

   /**Get the number of nodes that are local to this processor (includes nodes
      that are shared by other processors).
      @param numNodes Output. Number of local nodes.
      @return error-code 0 if successful
   */
    int getNumLocalNodes(int& numNodes);

   /**Get a list of the nodeIDs that are local to this processor (includes nodes
      that are shared by other processors).
      @param numNodes Output. Same as the value output by 'getNumLocalNodes'.
      @param nodeIDs Caller-allocated array, contents to be filled by this
      function.
      @param lenNodeIDs Input. Length of the caller-allocated nodeIDs array. If
      lenNodeIDs is less than numNodes, then only 'lenNodeIDs' nodeIDs are
      provided, of course. If lenNodeIDs is greater than numNodes, then only
      'numNodes' positions of the nodeIDs array are referenced.
      @return error-code 0 if successful
   */
    int getLocalNodeIDList(int& numNodes,
                           GlobalID* nodeIDs,
                           int lenNodeIDs);

    /** Pass nodal data for a specified field through to the solver. Example
      is geometric coordinates, etc.

    @param fieldID field identifier for the data to be passed. This is probably
       not one of the 'solution-field' identifiers. It should be an identifier
      that the solver is expecting and knows how to handle. This field id must
      previously have been provided to the fei implementation via the normal
      initFields method.

    @param numNodes number of nodes for which data is being provided
    @param nodeIDs List of length numNodes, giving node-identifiers for which
       data is being provided.
    @param nodeData List of length fieldSize*numNodes, where fieldSize is the
      size which was associated with fieldID when initFields was called. The
      data for nodeIDs[i] begins in position i*fieldSize of this array.
    */
    int putNodalFieldData(int fieldID,
                          int numNodes,
                          const GlobalID* nodeIDs,
                          const double* nodeData);

  private: //methods
    void basic_initializations();

    int inputRHS(int IDType,
                 int fieldID,
                 int numIDs,
                 const GlobalID* IDs,
                 const double* rhsEntries,
                 bool sumInto);

    int fillNodeset(int blockID) const;

    int aggregateSystem();

    int aggregateSystem_LinSysCore();

  private: //attributes

    std::vector<fei::SharedPtr<LibraryWrapper> > wrapper_;

    int nodeIDType_;
    int elemIDType_;
    int constraintIDType_;

    std::vector<fei::SharedPtr<fei::Factory> > factory_;
    bool createdFactory_;
    fei::SharedPtr<fei::VectorSpace> rowSpace_;
    fei::SharedPtr<fei::MatrixGraph> matGraph_;

    fei::SharedPtr<fei::Vector> x_;
    std::vector<fei::SharedPtr<fei::Vector> > b_;
    std::vector<fei::SharedPtr<fei::Matrix> > A_;
    fei::SharedPtr<fei::LinearSystem> linSys_;

    bool newData_;

    Data *soln_fei_matrix_;
    Data *soln_fei_vector_;

    MPI_Comm comm_;
    int masterRank_;
    int localProc_;
    int numProcs_;

    int numParams_;
    char** paramStrings_;

    std::vector<int> matrixIDs_;
    std::vector<int> rhsIDs_;

    std::vector<double> matScalars_;
    bool matScalarsSet_;
    std::vector<double> rhsScalars_;
    bool rhsScalarsSet_;

    int constraintID_;

    int index_soln_;
    int index_current_;
    int index_current_rhs_row_;

    int solveType_;
    int iterations_;

    bool setSolveTypeCalled_;
    bool initPhaseIsComplete_;

    bool aggregateSystemFormed_;
    int newMatrixDataLoaded_;

    int solveCounter_;

    double initTime_, loadTime_, solveTime_, solnReturnTime_;

    std::vector<int> iwork_;

    mutable std::set<int> nodeset_;
    mutable bool nodeset_filled_;
    mutable int nodeset_blockid_;

    mutable std::map<int,int> block_dof_per_elem_;
    bool any_blocks_have_elem_dof_;
  };//class FEI_Impl

}//namespace fei

#endif

