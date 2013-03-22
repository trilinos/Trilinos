#ifndef _SNL_FEI_Structure_hpp_
#define _SNL_FEI_Structure_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_fwd.hpp>
#include <fei_defs.h>
#include <fei_constants.hpp>
#include "fei_TemplateUtils.hpp"
#include <snl_fei_PointBlockMap.hpp>
#include <fei_EqnBuffer.hpp>
#include <fei_FieldDofMap.hpp>
#include <fei_CSRMat.hpp>
#include <fei_CSVec.hpp>

#include <fei_NodeCommMgr.hpp>
#include <fei_NodeDatabase.hpp>
#include <fei_NodeDescriptor.hpp>

#include <fei_Lookup.hpp>

/** The SNL_FEI_Structure class is a container for the data that makes up
  the structure of a finite element problem. It contains things like the
  FieldDatabase, NodeDatabase, ConnectivityTable, etc.<p> This class is the
  primary class that implements the general "FEI initialization phase", where
  the user passes all structure-defining data, which ends up being translated
  into a corresponding sparse matrix structure.

  The problem structure is initialized gradually by calling
  the 'init' functions in this class. They also exist in the FEI class, where
  they are documented.

  Once all structure-defining data has been passed, the process is completed
  by calling the initComplete function.
*/

class SNL_FEI_Structure : public Lookup {
 public:
  /** Constructor.
      @param comm MPI_Communicator
      @param name String by which this structure instance may be referred to.
      @param debugOutputLevel Determines whether this object will produce a file
      containing a log of debug output information.
      @param path If debugOutputLevel is not zero, the debug log file will be
      placed in the location specified by path.
  */
   SNL_FEI_Structure(MPI_Comm comm);

   /** Destructor. */
   virtual ~SNL_FEI_Structure();

   /** Set parameters on this object. Currently three parameters are recognized:
       "debugOutput 'path'" where 'path' is the path to the location where
       debug-log files will be produced.<br>
       "debugOutputOff" which will turn off debug-output.<br>
       "checkSharedNodes" which specifies that a 'sanity-check' should be done
       to make sure that the shared-node information is globally consistent,
       before any communication is attempted which depends on the consistency of
       that information.<br>
       "sharedNodeOwnership <string>" where <string> can be either
       LowNumberedProc or ProcWithLocalElem. LowNumberedProc specifies that the
       lowest-numbered sharing processor will be the owner of shared nodes,
       while ProcWithLocalElem specifies that ownership will be given to the 
       proc with a local element containing the shared node, if not all sharing
       procs have local elements containing the shared node.
   */
   int parameters(int numParams, const char*const* paramStrings);

   int initFields(int numFields, const int* fieldSizes, const int* fieldIDs,
                  const int* fieldTypes = NULL);

   int initElemBlock(GlobalID elemBlockID,
		     int numElements,
		     int numNodesPerElement,
		     const int* numFieldsPerNode,
		     const int* const* nodalFieldIDs,
		     int numElemDofFieldsPerElement,
		     const int* elemDofFieldIDs,
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

   // constraint relation initialization
   //- lagrange multiplier formulation
   int initCRMult(int numCRNodes,
		  const GlobalID* CRNodes,
		  const int *CRFields,
		  int& CRID);

   // - penalty function formulation
   int initCRPen(int numCRNodes,
		 const GlobalID* CRNodes, 
		 const int *CRFields,
		 int& CRID); 

   int initComplete(bool generateGraph = true);

   const std::vector<int>& getFieldIDs() const
     { return fieldIDs_; }

   /** implementation of Lookup::getFieldIDsPtr */
   const int* getFieldIDsPtr()
     {
       int len = fieldDatabase_->size();
       workarray_.resize(len*2);
       fei::copyToArrays<std::map<int,int> >(*fieldDatabase_, len,
						 &workarray_[0],
						 &workarray_[0]+len);
       return( &workarray_[0] );
     }

   /** implementation of Lookup::getFieldSizesPtr */
   const int* getFieldSizesPtr()
     {
       int len = fieldDatabase_->size();
       workarray_.resize(len*2);
       fei::copyToArrays<std::map<int,int> >(*fieldDatabase_, len,
						 &workarray_[0],
						 &workarray_[0]+len);
       return( &workarray_[0]+len );
     }

   /** implementation of Lookup::getNumFields */
   int getNumFields() { return( fieldDatabase_->size() ); };

   /** implementation of Lookup::getFieldSize */
   int getFieldSize(int fieldID)
     {
       std::map<int,int>::const_iterator
	 f_iter = fieldDatabase_->find(fieldID);

       return(f_iter != fieldDatabase_->end() ? (*f_iter).second : -1);
     }

   fei::FieldDofMap<int>& getFieldDofMap()
     { return fieldDofMap_; }

   /** implementation of Lookup::isInLocalElement */
   bool isInLocalElement(int nodeNumber);

   const int* getNumFieldsPerNode(GlobalID blockID);

   const int* const* getFieldIDsTable(GlobalID blockID);

   /** Given a nodeNumber/fieldID pair, return the associated first equation-
       number.
       @return eqnNumber
   */
   int getEqnNumber(int nodeNumber, int fieldID);

   /** Given a global equation number, return the processor on which that
       equation resides.
       @param eqn Global 0-based equation number.
       @return proc Owning processor. If eqn is out of range (less than 0 or
       greater than global-number-of-equations) then -1 is returned.
   */
   int getOwnerProcForEqn(int eqn);


     //////////////////////////////////////////////////////////////////////////
   //now the element-block functions

   int getNumElemBlocks() {return(blockIDs_.size());};
   const GlobalID* getElemBlockIDs() {return(&blockIDs_[0]);};

   void getElemBlockInfo(GlobalID blockID,
                         int& interleaveStrategy, int& lumpingStrategy,
                         int& numElemDOF, int& numElements,
                         int& numNodesPerElem, int& numEqnsPerElem);

   int addBlock(GlobalID blockID);

   int getBlockDescriptor(GlobalID blockID, BlockDescriptor*& block);

   /** Given an index, return the corresponding block-descriptor.
       @return 0 if successful, non-zero if index out of range, etc.
    */
   int getBlockDescriptor_index(int index, BlockDescriptor*& block);

   /** Given a blockID, return its index. Returns -1 if blockID not found. */
   int getIndexOfBlock(GlobalID blockID) const;
   /** Given an index, return a blockID. Returns -1 if index out of range. */
   int getBlockID(unsigned index) const
   {
     if (index < blockIDs_.size()) return blockIDs_[index];
     return -1;
   }

   int allocateBlockConnectivity(GlobalID blockID);
   void destroyConnectivityTables();

   ConnectivityTable& getBlockConnectivity(GlobalID blockID);

   void getScatterIndices_ID(GlobalID blockID, GlobalID elemID,
                             int interleaveStrategy,
                             int* scatterIndices);

   void getScatterIndices_index(int blockIndex, int elemIndex,
                             int interleaveStrategy,
                             int* scatterIndices);

   int getBlkScatterIndices_index(int blockIndex,
				  int elemIndex,
				  int* scatterIndices);

   void getScatterIndices_ID(GlobalID blockID, GlobalID elemID,
                             int interleaveStrategy,
                             int* scatterIndices,
			     int* blkScatterIndices,
			     int* blkSizes);

   void getScatterIndices_index(int blockIndex, int elemIndex,
                             int interleaveStrategy,
                             int* scatterIndices,
			     int* blkScatterIndices,
				int* blkSizes);


   /////////////////////////////////////////////////////////////////////////////
   //now the shared-node lookup functions from the Lookup interface.

   int getNumSharedNodes() {return(nodeCommMgr_->getNumSharedNodes());};

   const int* getSharedNodeNumbers() {
      return(&(nodeCommMgr_->getSharedNodeNumbers())[0]);
   };

   const int* getSharedNodeProcs(int nodeNumber) {
      int index = nodeCommMgr_->getSharedNodeIndex_num(nodeNumber);
      if (index < 0) return(NULL);
      return(&(nodeCommMgr_->getSharedNodeProcs(index))[0]);
   };

   int getNumSharingProcs(int nodeNumber) {
      int index = nodeCommMgr_->getSharedNodeIndex_num(nodeNumber);
      if (index < 0) return(-1);
      return(nodeCommMgr_->getSharedNodeProcs(index).size());
   };

   int getNumSubdomains(int nodeNumber) {
     const NodeDescriptor* node = NULL;
     int err = nodeDatabase_->getNodeWithNumber(nodeNumber, node);
     if (err != 0) return(-1);
     GlobalID nodeID = node->getGlobalNodeID();
     return(nodeCommMgr_->getSharedNodeNumSubdomains(nodeID));
   }

   int* getSubdomainList(int nodeNumber) {
     const NodeDescriptor* node = NULL;
     int err = nodeDatabase_->getNodeWithNumber(nodeNumber, node);
     if (err != 0) return(NULL);
     GlobalID nodeID = node->getGlobalNodeID();
     return(&(*(nodeCommMgr_->getSharedNodeSubdomainList(nodeID)))[0]);

   }

   int translateToReducedNodeNumber(int nodeNumber, int proc);

   /////////////////////////////////////////////////////////////////////////////

   /** implementation of Lookup::getAssociatedNodeNumber */
   int getAssociatedNodeNumber(int eqnNumber)
     {
       int eqn = translateFromReducedEqn(eqnNumber);
       int nodeNumber = nodeDatabase_->getAssociatedNodeNumber(eqn);
       int reducedNodeNumber = -1;
       if (nodeNumber >= 0) {
	 reducedNodeNumber = translateToReducedNodeNumber(nodeNumber, localProc_);
       }
       return( reducedNodeNumber );
     }

   /** implementation of Lookup::getAssociatedFieldID */
   int getAssociatedFieldID(int eqnNumber)
     {
       int eqn = translateFromReducedEqn(eqnNumber);
       return( nodeDatabase_->getAssociatedFieldID(eqn) );
     }


   /////////////////////////////////////////////////////////////////////////////
   //now the point-eqn to block-eqn queries...

   bool isExactlyBlkEqn(int ptEqn) {
      return(blkEqnMapper_->isExactlyBlkEqn(ptEqn));
   };

   int ptEqnToBlkEqn(int ptEqn) {
      return(blkEqnMapper_->eqnToBlkEqn(ptEqn));
   };
 
   int getOffsetIntoBlkEqn(int blkEqn, int ptEqn) {
      return(blkEqnMapper_->getBlkEqnOffset(blkEqn, ptEqn));
   };

   int getBlkEqnSize(int blkEqn) {
      return(blkEqnMapper_->getBlkEqnSize(blkEqn));
   }
   /////////////////////////////////////////////////////////////////////////////


   int getNumActiveNodes() {return(nodeDatabase_->getNodeIDs().size());}

   NodeDatabase& getNodeDatabase() { return( *nodeDatabase_ ); }

   std::map<GlobalID,int>& getActiveNodeIDList()
     { return( nodeDatabase_->getNodeIDs() ); }

   std::vector<int>& getGlobalNodeOffsets() {return(globalNodeOffsets_);}
   std::vector<int>& getGlobalEqnOffsets() {return(globalEqnOffsets_);}
   std::vector<int>& getGlobalBlkEqnOffsets() {return(globalBlkEqnOffsets_);}

   NodeCommMgr& getNodeCommMgr() {return(*nodeCommMgr_);}
   EqnCommMgr&  getEqnCommMgr()  {return(*eqnCommMgr_ );}

   void initializeEqnCommMgr();

   void getEqnInfo(int& numGlobalEqns, int& numLocalEqns,
                   int& localStartRow, int& localEndRow);

   int getEqnNumbers(GlobalID ID, int idType, int fieldID,
		     int& numEqns, int* eqnNumbers);

   int getEqnNumbers(int numIDs, const GlobalID* IDs,
		     int idType, int fieldID,
		     int& numEqns, int* eqnNumbers);

   void getEqnBlkInfo(int& numGlobalEqnBlks, int& numLocalEqnBlks,
                      int& localBlkOffset);

   snl_fei::PointBlockMap& getBlkEqnMapper() {return(*blkEqnMapper_);}

   void destroyMatIndices();

   int getNumMultConstRecords() {return(multCRs_.size());};

   std::map<GlobalID,snl_fei::Constraint<GlobalID>*>&
     getMultConstRecords()
     {return(multCRs_);};

   int getMultConstRecord(int CRID, snl_fei::Constraint<GlobalID>*& multCR)
   {
     std::map<int,snl_fei::Constraint<GlobalID>*>::iterator
       cr_iter = multCRs_.find(CRID);
     int returncode = -1;
     if (cr_iter != multCRs_.end()) {
       multCR = (*cr_iter).second;
       returncode = 0;
     }

     return( returncode );
   }

   int getNumPenConstRecords() {return(penCRs_.size());}
   std::map<GlobalID,snl_fei::Constraint<GlobalID>*>&
     getPenConstRecords()
     { return(penCRs_); }

   int getPenConstRecord(int CRID, snl_fei::Constraint<GlobalID>*& penCR)
   {
     std::map<int,snl_fei::Constraint<GlobalID>*>::iterator
       cr_iter = penCRs_.find(CRID);
     int returncode = -1;
     if (cr_iter != penCRs_.end()) {
       penCR = (*cr_iter).second;
       returncode = 0;
     }

     return( returncode );
   }

   void addSlaveVariable(SlaveVariable* svar) {slaveVars_->push_back(svar);}

   int calculateSlaveEqns(MPI_Comm comm);

   fei::FillableMat* getSlaveDependencies() {return(slaveMatrix_);}

   EqnBuffer* getSlaveEqns() { return(slaveEqns_); }

   int numSlaveEquations() { return(numSlvs_); }

   /** Given a "global" equation number, return true if it is a slave equation.
    */
   bool isSlaveEqn(int eqn);

   /** Given a "global" equation number, return the corresponding equation
       number in the "reduced" equation space. This is a trivial one-to-one
       mapping unless there are slave equations. If there are slave equations,
       then the number of 'reduced' equations is reduced by num-slave-eqns, and 
       any particular global equation number is reduced by the number of slave
       equations that are smaller than it.
       @return true if eqn is a slave-equation, false otherwise.
   */
   bool translateToReducedEqn(int eqn, int& reducedEqn);

   /** Given an EqnCommMgr object with global equation numbers, translate all
       of its indices (row-numbers and column-indices) to the "reduced" equation
       space.
   */
   int translateToReducedEqns(EqnCommMgr& eqnCommMgr);

   /** Given an EqnBuffer object with global equation numbers, translate all of
       its indices (row-numbers and column-indices) to the "reduced" equation
       space.
   */
   int translateToReducedEqns(EqnBuffer& eqnBuf);

   /** Given a ProcEqns object with global equation numbers, translate all of
       its indices (row-numbers and column-indices) to the "reduced" equation
       space.
   */
   int translateToReducedEqns(ProcEqns& procEqns);

   /** Given a CSRMat object with global equation numbers, translate all of its
       indices (row-numbers and column-indices) to the "reduced" equation
       space.
   */
   int translateMatToReducedEqns(fei::CSRMat& mat);

   /** Given an "Reduced" equation number, translate it to the "global"
       numbering, which is the numbering that includes slave equations. This
       is the inverse of the 'translateToReducedEqn' function above.
       @return eqn
   */
   int translateFromReducedEqn(int reducedEqn);

   /** Given a slave equation, fill a std::vector with the
       equation-numbers upon which the slave depends.
       @param slaveEqn
       @param masterEqns Output. NULL if slaveEqn is not a slave equation.
       @return error-code
   */
   int getMasterEqnNumbers(int slaveEqn, std::vector<int>*& masterEqns);

   /** Given a slave equation, fill a std::vector with the
       coefficients of the equations upon which the slave depends.
       @param slaveEqn
       @param masterCoefs Output. NULL if slaveEqn is not a slave equation.
       @return error-code
   */
   int getMasterEqnCoefs(int slaveEqn, std::vector<double>*& masterCoefs);

   /** Given a slave equation, provide the rhs-value associated with the master
       equation that the slave is defined by.
       @param slaveEqn
       @param rhsValue Output. Not referenced if slaveEqn is not a slave
       equation.
       @return error-code 0 if successful, -1 if slaveEqn is not a slave.
   */
   int getMasterEqnRHS(int slaveEqn, double& rhsValue);

   int getNumGlobalEqns() { return( numGlobalEqns_ ); }
   int getNumLocalEqns() { return( numLocalEqns_ ); }
   int getFirstLocalEqn() { return( localStartRow_ ); }
   int getLastLocalEqn() { return( localEndRow_ ); }

   int getFirstReducedEqn() { return( reducedStartRow_ ); }
   int getLastReducedEqn()  { return( reducedEndRow_   ); }

   int getNumGlobalEqnBlks() { return( numGlobalEqnBlks_ ); }
   int getNumLocalEqnBlks() { return( numLocalEqnBlks_ ); }
   int getNumLocalReducedEqnBlks() { return( numLocalReducedEqnBlks_ ); }
   int getGlobalMaxBlkSize() { return(globalMaxBlkSize_); }

   int getNumLocalReducedEqns() { return( numLocalReducedRows_ ); }

   int getMatrixRowLengths(std::vector<int>& rowLengths);
   int getMatrixStructure(int** colIndices, std::vector<int>& rowLengths);

   int getMatrixStructure(int** ptColIndices, std::vector<int>& ptRowLengths,
			  int** blkColIndices, int* blkIndices_1D,
			  std::vector<int>& blkRowLengths,
			  std::vector<int>& numPtRowsPerBlkRow);

   static int gatherSlaveEqns(MPI_Comm comm,
			      EqnCommMgr* eqnCommMgr,
			      EqnBuffer* slaveEqns);

   static int removeCouplings(EqnBuffer& eqnbuf, int& levelsOfCoupling);

   int calcTotalNumElemDOF();
   int calcNumMultCREqns();

   MPI_Comm getCommunicator() const { return( comm_ ); }

#ifdef FEI_HAVE_IOSFWD
   int setDbgOut(std::ostream& ostr, const char* path, const char* feiName);
#else
   int setDbgOut(ostream& ostr, const char* path, const char* feiName);
#endif

 private:

   /// Return a pointer to a NodeDescriptor, may return NULL.
   /// The caller is responsible for checking for a NULL pointer
   /// and doing the right thing.
   NodeDescriptor* findNode(GlobalID nodeID);
   /// Return a reference to a NodeDescriptor, will abort in debug
   /// with an error and seg-fault in optimized mode if a NodeDescriptor
   /// object is not found (i.e. it dereferences a NULL ptr).
   NodeDescriptor& findNodeDescriptor(GlobalID nodeID);

   int writeEqn2NodeMap();

   int getElemNodeDescriptors(int blockIndex, int elemIndex,
                              NodeDescriptor** nodes);

   int getNodeIndices_simple(NodeDescriptor** nodes, int numNodes,
			     int fieldID,
			     int* scatterIndices, int& offset);

   int getNodeIndices_simple(NodeDescriptor** nodes, int numNodes,
			     int fieldID,
			     int* scatterIndices, int& offset,
			     int* blkScatterIndices,
			     int* blkSizes, int& blkOffset);

   int getNodeMajorIndices(NodeDescriptor** nodes, int numNodes,
                           int** fieldIDs, int* fieldsPerNode,
                           int* scatterIndices, int& offset);

   int getNodeBlkIndices(NodeDescriptor** nodes, int numNodes,
			 int* scatterIndices, int& offset);

   int getNodeMajorIndices(NodeDescriptor** nodes, int numNodes,
                           int** fieldIDs, int* fieldsPerNode,
                           int* scatterIndices, int& offset,
			   int* blkScatterIndices,
			   int* blkSizes, int& blkOffset);

   int getNodeMajorIndices(NodeDescriptor** nodes, int numNodes,
                           std::vector<int>* fieldIDs,
                           std::vector<int>& fieldsPerNode,
                           std::vector<int>& scatterIndices);

   int getFieldMajorIndices(NodeDescriptor** nodes, int numNodes,
                            int** fieldIDs, int* fieldsPerNode,
                            int* scatterIndices, int& offset);

   int getFieldMajorIndices(NodeDescriptor** nodes, int numNodes,
                                          std::vector<int>* fieldIDs,
                                          std::vector<int>& fieldsPerNode,
                                          std::vector<int>& scatterIndices);

   void calcGlobalEqnInfo(int numLocallyOwnedNodes, 
			  int numLocalEqns,
			  int numLocalEqnBlks);

   int finalizeActiveNodes();
   int finalizeNodeCommMgr();
   bool activeNodesInitialized();

   int formMatrixStructure();

   int initElemBlockStructure();
   int initMultCRStructure();
   int initPenCRStructure();
   int createMatrixPosition(int row, int col, const char* callingFunction);
   int createMatrixPositions(int row, int numCols, int* cols,
			     const char* callingFunction);
   int createMatrixPositions(fei::CSRMat& mat);

   int createSymmEqnStructure(std::vector<int>& scatterIndices);
   int createBlkSymmEqnStructure(std::vector<int>& scatterIndices);

   int storeElementScatterIndices(std::vector<int>& scatterIndices);
   int storeElementScatterIndices_noSlaves(std::vector<int>& scatterIndices);
   int storeElementScatterBlkIndices_noSlaves(std::vector<int>& scatterIndices);

   void storeLocalNodeIndices(NodeDescriptor& iNode, int iField,
                                     NodeDescriptor& jNode, int jField);
   void storeNodalColumnIndices(int eqn, NodeDescriptor& node,
					       int fieldID);
   void storeNodalRowIndices(NodeDescriptor& node, int fieldID, int eqn);
   void storeNodalSendIndex(NodeDescriptor& node, int fieldID, int col);
   void storeNodalSendIndices(NodeDescriptor& iNode, int iField,
                                  NodeDescriptor& jNode, int jField);

   int assembleReducedStructure();

   bool nodalEqnsAllSlaves(const NodeDescriptor* node, std::vector<int>& slaveEqns);

   int initializeBlkEqnMapper();

   int setNumNodesAndEqnsPerBlock();

   void destroyBlockRoster();

#ifdef FEI_HAVE_IOSFWD
   std::ostream& dbgOut() { return( *dbgOStreamPtr_ ); }
#else
   ostream& dbgOut() { return( *dbgOStreamPtr_ ); }
#endif

   void addCR(int CRID,
	     snl_fei::Constraint<GlobalID>*& cr,
	     std::map<GlobalID,snl_fei::Constraint<GlobalID>* >& crDB);

   int setNodalEqnInfo();
   void setElemDOFEqnInfo();
   int setMultCREqnInfo();

   SNL_FEI_Structure(const SNL_FEI_Structure& src);
   SNL_FEI_Structure& operator=(const SNL_FEI_Structure& src);

   MPI_Comm comm_;

   int localProc_, masterProc_, numProcs_;

   std::vector<int> fieldIDs_;
   std::vector<int> fieldSizes_;
   std::map<int,int>* fieldDatabase_;
   fei::FieldDofMap<int> fieldDofMap_;
   std::vector<int> workarray_;

   std::vector<GlobalID> blockIDs_;
   std::vector<BlockDescriptor*> blocks_;
   std::vector<ConnectivityTable*> connTables_;

   NodeDatabase* nodeDatabase_;

   bool activeNodesInitialized_;

   std::vector<int> globalNodeOffsets_;
   std::vector<int> globalEqnOffsets_;
   std::vector<int> globalBlkEqnOffsets_;

   std::vector<SlaveVariable*>* slaveVars_;
   EqnBuffer* slaveEqns_;
   std::vector<int>* slvEqnNumbers_;
   int numSlvs_, lowestSlv_, highestSlv_;
   fei::FillableMat* slaveMatrix_;
   std::vector<int> globalNumNodesVanished_;
   std::vector<int> localVanishedNodeNumbers_;

   NodeCommMgr* nodeCommMgr_;
   EqnCommMgr* eqnCommMgr_;
   EqnCommMgr* slvCommMgr_;

   int numGlobalEqns_;
   int numLocalEqns_;
   int localStartRow_;
   int localEndRow_;

   int numLocalNodalEqns_;
   int numLocalElemDOF_;
   int numLocalMultCRs_;

   int reducedStartRow_, reducedEndRow_, numLocalReducedRows_;
   fei::FillableMat *Kid_, *Kdi_, *Kdd_;
   fei::CSRMat csrD, csrKid, csrKdi, csrKdd, tmpMat1_, tmpMat2_;
   int reducedEqnCounter_, reducedRHSCounter_;
   std::vector<int> rSlave_, cSlave_;
   std::vector<NodeDescriptor*> work_nodePtrs_;

   bool structureFinalized_;
   bool generateGraph_;

   fei::ctg_set<int>* sysMatIndices_;

   bool blockMatrix_;
   int numGlobalEqnBlks_;
   int numLocalEqnBlks_;
   int numLocalReducedEqnBlks_;
   int localBlkOffset_;
   int localReducedBlkOffset_;
   int globalMaxBlkSize_;

   int firstLocalNodeNumber_;
   int numGlobalNodes_;

   fei::ctg_set<int>* sysBlkMatIndices_;
   bool matIndicesDestroyed_;

   std::vector<int> workSpace_;

   snl_fei::PointBlockMap* blkEqnMapper_;

   std::map<GlobalID, snl_fei::Constraint<GlobalID>* > multCRs_;

   std::map<GlobalID, snl_fei::Constraint<GlobalID>* > penCRs_;

   bool checkSharedNodes_;

   std::string name_;

   int outputLevel_;
   bool debugOutput_;
   std::string dbgPath_;
#ifdef FEI_HAVE_IOSFWD
   std::ostream* dbgOStreamPtr_;
#else
   ostream* dbgOStreamPtr_;
#endif
   bool setDbgOutCalled_;
};

#endif

