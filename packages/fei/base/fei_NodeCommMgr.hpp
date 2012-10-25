#ifndef _NodeCommMgr_hpp_
#define _NodeCommMgr_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_fwd.hpp>

#include <fei_CommUtils.hpp>

/**
  NodeCommMgr (Node communication manager) is responsible for
  keeping track of nodes that require communication (shared nodes).
  
  All shared nodes are put into this class through the addSharedNodes function.
  When they've all been put in, initComplete() is called, which
  allocates a list of pointer-to-NodeDescriptors.

  Before initComplete is called, all nodes that appear in local element-
  connectivities should be passed to NodeCommMgr using the informLocal
  function.
  
  The communication that must then happen is this:
  
    For all nodes that we own, we must send their field-IDs and (global)
    equation numbers to the remote processors that share the node.
  
    For all nodes that we don't own, we need to receive the field-IDs
    and equation numbers.
*/

class NodeCommMgr : public fei::MessageHandler<int> {
 public:
   enum { STRICTLY_LOW_PROC, PROC_WITH_LOCAL_ELEM, CALLER_SPECIFIES };

   NodeCommMgr(MPI_Comm comm, const SNL_FEI_Structure& problemStructure, int sharedNodeOwnership=STRICTLY_LOW_PROC);
   virtual ~NodeCommMgr();

   size_t getNumSharedNodes() {return(sharedNodeIDs.size());}
   std::vector<GlobalID>& getLocalNodeIDs() {return(localNodeIDs);}
   std::vector<GlobalID>& getSharedNodeIDs() {return(sharedNodeIDs);}
   std::vector<int>& getSharedNodeNumbers() {return(sharedNodeNumbers);}

   int getSharedNodeIndex_num(int nodeNumber);

   int addSharedNodes(const GlobalID* nodeIDs, int numNodes,
		      const int* const* procs, const int* numProcs);

   int initComplete(NodeDatabase& nodeDB, bool safetyCheck);

   int informLocal(const NodeDescriptor& node);

   int exchangeEqnInfo();

   int getSharedNodeIndex(GlobalID nodeID);

   int getSharedNodeNumSubdomains(GlobalID nodeID);
   std::vector<int>* getSharedNodeSubdomainList(GlobalID nodeID);

   NodeDescriptor& getSharedNodeAtIndex(int index)
     {return(*(sharedNodes_[index]));}

   std::vector<int>& getSharedNodeProcs(int index) 
     {return(*(sharingProcs_[index]));};

   void setSharedOwnershipRule(int ownershipRule)
     { sharedNodeOwnership_ = ownershipRule; }

   std::vector<int>& getSendProcs();
   std::vector<int>& getRecvProcs();

   int getSendMessageLength(int destProc, int& messageLength);
   int getSendMessage(int destProc, std::vector<int>& message);
   int processRecvMessage(int srcProc, std::vector<int>& message);

 private:
   NodeCommMgr(const NodeCommMgr& src);
   NodeCommMgr& operator=(const NodeCommMgr& src);

   int allocateNodeDescriptorPtrs(NodeDatabase& nodeDB);

   int storeNodeProcs(int index, std::vector<std::vector<int>*>& procTable,
		      const int* procs, int numProcs);

   int checkSharedNodeInfo();

   int checkCommArrays(const char* whichCheck,
		       std::vector<int>& globalRemoteProcs,
		       std::vector<int>& globalNodesPerRemoteProc,
		       std::vector<int>& globalRemoteProcLengths,
		       std::vector<int>& nodesPerRemoteProc,
		       std::vector<int>& remoteProcs);

   void setNodeNumbersArray();

   void packLocalNodesAndData(int* data, int proc,
                             int numNodes, int len);
   void packRemoteNodesAndData(GlobalID* data, int proc,
                             int numNodes, int len);

   int adjustSharedOwnership();

   int createProcLists();

   int createProcList(std::vector<int>& itemsPerProc,
		       std::vector<int>& procs);

   int exchangeSharedRemoteFieldsBlks();

   int getGlobalMaxFieldsBlocks(int& maxFields, int& maxBlocks);

   int getGlobalMaxFieldsBlocksSubdomains();

   NodeDescriptor** sharedNodes_;
   bool sharedNodesAllocated_;

   int sharedNodeOwnership_;

   std::vector<GlobalID> localNodeIDs;
   std::vector<GlobalID> remoteNodeIDs;

   std::vector<GlobalID> sharedNodeIDs; //global node identifiers

   std::vector<std::vector<int> > sharedNodeSubdomains;
                                          //subdomains each shared node
                                          //appears in.
   std::vector<int> trivialSubdomainList;

   std::vector<std::vector<int>*> sharingProcs_;//table, i-th row is a list of procs
                                          //associated with i-th shared node

   std::vector<int> sharedNodeNumbers;

   std::vector<int> remoteOwnerProcs_, remoteSharingProcs_;
   std::vector<int> nodesPerOwnerProc_, nodesPerSharingProc_;

   MPI_Comm comm_;
   int numProcs_, localProc_;

   int maxFields_;
   int maxBlocks_;
   int maxSubdomains_;

   bool initCompleteCalled_;
   const SNL_FEI_Structure& probStruc;
};

#endif

