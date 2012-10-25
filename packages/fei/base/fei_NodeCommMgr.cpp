/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <fei_mpi.h>

#include <fei_defs.h>

#include <fei_TemplateUtils.hpp>
#include <fei_mpiTraits.hpp>
#include <fei_CommUtils.hpp>
#include <fei_NodeDescriptor.hpp>
#include <fei_NodeCommMgr.hpp>
#include <SNL_FEI_Structure.hpp>

#include <fei_NodeDatabase.hpp>

#undef fei_file
#define fei_file "fei_NodeCommMgr.cpp"
#include <fei_ErrMacros.hpp>

//------Constructor-------------------------------------------------------------
NodeCommMgr::NodeCommMgr(MPI_Comm comm, const SNL_FEI_Structure& problemStructure, int sharedNodeOwnership)
  : sharedNodes_(NULL),
    sharedNodesAllocated_(false),
    sharedNodeOwnership_(sharedNodeOwnership),
    localNodeIDs(),
    remoteNodeIDs(),
    sharedNodeIDs(),
    sharedNodeSubdomains(),
    trivialSubdomainList(1),
    sharingProcs_(),
    sharedNodeNumbers(),
    remoteOwnerProcs_(),
    remoteSharingProcs_(),
    nodesPerOwnerProc_(),
    nodesPerSharingProc_(),
    comm_(comm),
    numProcs_(1),
    localProc_(0),
    maxFields_(0),
    maxBlocks_(0),
    maxSubdomains_(0),
    initCompleteCalled_(false),
    probStruc(problemStructure)
{
  numProcs_ = fei::numProcs(comm_);
  localProc_= fei::localProc(comm_);
  trivialSubdomainList[0] = localProc_;
}

//-----Destructor---------------------------------------------------------------
NodeCommMgr::~NodeCommMgr() {

   for(unsigned i=0; i<sharedNodeIDs.size(); i++) {
      delete sharingProcs_[i];
   }

   delete [] sharedNodes_;
   sharedNodesAllocated_ = false;
}

//------------------------------------------------------------------------------
int NodeCommMgr::getSharedNodeIndex(GlobalID nodeID)
{
  return( fei::binarySearch(nodeID, &sharedNodeIDs[0], sharedNodeIDs.size()) );
}

//------------------------------------------------------------------------------
int NodeCommMgr::getSharedNodeNumSubdomains(GlobalID nodeID)
{
  int index = getSharedNodeIndex(nodeID);

  //If the node isn't one of our shared nodes, then return 1, signifying that
  //the node is in 1 subdomain, that being the local subdomain.
  if (index < 0) return(1);

  //Since the node is one of our shared nodes, it should have an entry in our
  //sharedNodeNumSubdomains array indicating how many subdomains it
  //appears in. So return this number.
  return(sharedNodeSubdomains[index].size());
}

//------------------------------------------------------------------------------
std::vector<int>* NodeCommMgr::getSharedNodeSubdomainList(GlobalID nodeID)
{
  int index = getSharedNodeIndex(nodeID);

  //If the node isn't one of our shared nodes, then return 1, signifying that
  //the node is in 1 subdomain, that being the local subdomain.
  if (index < 0) return( &trivialSubdomainList );

  //Since the node is one of our shared nodes, it should have an entry in our
  //sharedNodeSubdomains array.
  return( &(sharedNodeSubdomains[index]) );
}

//------------------------------------------------------------------------------
int NodeCommMgr::informLocal(const NodeDescriptor& node) {
//
//NodeCommMgr is being informed that 'node' is present in the local
//active node list.
//
//This means that either:
// 1. it is a locally-owned shared node, or
// 3. it isn't even a shared node (it's a purely local node), in which
//    case we'll do nothing.
//

   if (numProcs_ == 1) return(0);

   GlobalID nodeID = node.getGlobalNodeID();

   int sharedIndex = getSharedNodeIndex(nodeID);

   //if this node isn't a shared node, then simply return.
   if (sharedIndex < 0) return(0);

   //Since this node is present as a shared node, let's put nodeID in the
   //localNodeIDs list if it isn't already there.

   int index = fei::sortedListInsert(nodeID, localNodeIDs);

   //if index is -2, it means the localNodeIDs array had an allocation failure.
   if (index == -2) return(-2);

   return(0);
}

//------------------------------------------------------------------------------
int NodeCommMgr::getGlobalMaxFieldsBlocks(int& maxFields, int& maxBlocks)
{
  std::vector<int> localMax(2, 0), globalMax(2, 0);

  for(unsigned i=0; i<sharedNodeIDs.size(); i++) {
    int numFlds = sharedNodes_[i]->getNumFields();
    if (numFlds > localMax[0]) localMax[0] = numFlds;

    int numBlks = sharedNodes_[i]->getNumBlocks();
    if (numBlks > localMax[1]) localMax[1] = numBlks;
  }

  int err = fei::GlobalMax(comm_, localMax, globalMax);
  if (err != 0) return(err);

  maxFields = globalMax[0];
  maxBlocks = globalMax[1];

  return(0);
}

//------------------------------------------------------------------------------
int NodeCommMgr::getGlobalMaxFieldsBlocksSubdomains()
{
  std::vector<int> localMax(3, 0), globalMax(3, 0);

  for(unsigned i=0; i<sharedNodeIDs.size(); i++) {
    int numFlds = sharedNodes_[i]->getNumFields();
    if (numFlds > localMax[0]) localMax[0] = numFlds;

    int numBlks = sharedNodes_[i]->getNumBlocks();
    if (numBlks > localMax[1]) localMax[1] = numBlks;

    int numShrd = sharingProcs_[i]->size();
    if (numShrd > localMax[2]) localMax[2] = numShrd;
  }

  int err = fei::GlobalMax(comm_, localMax, globalMax);
  if (err != 0) return(err);

  maxFields_     = globalMax[0];
  maxBlocks_     = globalMax[1];
  maxSubdomains_ = globalMax[2];

  return(0);
}

//------------------------------------------------------------------------------
std::vector<int>& NodeCommMgr::getSendProcs()
{
  return( remoteSharingProcs_ );
}

//------------------------------------------------------------------------------
std::vector<int>& NodeCommMgr::getRecvProcs()
{
  return( remoteOwnerProcs_ );
}

//------------------------------------------------------------------------------
int NodeCommMgr::getSendMessageLength(int destProc, int& messageLength)
{
  std::vector<int>::iterator
   rs_iter = std::lower_bound(remoteSharingProcs_.begin(),
                              remoteSharingProcs_.end(), destProc);
  if (rs_iter == remoteSharingProcs_.end() || destProc != *rs_iter) {
    ERReturn(-1);
  }

  int idx = rs_iter - remoteSharingProcs_.begin();

  int len = 7+maxFields_*2 + maxBlocks_ + maxSubdomains_;
  messageLength = nodesPerSharingProc_[idx] * (len+1);
  return(0);
}

//------------------------------------------------------------------------------
int NodeCommMgr::getSendMessage(int destProc, std::vector<int>& message)
{
  std::vector<int>::iterator
   rs_iter = std::lower_bound(remoteSharingProcs_.begin(),
                              remoteSharingProcs_.end(), destProc);
  if (rs_iter == remoteSharingProcs_.end() || destProc != *rs_iter) {
    ERReturn(-1);
  }

  int idx = rs_iter - remoteSharingProcs_.begin();
  int len = 0;
  CHK_ERR( getSendMessageLength(destProc, len) );
  message.resize(len);

  packLocalNodesAndData(&message[0], destProc,
			nodesPerSharingProc_[idx], len);
  return(0);
}

//------------------------------------------------------------------------------
int NodeCommMgr::processRecvMessage(int srcProc, std::vector<int>& message)
{
  int idx = fei::binarySearch(srcProc, &remoteOwnerProcs_[0],
                                  remoteOwnerProcs_.size());
  int numNodes = nodesPerOwnerProc_[idx];
  int* msgPtr = &message[0];
  int offset = 0;

  for(int j=0; j<numNodes; j++) {
    int nIndex = fei::binarySearch(msgPtr[j], &sharedNodeIDs[0], sharedNodeIDs.size());
    if (nIndex < 0) return(-1);
    NodeDescriptor* node = sharedNodes_[nIndex];

    int nodeNum         = msgPtr[numNodes+offset++];
    int numFields       = msgPtr[numNodes+offset++];
    int numBlocks       = msgPtr[numNodes+offset++];
    int numSubdomains   = msgPtr[numNodes+offset++];

    node->setNodeNumber(nodeNum);
    node->setNumNodalDOF( msgPtr[numNodes+offset++]);
    node->setBlkEqnNumber(msgPtr[numNodes+offset++]);

    for(int fld=0; fld<numFields; fld++) {
      int fieldID = msgPtr[numNodes+offset++];
      int eqnNum  = msgPtr[numNodes+offset++];
      node->addField(fieldID);
      node->setFieldEqnNumber(fieldID, eqnNum);
    }

    for(int blk=0; blk<numBlocks; blk++) {
      int blk_idx = probStruc.getIndexOfBlock(msgPtr[numNodes+offset++]);
      //if blk_idx < 0 it means the incoming blockID doesn't exist on this proc
      if (blk_idx >= 0) {
        node->addBlockIndex(blk_idx);
      }
    }

    sharedNodeSubdomains[nIndex].resize(numSubdomains);
    for(int sd=0; sd<numSubdomains; sd++) {
      (sharedNodeSubdomains[nIndex])[sd] =
	msgPtr[numNodes+offset++];
    }
  }

  return(0);
}

//------------------------------------------------------------------------------
int NodeCommMgr::exchangeEqnInfo()
{
  //
  //This function will perform the communication necessary to:
  //
  //   1. For each locally owned node, send to each remote sharing proc:
  //      the node's nodeNumber, whether or not it appears in local elements,
  //      fieldID(s), the field-size(s), the first global equation numbers
  //      for those fields, and the processor-subdomains that contain the node.
  //   2. For each remotely owned node, receive the above information from
  //      the owners of these nodes.
  //
  //This is a collective function. All procs must enter it before any can
  //exit it.
  //
  //Most of this function is #ifdef'd according to whether FEI_SER is
  //defined.

#ifndef FEI_SER
   if (numProcs_ == 1) return(0);

   //each proc will find out a max. number of fields, blocks
   //and subdomains to expect per node.

   CHK_ERR( getGlobalMaxFieldsBlocksSubdomains() );

   CHK_ERR( fei::exchange(comm_, this) );

   setNodeNumbersArray();

#endif //#ifndef FEI_SER

   return(0);
}

//------------------------------------------------------------------------------
void NodeCommMgr::packLocalNodesAndData(int* data, 
                                       int proc, int numNodes, int len)
{
//This function packs up nodeIDs, as well as the list containing, for
//each node, the following:
//   node-number
//   numFields
//   numBlocks
//   numSubdomains
//   num-nodal-dof
//   blk-eqn-number
//     'numFields' pairs of (fieldID,eqnNumber)
//     subdomain list, length 'numSubdomains'
//
//Incoming parameter len is:
//  numNodes * (7 + maxFields*2 + maxBlocks + maxSubdomains),
//where maxFields is the maximum number of fields associated with any node,
//maxBlocks is the maximum number of blocks associated with any node, and
//maxSubdomains is the maximum number of subdomains containing any node.
//data is of length numNodes*(len+1).
//
//The above data will all be packed into the 'data' list, with nodeIDs 
//occupying the first numNodes positions, followed by the rest of the data.
//
   int nodeCounter = 0;
   int offset = 0;

   for(unsigned i=0; i<sharedNodeIDs.size(); i++) {
     if (sharedNodes_[i]->getOwnerProc() != localProc_) continue;

      NodeDescriptor* node = sharedNodes_[i];

      //is this local node associated with processor 'proc'?

      std::vector<int>& sProcs = *(sharingProcs_[i]);
      int index = fei::binarySearch(proc, &sProcs[0], sProcs.size());

      //if not, skip to the next iteration...
      if (index < 0) continue;

      if (nodeCounter >= numNodes) {
	fei::console_out() << "NodeCommMgr::packLocalNodesAndData: ERROR,"
	     << " nodeCounter >= numNodes." << FEI_ENDL;
      }

      data[nodeCounter++] = (int)(node->getGlobalNodeID());

      int nodeNum = node->getNodeNumber();
      int numFields = node->getNumFields();
      int numBlocks = node->getNumBlocks();
      const int* fieldIDsPtr = node->getFieldIDList();
      const int* fieldEqnNums = node->getFieldEqnNumbers();
      int blkEqnNumber = node->getBlkEqnNumber();

      const std::vector<unsigned>& nodeBlocks = node->getBlockIndexList();
      std::vector<int>& subdomains = sharedNodeSubdomains[i];

      data[numNodes+offset++] = nodeNum;
      data[numNodes+offset++] = numFields;
      data[numNodes+offset++] = numBlocks;
      data[numNodes+offset++] = subdomains.size();
      data[numNodes+offset++] = node->getNumNodalDOF();
      data[numNodes+offset++] = blkEqnNumber;

      for(int j=0; j<numFields; j++) {
	data[numNodes+offset++] = fieldIDsPtr[j];

	if (offset >= len) {
	  fei::console_out() << "NodeCommMgr::packLocalNodesAndData: ERROR,"
	       << " offset >= len." << FEI_ENDL;
	}

	data[numNodes+offset++] = fieldEqnNums[j];
      }

      for(int kk=0; kk<numBlocks; kk++) {
        GlobalID blkID = probStruc.getBlockID(nodeBlocks[kk]);
        data[numNodes+offset++] = blkID;
      }

      for(unsigned k=0; k<subdomains.size(); k++) {
	data[numNodes+offset++] = subdomains[k];
      }
   }
}

//------------------------------------------------------------------------------
void NodeCommMgr::packRemoteNodesAndData(GlobalID* data,
					 int proc, int numNodes, int len)
{
//
//This function packs up the nodeIDs owned by proc, as well as the list
//containing, for each node, the following:
//   residesLocally (0 or 1) indicating whether it appears in the local
//       processor's element domain.
//   numFields
//   numBlocks
//   numNodalDOF
//     'numFields' entries of (fieldID)
//     'numBlocks' entries of (block)
//
//Incoming parameter len is numNodes * (3 + maxFields + maxBlocks),
//where maxFields is the maximum number of fields associated with any node,
//and maxBlocks is the maximum number of blocks associated with any node.
//nodeIDs is of length numNodes, and
//data is of length numNodes*len.
//
//The above data will all be put in the 'data' list, with the nodeIDs
//occupying the first numNodes positions, followed by the rest of the data.
//
   int nodeCounter = 0;
   int offset = 0;

   for(unsigned i=0; i<sharedNodeIDs.size(); i++) {
      NodeDescriptor* node = sharedNodes_[i];

      int thisProc = node->getOwnerProc();
      if (thisProc != proc) continue;

      if (nodeCounter >= numNodes) {
         fei::console_out() << localProc_ << ": NodeCommMgr::packRemoteNodesAndData: ERROR,"
              << " nodeCounter >= numNodes: " << numNodes << FEI_ENDL;
      }

      data[nodeCounter++] = node->getGlobalNodeID();

      int numFields = node->getNumFields();
      int numBlocks = node->getNumBlocks();
      const int* fieldIDsPtr = node->getFieldIDList();

      const std::vector<unsigned>& nodeBlocks = node->getBlockIndexList();
      int lindex = fei::binarySearch(sharedNodeIDs[i], &localNodeIDs[0], localNodeIDs.size());

      data[numNodes+offset++] = (lindex >= 0) ? 1 : 0;
      data[numNodes+offset++] = (GlobalID)numFields;
      data[numNodes+offset++] = (GlobalID)numBlocks;
      data[numNodes+offset++] = (GlobalID)node->getNumNodalDOF();

      for(int j=0; j<numFields; j++) {
         if (offset >= len) {
            fei::console_out() << "NodeCommMgr::packRemoteNodesAndData: ERROR,"
                 << " offset >= len." << FEI_ENDL;
         }

         data[numNodes+offset++] = (GlobalID)fieldIDsPtr[j];
      }

      for(int k=0; k<numBlocks; k++) {
         if (offset >= len) {
            fei::console_out() << "NodeCommMgr::packRemoteNodesAndData: ERROR,"
                 << " offset >= len." << FEI_ENDL;
         }

         data[numNodes+offset++] = probStruc.getBlockID(nodeBlocks[k]);
      }
   }
}

//------------------------------------------------------------------------------
int NodeCommMgr::createProcList(std::vector<int>& itemsPerProc,
                                std::vector<int>& procs)
{
//
//This function looks through the itemsPerProc list and counts how many
//positions in this list are greater than 0. Then it creates a list of
//the indices of those positions. i.e., itemsPerProc is a list of how many
//items are to be sent to or recvd from each proc. When itemsPerProc is
//greater than 0, that proc is put in the sharingProcs list.
//
   int numProcs = 0;
   int len = itemsPerProc.size();

   for(int i=0; i<len; i++) {
      if (itemsPerProc[i] > 0) numProcs++;
   }

   procs.resize(numProcs);

   int offset = 0;

   for(int i=0; i<len; i++) {
      if (itemsPerProc[i] > 0) procs[offset++] = i;
   }
   return(0);
}

//------------------------------------------------------------------------------
int NodeCommMgr::getSharedNodeIndex_num(int nodeNumber)
{
  for(unsigned i=0; i<sharedNodeNumbers.size(); ++i) {
    if (sharedNodeNumbers[i] == nodeNumber) return(i);
  }

  return(-1);
}

//------------------------------------------------------------------------------
int NodeCommMgr::addSharedNodes( const GlobalID* nodeIDs,
				 int numNodes, 
				 const int* const* procs,
				 const int* numProcs )
{
  //
  //Store the incoming nodeIDs and proc-numbers in the sharedNodeIDs array and
  //sharingProcs_ table.
  //

  try {

  for(int i=0; i<numNodes; i++) {
    int insertPoint = -1;
    int index = fei::binarySearch(nodeIDs[i], sharedNodeIDs, insertPoint);
    if (index < 0) {
      sharingProcs_.insert(sharingProcs_.begin()+insertPoint, new std::vector<int>);

      sharedNodeIDs.insert(sharedNodeIDs.begin()+insertPoint, nodeIDs[i]);

      index = insertPoint;
    }
    
    int err = storeNodeProcs(index, sharingProcs_, procs[i], numProcs[i]);
    if (err != 0) return(err);
  }

  }
  catch(std::runtime_error& exc) {
    fei::console_out() << exc.what() << FEI_ENDL;
    ERReturn(-1);
  }

  return(0);
}

//------------------------------------------------------------------------------
int NodeCommMgr::allocateNodeDescriptorPtrs(NodeDatabase& nodeDB)
{
  //This function is called when all shared nodes have been added. We now
  //allocate a list of pointer-to-NodeDescriptor, of length 
  //sharedNodeIDs.length(), and fill that list with NodeDescriptor-pointers
  //from the node-database.

  if (sharedNodeIDs.size() == 0) return(0);

  if (sharedNodes_ != NULL) delete [] sharedNodes_;
  sharedNodes_ = new NodeDescriptor*[sharedNodeIDs.size()];
  if (sharedNodes_ == NULL) return(-1);

  for(unsigned i=0; i<sharedNodeIDs.size(); i++) {
    NodeDescriptor* node = NULL;
    int err = nodeDB.getNodeWithID(sharedNodeIDs[i], node);
    if (err != 0) return(-1);

    sharedNodes_[i] = node;
  }

  sharedNodesAllocated_ = true;
  return(0);
}

//------------------------------------------------------------------------------
int NodeCommMgr::initComplete(NodeDatabase& nodeDB, bool safetyCheck)
{
//
//This function is called when initialization is complete (i.e., when
//all sharedNodes have been added, allocatedNodeDescriptorPtrs() has been
//called, and informLocal() has been called for all nodes that appear in
//the local finite-element structure.
//
//The task of this function is to assign owner-procs to nodes.
//
//if 'safetyCheck' is true, a global consistency check of the shared node info
//will be performed before the communication is attempted.
//
//return value is 0 if successful, non-zero if an error was encountered
//
  int err = allocateNodeDescriptorPtrs(nodeDB);
  if (err != 0) return(err);

  //Run through the shared nodes, and for each one that has been
  //identified as local, assign its owner to be the lowest-numbered sharing
  //processor, which may or may not be localProc_.

  for(unsigned ii=0; ii<sharedNodeIDs.size(); ii++) {
    std::vector<int>& shProcs = *(sharingProcs_[ii]);

    //first, insert localProc_ in this node's list of sharing proc, since the
    //FEI's initSharedNodes function doesn't mandate that the local processor be
    //included in the list of sharing processors. (i.e., localProc_ may not be
    //in this list yet...)
    std::vector<int>::iterator sh_iter =
      std::lower_bound(shProcs.begin(), shProcs.end(), localProc_);
    if (sh_iter == shProcs.end() || localProc_ != *sh_iter) {
      shProcs.insert(sh_iter, localProc_);
    }

    int proc = shProcs[0];

    sharedNodes_[ii]->setOwnerProc(proc);
  }

  //One of the tasks of this object is to gather information on the number
  //of subdomains each shared node appears in. So one thing we'll do here is
  //size and zero the array that will hold that information.
  sharedNodeSubdomains.resize(sharedNodeIDs.size());

  for(unsigned i=0; i<sharedNodeSubdomains.size(); ++i) {
    sharedNodeSubdomains[i].resize(0);
  }

  //now add the local processor to the sharedNodeSubdomains for each node that
  //appears in our localNodeIDs list.
  for(unsigned i=0; i<sharedNodeIDs.size(); i++) {
    int index = fei::binarySearch(sharedNodeIDs[i], &localNodeIDs[0], localNodeIDs.size());
    if (index >= 0) {
      sharedNodeSubdomains[i].push_back(localProc_);
    }
  }

  if (sharedNodeOwnership_ == PROC_WITH_LOCAL_ELEM) {
    err = adjustSharedOwnership();
    if (err != 0) return(err);
  }

  err = createProcLists();

  if (safetyCheck) {
    err = checkSharedNodeInfo();
    if (err != 0) return(-1);
  }

  exchangeSharedRemoteFieldsBlks();

  initCompleteCalled_ = true;

  return(0);
}

//------------------------------------------------------------------------------
#undef _feiFunc_
#define _feiFunc_ "NodeCommMgr::checkSharedNodeInfo"
int NodeCommMgr::checkSharedNodeInfo()
{
  //This function's task is to "audit" the shared-node info. I.e., to make sure
  //that the info is globally symmetric (e.g., if the local processor thinks it
  //shares a node with another processor, does that other processor also think
  //it shares a node with the local proc?).
  //If this function finds that the shared-node info is consistent/correct, then
  //the return-value is 0. If the shared-node info is found to be wrong, then
  //one or more messages will be written to stderr, and the return-value is -1.
  //
  //This is a collective function, which is relatively expensive. It does a
  //few global reductions...
  //

  if (numProcs_==1) return(0);

  //Make sure that the processors we think are "remote owner"
  //procs, think we are a "remote sharing" proc, and vice-versa.

  std::vector<int> globalOwnerProcs, globalSharingProcs;
  std::vector<int> recvOwnerLengths, recvSharingLengths;

  std::vector<int> globalNodesPerOwnerProcs, globalNodesPerSharingProcs;
  std::vector<int> recvNodesPerOwnerLengths, recvNodesPerSharingLengths;

  //First, gather up each processor's list of remote procs and nodes-per-proc
  //onto all other processors...

  CHK_ERR( fei::Allgatherv(comm_, remoteOwnerProcs_,
				 recvOwnerLengths, globalOwnerProcs) );

  CHK_ERR( fei::Allgatherv(comm_, nodesPerOwnerProc_,
				 recvNodesPerOwnerLengths,
				 globalNodesPerOwnerProcs) );

  CHK_ERR( fei::Allgatherv(comm_, remoteSharingProcs_,
				 recvSharingLengths, globalSharingProcs) );

  CHK_ERR( fei::Allgatherv(comm_, nodesPerSharingProc_,
				 recvNodesPerSharingLengths,
				 globalNodesPerSharingProcs) );

  //Now check the consistency of the global "owners" data against local "sharing"
  //data.
  int err =  checkCommArrays( "owners",
			      globalOwnerProcs, globalNodesPerOwnerProcs,
			      recvOwnerLengths,
			      nodesPerSharingProc_, remoteSharingProcs_ );

  //Now check the consistency of the global "sharing" data against local "owners"
  //data.
  err +=  checkCommArrays( "sharing",
			   globalSharingProcs, globalNodesPerSharingProcs,
			   recvSharingLengths,
			   nodesPerOwnerProc_, remoteOwnerProcs_ );

  int globalErr = 0;

  CHK_ERR( fei::GlobalSum(comm_, err, globalErr) );

  return(globalErr);
}

//------------------------------------------------------------------------------
int NodeCommMgr::checkCommArrays(const char* whichCheck,
				 std::vector<int>& globalRemoteProcs,
				 std::vector<int>& globalNodesPerRemoteProc,
				 std::vector<int>& globalRemoteProcLengths,
				 std::vector<int>& nodesPerRemoteProc,
				 std::vector<int>& remoteProcs)
{
  int offset = 0;

  for(int i=0; i<numProcs_; i++) {
    int length = globalRemoteProcLengths[i];

    if (i==localProc_) { offset += length; continue; }

    for(int j=0; j<length; j++) {
      if (globalRemoteProcs[offset+j] == localProc_) {
	//proc i says that we (localProc_) own nodes that it shares.
	int numShared = globalNodesPerRemoteProc[offset+j];

	int index = fei::binarySearch(i, &remoteProcs[0], remoteProcs.size());
	if (index < 0) {
	  //we don't think proc i shares any nodes that we own.
	  fei::console_out() << "FEI NodeCommMgr::checkSharedNodeInfo "<<whichCheck
	       << " ERROR. Local proc (" << localProc_ 
	       << ") doesn't share nodes with proc " << i << " but proc " << i
	       << " thinks it shares nodes with proc " << localProc_ << FEI_ENDL;
	  return(-1);
	}

	//We think that we own nodesPerRemoteProc[index] nodes that proc i
	//shares.
	int numWeThinkWeShare = nodesPerRemoteProc[index];
	if (numWeThinkWeShare != numShared) {
	  fei::console_out() << "FEI NodeCommMgr::checkSharedNodeInfo "<<whichCheck
	       << " ERROR. Local proc (" << localProc_ << ") thinks it shares "
	       << numWeThinkWeShare << " nodes with proc " << i << ", but proc " 
	       << i << " thinks it shares " << numShared << " nodes with proc "
	       << localProc_ << "." << FEI_ENDL;
	  return(-1);
	}
      }
    }

    offset += length;
  }

  return(0);
}

//------------------------------------------------------------------------------
int NodeCommMgr::adjustSharedOwnership()
{
  //For each shared node that has not been identified as local, assign its
  //owner to be the next lowest-numbered sharing proc. (Each list of sharing
  //procs is sorted by processor-number, so we just assign the owner to be the
  //next one in the list.)
  //
  //If a node is not local, and localProc_ is the lowest sharing proc, then we
  //also need to flag that node as remote and tell other processors that we
  //don't own it.
  //
  remoteNodeIDs.resize(0);
  int err;
  for(unsigned i=0; i<sharedNodeIDs.size(); i++) {
    GlobalID nodeID = sharedNodeIDs[i];

    std::vector<int>& shProcs = *(sharingProcs_[i]);

    if (fei::binarySearch(nodeID, &localNodeIDs[0], localNodeIDs.size()) >= 0) continue;

    int proc = shProcs[0];

    if (proc == localProc_) {
      sharedNodes_[i]->setOwnerProc(shProcs[1]);
      err = fei::sortedListInsert(nodeID, remoteNodeIDs);
      if (err == -2) return(err);
    }
  }

  //Now we need to let the other processors know that the remote nodes 
  //aren't owned by us. This is going to require some communication. We'll
  //gather the nodeIDs onto all processors, after which each processor
  //will reset the owner proc for that node. (Later, as an optimization, I'll
  //do this without all-to-all communication.)

  std::vector<GlobalID> allRemoteNodeIDs;
  std::vector<int> numPerProc;

  err = fei::Allgatherv(comm_, remoteNodeIDs, numPerProc, allRemoteNodeIDs);
  if (err != 0) return(-1);

  //Now we need to run through the global list of 'special' nodes, and for the ones
  //that we do own (appear locally), add them to a new list that will be once again
  //globally gathered. That new list will then be used by each processor in setting
  //the nodes' real owners.

  //we'll keep the 'new' list in remoteNodeIDs.
  remoteNodeIDs.resize(0);

  int offset = 0;
  for(unsigned i=0; i<numPerProc.size(); i++) {
    for(int j=0; j<numPerProc[i]; j++) {

      //skip the nodes that we sent, we already know we don't own those.
      if ((int)i==localProc_) {offset++; continue;}

      GlobalID nodeID = allRemoteNodeIDs[offset++];
      int index = getSharedNodeIndex(nodeID);

      //if it's not even one of our shared nodes, then continue.
      if (index < 0) continue;

      if (fei::binarySearch(nodeID, &localNodeIDs[0], localNodeIDs.size()) >= 0) {
        fei::sortedListInsert(nodeID, remoteNodeIDs);
      }
    }
  }

  //now re-gather the remoteNodeIDs list to all processors. This time, we should only
  //receive nodeIDs from processors that can be valid owners. i.e., processors that
  //have those nodes in at least one local element.
  err = fei::Allgatherv(comm_, remoteNodeIDs, numPerProc, allRemoteNodeIDs);
  if (err != 0) return(-1);

  //Now we run the 'allRemoteNodeIDs' list for the last time, setting the owner-proc
  //for each node. We'll run the list from the back to the front so that if multiple
  //processors are possible owners, the lowest-numbered one will be the last one that
  //get's set.
  offset = allRemoteNodeIDs.size()-1;
  for(int i=(int)numPerProc.size()-1; i>=0; i--) {
    for(int j=0; j<numPerProc[i]; j++) {
      GlobalID nodeID = allRemoteNodeIDs[offset--];
      int index = getSharedNodeIndex(nodeID);

      if (index < 0) continue;

      sharedNodes_[index]->setOwnerProc(i);
    }
  }

  return(0);
}

//------------------------------------------------------------------------------
void NodeCommMgr::setNodeNumbersArray()
{
  sharedNodeNumbers.resize(sharedNodeIDs.size());

  for(unsigned i=0; i<sharedNodeIDs.size(); i++) {
    sharedNodeNumbers[i] = sharedNodes_[i]->getNodeNumber();
  }
}

//------------------------------------------------------------------------------
int NodeCommMgr::createProcLists()
{
  std::vector<int> localNodesPerProc(numProcs_, 0);
  std::vector<int> remoteNodesPerProc(numProcs_, 0);

  //first, figure out how many locally-owned nodes each remote processor is
  //associated with, and how many remotely-owned nodes we'll be recv'ing info
  //about from each remote processor.

  for(unsigned i=0; i<sharedNodeIDs.size(); i++) {
    int proc = sharedNodes_[i]->getOwnerProc();

    if (proc != localProc_) {
      remoteNodesPerProc[proc]++;
    }
    else {
      std::vector<int>& shProcs = *(sharingProcs_[i]);
      for(unsigned j=0; j<shProcs.size(); j++) {
	int sproc = shProcs[j];

	if (sproc != localProc_) {
	  localNodesPerProc[sproc]++;
	}
      }
    }
  }

  //now create condensed lists of remote owner procs, and
  //remote sharing procs.
  int err = createProcList(remoteNodesPerProc, remoteOwnerProcs_);
  if (err != 0) return(err);

  err = createProcList(localNodesPerProc, remoteSharingProcs_);
  if (err != 0) return(err);


  nodesPerOwnerProc_.resize(remoteOwnerProcs_.size());

  nodesPerSharingProc_.resize(remoteSharingProcs_.size());

  int offset = 0;
  for(int i=0; i<numProcs_; i++) {
    if (remoteNodesPerProc[i] > 0) 
      nodesPerOwnerProc_[offset++] = remoteNodesPerProc[i];
  }

  offset = 0;
  for(int i=0; i<numProcs_; i++) {
    if (localNodesPerProc[i] > 0) 
      nodesPerSharingProc_[offset++] = localNodesPerProc[i];
  }

  return(0);
}

//------------------------------------------------------------------------------
int NodeCommMgr::exchangeSharedRemoteFieldsBlks()
{
  //first each proc will find out a max. number of fields and blocks to expect
  //per node.

  //most of this function is #ifdef'd according to whether FEI_SER is
  //defined.
#ifndef FEI_SER
  int maxFields, maxBlocks;
  int err = getGlobalMaxFieldsBlocks(maxFields, maxBlocks);
  if (err) return(-1);

  //now we can allocate lists to recv into and launch the irecv's.
  //from each processor, we'll recv a list of length:
  //            num-nodes*(4+ maxFields + maxBlocks)

  int len = 4 + maxFields + maxBlocks;

  GlobalID** recvData = NULL;
  MPI_Request* recvDataReqs = NULL;

  unsigned i, numRecvProcs = remoteSharingProcs_.size();

  if (numRecvProcs > 0) {
    recvData = new GlobalID*[numRecvProcs];
    recvDataReqs = new MPI_Request[numRecvProcs];
  }

  int dataTag = 19904;

  int numRcvStarted = 0;
  for(i=0; i<remoteSharingProcs_.size(); i++) {
    int numRecvNodes = nodesPerSharingProc_[i];
    recvData[i] = new GlobalID[numRecvNodes*(len+1)];
    MPI_Irecv(recvData[i], numRecvNodes*(len+1),
	      fei::mpiTraits<GlobalID>::mpi_type(),
	      remoteSharingProcs_[i], dataTag, comm_, &recvDataReqs[i]);
    numRcvStarted++;
  }

  //next, send all outgoing messages.

  fei::Barrier(comm_);

  for(i=0; i<remoteOwnerProcs_.size(); i++) {
    int numSendNodes = nodesPerOwnerProc_[i];

    std::vector<GlobalID> sendData(numSendNodes*(len+1), 0);

    packRemoteNodesAndData(&sendData[0], remoteOwnerProcs_[i],
                           numSendNodes, numSendNodes*len);

    MPI_Send(&sendData[0], sendData.size(),
	      fei::mpiTraits<GlobalID>::mpi_type(),
	     remoteOwnerProcs_[i], dataTag, comm_);
  }

  //finally, complete the irecvs and put away the node field info.
  int numCompleted = 0;
  for(i=0; i<remoteSharingProcs_.size(); i++) {
    MPI_Status status;
    int index = i;
    MPI_Wait(&recvDataReqs[index], &status);
    numCompleted++;
    int remoteProc = status.MPI_SOURCE;

    int offset = 0;
    int numNodes = nodesPerSharingProc_[index];

    for(int j=0; j<numNodes; j++) {
      int nIndex = fei::binarySearch(recvData[index][j], &sharedNodeIDs[0], sharedNodeIDs.size());
      if (nIndex < 0) {
	fei::console_out() << "NodeCommMgr::exchangeSharedRemote...: error, unknown nodeID "
	     << (int)recvData[index][j] << ", " << j
	     << "th node recvd from proc "
	     <<remoteSharingProcs_[index]
	     << ". Probably a communication mis-match, we expected " 
	     << numNodes
	     << " nodes from that proc, but recvd less than that." << FEI_ENDL;
	std::abort();
      }

      int residesRemotely = (int)recvData[index][numNodes+offset++];

      if (residesRemotely) {
        std::vector<int>& snSubd = sharedNodeSubdomains[nIndex];
        std::vector<int>::iterator sn_iter =
          std::lower_bound(snSubd.begin(), snSubd.end(), remoteProc);
        if (sn_iter == snSubd.end() || remoteProc != *sn_iter) {
          snSubd.insert(sn_iter, remoteProc);
        }
      }
      int numFields       = (int)recvData[index][numNodes+offset++];
      int numBlocks       = (int)recvData[index][numNodes+offset++];
      sharedNodes_[nIndex]->
 	     setNumNodalDOF((int)recvData[index][numNodes+offset++]);

      for(int fld=0; fld<numFields; fld++) {
        int fieldID = (int)recvData[index][numNodes+offset++];

        sharedNodes_[nIndex]->addField(fieldID);
      }

      for(int blk=0; blk<numBlocks; blk++) {
        int blk_idx = probStruc.getIndexOfBlock(recvData[index][numNodes+offset++]);
        //if blk_idx < 0 it means the incoming blockID doesn't exist on this proc
        if (blk_idx >= 0) {
          sharedNodes_[nIndex]->addBlockIndex(blk_idx);
        }
      }
    }
  }

  if (numRcvStarted != numCompleted) {
    fei::console_out() << "NodeCommMgr::exchangeSharedRemote...: recv-send mismatch;"
         << " numRcvStarted: " << numRcvStarted << ", numCompleted: "
         << numCompleted << FEI_ENDL;
    std::abort();
  }

  for(i=0; i<numRecvProcs; i++) {
    delete [] recvData[i];
  }

  delete [] recvData;
  delete [] recvDataReqs;

#endif //#ifndef FEI_SER

  return(0);
}

//------------------------------------------------------------------------------
int NodeCommMgr::storeNodeProcs(int index,
				std::vector<std::vector<int>*>& procTable,
				const int* procs, int numProcs)
{
//Private NodeCommMgr function.
//
//This function stores 'procs' in row 'index' of procTable, maintaining order
//in that row.
//
  std::vector<int>& row_index = *(procTable[index]);
  for(int i=0; i<numProcs; i++) {
    std::vector<int>::iterator r_iter =
      std::lower_bound(row_index.begin(), row_index.end(), procs[i]);
    if (r_iter == row_index.end() || procs[i] != *r_iter) {
      row_index.insert(r_iter, procs[i]);
    }
  }

  return(0);
}

