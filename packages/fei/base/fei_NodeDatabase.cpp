/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_NodeDatabase.hpp>

#include <fei_CommUtils.hpp>
#include <fei_NodeDescriptor.hpp>
#include <fei_NodeCommMgr.hpp>

#include <fei_macros.hpp>
#include <fei_TemplateUtils.hpp>
#undef fei_file
#define fei_file "fei_NodeDatabase.cpp"

#include <fei_ErrMacros.hpp>

//------------------------------------------------------------------------------
NodeDatabase::NodeDatabase(std::map<int,int>* fieldDatabase,
                           NodeCommMgr* nodeCommMgr)
  : nodePtrs_(),
    eqnNumbers_(0), eqnNodeIndices_(),
    nodeIDs_(),
    nodeNumbers_(),
    synchronized_(false),
    need_to_alloc_and_sync_(true),
    fieldDB_(fieldDatabase),
    nodeCommMgr_(nodeCommMgr),
    numLocalNodes_(0),
    firstLocalNodeNumber_(-1), lastLocalNodeNumber_(-1),
    nodePool_()
{
}

//------------------------------------------------------------------------------
NodeDatabase::~NodeDatabase()
{
  deleteMemory();
}

//------------------------------------------------------------------------------
void NodeDatabase::deleteMemory()
{
  for(size_t i=0; i<nodePtrs_.size(); ++i) {
    nodePool_.destroy(nodePtrs_[i]);
    nodePool_.deallocate(nodePtrs_[i], 1);
  }
}

//------------------------------------------------------------------------------
int NodeDatabase::getNodeWithID(GlobalID nodeID, const NodeDescriptor*& node) const
{
  int index = getIndexOfID(nodeID);
  if (index < 0) {
    //fei::console_out() << "FEI NodeDatabase: node " << (int)nodeID << " not found."<<FEI_ENDL;
    return(-1);
  }

  node = nodePtrs_[index];
  return(0);
}

//------------------------------------------------------------------------------
int NodeDatabase::getNodeWithID(GlobalID nodeID, NodeDescriptor*& node)
{
  int index = getIndexOfID(nodeID);
  if (index < 0) {
    //fei::console_out() << "FEI NodeDatabase: node " << (int)nodeID << " not found."<<FEI_ENDL;
    return(-1);
  }

  node = nodePtrs_[index];
  return(0);
}

//------------------------------------------------------------------------------
int NodeDatabase::getNodeWithNumber(int nodeNumber, const NodeDescriptor*& node) const
{
  if (!synchronized_) ERReturn(-1);

  std::map<int,int>::const_iterator iter = nodeNumbers_.find(nodeNumber);
  if (iter == nodeNumbers_.end()) {
    // The node wasn't found, return a NULL ptr.
    node = NULL;
    // Indicate that the node is NULL.
    return -1;
  }

  int index = iter->second;
  node = nodePtrs_[index];

  return(0);
}

//------------------------------------------------------------------------------
int NodeDatabase::getNodeWithEqn(int eqnNumber, const NodeDescriptor*& node) const
{
  int insertPoint = -1;
  int index = fei::binarySearch(eqnNumber, eqnNumbers_, insertPoint);

  if (index >= 0) {
    node = nodePtrs_[eqnNodeIndices_[index]];
  }
  else if (insertPoint > 0) {
    node = nodePtrs_[eqnNodeIndices_[insertPoint-1]];
  }
  else {
    //We only reach this line if insertPoint==0, which means the specified
    //eqnNumber lies below the equation-number for the first node in our list.
    //In other words, we don't have the requested node.
    return(-1);
  }

  //Now let's make sure that the specified equation is less than or equal to
  //this node's last equation.

  int numNodeFields = node->getNumFields();
  if (numNodeFields <= 0) return(-1);

  int lastFieldOnNode = node->getFieldIDList()[numNodeFields-1];

  int lastFieldSize = (*fieldDB_)[lastFieldOnNode];

  int lastEqnOnNode = node->getFieldEqnNumbers()[numNodeFields-1] +
            lastFieldSize - 1;

  if (eqnNumber <= lastEqnOnNode) return(0);

  return(-1);
}

//------------------------------------------------------------------------------
void NodeDatabase::getNodeAtIndex(int i, const NodeDescriptor*& node) const
{
  int nnodes = nodePtrs_.size();
  if (i>=0 && i < nnodes) {
    node = nodePtrs_[i];
  }
  else {
    node = NULL;
  }
}

//------------------------------------------------------------------------------
void NodeDatabase::getNodeAtIndex(int i, NodeDescriptor*& node)
{
  int nnodes = nodePtrs_.size();
  if (i>=0 && i < nnodes) {
    node = nodePtrs_[i];
  }
  else {
    node = NULL;
  }
}

//------------------------------------------------------------------------------
int NodeDatabase::countLocalNodalEqns(int localRank)
{
  int numEqns = 0;

  for(size_t i=0; i<nodePtrs_.size(); i++) {
    NodeDescriptor* node = nodePtrs_[i];

    if (node->getOwnerProc() == localRank) {
      int numFields = node->getNumFields();
      const int* fieldIDList = node->getFieldIDList();

      for(int j=0; j<numFields; j++) {
        int numParams =        (*fieldDB_)[fieldIDList[j]];

        numEqns += numParams;
      }
    }
  }

  return(numEqns);
}

//------------------------------------------------------------------------------
int NodeDatabase::countLocalNodeDescriptors(int localRank)
{
  int numLocal = 0;
  for(size_t i=0; i<nodePtrs_.size(); i++) {
    if (nodePtrs_[i]->getOwnerProc() == localRank) numLocal++;
  }

  return(numLocal);
}

//------------------------------------------------------------------------------
int NodeDatabase::getIndexOfID(GlobalID nodeID) const
{
  std::map<GlobalID,int>::const_iterator
    iter = nodeIDs_.find(nodeID);

  if (iter == nodeIDs_.end()) return(-1);

  return( iter->second );
}

//------------------------------------------------------------------------------
int NodeDatabase::initNodeID(GlobalID nodeID)
{
  static NodeDescriptor dummyNode;

  int index = nodeIDs_.size();
  std::map<GlobalID,int>::iterator
    iter = nodeIDs_.lower_bound(nodeID);

  if (iter == nodeIDs_.end() || iter->first != nodeID) {
    nodeIDs_.insert(iter, std::make_pair(nodeID,index));

    NodeDescriptor* nodePtr = nodePool_.allocate(1);
    nodePool_.construct(nodePtr, dummyNode);

    nodePtr->setGlobalNodeID(nodeID);

    nodePtrs_.push_back(nodePtr);

    need_to_alloc_and_sync_ = true;
  }

  return(0);
}

//------------------------------------------------------------------------------
int NodeDatabase::initNodeIDs(GlobalID* nodeIDs, int numNodes)
{
  static NodeDescriptor dummyNode;

  for(int i=0; i<numNodes; i++) {
    initNodeID(nodeIDs[i]);
  }

  return(0);
}

//------------------------------------------------------------------------------
int NodeDatabase::synchronize(int firstLocalNodeNumber,
                              int firstLocalEqn,
                              int localRank,
                              MPI_Comm comm)
{
  eqnNumbers_.reserve(nodePtrs_.size());
  eqnNodeIndices_.reserve(nodePtrs_.size());

  eqnNumbers_.resize(0);
  eqnNodeIndices_.resize(0);

  firstLocalNodeNumber_ = firstLocalNodeNumber;
  int nodeNumber = firstLocalNodeNumber;
  int numEqns = 0;
  
  nodeNumbers_.clear();

  numLocalNodes_ = 0;
  std::map<GlobalID,int>::iterator
    iter = nodeIDs_.begin(), iter_end = nodeIDs_.end();

  for(; iter!=iter_end; ++iter) {
    int i = iter->second;
    NodeDescriptor* node = NULL;
    getNodeAtIndex(i, node);
    if (node==NULL) continue;

    int numFields = node->getNumFields();
    const int* fieldIDList = node->getFieldIDList();

    int numNodalDOF = 0;
    int firstEqnNumber, eqnNumber;

    for(int j=0; j<numFields; j++) {
      int numFieldParams = (*fieldDB_)[fieldIDList[j]];
      numNodalDOF += numFieldParams;

      if (node->getOwnerProc() == localRank) {
        eqnNumber = firstLocalEqn + numEqns;
        if (j==0) firstEqnNumber = eqnNumber;

        numEqns += numFieldParams;

        node->setFieldEqnNumber(fieldIDList[j], eqnNumber);
      }
    }

    if (node->getOwnerProc() == localRank) {
      node->setNodeNumber(nodeNumber++);
      numLocalNodes_++;

      int insertPoint = fei::sortedListInsert(firstEqnNumber, eqnNumbers_);
      if (insertPoint == -2) ERReturn(-2);
      if (insertPoint >= 0) eqnNodeIndices_.insert(eqnNodeIndices_.begin()+insertPoint, i);
    }

    node->setNumNodalDOF(numNodalDOF);

    int thisNodeNumber = node->getNodeNumber();
    nodeNumbers_.insert(std::make_pair(thisNodeNumber, i));
  }

  lastLocalNodeNumber_ = nodeNumber - 1;

  //  next, have the node comm manager get the field IDs and 
  //  equation numbers for all of the nodes that we know about but don't
  //  own. i.e., the remotely-owned shared nodes.
  //  Again, we'll get the nodeNumber info for those nodes while we're at it.
  //
  CHK_ERR( nodeCommMgr_->exchangeEqnInfo() );

  //Now finish up by inserting equation-numbers for shared nodes into our
  //eqnNumbers_ and eqnNodeIndices_ lists, for future lookups...

  int numSharedNodes = nodeCommMgr_->getNumSharedNodes();
  for(int i=0; i<numSharedNodes; i++) {
    NodeDescriptor& node = nodeCommMgr_->getSharedNodeAtIndex(i);
    GlobalID nodeID = node.getGlobalNodeID();
    int index = getIndexOfID(nodeID);
    int nDOF = node.getNumNodalDOF();
    if (nDOF <= 0) {
      continue;
      //FEI_COUT << "localRank " << localRank << ", node "<<nodeID<<" has nDOF=" << nDOF<<FEI_ENDL;
      //ERReturn(-1);
    }
    int firstEqn = node.getFieldEqnNumbers()[0];
    int insertPoint = fei::sortedListInsert(firstEqn, eqnNumbers_);
    if (insertPoint == -2) ERReturn(-2);
    if (insertPoint >= 0) eqnNodeIndices_.insert(eqnNodeIndices_.begin()+insertPoint, index);

    int thisNodeNumber = node.getNodeNumber();
    nodeNumbers_.insert(std::make_pair(thisNodeNumber, index));
  }

  synchronized_ = true;
  need_to_alloc_and_sync_ = false;

  return(0);
}

//------------------------------------------------------------------------------
int NodeDatabase::getAssociatedNodeNumber(int eqnNumber)
{
  int insertPoint = -1;
  int index = fei::binarySearch(eqnNumber, eqnNumbers_, insertPoint);

  if (index >= 0) {
    return( nodePtrs_[eqnNodeIndices_[index]]->getNodeNumber() );
  }

  if (insertPoint > 0) {
    NodeDescriptor& node = *(nodePtrs_[eqnNodeIndices_[insertPoint-1]]);
    const int* fieldEqnNumbers = node.getFieldEqnNumbers();
    const int* fieldIDList = node.getFieldIDList();
    int numFields = node.getNumFields();

    int lastEqn = fieldEqnNumbers[numFields-1];

    int fieldSize = -1;
    std::map<int,int>::const_iterator f_iter
      = fieldDB_->find(fieldIDList[numFields-1]);
    if (f_iter == fieldDB_->end()) ERReturn(-1);
    fieldSize = (*f_iter).second;

    //if eqnNumber is inside the range of eqn-numbers assocated with this node,
    //then return this node's node-number
    if (eqnNumber >= fieldEqnNumbers[0] && (lastEqn+fieldSize - 1) >= eqnNumber) {
      return( node.getNodeNumber() );
    }
  }

  //if we get to here, then we don't know how to locate the node with this
  //eqn-number...
  return(-1);
}

//------------------------------------------------------------------------------
int NodeDatabase::getAssociatedFieldID(int eqnNumber)
{
  int insertPoint = -1;
  int index = fei::binarySearch(eqnNumber, eqnNumbers_, insertPoint);

  int index2 = index;
  if (index2 < 0) index2 = insertPoint-1;

  if (index2 < 0) ERReturn(-1);

  NodeDescriptor& node = *(nodePtrs_[eqnNodeIndices_[index2]]);

  const int* fieldEqnNumbers = node.getFieldEqnNumbers();
  const int* fieldIDList = node.getFieldIDList();
  int numFields = node.getNumFields();

  int lastEqn = fieldEqnNumbers[numFields-1];

  int fieldSize = -1;
  std::map<int,int>::const_iterator f_iter
    = fieldDB_->find(fieldIDList[numFields-1]);
  if (f_iter == fieldDB_->end()) ERReturn(-1);
  fieldSize = (*f_iter).second;

  //if eqnNumber is outside the range of eqn-numbers that are associated with
  //this node, then we're in trouble...
  if (eqnNumber < fieldEqnNumbers[0] || eqnNumber > lastEqn+fieldSize) {
    ERReturn(-1);
  }

  //ok, so we're ready to figure out which fieldID eqnNumber is associated with.
  for(int i=0; i<numFields-1; i++) {
    if (eqnNumber >= fieldEqnNumbers[i] && eqnNumber < fieldEqnNumbers[i+1]) {
      return(fieldIDList[i]);
    }
  }

  //if we get to here, then eqnNumber is associated with the last fieldID on the
  //node.
  return(fieldIDList[numFields-1]);
}
