#ifndef _NodeDatabase_hpp_
#define _NodeDatabase_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_macros.hpp"
#include "fei_fwd.hpp"
#include "fei_defs.h"
#include "fei_NodeDescriptor.hpp"
#include "fei_Pool_alloc.hpp"
#include "fei_mpi.h"

#include <map>

/** Container that holds NodeDescriptors, and is able to reference them by
 global identifiers, or by nodeNumbers or eqnNumbers.

These are the 3 fundamental ways to refer to a node: 
<ul>
<li> nodeID -- global node identifier. nodeID's are provided by the calling
finite-element application. nodeID's are not necessarily contiguous,
nor do they necessarily start at 0 or 1 or anything else. The only requirement
is that they are globally unique.
<li> nodeNumber -- node-numbers are assigned by the FEI implementation. They are
globally 0-based and contiguous. Each processor will have a first and last
node-number, and all nodes lying in that range may be assumed to be locally
owned. However, processors usually share nodes that are not locally owned, and
thus those nodes will have node-numbers that lie outside the 'local range'.
<li> eqnNumber -- equation-numbers are also assigned by the FEI implementation,
and are also globally 0-based and contiguous. Multiple eqnNumbers may refer to
a single node, whereas each node only has one nodeNumber and one nodeID.
</ul>

The goal of this class is to be able to provide a node's descriptor as fast as
possible, given a nodeID, nodeNumber, or eqnNumber. Binary searches are used wherever
possible.

The return-value of all functions in this class (except trivial query/accessor
functions) is an error-code. If the function is successful, the error-code is 
0. If a node-not-found error occurs, -1 is returned. If an allocation fails, 
-2 is returned.

This class is intended to be used by an FEI implementation, with usage
proceeding in 3 distinct 'phases':
<ul>
<li>nodeID initialization. After construction, NodeDatabase can accept an
arbitrary number of calls to the initNodeID function.
<li>NodeDescriptor initialization. After all nodeIDs have been 'pushed' into
the NodeDatabase using initNodeID, the function allocateNodeDescriptors()
should be called. After this, some of the node-query functions are available 
to serve out NodeDescriptors on which the calling code can set things like
fieldIDs with associated equation-numbers, owning processors, etc. After the 
NodeDescriptors have been initialized, NodeDatabase::synchronize should be 
called.
<li>general node queries. After synchronize has been called, the user is
able to request NodeDescriptors by supplying a nodeNumber or an eqnNumber.
</ul>
*/

class NodeDatabase {
 public:
  /** Constructor. */
  NodeDatabase(std::map<int,int>* fieldDatabase,
	       NodeCommMgr* nodeCommMgr);

  /** Destructor. */
  virtual ~NodeDatabase();

  /** Obtain number-of-node-descriptors (in function's return-value). Note that
      this remains 0 until after allocateNodeDescriptors() is called.
   */
  int getNumNodeDescriptors() const { return( nodePtrs_.size() ); };

  /** Obtain the set of nodeIDs. This is available anytime, and doesn't
      necessarily imply that there is a corresponding array of NodeDescriptors.
  */
  std::map<GlobalID,int>& getNodeIDs() { return( nodeIDs_ ); };

  /** Given a nodeID, return the corresponding node-descriptor. This function is
      only available after allocateNodeDescriptors() has been called 
      (returns -1 if allocateNodeDescriptors() hasn't been called).
      @param nodeID Input. nodeID for which a NodeDescriptor is required.
      @param node Output. NodeDescriptor corresponding to node 'nodeID'. If no
      corresponding node is found, or if allocatedNodeDescriptors() hasn't been
      called yet, node is not referenced.
      @return error-code 0 if successful, -1 if nodeID not found.
  */
  int getNodeWithID(GlobalID nodeID, const NodeDescriptor*& node) const;
  int getNodeWithID(GlobalID nodeID, NodeDescriptor*& node);

  /** Given a nodeNumber, return the corresponding node-descriptor.
      @param nodeNumber Input. nodeNumber for which a NodeDescriptor is
      required. This function returns -1 (and doesn't reference 'node') if
      synchronize() hasn't been called yet.
      @param node Output. NodeDescriptor corresponding to node 'nodeNumber'. 
      If no corresponding node is found, node is not referenced.
      @return error-code 0 if successful, -1 if node with nodeNumber not 
      present.
  */
  int getNodeWithNumber(int nodeNumber, const NodeDescriptor*& node) const;

  /** Given an equation-number, return the corresponding node-descriptor.
      @param eqnNumber Input. eqnNumber for which a NodeDescriptor is required.
      @param node Output. NodeDescriptor corresponding to eqnNumber. If no
      corresponding node is found, node is not referenced.
      @return error-code 0 if successful, -1 if node with eqnNumber not present.
  */
  int getNodeWithEqn(int eqnNumber, const NodeDescriptor*& node) const;

  /** Given an index i, return the i-th node-descriptor.
      @param i Input. Offset of requested NodeDescriptor.
      @param node Output. i-th NodeDescriptor.
      @return error-code. 0 if successful. -1 if i is out of bounds.
  */
  void getNodeAtIndex(int i, const NodeDescriptor*& node) const;
  void getNodeAtIndex(int i, NodeDescriptor*& node);

  /** Run through the locally-owned NodeDescriptors and count the number of
      nodal equations. (Returns 0 if the node-descriptors haven't
      been allocated or initialized yet.) The result is produced in the
      function's return-value.
  */
  int countLocalNodalEqns(int localRank);

  /** Return the number of internal node-descriptors that are
      locally-owned. This function runs through all currently allocated
      NodeDescriptors, counting the ones for which ownerProc()==localRank.
      The result is produced in the return-value. (Obviously returns 0 if the
      node-descriptors haven't been allocated or initialized yet.)
  */
  int countLocalNodeDescriptors(int localRank);

  /** Given a nodeID, return (in function's return-value) the index of that
      nodeID.
      @param nodeID Input
      @return offset if nodeID is found, -1 if nodeID not present.
   */
  int getIndexOfID(GlobalID nodeID) const;

  /** Initialization, add node with global identifier 'nodeID'. Only available
      before 'allocatedNodeDescriptors' has been called. Note that nodeIDs can
      be initialized multiple times using this function, but nodeID is only
      stored once in the internal array.
      @param nodeID Input.
      @return error-code 0 if successful. -2 if allocation failed.
  */
  int initNodeID(GlobalID nodeID);

  /** Initialization, add nodes with global identifiers 'nodeIDs'. Only
      available before 'allocatedNodeDescriptors' has been called. Note that
      nodeIDs can be initialized multiple times using this function, but each
      nodeID is only stored once in the internal array.
      @param nodeIDs Input.
      @return error-code 0 if successful. -2 if allocation failed.
  */
  int initNodeIDs(GlobalID* nodeIDs, int numNodes);

  /** Signal that node-descriptor initialization is complete. In other words,
      that the calling code has set fieldIDs, equation-numbers and owner-proc
      information on each node-descriptor. At this point, this class will run
      through the array of node-descriptors and set the nodeNumber for each
      node, as well as setting the cross-referencing info required to look up 
      nodes by equation-number.
      @param firstLocalNodeNumber Input. This will be the starting point for
      assigning node-numbers within this function.
      @param firstLocalEqn Input. This will be the starting point for assigning
      equation-numbers within this function.
      @param localRank Input. MPI rank of the local processor.
      @return error-code 0 if successful
  */
  int synchronize(int firstLocalNodeNumber,
		  int firstLocalEqn,
		  int localRank,
		  MPI_Comm comm);

  /** Given an equation-number, return the associated nodeNumber. i.e., return
      the nodeNumber of the node that contains eqnNumber.
  */
  int getAssociatedNodeNumber(int eqnNumber);

  /** Given an equation-number, return the associated fieldID. i.e., return
      the fieldID of the nodal solution field that contains eqnNumber.
  */
  int getAssociatedFieldID(int eqnNumber);

  /** Query whether synchronize() has been called. */
  bool isSynchronized() { return( synchronized_ ); };

 private:
  NodeDatabase(const NodeDatabase& src);
  NodeDatabase& operator=(const NodeDatabase& src);

  void deleteMemory();

  std::vector<NodeDescriptor*> nodePtrs_;

  std::vector<int> eqnNumbers_;  //eqnNumbers_ will be a sorted list of the
                                  //first global equation number at each node
                                  //in nodePtrs_.
                                  //the relationship between eqnNumbers_ and
  std::vector<int> eqnNodeIndices_;  //eqnNodeIndices_ is like this:
                                  //if eqn == eqnNumbers_[i], then
                                  //  nodePtrs_[eqnNodeIndices_[i]] points to
                                  //  the node with 'eqn'

  std::map<GlobalID,int> nodeIDs_; //nodeIDs_ maps node-ID to an index into
                                //the nodePtrs_ array of NodeDescriptors.

  std::map<int,int> nodeNumbers_;

  bool synchronized_;
  bool need_to_alloc_and_sync_;

  std::map<int,int>* fieldDB_;
  NodeCommMgr* nodeCommMgr_;

  int numLocalNodes_;
  int firstLocalNodeNumber_, lastLocalNodeNumber_;

  fei_Pool_alloc<NodeDescriptor> nodePool_;
};

#endif
