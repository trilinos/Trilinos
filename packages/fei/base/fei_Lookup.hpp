#ifndef _fei_Lookup_hpp_
#define _fei_Lookup_hpp_

#include <fei_defs.h>

/**  
  This interface is intended to be used by a LinearSystemCore implementation
  or a FiniteElementData implementation to look up various
  information about the structure of the finite-element problem being assembled
  into LinearSystemCore by a FEI implementation.
  
  Background:
    The finite element problem consists of a set of element-blocks, each of
    which contains a set of elements. Each element has a list of connected
    nodes, and each node has a set of fields. Each field consists of 1 to
    several scalar quantities. Each of those scalar quantities corresponds
    to an equation in the linear system. Exception: some fields may not be
    solution-fields. This is indicated by a negative fieldID. There are no
    equations corresponding to fields with negative fieldIDs. Data that is
    passed in associated with negative fieldIDs is probably coordinate or
    nullspace data, or other data passed from the application to the solver
    (note that this non-solution-field data won't appear in element-matrices,
    it will be passed separately through a special function for supplying
    nodal data).
  
    elem-block IDs and field IDs are application-provided numbers, and no
    assumption may be made regarding their order, contiguousness, etc.
  
    Equation numbers are assigned to node/field pairs by the FEI implementation,
    and are also globally unique, and non-shared. Each equation resides on
    exactly one processor. Equation numbers are 0-based.
  
  
  NOTES: 1. functions that return an equation number, or a size (e.g., 
         num-equations, num-fields, etc.) may indicate an error or 'not-found'
         condition by returning a negative number. Functions that return a
         pointer may indicate an error by returning NULL.
*/

class Lookup {
 public:
  /** Destructor. */
   virtual ~Lookup(){};

  /** Get the (global) number of fields defined for the problem. A field may
    consist of 1 to several scalar quantities. Examples include pressure,
    displacement (a 3-vector in 3D), etc.
 */
   virtual int getNumFields() = 0;

  /** Given a fieldID, obtain the associated fieldSize (number of scalar components)
     @param fieldID integer identifier for a field
  */
   virtual int getFieldSize(int fieldID) = 0;


  /**
    Return a pointer to a list (of length numFields) of the fieldIDs 
  */
   virtual const int* getFieldIDsPtr() = 0;

  /**
    Return a pointer to a list (of length numFields) of the fieldSizes 
  */
   virtual const int* getFieldSizesPtr() = 0;


  /** Return the number of element-blocks in the (local) finite-element
    problem.
  */
   virtual int getNumElemBlocks() = 0;


  /** Return a pointer to the list (of length numElemBlocks) containing the
   element-block identifiers for the (local) finite-element problem.
  */
   virtual const GlobalID* getElemBlockIDs() = 0;


  /** Given a blockID, provide several pieces of element-block information.
     @param interleaveStrategy element-equation ordering: 0 => node-major,
                       1 => field-major
   @param lumpingStrategy element-matrices may be lumped if they're mass 
          matrices, 0 => not lumped, 1 => lumped
   @param numElemDOF number of element-dof at each element in this block
   @param numElements number of elements in this block
   @param numNodesPerElem number of nodes connected to each element in this 
       block
   @param numEqnsPerElem number of scalar equations at each element in this
        block
  */
   virtual void getElemBlockInfo(GlobalID blockID,
                         int& interleaveStrategy, int& lumpingStrategy,
                         int& numElemDOF, int& numElements,
                         int& numNodesPerElem, int& numEqnsPerElem) = 0;


   /** Given a blockID, return a pointer to a list (of length numNodesPerElem)
      of numFieldsPerNode.
     @param blockID identifier of the elem-block in question
   */

   virtual const int* getNumFieldsPerNode(GlobalID blockID) = 0;


   /** Given a blockID, return a pointer to a table,
     (num-rows == numNodesPerElem, row-length[i] == fieldsPerNode[i])
     containing the fieldIDs at each node of elements in that element-block.
     @param blockID identifier of the elem-block in question
   */
   virtual const int* const* getFieldIDsTable(GlobalID blockID) = 0;


   /** Given a nodeNumber/fieldID pair, this function returns the first global
     (0-based) equation number associated with that nodeNumber/fieldID pair.
      @param nodeNumber
      @param fieldID
   */
   virtual int getEqnNumber(int nodeNumber, int fieldID) = 0;


   /** Given a global (0-based) equation number, return the node-number.
       @param eqnNumber
   */
   virtual int getAssociatedNodeNumber(int eqnNumber) = 0;


   /** Given a global (0-based) equation number, return the fieldID.
       @param eqnNumber
   */
   virtual int getAssociatedFieldID(int eqnNumber) = 0;


   /** Given a nodeNumber, determine whether that node is connected to a local
       element.
       @param nodeNumber
       @return true if nodeNumber is connected to a local element, false
       otherwise.
   */
   virtual bool isInLocalElement(int nodeNumber) = 0;

   /** Given a nodeNumber, return the number of subdomains that contain this
       node. subdomains correspond to processors. The number of subdomains that
       contain a node does not always equal the number of processors that share
       a node. There are two kinds of "sharing" -- the "normal" kind, where a
       node is shared because it is connected to elements that reside on more
       than one processor, and the "wierd" kind where nodes are considered 
       shared simply because of being in cross-processor constraints. This
       function describes how many processors share this node in the "normal"
       sense. Thus, in general, this relationship holds:
       getNumSubdomains(nodeNum) <= getNumSharingProcs(nodeNum)
       @param nodeNumber
       @return numSubdomains
   */
   virtual int getNumSubdomains(int nodeNumber) = 0;

   /** Given a nodeNumber, return a list of the subdomains that contain this
       node. This is the list of subdomains that's counted by the method
       'getNumSubdomains' above.
       @param nodeNumber
       @return pointer to list of subdomains (processor ranks). NULL if the
       specified nodeNumber is not found.
   */
   virtual int* getSubdomainList(int nodeNumber) = 0;

   /** Return the number of local nodes that are shared by multiple processors
    */
   virtual int getNumSharedNodes() = 0;

   /** Return a pointer to the list of shared nodeNumbers 
    */
   virtual const int* getSharedNodeNumbers() = 0;

   /** Given a shared nodeNumber, return a pointer to the list of sharing procs.
     @param nodeNumber The subject of the query. Function returns NULL if 
        'nodeNumber' is not a shared node.
   */
   virtual const int* getSharedNodeProcs(int nodeNumber) = 0;

   /** Given a shared nodeNumber, return the number of processors that share it.
     @param nodeNumber Function returns -1 if 'nodeNumber' is not a shared
         node.
   */
   virtual int getNumSharingProcs(int nodeNumber) = 0;

   //- given a blk-eqn and a pt-eqn, return the pt-eqn's offset into the blk-eqn
   //     (i.e., distance from the 'beginning' of the blk-eqn)
   //- given a blk-eqn, return the 'size', or number of pt-eqns corresponding
   //     to it.
   /** Query whether a pt-eqn corresponds exactly to a blk-eqn. in other words,
       is pt-eqn the first point equation in a block-equation.
      @param ptEqn
   */
   virtual bool isExactlyBlkEqn(int ptEqn) = 0;

   /** Given a pt-eqn, return the corresponding blk-eqn.
     @param ptEqn
   */
   virtual int ptEqnToBlkEqn(int ptEqn) = 0;

   /** Given a blk-eqn and a pt-eqn, return the pt-eqn's offset into the blk-eqn
       (i.e., distance from the 'beginning' of the blk-eqn)
     @param blkEqn
     @param ptEqn
   */
   virtual int getOffsetIntoBlkEqn(int blkEqn, int ptEqn) = 0;

   /** Given a blk-eqn, return the 'size', or number of pt-eqns corresponding
          to it.
     @param blkEqn
   */
   virtual int getBlkEqnSize(int blkEqn) = 0;
};

#endif

