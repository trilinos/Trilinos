/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_iostream.hpp>
#include <fei_defs.h>

#include <fei_BlockDescriptor.hpp>

//====Constructor===============================================================
BlockDescriptor::BlockDescriptor()
 : blockID_(-1),
   numNodesPerElement_(0),
   numFieldsPerNode_(NULL),
   nodalFieldIDs_(NULL),
   fieldIDsAllocated_(false),
   numDistinctFields_(0),
   elemDofFieldIDs_(),
   interleaveStrategy_(0),
   lumpingStrategy_(0),
   numElements_(0),
   numElemDOFPerElement_(0),
   elemDOFEqnNumbers_(),
   numEqnsPerElement_(0),
   numBlkEqnsPerElement_(0),
   numActiveNodes_(0),
   totalNumEqns_(0)
{
   //There's nothing else for this constructor to do.
}

//====Destructor================================================================
BlockDescriptor::~BlockDescriptor() {

   destroyFieldArrays();

   numElemDOFPerElement_ = 0;
}

//==============================================================================
void BlockDescriptor::destroyFieldArrays() {
   if (numNodesPerElement_ == 0) return;

   for(int i=0; i<numNodesPerElement_; i++) {
      delete [] nodalFieldIDs_[i];
   }

   delete [] nodalFieldIDs_;
   nodalFieldIDs_ = NULL;
   delete [] numFieldsPerNode_;
   numFieldsPerNode_ = NULL;
   numNodesPerElement_ = 0;
}

//==============================================================================
int BlockDescriptor::setNumNodesPerElement(int numNodes)
{
  if (numNodes < 1) {
    return(-1);
  }

  destroyFieldArrays();

  numNodesPerElement_ = numNodes;

  numFieldsPerNode_ = new int[numNodesPerElement_];

  for(int i=0; i<numNodesPerElement_; i++) {
    numFieldsPerNode_[i] = 0;
  }
  return(0);
}

//==============================================================================
int* BlockDescriptor::fieldsPerNodePtr() {

   return(numFieldsPerNode_);
}

//==============================================================================
int BlockDescriptor::allocateFieldIDsTable()
{
  nodalFieldIDs_ = new int*[numNodesPerElement_];
  bool rowsAllZeroLength = true;

  for(int i=0; i<numNodesPerElement_; i++) {
    if (numFieldsPerNode_[i] > 0) {
      nodalFieldIDs_[i] = new int[numFieldsPerNode_[i]];
      rowsAllZeroLength = false;
    }
    else nodalFieldIDs_[i] = NULL;
  }

  if (rowsAllZeroLength || numNodesPerElement_ == 0) {
    fei::console_out() << "BlockDescriptor::allocateFieldIDsTable: ERROR, all rows of"
	 << " fieldIDs table have zero length. Set fieldsPerNode entries"
	 << " first." << FEI_ENDL;
    return(-1);
  }

  fieldIDsAllocated_ = true;
  return(0);
}

//==============================================================================
bool BlockDescriptor::containsField(int fieldID) {
//
//This function will mostly be called by the BASE_FEI function for
//getting solutions to return to the user.
//
//For cases where each of the nodes in an element have the same fields,
//this function will be quite fast.
//
//It will be slow for cases where there are quite a few nodes per element
//and the different nodes have different solution fields. (i.e., the search
//below has to step through most of the fieldIDs table before finding the
//fieldID in question.
//
//In general though, this function won't be called if the fieldID isn't
//associated with ANY node in this block, because the calling code can first
//query the node to find out if IT is associated with this block. And if the
//node is associated with this block, then the node's fields usually will be
//also, unless the node lies on a block boundary and 'fieldID' is only in
//the other block.
//
   for(int i=0; i<numNodesPerElement_; i++) {
      for(int j=0; j<numFieldsPerNode_[i]; j++) {
         if (nodalFieldIDs_[i][j] == fieldID) return(true);
      }
   }

   return(false);
}

//==============================================================================
int BlockDescriptor::setElemDofFieldIDs(int numFields, const int* fieldIDs)
{
   if (numFields <= 0) {
      elemDOFEqnNumbers_.resize(0);
      return(0);
   }

   for(int i=0; i<numFields; i++) elemDofFieldIDs_.push_back(fieldIDs[i]);

   elemDOFEqnNumbers_.resize(numElements_);

   for(int j=0; j<numElements_; j++) {
      elemDOFEqnNumbers_[j] = -1;
   }

   return(0);
}
