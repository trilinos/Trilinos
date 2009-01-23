/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_defs.h>

#include <fei_NodeDescriptor.hpp>

//======Constructor=============================================================
NodeDescriptor::NodeDescriptor()
 : nodeID_((GlobalID)-1),
   nodeNumber_(-1),
   numNodalDOF_(0),
   fieldIDList_(NULL),
   fieldEqnNumbers_(NULL),
   numFields_(0),
   blkEqnNumber_(0),
   ownerProc_(-1),
   blockList_(NULL),
   numBlocks_(0)
{
   //There's nothing for this constructor to do, apart from the
   //above initializations.
}

//======Destructor==============================================================
NodeDescriptor::~NodeDescriptor() {
  delete [] fieldIDList_;
  delete [] fieldEqnNumbers_;
  numFields_ = 0;

  delete [] blockList_;
  numBlocks_ = 0;
}

//==============================================================================
void NodeDescriptor::addField(int fieldID) {
//
//Add a field identifier to this node, ONLY if that field identifier
//is not already present.
//
//If fieldID is added, lengthen the corresponding list for equation numbers.
//

   int tmp = numFields_;
   int allocLen = numFields_;
   int index = snl_fei::sortedListInsert(fieldID, fieldIDList_, numFields_,
                                         allocLen);

   //index is the position at which fieldID was inserted, or found

   //if tmp < numFields_ then fieldID wasn't already present
   if (tmp < numFields_) {
      //
      //if the length of fieldIDList_ changed, let's lengthen the 
      //fieldEqnNumbers_ list.
      //fieldEqnNumbers_ will have an empty position 'index', which we'll set to
      //-99 for now. The calling code (BASE_FEI) will set the fieldEqnNumber for
      //this fieldID using setFieldEqnNumber(...).
      //

      allocLen = numFields_ - 1;
      snl_fei::listInsert(-99, index, fieldEqnNumbers_, tmp, allocLen);
   }
}

//==============================================================================
void NodeDescriptor::setFieldEqnNumber(int fieldID, int eqn) {
//
//Set the equation number corresponding to fieldID. fieldID must
//already have been added to this node using the addField function.
//If it was already added, then the fieldEqnNumbers_ list was lengthened
//appropriately, with an empty spot left for this eqn number.
//
   int insert = -1;
   int index = snl_fei::binarySearch(fieldID, fieldIDList_,
					    numFields_, insert);

   if (index < 0) {
      return;
   }

   fieldEqnNumbers_[index] = eqn;
}

//==============================================================================
bool NodeDescriptor::getFieldEqnNumber(int fieldID, int& eqnNumber)
{
   int insert = -1;
   int index = snl_fei::binarySearch(fieldID, fieldIDList_,
					    numFields_, insert);

   if (index < 0) {
      return(false);
   }

   eqnNumber = fieldEqnNumbers_[index];
   return(true);
}

//==============================================================================
bool NodeDescriptor::containedInBlock(GlobalID blk)
{
  //return true if this node is contained in element-block 'blk'.

   int insert;
   int index = snl_fei::binarySearch(blk, blockList_, numBlocks_, insert);
   if (index >= 0) return(true);
   else return(false);
}

