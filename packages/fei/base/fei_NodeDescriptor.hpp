/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/

#ifndef _NodeDescriptor_hpp_
#define _NodeDescriptor_hpp_


#include <fei_macros.hpp>

#include <fei_ArrayUtils.hpp>
#include <iostream>
/**
  The NodeDescriptor class holds the information that the FEI implementation
  needs to know about the nodes in the finite element problem:
  
     Global node identifier
     number of nodal degrees-of-freedom
     list of associated field identifiers, with their (global) equation numbers
     (global) blk-eqn number
     which processor is this node's owner
     list of (local) blocks that contain this node
  
  Note: 'block' is used in two contexts here. There are element-blocks, and
    there are block-equations. Element-blocks refer to the blocks of elements
    in the finite-element problem, with all elements in a block containing the
    same number of solution fields per node, etc. Block-equations refer to the
    small dense sub-blocks of a block-entry sparse matrix. Each node is 
    associated with a number of element-blocks, and each node has exactly one
    associated global 0-based block-equation number.
*/

class NodeDescriptor {
 public:
   NodeDescriptor();

   NodeDescriptor(const NodeDescriptor& src)
    : nodeID_(src.nodeID_), nodeNumber_(src.nodeNumber_),
      numNodalDOF_(0), fieldIDList_(NULL),
      fieldEqnNumbers_(NULL), numFields_(0), blkEqnNumber_(0),
      ownerProc_(src.ownerProc_), blockList_()
   {}

   virtual ~NodeDescriptor();

   GlobalID getGlobalNodeID() const {return(nodeID_);}
   void setGlobalNodeID(GlobalID node) {nodeID_ = node;}

   int getNodeNumber() const {return(nodeNumber_);}
   void setNodeNumber(int nn) {nodeNumber_ = nn;};

   int getBlkEqnNumber() const {return(blkEqnNumber_);}
   void setBlkEqnNumber(int blkEqn) {blkEqnNumber_ = blkEqn;}

   int getNumNodalDOF() const {return(numNodalDOF_);}
   void setNumNodalDOF(int dof) {numNodalDOF_ = dof;}

   void addField(int fieldID);
   void setFieldEqnNumber(int fieldID, int eqn);

   int getNumFields() const {return(numFields_);}
   const int* getFieldIDList() const {return(fieldIDList_);}
   const int* getFieldEqnNumbers() const {return(fieldEqnNumbers_);}

   /** Given a fieldID, return the first equation number associated with that
       field at this node.
       @param fieldID
       @param eqnNumber
       @return false if fieldID is not present at this node
   */
   bool getFieldEqnNumber(int fieldID, int& eqnNumber) const;

   /** Given an equation-number, get the associated fieldID and offset-into-field.
     throws an exception if the equation-number is not associated with this node.
   */
   void getFieldID(int eqnNumber, int& fieldID, int& offset_into_field) const;

   bool operator==(const NodeDescriptor& nd) const
     { return( nodeID_ == nd.nodeID_ ); }

   bool operator!=(const NodeDescriptor& nd) const
     { return( nodeID_ != nd.nodeID_ ); }

   bool operator<(const NodeDescriptor& nd) const
     { return( nodeID_ < nd.nodeID_ ); }

   bool operator>(const NodeDescriptor& nd) const
     { return( nodeID_ > nd.nodeID_ ); }

   int getOwnerProc() const {return(ownerProc_);}
   void setOwnerProc(int proc) {ownerProc_ = proc;}

   void addBlockIndex(unsigned blk_idx)
     { fei::sortedListInsert(blk_idx, blockList_); }

   size_t getNumBlocks() const {return blockList_.size();}
   const std::vector<unsigned>& getBlockIndexList() const {return(blockList_);}
   bool hasBlockIndex(unsigned blk_idx) const;

 private:
   NodeDescriptor& operator=(const NodeDescriptor& src);

   void allocFieldLists();
   void allocBlockList();

   GlobalID nodeID_;

   int nodeNumber_;

   int numNodalDOF_;      //total number of nodal degrees-of-freedom

   int* fieldIDList_;     //list of field identifiers
   int* fieldEqnNumbers_; //list of starting (global) equation numbers.
                          //fields can consist of more than one scalar (and
                          //have more than one associated equation), this
                          //is the first equation number
   int numFields_;

   int blkEqnNumber_;

   int ownerProc_;        //processor that owns the equations for this node

   std::vector<unsigned> blockList_;  //indexes of blocks that contain this node
};

#endif

