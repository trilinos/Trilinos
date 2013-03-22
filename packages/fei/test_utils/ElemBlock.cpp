/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_defs.h>
#include <test_utils/ElemBlock.hpp>
#include <cstdlib>

//==============================================================================
ElemBlock::ElemBlock()
 : blockID_(0),
   numElements_(0),
   numNodesPerElement_(0),
   numFieldsPerNode_(NULL),
   nodalFieldIDs_(NULL),
   elemIDs_(NULL),
   elemConn_(NULL),
   numStiffRows_(0),
   elemFormat_(0),
   elemStiff_(NULL),
   elemLoad_(NULL),
   numElemDOF_(0),
   elemDOFFieldIDs_(NULL),
   interleaveStrategy_(0),
   lumpingStrategy_(0)
{
}

//==============================================================================
ElemBlock::~ElemBlock() {
   deleteMemory();
}

//==============================================================================
void ElemBlock::deleteMemory() {
   for(int i=0; i<numElements_; i++) {
      for(int j=0; j<numStiffRows_; j++) {
         delete [] elemStiff_[i][j];
      }
      delete [] elemStiff_[i];
      delete [] elemLoad_[i];
      delete [] elemConn_[i];
   }

   for(int j=0; j<numNodesPerElement_; j++) {
      delete [] nodalFieldIDs_[j];
   }
   delete [] nodalFieldIDs_;
   delete [] numFieldsPerNode_;

   delete [] elemStiff_;
   delete [] elemLoad_;
   delete [] elemConn_;
   delete [] elemIDs_;

   if (numElemDOF_ > 0) {
      delete [] elemDOFFieldIDs_;
      numElemDOF_ = 0;
   }

   numElements_ = 0;
   numNodesPerElement_ = 0;
}

