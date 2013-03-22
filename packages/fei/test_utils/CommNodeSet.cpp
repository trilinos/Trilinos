/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_defs.h>
#include <test_utils/CommNodeSet.hpp>
#include <cstdlib>

//==============================================================================
CommNodeSet::CommNodeSet()
 : numNodes_(0),
   nodeIDs_(NULL),
   procs_(NULL),
   procsPerNode_(NULL)
{
}

//==============================================================================
CommNodeSet::~CommNodeSet() {
   deleteMemory();
}

//==============================================================================
void CommNodeSet::deleteMemory() {
   for(int i=0; i<numNodes_; i++) {
      delete [] procs_[i];
   }

   delete [] procs_;
   delete [] procsPerNode_;
   delete [] nodeIDs_;
   numNodes_ = 0;
}

