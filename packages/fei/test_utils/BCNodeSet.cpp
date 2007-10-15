/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_defs.h>
#include <test_utils/BCNodeSet.hpp>
#include <cstdlib>

//==============================================================================
BCNodeSet::BCNodeSet()
 : numNodes_(0),
   nodeIDs_(NULL),
   fieldID_(-1),
   alpha_(NULL),
   beta_(NULL),
   gamma_(NULL)
{
}

//==============================================================================
BCNodeSet::~BCNodeSet() {
   deleteMemory();
}

//==============================================================================
void BCNodeSet::deleteMemory() {
   for(int i=0; i<numNodes_; i++) {
      delete [] alpha_[i];
      delete [] beta_[i];
      delete [] gamma_[i];
   }
   delete [] alpha_;
   delete [] beta_;
   delete [] gamma_;

   delete [] nodeIDs_;

   numNodes_ = 0;
}

