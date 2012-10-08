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
   offsetsIntoField_(NULL),
   prescribed_values_(NULL)
{
}

//==============================================================================
BCNodeSet::~BCNodeSet() {
   deleteMemory();
}

//==============================================================================
void BCNodeSet::deleteMemory() {
   delete [] offsetsIntoField_;
   delete [] prescribed_values_;

   delete [] nodeIDs_;

   numNodes_ = 0;
}

