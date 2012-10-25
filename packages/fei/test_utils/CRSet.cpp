/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_defs.h>
#include <test_utils/CRSet.hpp>
#include <cstdlib>

//==============================================================================
CRSet::CRSet()
 : crID_(-1),
   numNodes_(0),
   nodeIDs_(NULL),
   fieldIDs_(NULL),
   weights_(NULL),
   values_(NULL),
   penValues_(NULL)
{
}

//==============================================================================
CRSet::~CRSet() {
   deleteMemory();
}

//==============================================================================
void CRSet::deleteMemory() {

   for(int j=0; j<1; j++) {
      delete [] nodeIDs_[j];
   }

   delete [] weights_;
   delete [] values_;
   delete [] penValues_;
   delete [] fieldIDs_;
   delete [] nodeIDs_;

   numNodes_ = 0;
}

