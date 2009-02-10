/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <fei_PatternDescriptor.hpp>

//------------------------------------------------------------------------------
PatternDescriptor::PatternDescriptor()
 : patternID_(-1),
   fieldsPerRow_(),
   rowFieldIDs_(NULL),
   fieldsPerCol_(),
   colFieldIDs_(NULL),
   interleaveStrategy_(-1)
{
}

//------------------------------------------------------------------------------
PatternDescriptor::~PatternDescriptor()
{
   if (fieldsPerRow_.size() > 0) delete [] rowFieldIDs_;
   if (fieldsPerCol_.size() > 0) delete [] colFieldIDs_;
}

//------------------------------------------------------------------------------
int PatternDescriptor::setNumRowIDs(int num)
{
//This function is responsible for adjusting the size of this class' internal
//data structures. If a resize or an allocation fails, return -1. Otherwise,
//return 0.
//
   if (fieldsPerRow_.size() > 0) {
      delete [] rowFieldIDs_;
   }

   fieldsPerRow_.resize(num);

   if (num == 0) return(0);

   rowFieldIDs_ = new std::vector<int>[num];

   return(0);
}

//------------------------------------------------------------------------------
int PatternDescriptor::setNumColIDsPerRow(int num)
{
//This function is responsible for adjusting the size of this class' internal
//data structures. If a resize or an allocation fails, return -1. Otherwise,
//return 0.
//
   if (fieldsPerCol_.size() > 0) {
      delete [] colFieldIDs_;
   }

   fieldsPerCol_.resize(num);

   if (num == 0) return(0);

   colFieldIDs_ = new std::vector<int>[num];

   return(0);
}

