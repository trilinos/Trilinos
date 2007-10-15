/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_iostream.hpp>
#include <fei_defs.h>

#include <fei_BCRecord.hpp>

//==========================================================================
BCRecord::BCRecord()
 : idType_(0), nodeID_(-1), fieldID_(-999), fieldSize_(0), coef_(NULL)
{
}

//==============================================================================
int BCRecord::init(GlobalID nodeID, int fieldID, int fieldSize, double* coefPtr)
{
  if (fieldSize <= 0) {
    FEI_CERR << "fei: BCRecord ERROR, fieldSize <= 0." << FEI_ENDL;
    return(-1);
  }
  nodeID_ = nodeID;
  fieldID_ = fieldID;
  fieldSize_ = fieldSize;
  coef_ = coefPtr;
  return(0);
}

//==============================================================================
void BCRecord::setDoubleList(int bcComponent, const double* input)
{
  int offset = bcComponent*fieldSize_;
  for(int i=0; i<fieldSize_; i++) {
    coef_[offset+i] = input[i];
  }
}

