/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_Record.hpp>

#include <fei_FieldMask.hpp>

fei::Record::Record()
  : isInLocalSubdomain_(false),
    ID_(-1),
    number_(-1),
    fieldMask_(NULL),
    offsetIntoEqnNumbers_(0),
    ownerProc_(-1),
    hasSlaveDof_(false)
{
}

fei::Record::~Record()
{
}

int fei::Record::deepCopy(const Record& rcd)
{
  ID_ = rcd.ID_;
  number_ = rcd.number_;
  fieldMask_ = rcd.fieldMask_;
  offsetIntoEqnNumbers_ = rcd.offsetIntoEqnNumbers_;
  ownerProc_ = rcd.ownerProc_;
  hasSlaveDof_ = rcd.hasSlaveDof_;

  return(0);
}

