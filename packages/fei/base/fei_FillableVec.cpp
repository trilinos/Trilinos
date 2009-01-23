/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_FillableVec.hpp"

namespace fei {

FillableVec::FillableVec()
 : vecdata_()
{
}

FillableVec::~FillableVec()
{
}

void
FillableVec::addEntries(unsigned numEntries,
                  const double* coefs,
                  const int* indices)
{
  for(unsigned i=0; i<numEntries; ++i) {
    addEntry(indices[i], coefs[i]);
  }
}

void
FillableVec::putEntries(unsigned numEntries,
                  const double* coefs,
                  const int* indices)
{
  for(unsigned i=0; i<numEntries; ++i) {
    putEntry(indices[i], coefs[i]);
  }
}

void
FillableVec::setValues(double value)
{
  if (size() > 0) {
    feipoolmap::iterator iter = vecdata_.begin(), iter_end = vecdata_.end();
    for(; iter != iter_end; ++iter) {
      iter->second = value;
    }
  }
}

void
FillableVec::clear()
{
  vecdata_.clear();
}

unsigned
FillableVec::size() const
{
  return vecdata_.size();
}

bool
FillableVec::hasEntry(int index) const
{
  feipoolmap::const_iterator iter = vecdata_.find(index);
  return iter != vecdata_.end();
}

double
FillableVec::getEntry(int index) const
{
  feipoolmap::const_iterator iter = vecdata_.find(index);
  if (iter == vecdata_.end()) {
    throw std::runtime_error("FillableVec::getEntry: index not found.");
  }

  return iter->second;
}

void
FillableVec::removeEntry(int index)
{
  feipoolmap::iterator iter = vecdata_.find(index);
  if (iter != vecdata_.end()) {
    vecdata_.erase(iter);
  }
}

}//namespace fei

