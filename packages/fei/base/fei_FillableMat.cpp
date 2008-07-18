/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_FillableMat.hpp"
#include "fei_Exception.hpp"

namespace fei {

FillableMat::FillableMat()
 : matdata_(),
   vecpool_()
{
}

FillableMat::~FillableMat()
{
  feipoolmat::iterator
    iter = matdata_.begin(), end = matdata_.end();
  for(; iter!=end; ++iter) {
    vecpool_.destroy(iter->second);
    vecpool_.deallocate(iter->second, 1);
  }
}

void
FillableMat::zero()
{
  feipoolmat::iterator
    iter = matdata_.begin(), iter_end = matdata_.end();

  for(; iter != iter_end; ++iter) {
    iter->second->zero();
  }
}

void
FillableMat::sumInCoef(int row, int col, double coef)
{
  static FillableVec dummy;

  feipoolmat::iterator iter = matdata_.lower_bound(row);

  if (iter == matdata_.end() || iter->first != row) {
    FillableVec* vptr = vecpool_.allocate(1);
    vecpool_.construct(vptr, dummy);
    iter = matdata_.insert(iter, std::make_pair(row, vptr));
  }

  iter->second->addEntry(col, coef);
}

void
FillableMat::putCoef(int row, int col, double coef)
{
  static FillableVec dummy;
  feipoolmat::iterator iter = matdata_.lower_bound(row);

  if (iter == matdata_.end() || iter->first != row) {
    FillableVec* vptr = vecpool_.allocate(1);
    vecpool_.construct(vptr, dummy);
    iter = matdata_.insert(iter, std::make_pair(row, vptr));
  }

  iter->second->putEntry(col, coef);
}

void
FillableMat::sumInRow(int row, const int* cols, const double* coefs,
                      unsigned len)
{
  static FillableVec dummy;
  feipoolmat::iterator iter = matdata_.lower_bound(row);

  if (iter == matdata_.end() || iter->first != row) {
    FillableVec* vptr = vecpool_.allocate(1);
    vecpool_.construct(vptr, dummy);
    iter = matdata_.insert(iter, std::make_pair(row, vptr));
  }

  for(unsigned i=0; i<len; ++i) {
    iter->second->addEntry(cols[i], coefs[i]);
  }
}

void
FillableMat::putRow(int row, const int* cols, const double* coefs,
                    unsigned len)
{
  static FillableVec dummy;
  feipoolmat::iterator iter = matdata_.lower_bound(row);

  if (iter == matdata_.end() || iter->first != row) {
    FillableVec* vptr = vecpool_.allocate(1);
    vecpool_.construct(vptr, dummy);
    iter = matdata_.insert(iter, std::make_pair(row, vptr));
  }

  for(unsigned i=0; i<len; ++i) {
    iter->second->putEntry(cols[i], coefs[i]);
  }
}

unsigned
FillableMat::getNumRows() const
{
  return matdata_.size();
}

bool
FillableMat::hasRow(int row) const
{
  feipoolmat::const_iterator iter = matdata_.find(row);
  return iter != matdata_.end();
}

FillableVec*
FillableMat::getRow(int row)
{
  feipoolmat::iterator iter = matdata_.find(row);

  if (iter == matdata_.end()) {
    throw fei::Exception("fei::FillableMat: row not found.");
  }

  return iter->second;
}

void
FillableMat::clear()
{
  feipoolmat::iterator
    iter = matdata_.begin(), end = matdata_.end();
  for(; iter!=end; ++iter) {
    vecpool_.destroy(iter->second);
    vecpool_.deallocate(iter->second, 1);
  }

  matdata_.clear();
}

}//namespace fei

