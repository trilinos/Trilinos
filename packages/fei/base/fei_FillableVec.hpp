/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_FillableVec_hpp_
#define _fei_FillableVec_hpp_

#include "fei_macros.hpp"
#include "fei_Pool_alloc.hpp"
#include <map>

namespace fei {

class FillableVec {
 public:
  FillableVec();
  virtual ~FillableVec();

  void addEntry(int index, double coef);

  void addEntries(unsigned numEntries,
                  const double* coefs,
                  const int* indices);

  void putEntry(int index, double coef);

  void putEntries(unsigned numEntries,
                  const double* coefs,
                  const int* indices);

  void setValues(double value);

  unsigned size() const;

  bool hasEntry(int index) const;

  /** Return coef for index if index is present, otherwise throw.
  */
  double getEntry(int index) const;

  void removeEntry(int index);

  void clear();

  typedef
    std::map<int,double,std::less<int>,
     fei_Pool_alloc<std::pair<const int,int> > > feipoolmap;

  typedef feipoolmap::iterator iterator;
  typedef feipoolmap::const_iterator const_iterator;

  iterator begin() {return vecdata_.begin();}
  iterator end() {return vecdata_.end();}

  const_iterator begin() const {return vecdata_.begin();}
  const_iterator end() const {return vecdata_.end();}

  bool operator==(const FillableVec& rhs) const;

  bool operator!=(const FillableVec& rhs) const;

 private:
  feipoolmap vecdata_;
};//class FillableVec

inline void
FillableVec::addEntry(int index, double coef)
{
  feipoolmap::iterator iter = vecdata_.lower_bound(index);
  if (iter == vecdata_.end() || iter->first != index) {
    vecdata_.insert(iter, std::make_pair(index, coef));
  }
  else {
    iter->second += coef;
  }
}

inline void
FillableVec::putEntry(int index, double coef)
{
  feipoolmap::iterator iter = vecdata_.lower_bound(index);
  if (iter == vecdata_.end() || iter->first != index) {
    vecdata_.insert(iter, std::make_pair(index, coef));
  }
  else {
    iter->second = coef;
  }
}

inline
bool
FillableVec::operator==(const FillableVec& rhs) const
{
  return vecdata_ == rhs.vecdata_;
}

inline
bool
FillableVec::operator!=(const FillableVec& rhs) const
{
  return vecdata_ != rhs.vecdata_;
}

}//namespace fei

#endif

