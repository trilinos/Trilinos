/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_FillableMat_hpp_
#define _fei_FillableMat_hpp_

#include "fei_FillableVec.hpp"

namespace fei {

class FillableMat {
 public:
  FillableMat();
  virtual ~FillableMat();

  void zero();

  void sumInCoef(int row, int col, double coef);
  void putCoef(int row, int col, double coef);

  void sumInRow(int row, const int* cols, const double* coefs, unsigned len);
  void putRow(int row, const int* cols, const double* coefs, unsigned len);

  unsigned getNumRows() const;

  bool hasRow(int row) const;

  /** If row is not present, throw std::runtime_error. */
  FillableVec* getRow(int row);

  typedef std::map<int, FillableVec*, std::less<int>,
                fei_Pool_alloc<std::pair<const int,FillableVec*> > > feipoolmat;

  typedef feipoolmat::iterator iterator;
  typedef feipoolmat::const_iterator const_iterator;

  iterator begin() {return matdata_.begin();}
  iterator end() {return matdata_.end();}

  const_iterator begin() const {return matdata_.begin();}
  const_iterator end() const {return matdata_.end();}

  void clear();

 private:
  feipoolmat matdata_;
  fei_Pool_alloc<FillableVec> vecpool_;
}; //class FillableMat
}//namespace fei

#endif

