/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_CSVec.hpp"
#include <algorithm>

namespace fei {

CSVec::CSVec(const FillableVec& invec)
 : indices_(invec.size()),
   coefs_(invec.size())
{
  operator=(invec);
}

CSVec::CSVec(unsigned sz)
 : indices_(sz, 0),
   coefs_(sz, 0.0)
{
}

CSVec::~CSVec()
{
}

CSVec&
CSVec::operator=(const FillableVec& invec)
{
  indices_.resize(invec.size());
  coefs_.resize(invec.size());

  FillableVec::const_iterator iter = invec.begin(), iter_end = invec.end();

  unsigned i=0;
  for(; iter != iter_end; ++iter, ++i) {
    indices_[i] = iter->first;
    coefs_[i] = iter->second;
  }

  return *this;
}

void add_entry(CSVec& vec, int eqn, double coef)
{
  std::vector<int>& v_ind = vec.indices();
  std::vector<double>& v_coef = vec.coefs();

  std::vector<int>::iterator
    iter = std::lower_bound(v_ind.begin(), v_ind.end(), eqn);

  size_t offset = iter - v_ind.begin();

  if (iter == v_ind.end() || *iter != eqn) {
    v_ind.insert(iter, eqn);
    v_coef.insert(v_coef.begin()+offset, coef);
  }
  else {
    v_coef[offset] += coef;
  }
}

void put_entry(CSVec& vec, int eqn, double coef)
{
  std::vector<int>& v_ind = vec.indices();
  std::vector<double>& v_coef = vec.coefs();

  std::vector<int>::iterator
    iter = std::lower_bound(v_ind.begin(), v_ind.end(), eqn);

  size_t offset = iter - v_ind.begin();

  if (iter == v_ind.end() || *iter != eqn) {
    v_ind.insert(iter, eqn);
    v_coef.insert(v_coef.begin()+offset, coef);
  }
  else {
    v_coef[offset] = coef;
  }
}

}//namespace fei

