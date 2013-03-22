/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_CSVec.hpp>
#include <algorithm>

namespace fei {

CSVec::CSVec(unsigned sz)
 : indices_(sz, 0),
   coefs_(sz, 0.0)
{
}

CSVec::~CSVec()
{
}

CSVec&
CSVec::operator=(const CSVec& invec)
{
  indices_ = invec.indices_;
  coefs_ = invec.coefs_;

  return *this;
}

void add_entries(CSVec& vec, int num, const int* eqns, const double* coefs)
{
  for(int i=0; i<num; ++i) add_entry(vec, eqns[i], coefs[i]);
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

double get_entry(const CSVec& vec, int eqn)
{
  const std::vector<int>& v_ind = vec.indices();
  const std::vector<double>& v_coef = vec.coefs();

  if (vec.size() == 0) {
    throw std::runtime_error("get_entry error, CSVec is empty");
  }

  std::vector<int>::const_iterator
    iter = std::lower_bound(v_ind.begin(), v_ind.end(), eqn);

  if (iter == v_ind.end()) {
    throw std::runtime_error("get_entry error, entry not found.");
  }

  return v_coef[iter - v_ind.begin()];
}

void remove_entry(CSVec& vec, int eqn)
{
  std::vector<int>& v_ind = vec.indices();
  std::vector<double>& v_coef = vec.coefs();

  std::vector<int>::iterator
    iter = std::lower_bound(v_ind.begin(), v_ind.end(), eqn);

  if (iter != v_ind.end() && *iter == eqn) {
    size_t offset = iter - v_ind.begin();
    v_ind.erase(iter);

    std::vector<double>::iterator coef_iter = v_coef.begin()+offset;
    v_coef.erase(coef_iter);
  }
}

void set_values(CSVec& vec, double scalar)
{
  std::fill(vec.coefs().begin(), vec.coefs().end(), scalar);
}

void add_CSVec_CSVec(const CSVec& u, CSVec& v)
{
  const std::vector<int>& indices = u.indices();
  const std::vector<double>& coefs = u.coefs();

  for(size_t i=0; i<indices.size(); ++i) {
    add_entry(v, indices[i], coefs[i]);
  }
}

}//namespace fei

