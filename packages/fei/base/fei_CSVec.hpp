#ifndef _fei_CSVec_hpp_
#define _fei_CSVec_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_FillableVec.hpp>
#include <fei_SSVec.hpp>
#include <vector>

namespace fei {

/** 'Compressed Sparse Vector' stored as two std::vectors: a vector of ints for the
    indices and a vector of doubles for the coefficients.

   Non-member functions add_entry and put_entry maintain sortedness of the
   vector when inserting new entries.
*/
class CSVec {
 public:
  CSVec(const FillableVec& invec);
  CSVec(unsigned sz=0);
  virtual ~CSVec();

  CSVec& operator=(const FillableVec& invec);

  CSVec& operator=(const SSVec& invec);

  std::vector<int>& indices() {return indices_;}
  const std::vector<int>& indices() const {return indices_;}
  std::vector<double>& coefs() {return coefs_;}
  const std::vector<double>& coefs() const {return coefs_;}

  size_t size() const {return indices_.size();}

  void clear() { indices_.clear(); coefs_.clear(); }

 private:
  std::vector<int> indices_;
  std::vector<double> coefs_;
};//class CSVec

void add_entry(CSVec& vec, int eqn, double coef);

void put_entry(CSVec& vec, int eqn, double coef);

void remove_entry(CSVec& vec, int eqn);

void set_values(CSVec& vec, double scalar);

/** form v = v + u
*/
void add_CSVec_CSVec(const CSVec& u, CSVec& v);

}//namespace fei

#endif

