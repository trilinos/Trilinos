/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_FillableMat_hpp_
#define _fei_FillableMat_hpp_

#include <vector>
#include <fei_FillableVec.hpp>
#include <fei_SSMat.hpp>

namespace fei {

class FillableMat {
 public:
  FillableMat();
  virtual ~FillableMat();

  FillableMat& operator=(const FillableMat& src);

  FillableMat& operator=(const SSMat& src);

  void setValues(double value);

  void createPosition(int row, int col);

  void sumInCoef(int row, int col, double coef);
  void putCoef(int row, int col, double coef);

  void sumInRow(int row, const int* cols, const double* coefs, unsigned len);
  void putRow(int row, const int* cols, const double* coefs, unsigned len);

  unsigned getNumRows() const;

  bool hasRow(int row) const;

  /** If row is not present and if 'create_if_not_already_present' is false,
      then throw std::runtime_error. */
  FillableVec* getRow(int row, bool create_if_not_already_present=false);

  typedef std::map<int, FillableVec*, std::less<int>,
                fei_Pool_alloc<std::pair<const int,FillableVec*> > > feipoolmat;

  typedef feipoolmat::iterator iterator;
  typedef feipoolmat::const_iterator const_iterator;

  iterator begin() {return matdata_.begin();}
  iterator end() {return matdata_.end();}

  const_iterator begin() const {return matdata_.begin();}
  const_iterator end() const {return matdata_.end();}

  void clear();

  bool operator==(const FillableMat& rhs) const;

  bool operator!=(const FillableMat& rhs) const;

 private:
  feipoolmat matdata_;
  fei_Pool_alloc<FillableVec> vecpool_;
}; //class FillableMat

/** Return the number of nonzeros in mat. */
int count_nnz(const FillableMat& mat);

/** Fill a std::vector with the row-numbers from the given matrix. */
void get_row_numbers(const FillableMat& mat, std::vector<int>& rows);

}//namespace fei

#endif

