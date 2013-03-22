/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_FillableMat_hpp_
#define _fei_FillableMat_hpp_

#include <fei_CSVec.hpp>
#include <fei_Pool_alloc.hpp>
#include <fei_EqnBuffer.hpp>
#include <vector>
#include <map>

namespace fei {

class FillableMat {
 public:
  FillableMat();
  FillableMat(EqnBuffer& eqnbuf);
  virtual ~FillableMat();

  FillableMat& operator=(const FillableMat& src);

  void setValues(double value);

  void createPosition(int row, int col);

  void sumInCoef(int row, int col, double coef);
  void putCoef(int row, int col, double coef);

  void sumInRow(int row, const int* cols, const double* coefs, unsigned len);
  void putRow(int row, const int* cols, const double* coefs, unsigned len);

  unsigned getNumRows() const;

  bool hasRow(int row) const;

  const CSVec* getRow(int row) const;
  CSVec* create_or_getRow(int row);

  typedef std::map<int, CSVec*, std::less<int>,
                fei_Pool_alloc<std::pair<const int,CSVec*> > > feipoolmat;

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
  fei_Pool_alloc<CSVec> vecpool_;
}; //class FillableMat

/** Print the contents of a FillableMat to the given std::ostream. */
void print(std::ostream& os, const FillableMat& mat);

/** Return the number of nonzeros in mat. */
int count_nnz(const FillableMat& mat);

/** Fill a std::vector with the row-numbers from the given matrix. */
void get_row_numbers(const FillableMat& mat, std::vector<int>& rows);

}//namespace fei

#endif

