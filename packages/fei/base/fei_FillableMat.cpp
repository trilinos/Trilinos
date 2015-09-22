/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_FillableMat.hpp>
#include <fei_EqnBuffer.hpp>
#include <fei_CSVec.hpp>

namespace fei {

//-----------------------------------------------------------------
FillableMat::FillableMat()
 : matdata_(),
   vecpool_()
{
}

//-----------------------------------------------------------------
FillableMat::FillableMat(EqnBuffer& eqnbuf)
 : matdata_(),
   vecpool_()
{
  std::vector<int>& eqnNums = eqnbuf.eqnNumbers();
  int numEqns = eqnNums.size();
  std::vector<fei::CSVec*>& eqns = eqnbuf.eqns();

  for(int i=0; i<numEqns; ++i) {
    int row = eqnNums[i];
    fei::CSVec* row_vec = eqns[i];
    int rowlen = row_vec->size();
    int* indices = &(row_vec->indices()[0]);
    double* coefs = &(row_vec->coefs()[0]);

    for(int j=0; j<rowlen; ++j) {
      putCoef(row, indices[j], coefs[j]);
    }
  }
}

//-----------------------------------------------------------------
FillableMat::~FillableMat()
{
  feipoolmat::iterator
    iter = matdata_.begin(), iter_end = matdata_.end();
  for(; iter!=iter_end; ++iter) {
    vecpool_.destroy(iter->second);
    vecpool_.deallocate(iter->second, 1);
  }
}

//-----------------------------------------------------------------
FillableMat&
FillableMat::operator=(const FillableMat& src)
{
  clear();

  FillableMat::const_iterator
    s_iter = src.begin(),
    s_end = src.end();

  for(; s_iter != s_end; ++s_iter) {
    int row = s_iter->first;
    const CSVec* srow = s_iter->second;
    const std::vector<int>& s_ind = srow->indices();
    const std::vector<double>& s_coef = srow->coefs();

    for(size_t i=0; i<s_ind.size(); ++i) {
      int col = s_ind[i];
      double coef = s_coef[i];

      putCoef(row, col, coef);
    }
  }

  return *this;
}

//-----------------------------------------------------------------
void
FillableMat::setValues(double value)
{
  feipoolmat::iterator
    iter = matdata_.begin(), iter_end = matdata_.end();

  for(; iter != iter_end; ++iter) {
    set_values(*(iter->second), value);
  }
}

//-----------------------------------------------------------------
void
FillableMat::createPosition(int row, int col)
{
  sumInCoef(row, col, 0.0);
}

//-----------------------------------------------------------------
FillableMat::feipoolmat::iterator
insert_row(FillableMat::feipoolmat& matdata,
           FillableMat::feipoolmat::iterator iter,
           int row,
           fei_Pool_alloc<CSVec>& vecpool)
{
  static CSVec dummy;

  CSVec* vptr = vecpool.allocate(1);
  vecpool.construct(vptr, dummy);

  if (vptr->indices().capacity() == 0) {
    vptr->indices().reserve(16);
    vptr->coefs().reserve(16);
  }

  return matdata.insert(iter, std::make_pair(row, vptr));
}

//-----------------------------------------------------------------
void
FillableMat::sumInCoef(int row, int col, double coef)
{
  CSVec* rowvec = create_or_getRow(row);

  add_entry(*rowvec, col, coef);
}

//-----------------------------------------------------------------
void
FillableMat::putCoef(int row, int col, double coef)
{
  CSVec* rowvec = create_or_getRow(row);

  put_entry(*rowvec, col, coef);
}

//-----------------------------------------------------------------
void
FillableMat::sumInRow(int row, const int* cols, const double* coefs,
                      unsigned len)
{
  CSVec* rowvec = create_or_getRow(row);

  for(unsigned i=0; i<len; ++i) {
    add_entry(*rowvec, cols[i], coefs[i]);
  }
}

//-----------------------------------------------------------------
void
FillableMat::putRow(int row, const int* cols, const double* coefs,
                    unsigned len)
{
  CSVec* rowvec = create_or_getRow(row);

  for(unsigned i=0; i<len; ++i) {
    put_entry(*rowvec, cols[i], coefs[i]);
  }
}

//-----------------------------------------------------------------
unsigned
FillableMat::getNumRows() const
{
  return matdata_.size();
}

//-----------------------------------------------------------------
bool
FillableMat::hasRow(int row) const
{
  feipoolmat::const_iterator iter = matdata_.find(row);
  return iter != matdata_.end();
}

//-----------------------------------------------------------------
const CSVec*
FillableMat::getRow(int row) const
{
  feipoolmat::const_iterator iter = matdata_.lower_bound(row);

  if (iter == matdata_.end() || iter->first != row) {
    throw std::runtime_error("fei::FillableMat: row not found.");
  }

  return iter->second;
}

//-----------------------------------------------------------------
CSVec*
FillableMat::create_or_getRow(int row)
{
  feipoolmat::iterator iter = matdata_.lower_bound(row);

  if (iter == matdata_.end() || iter->first != row) {
    iter = insert_row(matdata_, iter, row, vecpool_);
  }

  return iter->second;
}

//-----------------------------------------------------------------
void
FillableMat::clear()
{
  feipoolmat::iterator
    iter = matdata_.begin(), iter_end = matdata_.end();
  for(; iter!=iter_end; ++iter) {
    vecpool_.destroy(iter->second);
    vecpool_.deallocate(iter->second, 1);
  }

  matdata_.clear();
}

//-----------------------------------------------------------------
bool
FillableMat::operator==(const FillableMat& rhs) const
{
  if (getNumRows() != rhs.getNumRows()) return false;

  FillableMat::const_iterator
    this_it = begin(),
    this_end = end();

  FillableMat::const_iterator
    rhs_it = rhs.begin(),
    rhs_end = rhs.end();

  for(; this_it != this_end; ++this_it, ++rhs_it) {
    int this_row = this_it->first;
    int rhs_row = rhs_it->first;
    if (this_row != rhs_row) return false;

    const CSVec* this_row_vec = this_it->second;
    const CSVec* rhs_row_vec = rhs_it->second;

    if (*this_row_vec != *rhs_row_vec) return false;
  }

  return true;
}

//-----------------------------------------------------------------
bool
FillableMat::operator!=(const FillableMat& rhs) const
{
  return !(*this == rhs);
}

//-----------------------------------------------------------------
void print(std::ostream& os, const FillableMat& mat)
{
  FillableMat::const_iterator
    irow = mat.begin(), irowend = mat.end();
  for(; irow!=irowend; ++irow) {
    int row = irow->first;
    const CSVec* vec = irow->second;
    const std::vector<int>& v_ind = vec->indices();
    const std::vector<double>& v_coef = vec->coefs();
    os << "row " << row << ": ";
    for(size_t i=0; i<v_ind.size(); ++i) {
      os << "("<<v_ind[i]<<","<<v_coef[i]<<") ";
    }
    os << std::endl;
  }
}

//-----------------------------------------------------------------
int count_nnz(const FillableMat& mat)
{
  int nnz = 0;

  FillableMat::const_iterator
    r_iter = mat.begin(),
    r_end = mat.end();

  for(; r_iter != r_end; ++r_iter) {
    CSVec* row = r_iter->second;
    nnz += row->size();
  }

  return nnz;
}

//-----------------------------------------------------------------
void get_row_numbers(const FillableMat& mat, std::vector<int>& rows)
{
  rows.resize(mat.getNumRows());

  FillableMat::const_iterator
    m_iter = mat.begin(),
    m_end = mat.end();

  size_t offset = 0;
  for(; m_iter!=m_end; ++m_iter) {
    rows[offset++] = m_iter->first;
  }
}

}//namespace fei

