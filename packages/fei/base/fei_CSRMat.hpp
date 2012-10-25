#ifndef _fei_CSRMat_hpp_
#define _fei_CSRMat_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_macros.hpp"
#include "fei_FillableMat.hpp"
#include "fei_SparseRowGraph.hpp"
#include "fei_CSVec.hpp"

namespace fei {

/** Compressed Sparse Row Matrix object.
*/
class CSRMat {
 public:
  CSRMat();
  CSRMat(const FillableMat& fmat);
  virtual ~CSRMat();

  SparseRowGraph& getGraph() {return srg_;}
  const SparseRowGraph& getGraph() const {return srg_;}
 
  std::vector<double>& getPackedCoefs() {return packedcoefs_;}
  const std::vector<double>& getPackedCoefs() const {return packedcoefs_;}

  unsigned getNumRows() const {return srg_.rowNumbers.size();}

  CSRMat& operator=(const FillableMat& src);

  CSRMat& operator+=(const CSRMat& src);

  bool operator==(const CSRMat& rhs) const;

  bool operator!=(const CSRMat& rhs) const;

 private:
  SparseRowGraph srg_;
  std::vector<double> packedcoefs_;
};//class CSRMat

/** form y = A*x */
void multiply_CSRMat_CSVec(const CSRMat& A, const CSVec& x, CSVec& y);

/** form y = A^T*x */
void multiply_trans_CSRMat_CSVec(const CSRMat& A, const CSVec& x, CSVec& y);

/** form C = A*B */
void multiply_CSRMat_CSRMat(const CSRMat& A, const CSRMat& B, CSRMat& C,
                            bool storeResultZeros=false);

/** form C = A^T*B */
void multiply_trans_CSRMat_CSRMat(const CSRMat& A, const CSRMat& B, CSRMat& C,
                                  bool storeResultZeros=false);

void add_CSRMat_to_FillableMat(const CSRMat& csrm, FillableMat& fm);

}//namespace fei

#endif

