/*
// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
*/

#ifndef MY_BETTER_OPERATOR_HPP
#define MY_BETTER_OPERATOR_HPP

#include "BelosConfigDefs.hpp"
#include "BelosOperator.hpp"
#include "MyMultiVec.hpp"
#include "Teuchos_BLAS.hpp"

//! Simple example of a user's defined Belos::Operator class.
/*!
 * This is a simple, single processor example of user's defined
 * Belos::Operator-derived class. The class is templated with ScalarType;
 * possible choices are, for example, "float", "double", or
 * "std::complex<double>".
 *
 * This class is based on the MyOperator class written by
 * Oscar Chinallato (ETHZ/ICOS) and Marzio Sala (ETHZ/COLAB)
 *
 * \author Christopher Baker (FSU/SCS,SNL/CSRI)
 *
 */
template <class ScalarType>
class MyBetterOperator : public Belos::Operator<ScalarType>
{

public:

  //! Constructor
  MyBetterOperator(const int nrows, const int *colptr,
                   const int nnz, const int *rowin, const ScalarType *vals)
  : _nr(nrows), _nnz(nnz), _cptr(nrows+1), _rind(nnz), _vals(nnz)
  {
    std::copy<const int*,IntIter>(colptr,colptr+nrows+1,_cptr.begin());
    std::copy<const int*,IntIter>(rowin,rowin+nnz,_rind.begin());
    std::copy<const ScalarType*,STIter>(vals,vals+nnz,_vals.begin());
  }

  //! Deconstructor
  ~MyBetterOperator()
  { }

  //! Applies the matrix to a multivector.
  void Apply(const Belos::MultiVec<ScalarType>& X,
             Belos::MultiVec<ScalarType>& Y,
             Belos::ETrans trans = Belos::NOTRANS) const
  {
    const MyMultiVec<ScalarType>* MyX;
    MyX = dynamic_cast<const MyMultiVec<ScalarType>*>(&X);
    assert (MyX != 0);

    MyMultiVec<ScalarType>* MyY;
    MyY = dynamic_cast<MyMultiVec<ScalarType>*>(&Y);
    assert (MyY != 0);

    // Initialize output std::vector to zero.
    MyY->MvInit( Teuchos::ScalarTraits<ScalarType>::zero() );

    assert (X.GetNumberVecs() == Y.GetNumberVecs());
    assert (X.GetGlobalLength() == Y.GetGlobalLength());

    int nv = X.GetNumberVecs();

    // Apply operator
    int IA1, IA2, ri;
    ScalarType aval;
    int i,j,v;
    for (j=0; j<_nr; j++) {
      IA1 = _cptr[j]-1;
      IA2 = _cptr[j+1]-1;
      for (i=IA1; i<IA2; i++) {
        ri = _rind[i]-1;
        aval = _vals[i];
        for (v=0; v<nv; v++) {
          (*MyY)[v][ri] += aval*(*MyX)[v][j];
        }
      }
    }
  }

  void Print( std::ostream& os ) {
    for (int j=0; j<_nr; j++) {
      int IA1 = _cptr[j]-1;
      int IA2 = _cptr[j+1]-1;
      for (int i=IA1; i<IA2; i++) {
        os << "("<<_rind[i]-1<<","<<j<<")\t"<<_vals[i]<< std::endl;
      }
    }
  }

  private:
  typedef typename std::vector<ScalarType>::iterator STIter;
  typedef std::vector<int>::iterator        IntIter;
  //! Number of rows and columns
  int _nr, _nnz;
  //! Column pointers
  std::vector<int> _cptr;
  //! Row indices
  std::vector<int> _rind;
  //! Values
  std::vector<ScalarType> _vals;
};

#endif //MY_BETTER_OPERATOR_HPP
