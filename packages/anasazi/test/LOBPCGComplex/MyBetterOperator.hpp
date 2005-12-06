#ifndef MY_OPERATOR_HPP
#define MY_OPERATOR_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziOperator.hpp"
#include "MyMultiVec.hpp"
#include <algorithm>

//! Simple example of a user's defined Anasazi::Operator class.
/*! 
 * This is a simple, single processor example of user's defined
 * Anasazi::Operator-derived class. The class is templated with ScalarType;
 * possible choices are, for example, "float", "double", or
 * "complex<double>".
 *
 * This class is based on the MyOperator class written by
 * Oscar Chinallato (ETHZ/ICOS) and Marzio Sala (ETHZ/COLAB)
 *
 * \author Christopher Baker (FSU/SCS,SNL/CSRI)
 *
 */
template <class ScalarType>
class MyBetterOperator : public Anasazi::Operator<ScalarType>
{

public:

  MyBetterOperator(const int nrows, const int *colptr,
                   const int nnz, const int *rowin, const ScalarType *vals)
  {
    int i;

    _nr = nrows;
    std::copy<const int*,IntIter>(colptr,colptr+nrows+1,_cptr.begin());
  }

  //! Dtor
  ~MyBetterOperator()
  { }

  //! Applies the matrix to a multivector.
  Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType>& X, 
                            Anasazi::MultiVec<ScalarType>& Y) const
  {
    const MyMultiVec<ScalarType>* MyX;
    MyX = dynamic_cast<const MyMultiVec<ScalarType>*>(&X); 
    assert (MyX != 0);
    
    MyMultiVec<ScalarType>* MyY;
    MyY = dynamic_cast<MyMultiVec<ScalarType>*>(&Y); 
    assert (MyY != 0);
    
    assert (X.GetNumberVecs() == Y.GetNumberVecs());
    assert (X.GetVecLength() == Y.GetVecLength());
   
    return(Anasazi::Ok);
  }

private:
  typedef typename std::vector<ScalarType>::iterator STIter;
  typedef std::vector<int>::iterator        IntIter;
  //! Number of rows and columns
  int _nr;
  //! Column pointers 
  std::vector<int> _cptr;
  //! Row indices
  std::vector<int> _rind;
  //! Values
  std::vector<ScalarType> _vals;
};

#endif //MY_OPERATOR_HPP
