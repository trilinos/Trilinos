#ifndef MY_OPERATOR_HPP
#define MY_OPERATOR_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziOperator.hpp"
#include "MyMultiVec.hpp"

//! Simple example of a user's defined Anasazi::Operator class.
/*! 
 * This is a simple, single processor example of user's defined
 * Anasazi::Operator-derived class. The class is templated with ScalarType;
 * possible choices are, for example, "float", "double", or
 * "complex<double>".
 *
 * This file can be easily extended to tackle more interesting cases.
 *
 * \author Oscar Chinallato (ETHZ/ICOS) and Marzio Sala (ETHZ/COLAB)
 *
 * \date Last modified on 01-Nov-05
 */
class MyOperator : public Anasazi::Operator<ScalarType>
{

public:
  
  /* Constructs a square matrix with \c NumRows rows and columns.
   * The matrix is tridiagonal, and the computational stencil is 
   * [-1, 2, -1]
   */
  MyOperator(const int NumRows) :
    NumRows_(NumRows)
  {
    l_ = -1.0;
    d_ =  2.0;
    u_ = -1.0;
  }

  // Constructor for tridiagonal matrix.
  MyOperator(const int NumRows, std::vector<ScalarType> ldu) :
    NumRows_(NumRows)
  {
    l_ = ldu[0];
    d_ = ldu[1];
    u_ = ldu[2];    
  }

  // Constructor for a diagonal matrix with variable entries.
  MyOperator(std::vector<ScalarType> diag) :
    NumRows_(diag.size())
  {
    int i;

    diag_.resize(diag.size());
    for(i=0; i<diag_.size(); ++i)
      diag_[i] = diag[i];
  }

  //! Dtor
  ~MyOperator()
  {}
  
  //! Applies the tridiagonal or diagonal matrix to a multivector.
  Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType>& X, 
			    Anasazi::MultiVec<ScalarType>& Y) const
  {
    const MyMultiVec* MyX;
    MyX = dynamic_cast<const MyMultiVec*>(&X); 
    assert (MyX != 0);
    
    MyMultiVec* MyY;
    MyY = dynamic_cast<MyMultiVec*>(&Y); 
    assert (MyY != 0);
    
    assert (X.GetNumberVecs() == Y.GetNumberVecs());
    assert (X.GetVecLength() == Y.GetVecLength());
   
    if (diag_.size() == 0)
    {
      // This is a tridiagonal matrix
      for (int v = 0 ; v < X.GetNumberVecs() ; ++v)
      {
        for (int i = 0 ; i < X.GetVecLength() ; ++i)
        {
          if (i == 0)
            (*MyY)[v][i] = (d_ * (*MyX)[v][i] + u_ * (*MyX)[v][i + 1]);
          else if (i == X.GetVecLength() - 1)
            (*MyY)[v][i] = (d_ * (*MyX)[v][i] + l_ * (*MyX)[v][i-1]);
          else
            (*MyY)[v][i] = (d_ * (*MyX)[v][i] + l_ * (*MyX)[v][i-1] + u_ * (*MyX)[v][i+1]);
        }
      }
    } 
    else
    {
      // This is a diagonal matrix
      for (int v = 0 ; v < X.GetNumberVecs() ; ++v)
      {
        for (int i = 0 ; i < X.GetVecLength() ; ++i)
        {
          (*MyY)[v][i] = diag_[i] * (*MyX)[v][i];
        }
      }      
    }
    return(Anasazi::Ok);
  }
  
private:
  //! Number of rows and columns
  int NumRows_;
  //! Elements on subdiagonal, diagonal, and superdiagonal.
  ScalarType l_, d_, u_;
  //! Elements on diagonal (for variable-diagonal case).
  std::vector<ScalarType> diag_;
};

namespace Anasazi 
{
  template <> 
  class OperatorTraits < double, MyMultiVec, MyOperator >
  {
  public:
    
    /*! \brief This method takes the Epetra_MultiVector \c x and
      applies the Epetra_Operator \c Op to it resulting in the Epetra_MultiVector \c y.
    */    
    static ReturnType Apply (const MyOperator& Op, const MyMultiVec& x, 
			      MyMultiVec& y )
    { 
      return (Op.Apply( x, y ));
    }
  };

} // namespace Anasazi
  
#endif //MY_OPERATOR_HPP
