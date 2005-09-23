#ifndef MY_OPERATOR_HPP
#define MY_OPERATOR_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziOperator.hpp"
#include "MyMultiVec.hpp"

class MyOperator : public Anasazi::Operator<ScalarType>
{

public:
  
  MyOperator(const int NumRows) :
    NumRows_(NumRows)
  {
    h_ = 1.0/NumRows_;
  }
  
  ~MyOperator()
  {}
  
  Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType>&X, 
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
    
    for (int v = 0 ; v < X.GetNumberVecs() ; ++v)
      {
        for (int i = 0 ; i < X.GetVecLength() ; ++i)
	  {
	    if (i == 0)
	      (*MyY)[v][i] = (2.0 * (*MyX)[v][i] - (*MyX)[v][i + 1])/h_/h_;
	    else if (i == X.GetVecLength() - 1)
	      (*MyY)[v][i] = (2.0 * (*MyX)[v][i] - (*MyX)[v][i-1])/h_/h_;
	    else
	      (*MyY)[v][i] = (2.0 * (*MyX)[v][i] - (*MyX)[v][i-1] - (*MyX)[v][i+1])/h_/h_;
	  }
      }
    
    return(Anasazi::Ok);
  }
  
private:
  int NumRows_;
  double h_;
};

#endif //MY_OPERATOR_HPP
