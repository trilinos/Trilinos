#ifndef MYOPERATOR_HPP
#define MYOPERATOR_HPP

#include "AnasaziConfigDefs.hpp"
//#include "AnasaziBlockDavidson.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_MultiVector.h"

#include "AnasaziMultiVec.hpp"
#include "AnasaziOperator.hpp"
#include "AnasaziOperatorTraits.hpp"
//#include "Teuchos_SerialDenseMatrix.hpp"
//#include "Teuchos_ScalarTraits.hpp"

#include "MyMultiVec.hpp"

#include <complex>

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


//     Anasazi::ReturnType Apply(const Epetra_MultiVector &X, 
//                               Epetra_MultiVector &Y) const
//     {
//       for (int v = 0 ; v < X.NumVectors() ; ++v)
//       {
//         for (int i = 0 ; i < X.MyLength() ; ++i)
//         {
//           if (i == 0)
//             Y[v][i] = (2.0 * X[v][i] - X[v][i + 1])/h_/h_;
//           else if (i == X.MyLength() - 1)
//             Y[v][i] = (2.0 * X[v][i] - X[v][i-1])/h_/h_;
//           else
//             Y[v][i] = (2.0 * X[v][i] - X[v][i-1] - X[v][i+1])/h_/h_;
//         }
//       }

//       return(Anasazi::Ok);
//     }

  private:
  int NumRows_;
  double h_;
};




namespace Anasazi {

template <>
class OperatorTraits<ScalarType, MyMultiVec, MyOperator>
{
  public:

    static Anasazi::ReturnType Apply(const MyOperator& Op, 
                                     const MyMultiVec& x, 
				     MyMultiVec& y) 
    {
      return(Op.Apply(x, y));
    }
};

} // namespace Anasazi


#endif //MYOPERATOR_HPP
