#ifndef ML_EXPRESSIONS_H
#define ML_EXPRESSIONS_H

#include "ml_include.h"
#include "ml_epetra.h"
#include <iostream>
#include "MLAPI_Space.h"
#include "MLAPI_DoubleVector.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Smoother.h"

namespace MLAPI {

Operator operator+(const Operator& A, const Operator& B)
{
  assert (A.DomainSpace() == B.DomainSpace());
  assert (A.RangeSpace() == B.RangeSpace());

  ML_Operator* ML_AplusB = ML_Operator_Create(GetMLComm());
  ML_Operator_Add(A.GetOperator(),B.GetOperator(),ML_AplusB,MatrixType,1.);
  Operator AplusB(A.DomainSpace(),A.RangeSpace(), ML_AplusB,true);
  return(AplusB);
}

Operator operator-(const Operator& A, const Operator& B)
{
  assert (A.DomainSpace() == B.DomainSpace());
  assert (A.RangeSpace() == B.RangeSpace());

  ML_Operator* ML_AplusB = ML_Operator_Create(GetMLComm());
  ML_Operator_Add(A.GetOperator(),B.GetOperator(),ML_AplusB,MatrixType,-1.);
  Operator AplusB(A.DomainSpace(),A.RangeSpace(), ML_AplusB,true);
  return(AplusB);
}

Operator operator*(const Operator& A, const Operator& B)
{
  // FIXME check on spaces
  ML_Operator* ML_AtimesB = ML_Operator_Create(GetMLComm());
  ML_2matmult(A.GetOperator(), B.GetOperator(), ML_AtimesB, MatrixType);
  Operator AtimesB(B.DomainSpace(),A.RangeSpace(), ML_AtimesB,true);
  return(AtimesB);
}

DoubleVector&
operator/(const Smoother& left, const DoubleVector& right) {
  return(left.ApplyInverse(right));
}

template<class Left, class Right>
class BaseObjectSum {
  const Left& left; const Right& right;
public:
  BaseObjectSum(const Left& lhs, const Right& rhs)
    : left(lhs), right(rhs) {}

  double operator() (size_t i) const {
    return(left(i) + right(i));
  }
  const Space& VectorSpace() const
  {
    return(left.VectorSpace());
  }
};

template<>
class BaseObjectSum<DoubleVector,double> {
  const DoubleVector& left; const double right;
public:
  BaseObjectSum(const DoubleVector& lhs, const double rhs)
    : left(lhs), right(rhs) {}

  double operator() (size_t i) const {
    return(left(i) + right);
  }
  const Space& VectorSpace() const
  {
    return(left.VectorSpace());
  }
};

template<class Left, class Right>
class BaseObjectDiff {
  const Left& left; const Right& right;
public:
  BaseObjectDiff(const Left& lhs, const Right& rhs)
    : left(lhs), right(rhs) {}

  double operator() (size_t i) const {
    return(left(i) - right(i));
  }
  const Space& VectorSpace() const
  {
    return(left.VectorSpace());
  }
};

template<>
class BaseObjectDiff<DoubleVector,double> {
  const DoubleVector& left; const double right;
public:
  BaseObjectDiff(const DoubleVector& lhs, const double rhs)
    : left(lhs), right(rhs) {}

  double operator() (size_t i) const {
    return(left(i) - right);
  }
  const Space& VectorSpace() const
  {
    return(left.VectorSpace());
  }
};

template<class Left, class Right>
class BaseObjectMult {
  const Left& left; const Right& right;
public:
  BaseObjectMult(const Left& lhs, const Right& rhs)
    : left(lhs), right(rhs) {}

  double operator() (size_t i) const {
    return left(i) * right(i);
  }

  const Left& GetLeft() const
  {
    return(left);
  }

  const Right& GetRight() const
  {
    return(right);
  }

  const Space& VectorSpace() const
  {
    return(left.VectorSpace());
  }
    
};

Operator operator-(const Operator& A, const BaseObjectMult<double,Operator>& B)
{
  assert (A.DomainSpace() == B.GetRight().DomainSpace());
  assert (A.RangeSpace() == B.GetRight().RangeSpace());

  ML_Operator* ML_AplusB = ML_Operator_Create(GetMLComm());
  ML_Operator_Add(A.GetOperator(),B.GetRight().GetOperator(),
                  ML_AplusB,MatrixType, - B.GetLeft());
  Operator AplusB(A.DomainSpace(),A.RangeSpace(), ML_AplusB,true);
  return(AplusB);
}

Operator operator+(const Operator& A, const BaseObjectMult<double,Operator>& B)
{
  assert (A.DomainSpace() == B.GetRight().DomainSpace());
  assert (A.RangeSpace() == B.GetRight().RangeSpace());

  ML_Operator* ML_AplusB = ML_Operator_Create(GetMLComm());
  ML_Operator_Add(A.GetOperator(),B.GetRight().GetOperator(),
                  ML_AplusB,MatrixType, B.GetLeft());
  Operator AplusB(A.DomainSpace(),A.RangeSpace(), ML_AplusB,true);
  return(AplusB);
}

Operator operator+(const BaseObjectMult<double,Operator>& A, 
                    const BaseObjectMult<double,Operator>& B)
{
  assert (A.GetRight().DomainSpace() == B.GetRight().DomainSpace());
  assert (A.GetRight().RangeSpace() == B.GetRight().RangeSpace());

  ML_Operator* ML_AplusB = ML_Operator_Create(GetMLComm());
  ML_Operator_Add2(A.GetRight().GetOperator(),B.GetRight().GetOperator(),
                   ML_AplusB,MatrixType, A.GetLeft(),B.GetLeft());
  Operator AplusB(A.GetRight().DomainSpace(), A.GetRight().RangeSpace(),
                  ML_AplusB,true);
  return(AplusB);
}

Operator operator-(const BaseObjectMult<double,Operator>& A, 
                    const BaseObjectMult<double,Operator>& B)
{
  assert (A.GetRight().DomainSpace() == B.GetRight().DomainSpace());
  assert (A.GetRight().RangeSpace() == B.GetRight().RangeSpace());

  ML_Operator* ML_AplusB = ML_Operator_Create(GetMLComm());
  ML_Operator_Add2(A.GetRight().GetOperator(),B.GetRight().GetOperator(),
                   ML_AplusB,MatrixType, A.GetLeft(),-B.GetLeft());
  Operator AplusB(A.GetRight().DomainSpace(), A.GetRight().RangeSpace(),
                  ML_AplusB,true);
  return(AplusB);
}

template<>
class BaseObjectMult<double,DoubleVector> {
  const double left; const DoubleVector& right;
public:
  BaseObjectMult(const double lhs, const DoubleVector& rhs)
    : left(lhs), right(rhs) {}

  double operator() (size_t i) const {
    return left * right(i);
  }

  const Space& VectorSpace() const
  {
    return(right.VectorSpace());
  }
};

template<>
class BaseObjectMult<DoubleVector,double> {
  const DoubleVector& left; const double right;
public:
  BaseObjectMult(const DoubleVector& lhs, const double rhs)
    : left(lhs), right(rhs) {}

  double operator() (size_t i) const {
    return(left(i) * right);
  }

  const Space& VectorSpace() const
  {
    return(left.VectorSpace());
  }
};

template<class Left, class Right>
class BaseObjectDiv {
  const Left& left; const Right& right;
public:
  BaseObjectDiv(const Left& lhs, const Right& rhs)
    : left(lhs), right(rhs) {}

  double operator() (size_t i) const {
    return left(i) / right(i);
  }

  const Space& VectorSpace() const
  {
    return(left.VectorSpace());
  }
};

template<>
class BaseObjectDiv<DoubleVector,double> {
  const DoubleVector& left; const double right;
public:
  BaseObjectDiv(const DoubleVector& lhs, const double rhs)
    : left(lhs), right(rhs) {}

  double operator() (size_t i) const {
    return left(i) / right;
  }

  const Space& VectorSpace() const
  {
    return(left.VectorSpace());
  }
};

template<>
class BaseObjectDiv<double,DoubleVector> {
  const double left;
  const DoubleVector& right;
public:
  BaseObjectDiv(const double lhs, const DoubleVector& rhs)
    : left(lhs), right(rhs) {}

  double operator() (size_t i) const {
    return left / right(i);
  }

  const Space& VectorSpace() const
  {
    return(right.VectorSpace());
  }
};

DoubleVector& 
operator* (const Operator& Op, const DoubleVector& V)
{
  return(Op.Apply(V));
}

double
operator* (const DoubleVector& Left, const DoubleVector& Right)
{
  return(Left.DotProduct(Right));
}

template<>
class BaseObjectMult<Operator,DoubleVector> {
  const Operator& left; const DoubleVector& right;
public:
  BaseObjectMult(const Operator& lhs, const DoubleVector& rhs)
    : left(lhs), right(rhs) {}

  double operator() (size_t i) const {
    ML_EXIT(0);
    return(0.0);
  }
};

// operator+ just stores references
inline BaseObjectSum<DoubleVector,double>
operator+(const DoubleVector& left, double right) {
  return(BaseObjectSum<DoubleVector,double>(left,right));
}

// operator+ just stores references
inline BaseObjectDiff<DoubleVector,double>
operator-(const DoubleVector& left, double right) {
  return(BaseObjectDiff<DoubleVector,double>(left,right));
}

inline BaseObjectMult<DoubleVector,double>
operator*(const DoubleVector& left, double right) {
  return(BaseObjectMult<DoubleVector,double>(left,right));
}

inline BaseObjectDiv<DoubleVector,double>
operator/(const DoubleVector& left, double right) {
  return(BaseObjectDiv<DoubleVector,double>(left,right));
}

inline BaseObjectDiv<double,DoubleVector>
operator/(double left, const DoubleVector& right) {
  return(BaseObjectDiv<double,DoubleVector>(left,right));
}


inline BaseObjectSum<DoubleVector,DoubleVector>
operator+(const DoubleVector& left, const DoubleVector& right) {
  return(BaseObjectSum<DoubleVector,DoubleVector>(left,right));
}

inline BaseObjectDiff<DoubleVector,DoubleVector>
operator-(const DoubleVector& left, const DoubleVector& right) {
  return(BaseObjectDiff<DoubleVector,DoubleVector>(left,right));
}

template<class Left, class Right>
inline BaseObjectMult<BaseObjectMult<Left,Right>,DoubleVector>
operator*(const BaseObjectMult<Left,Right>& left,const DoubleVector& right) {
  return(BaseObjectMult<BaseObjectMult<Left,Right>,DoubleVector>(left, right));
}

template<class Left, class Right>
inline BaseObjectMult<BaseObjectMult<Left,Right>,BaseObjectMult<Left,Right> >
operator*(const BaseObjectMult<Left,Right>& left,const BaseObjectMult<Left,Right>& right) {
  return(BaseObjectMult<BaseObjectMult<Left,Right>,BaseObjectMult<Left,Right> >(left, right));
}

template<class Left, class Right>
inline BaseObjectSum<BaseObjectSum<Left,Right>,DoubleVector>
operator+(const BaseObjectSum<Left,Right>& left,const DoubleVector& right) {
  return(BaseObjectSum<BaseObjectSum<Left,Right>,DoubleVector>(left, right));
}

template<class Left, class Right>
inline BaseObjectSum<BaseObjectMult<Left,Right>, BaseObjectMult<Left,Right> >
operator+(const BaseObjectMult<Left,Right>& left, const BaseObjectMult<Left,Right>& right) {
  return(BaseObjectSum<BaseObjectMult<Left,Right>,BaseObjectMult<Left,Right> >(left, right));
}

template<class Left, class Right>
inline BaseObjectSum<BaseObjectSum<Left,Right>, BaseObjectMult<Left,Right> >
operator+(const BaseObjectSum<Left,Right>& left, const BaseObjectMult<Left,Right>& right) {
  return(BaseObjectSum<BaseObjectSum<Left,Right>,BaseObjectMult<Left,Right> >(left, right));
}

template<class Left, class Right>
inline BaseObjectSum<BaseObjectSum<Left,Right>, BaseObjectSum<Left,Right> >
operator+(const BaseObjectSum<Left,Right>& left, const BaseObjectSum<Left,Right>& right) {
  return(BaseObjectSum<BaseObjectSum<Left,Right>,BaseObjectSum<Left,Right> >(left, right));
}

template<class Left, class Right>
inline BaseObjectSum<Left,Right>
operator+(const Left& left, const Right& right) {
  return(BaseObjectSum<Left,Right>(left, right));
}

template<class Left, class Right>
inline BaseObjectMult<Left,Right>
operator*(const Left& left, const Right& right) {
  return(BaseObjectMult<Left,Right>(left, right));
}

// =====
std::ostream& operator<< (std::ostream& os, const DoubleVector& v) 
{
  for (size_t i = 0 ; i < v.VectorSpace().NumMyElements() ; ++i)
    os << v(i) << ' ';
  os << std::endl;
  return(os);
}

std::ostream& operator<< (std::ostream& os, const Space& v) 
{
  os << std::endl;
  os << "MLAPI::Space" << std::endl;
  os << "NumMyElements() = " << v.NumMyElements() << std::endl;
  os << "NumGlobalElements() = " << v.NumGlobalElements() << std::endl;

  os << "ProcID\t\tLID\t\tGID" << std::endl;
  for (size_t i = 0 ; i < v.NumMyElements() ; ++i)
    os << 0 << "\t\t" << i << "\t\t" << v(i) << std::endl;
  os << std::endl;
  return(os);
}

std::ostream& operator<< (std::ostream& os, const Operator& Op) 
{
  int    *bindx;
  double *val;
  int    allocated, row_length;
  ML_Operator* matrix = Op.GetOperator();

  if (matrix->getrow == NULL) 
    throw("getrow not set");

  allocated = 100;
  bindx = (int    *)  ML_allocate(allocated*sizeof(int   ));
  val   = (double *)  ML_allocate(allocated*sizeof(double));

  int NumGlobalRows = Op.DomainSpace().NumGlobalElements();
  int NumGlobalCols = Op.RangeSpace().NumGlobalElements();

  for (int iproc = 0 ; iproc < NumProc() ; ++iproc) {

    if (iproc == 0) {
      os << "Operator `" << Op.Name() << "'" << endl;
      os << "ProcID\tGlobal Row\tGlobal Col\tValue" << endl;
      os << endl;
    }

    if (MyPID() == iproc) {

      for (int i = 0 ; i < matrix->getrow->Nrows; i++) {
        ML_get_matrix_row(matrix, 1, &i, &allocated, &bindx, &val,
                          &row_length, 0);
        for  (int j = 0; j < row_length; j++) {
          int GlobalRow = Op.DomainSpace()(i);
          int GlobalCol = Op.ColumnSpace()(bindx[j]);
          os << iproc << "\t" << GlobalRow << "\t" << GlobalCol << "\t" << val[j] << endl;
        }
      }
    }
#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  ML_free(val);
  ML_free(bindx);
  return (os);
}

} // namespace MLAPI
#endif // if ML_EXPRESSIONS_H
