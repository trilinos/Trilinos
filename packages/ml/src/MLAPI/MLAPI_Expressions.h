#ifndef ML_EXPRESSIONS_H
#define ML_EXPRESSIONS_H

#include "ml_config.h"
#include "ml_epetra.h"
#include <iostream>
#include "MLAPI_Space.h"
#include "MLAPI_Vector.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Smoother.h"

namespace MLAPI {

Vector&
operator/(const Smoother& left, const Vector& right) {
  return(left.ApplyInverse(right));
}

template<class Left, class Right>
class BaseObjectSum {
  const Left& left; const Right& right;
public:
  BaseObjectSum(const Left& lhs, const Right& rhs)
    : left(lhs), right(rhs) {}

  double operator[] (size_t i) const {
    return left[i] + right[i];
  }
};

template<>
class BaseObjectSum<Vector,double> {
  const Vector& left; const double right;
public:
  BaseObjectSum(const Vector& lhs, const double rhs)
    : left(lhs), right(rhs) {}

  double operator[] (size_t i) const {
    return(left[i] + right);
  }
};

template<class Left, class Right>
class BaseObjectDiff {
  const Left& left; const Right& right;
public:
  BaseObjectDiff(const Left& lhs, const Right& rhs)
    : left(lhs), right(rhs) {}

  double operator[] (size_t i) const {
    return(left[i] - right[i]);
  }
};

template<>
class BaseObjectDiff<Vector,double> {
  const Vector& left; const double right;
public:
  BaseObjectDiff(const Vector& lhs, const double rhs)
    : left(lhs), right(rhs) {}

  double operator[] (size_t i) const {
    return(left[i] - right);
  }
};

template<class Left, class Right>
class BaseObjectMult {
  const Left& left; const Right& right;
public:
  BaseObjectMult(const Left& lhs, const Right& rhs)
    : left(lhs), right(rhs) {}

  double operator[] (size_t i) const {
    return left[i] * right[i];
  }
};

template<>
class BaseObjectMult<Vector,double> {
  const Vector& left; const double right;
public:
  BaseObjectMult(const Vector& lhs, const double rhs)
    : left(lhs), right(rhs) {}

  double operator[] (size_t i) const {
    return left[i] * right;
  }
};

template<class Left, class Right>
class BaseObjectDiv {
  const Left& left; const Right& right;
public:
  BaseObjectDiv(const Left& lhs, const Right& rhs)
    : left(lhs), right(rhs) {}

  double operator[] (size_t i) const {
    return left[i] / right[i];
  }
};

template<>
class BaseObjectDiv<Vector,double> {
  const Vector& left; const double right;
public:
  BaseObjectDiv(const Vector& lhs, const double rhs)
    : left(lhs), right(rhs) {}

  double operator[] (size_t i) const {
    return left[i] / right;
  }
};

Vector& 
operator* (const Operator& Op, const Vector& V)
{
  return(Op.Apply(V));
}

double
operator* (const Vector& Left, const Vector& Right)
{
  return(Left.DotProduct(Right));
}

template<>
class BaseObjectMult<Operator,Vector> {
  const Operator& left; const Vector& right;
public:
  BaseObjectMult(const Operator& lhs, const Vector& rhs)
    : left(lhs), right(rhs) {}

  double operator[] (size_t i) const {
    std::cout << "IN MATVEC" << std::endl;
    ML_EXIT(0);
    return(0.0);
  }
};

// operator+ just stores references
inline BaseObjectSum<Vector,double>
operator+(const Vector& left, double right) {
  return(BaseObjectSum<Vector,double>(left,right));
}

// operator+ just stores references
inline BaseObjectDiff<Vector,double>
operator-(const Vector& left, double right) {
  return(BaseObjectDiff<Vector,double>(left,right));
}

inline BaseObjectMult<Vector,double>
operator*(const Vector& left, double right) {
  return(BaseObjectMult<Vector,double>(left,right));
}

inline BaseObjectDiv<Vector,double>
operator/(const Vector& left, double right) {
  return(BaseObjectDiv<Vector,double>(left,right));
}


inline BaseObjectSum<Vector,Vector>
operator+(const Vector& left, const Vector& right) {
  return(BaseObjectSum<Vector,Vector>(left,right));
}

inline BaseObjectDiff<Vector,Vector>
operator-(const Vector& left, const Vector& right) {
  return(BaseObjectDiff<Vector,Vector>(left,right));
}

template<class Left, class Right>
inline BaseObjectMult<BaseObjectMult<Left,Right>,Vector>
operator*(const BaseObjectMult<Left,Right>& left,const Vector& right) {
  return(BaseObjectMult<BaseObjectMult<Left,Right>,Vector>(left, right));
}

template<class Left, class Right>
inline BaseObjectMult<BaseObjectMult<Left,Right>,BaseObjectMult<Left,Right> >
operator*(const BaseObjectMult<Left,Right>& left,const BaseObjectMult<Left,Right>& right) {
  return(BaseObjectMult<BaseObjectMult<Left,Right>,BaseObjectMult<Left,Right> >(left, right));
}

template<class Left, class Right>
inline BaseObjectSum<BaseObjectSum<Left,Right>,Vector>
operator+(const BaseObjectSum<Left,Right>& left,const Vector& right) {
  return(BaseObjectSum<BaseObjectSum<Left,Right>,Vector>(left, right));
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
std::ostream& operator<< (std::ostream& os, const Vector& v) 
{
  for (size_t i = 0 ; i < v.VectorSpace().NumMyElements() ; ++i)
    os << v[i] << ' ';
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
    os << 0 << "\t\t" << i << "\t\t" << v.GID(i) << std::endl;
  os << std::endl;
  return(os);
}

} // namespace MLAPI
#endif // if ML_EXPRESSIONS_H
