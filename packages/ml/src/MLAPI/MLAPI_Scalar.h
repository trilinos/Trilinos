#ifndef ML_SCALAR_H
#define ML_SCALAR_H
#include "ml_config.h"

namespace MLAPI {

class Scalar {

  double value_;
public:
  Scalar(double value = 0.0) :
    value_(value) {}

  Scalar& operator=(double rhs) {
    value_ = rhs;
  }

  inline const double& operator[](int i) const {
    return value_;
  }

  inline double& operator[] (int i) {
    return(value_);
  }

  double Value() const {
    return(value_);
  }
}; //

} // namespace MLAPI

#endif // if ML_SCALAR_H
