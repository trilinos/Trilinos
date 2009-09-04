#ifndef TPETRA_NODEOPS_HPP
#define TPETRA_NODEOPS_HPP

#ifndef KERNEL_PREFIX 
  #define KERNEL_PREFIX
#endif

namespace Tpetra {

  template <class Scalar>
  struct InitOp {
    Scalar val;
    Scalar *x;
    inline KERNEL_PREFIX void execute(size_t i) {
      x[i] = val;
    }
  };

  template <class Scalar>
  struct ScaleOp {
    Scalar val;
    Scalar *x;
    inline KERNEL_PREFIX void execute(size_t i) {
      x[i] *= val;
    }
  };

}

#endif
