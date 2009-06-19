#ifndef TEST_OPS_HPP_
#define TEST_OPS_HPP_

#ifndef KERNEL_PREFIX 
  #define KERNEL_PREFIX
#endif

template <class Scalar, class Node>
struct InitOp {
  typename Node::template buffer<Scalar>::buffer_t x;
  int n;
  inline KERNEL_PREFIX void execute(int i) {
    x[i] = 1;
  }
};

template <class Scalar, class Node>
struct SumOp {
  typedef Scalar ReductionType;
  typedef Node   NodeType;

  typename Node::template buffer<const Scalar>::buffer_t x;

  ReductionType result;

  inline SumOp() {
    result = identity();
  }

  static inline KERNEL_PREFIX ReductionType identity() {
    return (Scalar)0;
  }

  inline KERNEL_PREFIX ReductionType reduce(ReductionType x, ReductionType y) 
  {
    return x+y;
  }

  inline KERNEL_PREFIX Scalar generate(int i) {
    return x[i];
  }
};

#endif
