#ifndef TEST_OPS_HPP_
#define TEST_OPS_HPP_

#ifndef KERNEL_PREFIX 
  #define KERNEL_PREFIX
#endif

template <class Scalar, class Node>
struct InitOp {
  typename Node::template buffer<Scalar>::buffer_t x;
  unsigned int n;
  inline KERNEL_PREFIX void execute(int l, int h)
  {
    if (n < h) h = n;
    for (int ii=l; ii<h; ++ii) {
      x[ii] = ii;
    }
  }
};

template <class Scalar, class Node>
struct SumOp {
  typedef Scalar ReductionType;
  typedef Node   NodeType;

  typename Node::template buffer<const Scalar>::buffer_t x;

  ReductionType result;

  inline DotOp() {
    result = identity();
  }

  static inline KERNEL_PREFIX ReductionType identity() {
    return (Scalar)0;
  }

  inline KERNEL_PREFIX ReductionType reduce(ReductionType x, ReductionType y) 
  {
    return x+y;
  }

  inline KERNEL_PREFIX Scalar generate(int i)
  {
    return x[i];
  }
};

#endif
