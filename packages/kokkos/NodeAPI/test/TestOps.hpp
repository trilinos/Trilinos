#ifndef TEST_OPS_HPP_
#define TEST_OPS_HPP_

#ifndef KERNEL_PREFIX 
  #define KERNEL_PREFIX
#endif

template <class Scalar, class Node>
struct InitOp {
  typename Node::template buffer<Scalar>::buffer_t x;
  inline KERNEL_PREFIX void execute(int i) {
    x[i] = 1;
  }
};

template <class Node>
struct NullOp {
  typedef int ReductionType;
  static inline KERNEL_PREFIX int identity() {return 0;}
  static inline KERNEL_PREFIX int reduce(int x, int y) {return x+y;}
  static inline KERNEL_PREFIX int generate(int i) {return 0;}
};

template <class Scalar, class Node>
struct SumOp {
  typedef Scalar ReductionType;
  typedef Node   NodeType;

  typename Node::template buffer<const Scalar>::buffer_t x;

  static inline KERNEL_PREFIX ReductionType identity() {
    return (Scalar)0;
  }

  static inline KERNEL_PREFIX ReductionType reduce(ReductionType x, ReductionType y) 
  {
    return x+y;
  }

  inline KERNEL_PREFIX Scalar generate(int i) {
    return x[i];
  }
};

#endif
