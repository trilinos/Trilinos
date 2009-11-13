#ifndef TEST_OPS_HPP_
#define TEST_OPS_HPP_

#ifndef KERNEL_PREFIX 
  #define KERNEL_PREFIX
#endif

template <typename Scalar>
struct InitOp {
  Scalar *x;
  inline KERNEL_PREFIX void execute(int i) {
    x[i] = 1;
  }
};

struct NullOp {
  typedef int ReductionType;
  static inline KERNEL_PREFIX int identity() {return 0;}
  static inline KERNEL_PREFIX int reduce(int x, int y) {return x+y;}
  static inline KERNEL_PREFIX int generate(int i) {return 0;}
};

template <typename Scalar>
struct SumOp {
  typedef Scalar ReductionType;

  const Scalar *x;

  static inline KERNEL_PREFIX ReductionType identity() {
    return (Scalar)0;
  }

  static inline KERNEL_PREFIX ReductionType reduce(ReductionType x, ReductionType y) {
    return x+y;
  }

  inline KERNEL_PREFIX Scalar generate(int i) {
    return x[i];
  }
};

#endif
