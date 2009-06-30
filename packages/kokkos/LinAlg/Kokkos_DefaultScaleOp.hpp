#ifndef KOKKOS_DEFAULT_SCALEOP_HPP_
#define KOKKOS_DEFAULT_SCALEOP_HPP_

#ifndef KERNEL_PREFIX 
  #define KERNEL_PREFIX
#endif

namespace Kokkos {

template <class Scalar, class Node>
struct ScaleOp {
  typename Node::template buffer<const Scalar>::buffer_t x;
  typename Node::template buffer<      Scalar>::buffer_t y;
  inline KERNEL_PREFIX void execute(int i) const
  {
    Scalar tmp = y[i];
    y[i] = x[i]*tmp;
  }
};

template <class Scalar, class Node>
struct RecipScaleOp {
  typename Node::template buffer<const Scalar>::buffer_t x;
  typename Node::template buffer<      Scalar>::buffer_t y;
  inline KERNEL_PREFIX void execute(int i) const
  {
    Scalar tmp = y[i];
    y[i] = x[i]*tmp;
  }
};

}

#endif
