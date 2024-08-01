// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Kokkos_Core.hpp"

#include "Teuchos_Assert.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Phalanx_KokkosViewOfViews.hpp"

template<typename T,int dim>
class Vector {
  T values_[dim];
public:

  template<typename INDEX_I>
  KOKKOS_INLINE_FUNCTION
  T& operator[](const INDEX_I& index){return values_[index];}

  template<typename INDEX_I>
  KOKKOS_INLINE_FUNCTION
  const T& operator[](const INDEX_I& index)const{return values_[index];}

  template<typename INDEX_I>
  KOKKOS_INLINE_FUNCTION
  volatile T& operator[](const INDEX_I& index)volatile{return values_[index];}

  template<typename INDEX_I>
  KOKKOS_INLINE_FUNCTION
  const volatile T& operator[](const INDEX_I& index)const volatile{return values_[index];}
};

class MyClass {
  Kokkos::View<double*> a_;
  double b_[3];
  Kokkos::View<double*> c_;
  Vector<double,3> d_;
  // To test for cuda warnings when MyClass is lambda captured to
  // device
  PHX::ViewOfViews<1,Kokkos::View<double*>> e_;

public:
  MyClass() :
    a_("a",3),
    c_("c",3)
  {
    Kokkos::deep_copy(a_,1.0);
    b_[0] = 1.0;
    b_[1] = 2.0;
    b_[2] = 3.0;
    Kokkos::deep_copy(c_,2.0);
    d_[0] = 1.0;
    d_[1] = 2.0;
    d_[2] = 3.0;
  }

  void KOKKOS_FUNCTION checkInternalMethod1() const
  { this->callInternalMethod1(); }

  void KOKKOS_FUNCTION
  callInternalMethod1() const
  {
    Kokkos::printf("b_[0]=%f\n",b_[0]);
    Kokkos::printf("b_[1]=%f\n",b_[1]);
    Kokkos::printf("b_[2]=%f\n",b_[2]);
    a_(0)=b_[0];
    a_(1)=b_[1];
    a_(2)=b_[2];
  }

  void KOKKOS_FUNCTION checkInternalMethod2() const
  { this->callInternalMethod2(); }

  void KOKKOS_FUNCTION
  callInternalMethod2() const
  {
    a_(0)=c_(0);
    a_(1)=c_(1);
    a_(2)=c_(2);
  }

  void KOKKOS_FUNCTION checkInternalMethod3() const
  { this->callInternalMethod3(); }

  void KOKKOS_FUNCTION
  callInternalMethod3() const
  {
    a_(0)=d_[0];
    a_(1)=d_[1];
    a_(2)=d_[2];
  }
};

TEUCHOS_UNIT_TEST(KokkosClassOnDevice, One)
{
  MyClass my_class;

  Kokkos::parallel_for("test 1",1,KOKKOS_LAMBDA (const int ) {
      my_class.checkInternalMethod1();
  });

  Kokkos::parallel_for("test 2",1,KOKKOS_LAMBDA (const int ) {
      my_class.checkInternalMethod2();
  });

  Kokkos::parallel_for("test 3",1,KOKKOS_LAMBDA (const int ) {
      my_class.checkInternalMethod3();
  });
  
  Kokkos::fence();
}
