#ifndef KOKKOS_MULTIVECTOR_H_
#define KOKKOS_MULTIVECTOR_H_

#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>
#include <ctime>

#include <iostream>
#include <stdexcept> // For some recently added error handling

namespace Kokkos {

template<typename Scalar, class device>
struct MultiVectorDynamic{
#ifdef KOKKOS_USE_CUSPARSE
  typedef typename Kokkos::LayoutLeft layout;
#else
#ifdef KOKKOS_USE_MKL
  typedef typename Kokkos::LayoutRight layout;
#else
  typedef typename device::array_layout layout;
#endif
#endif
  typedef typename Kokkos::View<Scalar**  , layout, device>  type ;
  typedef typename Kokkos::View<const Scalar**  , layout, device>  const_type ;
  typedef typename Kokkos::View<const Scalar**  , layout, device, Kokkos::MemoryRandomAccess>  random_read_type ;
  MultiVectorDynamic() {}
  ~MultiVectorDynamic() {}
};

template<typename Scalar, class device, int n>
struct MultiVectorStatic{
  typedef Scalar scalar;
  typedef typename device::array_layout layout;
  typedef typename Kokkos::View<Scalar*[n]  , layout, device>  type ;
  typedef typename Kokkos::View<const Scalar*[n]  , layout, device>  const_type ;
  typedef typename Kokkos::View<const Scalar*[n]  , layout, device, Kokkos::MemoryRandomAccess>  random_read_type ;
  MultiVectorStatic() {}
  ~MultiVectorStatic() {}
};

/// \brief Functor for R.scale(array of alphas, MV X).
///
/// R(i,j) = alphas[j] * X(i,j), subject to the usual BLAS rules if
/// any of the alphas[j] coefficients are zero.
template<class RVector, class aVector, class XVector>
struct MV_MulScalarFunctor
{
  typedef typename XVector::execution_space execution_space;
  typedef typename XVector::size_type       size_type;

  RVector m_r;
  typename XVector::const_type m_x;
  typename aVector::const_type m_a;
  size_type n;

  MV_MulScalarFunctor () : n (1) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
    for (size_type k = 0; k < n; ++k) {
      if (m_a[k] == Kokkos::Details::ArithTraits<typename aVector::non_const_value_type>::zero ()) {
        m_r(i,k) = Kokkos::Details::ArithTraits<typename RVector::non_const_value_type>::zero ();
      } else {
        m_r(i,k) = m_a[k] * m_x(i,k);
      }
    }
  }
};

/// \brief Functor for X.scale(array of alphas).
///
/// X(i,j) *= alphas[j], subject to the usual BLAS rules if any of the
/// alphas[j] coefficients are zero.
template<class aVector, class XVector>
struct MV_MulScalarFunctorSelf
{
  typedef typename XVector::execution_space execution_space;
  typedef typename XVector::size_type       size_type;

  XVector m_x;
  typename aVector::const_type m_a;
  size_type n;

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
    for (size_type k = 0; k < n; ++k) {
      if (m_a[k] == Kokkos::Details::ArithTraits<typename aVector::non_const_value_type>::zero ()) {
        m_x(i,k) = Kokkos::Details::ArithTraits<typename XVector::non_const_value_type>::zero ();
      } else {
        m_x(i,k) *= m_a[k];
      }
    }
  }
};


//! Function for R.scale (array a, MV X) or X.scale (array a).
template<class RVector, class DataType, class Layout, class Device, class MemoryManagement, class Specialisation, class XVector>
RVector
MV_MulScalar (const RVector& r,
              const Kokkos::View<DataType,Layout,Device,MemoryManagement,Specialisation>& a,
              const XVector& x)
{
  typedef typename Kokkos::View<DataType,Layout,Device,MemoryManagement> aVector;
  if (r == x) {
    MV_MulScalarFunctorSelf<aVector,XVector> op ;
    op.m_x = x ;
    op.m_a = a ;
    op.n = x.dimension(1);
    Kokkos::parallel_for (x.dimension (0), op);
    return r;
  }

  MV_MulScalarFunctor<RVector,aVector,XVector> op ;
  op.m_r = r ;
  op.m_x = x ;
  op.m_a = a ;
  op.n = x.dimension(1);
  Kokkos::parallel_for( x.dimension(0) , op );
  return r;
}

/// \brief Functor for R.scale(value alpha, MV X).
///
/// R(i,j) = alpha * X(i,j), subject to the usual BLAS rules if alpha
/// is zero.
template<class RVector, class XVector>
struct MV_MulScalarFunctor<RVector, typename RVector::value_type, XVector>
{
  typedef typename XVector::execution_space execution_space;
  typedef typename XVector::size_type       size_type;

  RVector m_r;
  typename XVector::const_type m_x;
  typename RVector::value_type m_a;
  size_type n;

  MV_MulScalarFunctor () : n (1) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
    for (size_type k = 0; k < n; ++k) {
      if (m_a == Kokkos::Details::ArithTraits<typename RVector::value_type>::zero ()) {
        m_r(i,k) = Kokkos::Details::ArithTraits<typename XVector::non_const_value_type>::zero ();
      } else {
        m_r(i,k) = m_a * m_x(i,k);
      }
    }
  }
};

/// \brief Functor for X.scale(value alpha).
///
/// R(i,j) *= alpha, subject to the usual BLAS rules if alpha is zero.
template<class XVector>
struct MV_MulScalarFunctorSelf<typename XVector::non_const_value_type,XVector>
{
  typedef typename XVector::execution_space execution_space;
  typedef typename XVector::size_type       size_type;

  XVector m_x;
  typename XVector::non_const_value_type m_a;
  size_type n;

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
    for (size_type k = 0; k < n; ++k) {
      if (m_a == Kokkos::Details::ArithTraits<typename XVector::non_const_value_type>::zero ()) {
        m_x(i,k) = Kokkos::Details::ArithTraits<typename XVector::non_const_value_type>::zero ();
      } else {
        m_x(i,k) *= m_a;
      }
    }
  }
};


//! Function for R.scale (value a, MV X) or X.scale (value a).
template<class RVector, class XVector>
RVector
MV_MulScalar (const RVector& r,
              const typename XVector::non_const_value_type& a,
              const XVector& x)
{
  /*if(r.dimension_1()==1) {
    typedef View<typename RVector::value_type*,typename RVector::execution_space> RVector1D;
    typedef View<typename XVector::const_value_type*,typename XVector::execution_space> XVector1D;

    RVector1D r_1d = Kokkos::subview( r , ALL(),0 );
    XVector1D x_1d = Kokkos::subview( x , ALL(),0 );
    return V_MulScalar(r_1d,a,x_1d);
  }*/

  if (r == x) {
    MV_MulScalarFunctorSelf<typename XVector::non_const_value_type, XVector> op ;
    op.m_x = x;
    op.m_a = a;
    op.n = x.dimension (1);
    Kokkos::parallel_for (x.dimension (0), op);
    return r;
  }

  MV_MulScalarFunctor<RVector,typename XVector::non_const_value_type,XVector> op;
  op.m_r = r;
  op.m_x = x;
  op.m_a = a;
  op.n = x.dimension (1);
  Kokkos::parallel_for (x.dimension (0), op);
  return r;
}

/*------------------------------------------------------------------------------------------
 *-------------------------- Reciprocal element wise: y[i] = 1/x[i] ------------------------
 *------------------------------------------------------------------------------------------*/
template<class RVector, class XVector>
struct MV_ReciprocalFunctor
{
  typedef typename XVector::execution_space        execution_space;
  typedef typename XVector::size_type            size_type;

  RVector m_r;
  typename XVector::const_type m_x ;

  const size_type m_n;
  MV_ReciprocalFunctor(RVector r, XVector x, size_type n):m_r(r),m_x(x),m_n(n) {}
  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for(size_type k=0;k<m_n;k++)
     m_r(i,k) = Kokkos::Details::ArithTraits<typename XVector::non_const_value_type>::one() / m_x(i,k);
  }
};

template<class XVector>
struct MV_ReciprocalSelfFunctor
{
  typedef typename XVector::execution_space        execution_space;
  typedef typename XVector::size_type            size_type;

  XVector m_x ;

  const size_type m_n;
  MV_ReciprocalSelfFunctor(XVector x, size_type n):m_x(x),m_n(n) {}
  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for(size_type k=0;k<m_n;k++)
     m_x(i,k) = Kokkos::Details::ArithTraits<typename XVector::non_const_value_type>::one() / m_x(i,k);
  }
};

template<class RVector, class XVector>
RVector MV_Reciprocal( const RVector & r, const XVector & x)
{
  // TODO: Add error check (didn't link for some reason?)
  /*if(r.dimension_0() != x.dimension_0())
    Kokkos::Impl::throw_runtime_exception("Kokkos::MV_Reciprocal -- dimension(0) of r and x don't match");
  if(r.dimension_1() != x.dimension_1())
    Kokkos::Impl::throw_runtime_exception("Kokkos::MV_Reciprocal -- dimension(1) of r and x don't match");*/

  //TODO: Get 1D version done
  /*if(r.dimension_1()==1) {
    typedef View<typename RVector::value_type*,typename RVector::execution_space> RVector1D;
    typedef View<typename XVector::const_value_type*,typename XVector::execution_space> XVector1D;

    RVector1D r_1d = Kokkos::subview( r , ALL(),0 );
    XVector1D x_1d = Kokkos::subview( x , ALL(),0 );
    return V_MulScalar(r_1d,a,x_1d);
  }*/
  if(r==x) {
    MV_ReciprocalSelfFunctor<XVector> op(x,x.dimension_1()) ;
    Kokkos::parallel_for( x.dimension_0() , op );
    return r;
  }

  MV_ReciprocalFunctor<RVector,XVector> op(r,x,x.dimension_1()) ;
  Kokkos::parallel_for( x.dimension_0() , op );
  return r;
}

/*------------------------------------------------------------------------------------------
 *------------------- Reciprocal element wise with threshold: x[i] = 1/x[i] ----------------
 *------------------------------------------------------------------------------------------*/
template<class XVector>
struct MV_ReciprocalThresholdSelfFunctor
{
  typedef typename XVector::execution_space           execution_space;
  typedef typename XVector::size_type               size_type;
  typedef typename XVector::non_const_value_type   value_type;
  typedef Kokkos::Details::ArithTraits<value_type>        KAT;
  typedef typename KAT::mag_type                     mag_type;

  const XVector    m_x;
  const value_type m_min_val;
  const mag_type   m_min_val_mag;
  const size_type  m_n;

  MV_ReciprocalThresholdSelfFunctor(const XVector& x, const value_type& min_val, const size_type n) :
    m_x(x), m_min_val(min_val), m_min_val_mag(KAT::abs(min_val)), m_n(n) {}
  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
    for(size_type k=0;k<m_n;k++) {
      if (KAT::abs(m_x(i,k)) < m_min_val_mag)
        m_x(i,k) = m_min_val;
      else
        m_x(i,k) = KAT::one() / m_x(i,k);
    }
  }
};

template<class XVector>
XVector MV_ReciprocalThreshold( const XVector & x, const typename XVector::non_const_value_type& min_val )
{
  MV_ReciprocalThresholdSelfFunctor<XVector> op(x,min_val,x.dimension_1()) ;
  Kokkos::parallel_for( x.dimension_0() , op );
  return x;
}

/*------------------------------------------------------------------------------------------
 *-------------------------- Abs element wise: y[i] = abs(x[i]) ------------------------
 *------------------------------------------------------------------------------------------*/
template<class RVector, class XVector>
struct MV_AbsFunctor
{
  typedef typename XVector::execution_space        execution_space;
  typedef typename XVector::size_type            size_type;

  RVector m_r;
  typename XVector::const_type m_x ;

  const size_type m_n;
  MV_AbsFunctor(RVector r, XVector x, size_type n):m_r(r),m_x(x),m_n(n) {}
  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for(size_type k=0;k<m_n;k++)
     m_r(i,k) = Kokkos::Details::ArithTraits<typename XVector::non_const_value_type>::abs(m_x(i,k));
  }
};

template<class XVector>
struct MV_AbsSelfFunctor
{
  typedef typename XVector::execution_space        execution_space;
  typedef typename XVector::size_type            size_type;

  XVector m_x ;

  const size_type m_n;
  MV_AbsSelfFunctor(XVector x, size_type n):m_x(x),m_n(n) {}
  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for(size_type k=0;k<m_n;k++)
     m_x(i,k) = Kokkos::Details::ArithTraits<typename XVector::non_const_value_type>::abs(m_x(i,k));
  }
};

template<class RVector, class XVector>
RVector MV_Abs( const RVector & r, const XVector & x)
{
  // TODO: Add error check (didn't link for some reason?)
  /*if(r.dimension_0() != x.dimension_0())
    Kokkos::Impl::throw_runtime_exception("Kokkos::MV_Abs -- dimension(0) of r and x don't match");
  if(r.dimension_1() != x.dimension_1())
    Kokkos::Impl::throw_runtime_exception("Kokkos::MV_Abs -- dimension(1) of r and x don't match");*/

  //TODO: Get 1D version done
  /*if(r.dimension_1()==1) {
    typedef View<typename RVector::value_type*,typename RVector::execution_space> RVector1D;
    typedef View<typename XVector::const_value_type*,typename XVector::execution_space> XVector1D;

    RVector1D r_1d = Kokkos::subview( r , ALL(),0 );
    XVector1D x_1d = Kokkos::subview( x , ALL(),0 );
    return V_Abs(r_1d,x_1d);
  }*/
  if(r==x) {
    MV_AbsSelfFunctor<XVector> op(x,x.dimension_1()) ;
    Kokkos::parallel_for( x.dimension_0() , op );
    return r;
  }

  MV_AbsFunctor<RVector,XVector> op(r,x,x.dimension_1()) ;
  Kokkos::parallel_for( x.dimension_0() , op );
  return r;
}

/// \brief Functor for MultiVector::elementWiseMultiply.
///
/// C(i,j) = c * C(i,j) + ab * A(i) * B(i,j), subject to the usual BLAS update rules.
template<class CVector, class AVector, class BVector>
struct MV_ElementWiseMultiplyFunctor
{
  typedef typename CVector::execution_space        execution_space;
  typedef typename CVector::size_type            size_type;

  typename CVector::const_value_type m_c;
  CVector m_C;
  typename AVector::const_value_type m_ab;
  typename AVector::const_type m_A ;
  typename BVector::const_type m_B ;

  const size_type m_n;
  MV_ElementWiseMultiplyFunctor(
      typename CVector::const_value_type c,
      CVector C,
      typename AVector::const_value_type ab,
      typename AVector::const_type A,
      typename BVector::const_type B,
      const size_type n):
      m_c(c),m_C(C),m_ab(ab),m_A(A),m_B(B),m_n(n)
      {}

  KOKKOS_INLINE_FUNCTION void
  operator () (const size_type i) const
  {
    if (m_c == Kokkos::Details::ArithTraits<typename CVector::non_const_value_type>::zero ()) {
      if (m_ab == Kokkos::Details::ArithTraits<typename AVector::non_const_value_type>::zero ()) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
        for (size_type k = 0; k < m_n; ++k) {
          m_C(i,k) = Kokkos::Details::ArithTraits<typename CVector::non_const_value_type>::zero ();
        }
      }
      else { // m_ab != 0, but m_c == 0
        // BLAS update rules say that if m_c == 0, we must overwrite m_C.
        // This matters only if m_C has entries that are Inf or NaN.
        typename AVector::const_value_type Ai = m_A(i);
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
        for (size_type k = 0; k < m_n; ++k) {
          m_C(i,k) = m_ab * Ai * m_B(i,k);
        }
      }
    }
    else { // m_c != 0
      if (m_ab == Kokkos::Details::ArithTraits<typename AVector::non_const_value_type>::zero ()) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
        for (size_type k = 0; k < m_n; ++k) {
          m_C(i,k) = m_c * m_C(i,k);
        }
      }
      else { // m_ab != 0, and m_c != 0
        typename AVector::const_value_type Ai = m_A(i);
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
        for (size_type k = 0; k < m_n; ++k) {
          m_C(i,k) = m_c * m_C(i,k) + m_ab * Ai * m_B(i,k);
        }
      }
    }
  }
};


template<class CVector, class AVector, class BVector>
CVector MV_ElementWiseMultiply(
      typename CVector::const_value_type c,
      CVector C,
      typename AVector::const_value_type ab,
      AVector A,
      BVector B
    )
{
  // TODO: Add error check (didn't link for some reason?)
  /*if(r.dimension_0() != x.dimension_0())
    Kokkos::Impl::throw_runtime_exception("Kokkos::MV_ElementWiseMultiply -- dimension(0) of r and x don't match");
  if(r.dimension_1() != x.dimension_1())
    Kokkos::Impl::throw_runtime_exception("Kokkos::MV_ElementWiseMultiply -- dimension(1) of r and x don't match");*/

  //TODO: Get 1D version done
  /*if(r.dimension_1()==1) {
    typedef View<typename RVector::value_type*,typename RVector::execution_space> RVector1D;
    typedef View<typename XVector::const_value_type*,typename XVector::execution_space> XVector1D;

    RVector1D r_1d = Kokkos::subview( r , ALL(),0 );
    XVector1D x_1d = Kokkos::subview( x , ALL(),0 );
    return V_ElementWiseMultiply(r_1d,x_1d);
  }*/

  MV_ElementWiseMultiplyFunctor<CVector,AVector,BVector> op(c,C,ab,A,B,C.dimension_1()) ;
  Kokkos::parallel_for( C.dimension_0() , op );
  return C;
}

/*------------------------------------------------------------------------------------------
 *-------------------------- Vector Add: r = a*x + b*y -------------------------------------
 *------------------------------------------------------------------------------------------*/

/* Variants of Functors with a and b being vectors. */

//Unroll for n<=16
template<class RVector,class aVector, class XVector, class bVector, class YVector, int scalar_x, int scalar_y,int UNROLL>
struct MV_AddUnrollFunctor
{
  typedef typename RVector::execution_space        execution_space;
  typedef typename RVector::size_type            size_type;

  RVector   m_r ;
  XVector  m_x ;
  YVector   m_y ;
  aVector m_a;
  bVector m_b;
  size_type n;
  size_type start;

  MV_AddUnrollFunctor() {n=UNROLL;}
  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i ) const
  {
        if((scalar_x==1)&&(scalar_y==1)){
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
    for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = m_x(i,k) + m_y(i,k);
        }
        if((scalar_x==1)&&(scalar_y==-1)){
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = m_x(i,k) - m_y(i,k);
        }
        if((scalar_x==-1)&&(scalar_y==-1)){
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = -m_x(i,k) - m_y(i,k);
        }
        if((scalar_x==-1)&&(scalar_y==1)){
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = -m_x(i,k) + m_y(i,k);
        }
        if((scalar_x==2)&&(scalar_y==1)){
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = m_a(k)*m_x(i,k) + m_y(i,k);
        }
        if((scalar_x==2)&&(scalar_y==-1)){
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = m_a(k)*m_x(i,k) - m_y(i,k);
        }
        if((scalar_x==1)&&(scalar_y==2)){
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = m_x(i,k) + m_b(k)*m_y(i,k);
        }
        if((scalar_x==-1)&&(scalar_y==2)){
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = -m_x(i,k) + m_b(k)*m_y(i,k);
        }
        if((scalar_x==2)&&(scalar_y==2)){
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = m_a(k)*m_x(i,k) + m_b(k)*m_y(i,k);
        }
  }
};

template<class RVector,class aVector, class XVector, class bVector, class YVector, int scalar_x, int scalar_y>
struct MV_AddVectorFunctor
{
  typedef typename RVector::execution_space        execution_space;
  typedef typename RVector::size_type            size_type;

  RVector   m_r ;
  XVector  m_x ;
  YVector   m_y ;
  aVector m_a;
  bVector m_b;
  size_type n;

  MV_AddVectorFunctor() {n=1;}
  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i ) const
  {
        if((scalar_x==1)&&(scalar_y==1))
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for(size_type k=0;k<n;k++)
            m_r(i,k) = m_x(i,k) + m_y(i,k);
        if((scalar_x==1)&&(scalar_y==-1))
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for(size_type k=0;k<n;k++)
            m_r(i,k) = m_x(i,k) - m_y(i,k);
        if((scalar_x==-1)&&(scalar_y==-1))
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for(size_type k=0;k<n;k++)
            m_r(i,k) = -m_x(i,k) - m_y(i,k);
        if((scalar_x==-1)&&(scalar_y==1))
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for(size_type k=0;k<n;k++)
            m_r(i,k) = -m_x(i,k) + m_y(i,k);
        if((scalar_x==2)&&(scalar_y==1))
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for(size_type k=0;k<n;k++)
            m_r(i,k) = m_a(k)*m_x(i,k) + m_y(i,k);
        if((scalar_x==2)&&(scalar_y==-1))
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for(size_type k=0;k<n;k++)
            m_r(i,k) = m_a(k)*m_x(i,k) - m_y(i,k);
        if((scalar_x==1)&&(scalar_y==2))
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for(size_type k=0;k<n;k++)
            m_r(i,k) = m_x(i,k) + m_b(k)*m_y(i,k);
        if((scalar_x==-1)&&(scalar_y==2))
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for(size_type k=0;k<n;k++)
            m_r(i,k) = -m_x(i,k) + m_b(k)*m_y(i,k);
        if((scalar_x==2)&&(scalar_y==2))
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for(size_type k=0;k<n;k++)
            m_r(i,k) = m_a(k)*m_x(i,k) + m_b(k)*m_y(i,k);

  }
};

/* Variants of Functors with a and b being scalars. */

template<class RVector, class XVector, class YVector, int scalar_x, int scalar_y,int UNROLL>
struct MV_AddUnrollFunctor<RVector,typename XVector::non_const_value_type, XVector, typename YVector::non_const_value_type,YVector,scalar_x,scalar_y,UNROLL>
{
  typedef typename RVector::execution_space        execution_space;
  typedef typename RVector::size_type            size_type;

  RVector   m_r ;
  XVector  m_x ;
  YVector   m_y ;
  typename XVector::non_const_value_type m_a;
  typename YVector::non_const_value_type m_b;
  size_type n;
  size_type start;

  MV_AddUnrollFunctor() {n=UNROLL;}
  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i ) const
  {
  if((scalar_x==1)&&(scalar_y==1)){
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
    for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = m_x(i,k) + m_y(i,k);
  }
  if((scalar_x==1)&&(scalar_y==-1)){
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
    for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = m_x(i,k) - m_y(i,k);
  }
  if((scalar_x==-1)&&(scalar_y==-1)){
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = -m_x(i,k) - m_y(i,k);
  }
  if((scalar_x==-1)&&(scalar_y==1)){
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = -m_x(i,k) + m_y(i,k);
  }
  if((scalar_x==2)&&(scalar_y==1)){
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = m_a*m_x(i,k) + m_y(i,k);
  }
  if((scalar_x==2)&&(scalar_y==-1)){
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = m_a*m_x(i,k) - m_y(i,k);
  }
  if((scalar_x==1)&&(scalar_y==2)){
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = m_x(i,k) + m_b*m_y(i,k);
  }
  if((scalar_x==-1)&&(scalar_y==2)){
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = -m_x(i,k) + m_b*m_y(i,k);
  }
  if((scalar_x==2)&&(scalar_y==2)){
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = m_a*m_x(i,k) + m_b*m_y(i,k);
  }
  }
};

template<class RVector, class XVector, class YVector, int scalar_x, int scalar_y>
struct MV_AddVectorFunctor<RVector,typename XVector::non_const_value_type, XVector, typename YVector::non_const_value_type,YVector,scalar_x,scalar_y>
{
  typedef typename RVector::execution_space        execution_space;
  typedef typename RVector::size_type            size_type;

  RVector   m_r ;
  XVector  m_x ;
  YVector   m_y ;
  typename XVector::non_const_value_type m_a;
  typename YVector::non_const_value_type m_b;
  size_type n;

  MV_AddVectorFunctor() {n=1;}
  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i ) const
  {
  if((scalar_x==1)&&(scalar_y==1))
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for(size_type k=0;k<n;k++)
      m_r(i,k) = m_x(i,k) + m_y(i,k);
  if((scalar_x==1)&&(scalar_y==-1))
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for(size_type k=0;k<n;k++)
      m_r(i,k) = m_x(i,k) - m_y(i,k);
  if((scalar_x==-1)&&(scalar_y==-1))
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for(size_type k=0;k<n;k++)
      m_r(i,k) = -m_x(i,k) - m_y(i,k);
  if((scalar_x==-1)&&(scalar_y==1))
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for(size_type k=0;k<n;k++)
      m_r(i,k) = -m_x(i,k) + m_y(i,k);
  if((scalar_x==2)&&(scalar_y==1))
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for(size_type k=0;k<n;k++)
      m_r(i,k) = m_a*m_x(i,k) + m_y(i,k);
  if((scalar_x==2)&&(scalar_y==-1))
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for(size_type k=0;k<n;k++)
      m_r(i,k) = m_a*m_x(i,k) - m_y(i,k);
  if((scalar_x==1)&&(scalar_y==2))
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for(size_type k=0;k<n;k++)
      m_r(i,k) = m_x(i,k) + m_b*m_y(i,k);
  if((scalar_x==-1)&&(scalar_y==2))
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for(size_type k=0;k<n;k++)
      m_r(i,k) = -m_x(i,k) + m_b*m_y(i,k);
  if((scalar_x==2)&&(scalar_y==2))
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for(size_type k=0;k<n;k++)
      m_r(i,k) = m_a*m_x(i,k) + m_b*m_y(i,k);

  }
};

template<class RVector,class aVector, class XVector, class bVector, class YVector,int UNROLL>
RVector MV_AddUnroll( const RVector & r,const aVector &av,const XVector & x,
                const bVector &bv, const YVector & y, int n,
                int a=2,int b=2)
{
   if(a==1&&b==1) {
     MV_AddUnrollFunctor<RVector,aVector,XVector,bVector,YVector,1,1,UNROLL> op ;
     op.m_r = r ;
     op.m_x = x ;
     op.m_y = y ;
     op.m_a = av ;
     op.m_b = bv ;
     op.n = x.dimension(1);
     Kokkos::parallel_for( n , op );
     return r;
   }
   if(a==1&&b==-1) {
     MV_AddUnrollFunctor<RVector,aVector,XVector,bVector,YVector,1,-1,UNROLL> op ;
     op.m_r = r ;
     op.m_x = x ;
     op.m_y = y ;
     op.m_a = av ;
     op.m_b = bv ;
     op.n = x.dimension(1);
     Kokkos::parallel_for( n , op );
     return r;
   }
   if(a==-1&&b==1) {
     MV_AddUnrollFunctor<RVector,aVector,XVector,bVector,YVector,-1,1,UNROLL> op ;
     op.m_r = r ;
     op.m_x = x ;
     op.m_y = y ;
     op.m_a = av ;
     op.m_b = bv ;
     op.n = x.dimension(1);
     Kokkos::parallel_for( n , op );
     return r;
   }
   if(a==-1&&b==-1) {
     MV_AddUnrollFunctor<RVector,aVector,XVector,bVector,YVector,-1,-1,UNROLL> op ;
     op.m_r = r ;
     op.m_x = x ;
     op.m_y = y ;
     op.m_a = av ;
     op.m_b = bv ;
     op.n = x.dimension(1);
     Kokkos::parallel_for( n , op );
     return r;
   }
   if(a*a!=1&&b==1) {
     MV_AddUnrollFunctor<RVector,aVector,XVector,bVector,YVector,2,1,UNROLL> op ;
     op.m_r = r ;
     op.m_x = x ;
     op.m_y = y ;
     op.m_a = av ;
     op.m_b = bv ;
     op.n = x.dimension(1);
     Kokkos::parallel_for( n , op );
     return r;
   }
   if(a*a!=1&&b==-1) {
     MV_AddUnrollFunctor<RVector,aVector,XVector,bVector,YVector,2,-1,UNROLL> op ;
     op.m_r = r ;
     op.m_x = x ;
     op.m_y = y ;
     op.m_a = av ;
     op.m_b = bv ;
     op.n = x.dimension(1);
     Kokkos::parallel_for( n , op );
     return r;
   }
   if(a==1&&b*b!=1) {
     MV_AddUnrollFunctor<RVector,aVector,XVector,bVector,YVector,1,2,UNROLL> op ;
     op.m_r = r ;
     op.m_x = x ;
     op.m_y = y ;
     op.m_a = av ;
     op.m_b = bv ;
     op.n = x.dimension(1);
     Kokkos::parallel_for( n , op );
     return r;
   }
   if(a==-1&&b*b!=1) {
     MV_AddUnrollFunctor<RVector,aVector,XVector,bVector,YVector,-1,2,UNROLL> op ;
     op.m_r = r ;
     op.m_x = x ;
     op.m_y = y ;
     op.m_a = av ;
     op.m_b = bv ;
     op.n = x.dimension(1);
     Kokkos::parallel_for( n , op );
     return r;
   }
   MV_AddUnrollFunctor<RVector,aVector,XVector,bVector,YVector,2,2,UNROLL> op ;
   op.m_r = r ;
   op.m_x = x ;
   op.m_y = y ;
   op.m_a = av ;
   op.m_b = bv ;
   op.n = x.dimension(1);
   Kokkos::parallel_for( n , op );

   return r;
}

template<class RVector,class aVector, class XVector, class bVector, class YVector>
RVector MV_AddUnroll( const RVector & r,const aVector &av,const XVector & x,
                const bVector &bv, const YVector & y, int n,
                int a=2,int b=2)
{
        switch (x.dimension(1)){
      case 1: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 1>( r,av,x,bv,y,n,a,b);
                  break;
      case 2: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 2>( r,av,x,bv,y,n,a,b);
                  break;
      case 3: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 3>( r,av,x,bv,y,n,a,b);
                  break;
      case 4: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 4>( r,av,x,bv,y,n,a,b);
                  break;
      case 5: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 5>( r,av,x,bv,y,n,a,b);
                  break;
      case 6: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 6>( r,av,x,bv,y,n,a,b);
                  break;
      case 7: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 7>( r,av,x,bv,y,n,a,b);
                  break;
      case 8: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 8>( r,av,x,bv,y,n,a,b);
                  break;
      case 9: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 9>( r,av,x,bv,y,n,a,b);
                  break;
      case 10: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 10>( r,av,x,bv,y,n,a,b);
                  break;
      case 11: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 11>( r,av,x,bv,y,n,a,b);
                  break;
      case 12: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 12>( r,av,x,bv,y,n,a,b);
                  break;
      case 13: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 13>( r,av,x,bv,y,n,a,b);
                  break;
      case 14: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 14>( r,av,x,bv,y,n,a,b);
                  break;
      case 15: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 15>( r,av,x,bv,y,n,a,b);
                  break;
      case 16: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 16>( r,av,x,bv,y,n,a,b);
                  break;
        }
        return r;
}


template<class RVector,class aVector, class XVector, class bVector, class YVector>
RVector
MV_AddVector (const RVector& r,
              const aVector& av,
              const XVector& x,
              const bVector& bv,
              const YVector & y,
              const int n,
              const int a = 2,
              const int b = 2)
{
   if(a==1&&b==1) {
     MV_AddVectorFunctor<RVector,aVector,XVector,bVector,YVector,1,1> op;
     op.m_r = r;
     op.m_x = x;
     op.m_y = y;
     op.m_a = av;
     op.m_b = bv;
     op.n = x.dimension(1);
     Kokkos::parallel_for( n , op );
     return r;
   }
   if(a==1&&b==-1) {
     MV_AddVectorFunctor<RVector,aVector,XVector,bVector,YVector,1,-1> op;
     op.m_r = r;
     op.m_x = x;
     op.m_y = y;
     op.m_a = av;
     op.m_b = bv;
     op.n = x.dimension(1);
     Kokkos::parallel_for( n , op );
     return r;
   }
   if(a==-1&&b==1) {
     MV_AddVectorFunctor<RVector,aVector,XVector,bVector,YVector,-1,1> op;
     op.m_r = r;
     op.m_x = x;
     op.m_y = y;
     op.m_a = av;
     op.m_b = bv;
     op.n = x.dimension(1);
     Kokkos::parallel_for( n , op );
     return r;
   }
   if(a==-1&&b==-1) {
     MV_AddVectorFunctor<RVector,aVector,XVector,bVector,YVector,-1,-1> op;
     op.m_r = r;
     op.m_x = x;
     op.m_y = y;
     op.m_a = av;
     op.m_b = bv;
     op.n = x.dimension(1);
     Kokkos::parallel_for( n , op );
     return r;
   }
   if(a*a!=1&&b==1) {
     MV_AddVectorFunctor<RVector,aVector,XVector,bVector,YVector,2,1> op;
     op.m_r = r;
     op.m_x = x;
     op.m_y = y;
     op.m_a = av;
     op.m_b = bv;
     op.n = x.dimension(1);
     Kokkos::parallel_for( n , op );
     return r;
   }
   if(a*a!=1&&b==-1) {
     MV_AddVectorFunctor<RVector,aVector,XVector,bVector,YVector,2,-1> op;
     op.m_r = r;
     op.m_x = x;
     op.m_y = y;
     op.m_a = av;
     op.m_b = bv;
     op.n = x.dimension(1);
     Kokkos::parallel_for( n , op );
     return r;
   }
   if(a==1&&b*b!=1) {
     MV_AddVectorFunctor<RVector,aVector,XVector,bVector,YVector,1,2> op;
     op.m_r = r;
     op.m_x = x;
     op.m_y = y;
     op.m_a = av;
     op.m_b = bv;
     op.n = x.dimension(1);
     Kokkos::parallel_for( n , op );
     return r;
   }
   if(a==-1&&b*b!=1) {
     MV_AddVectorFunctor<RVector,aVector,XVector,bVector,YVector,-1,2> op;
     op.m_r = r;
     op.m_x = x;
     op.m_y = y;
     op.m_a = av;
     op.m_b = bv;
     op.n = x.dimension(1);
     Kokkos::parallel_for( n , op );
     return r;
   }
   MV_AddVectorFunctor<RVector,aVector,XVector,bVector,YVector,2,2> op;
   op.m_r = r;
   op.m_x = x;
   op.m_y = y;
   op.m_a = av;
   op.m_b = bv;
   op.n = x.dimension(1);
   Kokkos::parallel_for( n , op );

   return r;
}

template<class RVector,class aVector, class XVector, class bVector, class YVector>
RVector MV_AddSpecialise( const RVector & r,const aVector &av,const XVector & x,
                const bVector &bv, const YVector & y, unsigned int n,
                int a=2,int b=2)
{

        if(x.dimension(1)>16)
                return MV_AddVector( r,av,x,bv,y,a,b);

        if(x.dimension_1()==1) {
    typedef View<typename RVector::value_type*,typename RVector::execution_space> RVector1D;
    typedef View<typename XVector::const_value_type*,typename XVector::execution_space> XVector1D;
    typedef View<typename YVector::const_value_type*,typename YVector::execution_space> YVector1D;

    RVector1D r_1d = Kokkos::subview( r , ALL(),0 );
    XVector1D x_1d = Kokkos::subview( x , ALL(),0 );
    YVector1D y_1d = Kokkos::subview( y , ALL(),0 );

    V_Add(r_1d,av,x_1d,bv,y_1d,n);
    return r;
  } else
        return MV_AddUnroll( r,av,x,bv,y,a,b);
}

template<class RVector,class aVector, class XVector, class bVector, class YVector>
RVector
MV_Add (const RVector& r,
        const aVector& av,
        const XVector& x,
        const bVector& bv,
        const YVector& y,
        int n = -1)
{
  if (n == -1) {
    n = x.dimension_0 ();
  }
  if (x.dimension (1) > 16) {
    return MV_AddVector (r, av, x, bv, y, n, 2, 2);
  }

  if (x.dimension_1 () == 1) {
    typedef View<typename RVector::value_type*, typename RVector::execution_space> RVector1D;
    typedef View<typename XVector::const_value_type*, typename XVector::execution_space> XVector1D;
    typedef View<typename YVector::const_value_type*, typename YVector::execution_space> YVector1D;

    RVector1D r_1d = subview (r, ALL(), 0);
    XVector1D x_1d = subview (x, ALL(), 0);
    YVector1D y_1d = subview (y, ALL(), 0);

    V_Add (r_1d, av, x_1d, bv, y_1d, n);
    return r;
  } else {
    return MV_AddUnroll (r, av, x, bv, y, n, 2, 2);
  }
}

template<class RVector,class XVector,class YVector>
RVector
MV_Add (const RVector& r,
        const XVector& x,
        const YVector& y,
        int n = -1)
{
  if (n == -1) {
    n = x.dimension_0 ();
  }
  if (x.dimension_1 () == 1) {
    typedef View<typename RVector::value_type*, typename RVector::execution_space> RVector1D;
    typedef View<typename XVector::const_value_type*, typename XVector::execution_space> XVector1D;
    typedef View<typename YVector::const_value_type*, typename YVector::execution_space> YVector1D;

    RVector1D r_1d = subview (r , ALL(), 0);
    XVector1D x_1d = subview (x , ALL(), 0);
    YVector1D y_1d = subview (y , ALL(), 0);

    V_Add (r_1d, x_1d, y_1d, n);
    return r;
  } else {
    typename XVector::const_value_type a =
      Kokkos::Details::ArithTraits<typename XVector::non_const_value_type>::one ();
    return MV_AddSpecialise (r, a, x, a, y, n, 1, 1);
  }
}

template<class RVector,class XVector,class bVector, class YVector>
RVector MV_Add (const RVector& r,
                const XVector& x,
                const bVector& bv,
                const YVector& y,
                int n = -1)
{
  if (n == -1) {
    n = x.dimension_0 ();
  }
  if (x.dimension_1 () == 1) {
    typedef View<typename RVector::value_type*, typename RVector::execution_space> RVector1D;
    typedef View<typename XVector::const_value_type*, typename XVector::execution_space> XVector1D;
    typedef View<typename YVector::const_value_type*, typename YVector::execution_space> YVector1D;

    RVector1D r_1d = subview (r, ALL (), 0);
    XVector1D x_1d = subview (x, ALL (), 0);
    YVector1D y_1d = subview (y, ALL (), 0);

    V_Add (r_1d, x_1d, bv, y_1d, n);
    return r;
  } else {
    MV_AddSpecialise (r, bv, x, bv, y, n, 1, 2);
  }
}


template<class XVector,class YVector>
struct MV_DotProduct_Right_FunctorVector
{
  typedef typename XVector::execution_space         execution_space;
  typedef typename XVector::size_type             size_type;
  typedef typename XVector::non_const_value_type          xvalue_type;
  typedef Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef typename IPT::dot_type               value_type[];
  size_type value_count;


  typedef typename XVector::const_type        x_const_type;
  typedef typename YVector::const_type        y_const_type;
  x_const_type  m_x ;
  y_const_type  m_y ;

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i, value_type sum ) const
  {
    const size_type numVecs=value_count;

#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for(size_type k=0;k<numVecs;k++)
      sum[k]+=IPT::dot( m_x(i,k), m_y(i,k) );  // m_x(i,k) * m_y(i,k)
  }
  KOKKOS_INLINE_FUNCTION void init( value_type update) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for(size_type k=0;k<numVecs;k++)
      update[k] = 0;
  }
  KOKKOS_INLINE_FUNCTION void join( volatile value_type  update ,
                                    const volatile value_type  source ) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for(size_type k=0;k<numVecs;k++){
      update[k] += source[k];
    }
  }
};


// Implementation detail of Tpetra::MultiVector::dot, when both
// MultiVectors in the dot product have constant stride.  Compute the
// dot product of the local part of each corresponding vector (column)
// in two MultiVectors.
template<class MultiVecViewType>
struct MultiVecDotFunctor {
  typedef typename MultiVecViewType::execution_space execution_space;
  typedef typename MultiVecViewType::size_type size_type;
  typedef typename MultiVecViewType::value_type mv_value_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<mv_value_type> IPT;
  typedef typename IPT::dot_type value_type[];

  typedef MultiVecViewType mv_view_type;
  typedef typename MultiVecViewType::const_type mv_const_view_type;
  typedef Kokkos::View<typename IPT::dot_type*, execution_space> dot_view_type;

  mv_const_view_type X_, Y_;
  dot_view_type dots_;
  // Kokkos::parallel_reduce wants this, so it needs to be public.
  size_type value_count;

  MultiVecDotFunctor (const mv_const_view_type& X,
                      const mv_const_view_type& Y,
                      const dot_view_type& dots) :
    X_ (X), Y_ (Y), dots_ (dots), value_count (X.dimension_1 ())
  {
    if (value_count != dots.dimension_0 ()) {
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
      Kokkos::abort("Kokkos::MultiVecDotFunctor: value_count does not match the length of 'dots'");
#else
      std::ostringstream os;
      os << "Kokkos::MultiVecDotFunctor: value_count does not match the length "
        "of 'dots'.  X is " << X.dimension_0 () << " x " << X.dimension_1 () <<
        ", Y is " << Y.dimension_0 () << " x " << Y.dimension_1 () << ", "
        "dots has length " << dots.dimension_0 () << ", and value_count = " <<
        value_count << ".";
      throw std::invalid_argument (os.str ());
#endif
    }
  }

  KOKKOS_INLINE_FUNCTION void
  operator() (const size_type i, value_type sum) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      sum[k] += IPT::dot (X_(i,k), Y_(i,k));
    }
  }

  KOKKOS_INLINE_FUNCTION void
  init (value_type update) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      update[k] = Kokkos::Details::ArithTraits<typename IPT::dot_type>::zero ();
    }
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type update,
        const volatile value_type source) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      update[k] += source[k];
    }
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void
  final (const value_type dst) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      dots_(k) = dst[k];
    }
  }
};


// Implementation detail of Tpetra::MultiVector::norm2, when the
// MultiVector has constant stride.  Compute the square of the
// two-norm of each column of a multivector.
template<class MultiVecViewType, class NormsViewType>
struct MultiVecNorm2SquaredFunctor {
  typedef typename MultiVecViewType::execution_space execution_space;
  typedef typename MultiVecViewType::size_type size_type;
  typedef typename MultiVecViewType::value_type mv_value_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<mv_value_type> IPT;
  typedef typename IPT::mag_type value_type[];

  typedef MultiVecViewType mv_view_type;
  typedef typename MultiVecViewType::const_type mv_const_view_type;
  typedef NormsViewType norms_view_type;

  mv_const_view_type X_;
  norms_view_type norms_;
  // Kokkos::parallel_reduce wants this, so it needs to be public.
  size_type value_count;

  MultiVecNorm2SquaredFunctor (const mv_const_view_type& X,
                               const norms_view_type& norms) :
    X_ (X), norms_ (norms), value_count (X.dimension_1 ())
  {
    if (value_count != norms.dimension_0 ()) {
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
      Kokkos::abort("Kokkos::MultiVecNorm2SquaredFunctor: value_count does not match the length of 'norms'");
#else
      std::ostringstream os;
      os << "Kokkos::MultiVecNorm2SquaredFunctor: value_count does not match "
        "the length of 'norms'.  X is " << X.dimension_0 () << " x " <<
        X.dimension_1 () << ", norms has length " << norms.dimension_0 () <<
        ", and value_count = " << value_count << ".";
      throw std::invalid_argument (os.str ());
#endif
    }
  }

  KOKKOS_INLINE_FUNCTION void
  operator() (const size_type i, value_type sum) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      const typename IPT::mag_type tmp = IPT::norm (X_(i,k));
      sum[k] += tmp * tmp;
    }
  }

  KOKKOS_INLINE_FUNCTION void
  init (value_type update) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      update[k] = Kokkos::Details::ArithTraits<typename IPT::mag_type>::zero ();
    }
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type update,
        const volatile value_type source) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      update[k] += source[k];
    }
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void
  final (const value_type dst) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      norms_(k) = dst[k];
    }
  }
};


// Implementation detail of Tpetra::MultiVector::norm1, when the
// MultiVector has constant stride.  Compute the one-norm of each
// column of a multivector.
template<class MultiVecViewType, class NormsViewType>
struct MultiVecNorm1Functor {
  typedef typename MultiVecViewType::execution_space execution_space;
  typedef typename MultiVecViewType::size_type size_type;
  typedef typename MultiVecViewType::value_type mv_value_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<mv_value_type> IPT;
  typedef typename IPT::mag_type value_type[];

  typedef MultiVecViewType mv_view_type;
  typedef typename MultiVecViewType::const_type mv_const_view_type;
  typedef NormsViewType norms_view_type;

  mv_const_view_type X_;
  norms_view_type norms_;
  // Kokkos::parallel_reduce wants this, so it needs to be public.
  size_type value_count;

  MultiVecNorm1Functor (const mv_const_view_type& X,
                        const norms_view_type& norms) :
    X_ (X), norms_ (norms), value_count (X.dimension_1 ())
  {
    if (value_count != norms.dimension_0 ()) {
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
      Kokkos::abort("Kokkos::MultiVecNorm1Functor: value_count does not match the length of 'norms'");
#else
      std::ostringstream os;
      os << "Kokkos::MultiVecNorm1Functor: value_count does not match the "
         << "length of 'norms'.  X is " << X.dimension_0 () << " x "
         << X.dimension_1 () << ", norms has length " << norms.dimension_0 ()
         << ", and value_count = " << value_count << ".";
      throw std::invalid_argument (os.str ());
#endif
    }
  }

  KOKKOS_INLINE_FUNCTION void
  operator() (const size_type i, value_type sum) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      sum[k] += IPT::norm (X_(i,k)); // absolute value
    }
  }

  KOKKOS_INLINE_FUNCTION void
  init (value_type update) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      update[k] = Kokkos::Details::ArithTraits<typename IPT::mag_type>::zero ();
    }
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type update,
        const volatile value_type source) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      update[k] += source[k];
    }
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void
  final (const value_type dst) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      norms_(k) = dst[k];
    }
  }
};


// Implementation detail of Tpetra::MultiVector::normInf, when the
// MultiVector has constant stride.  Compute the infinity-norm of each
// column of a multivector.
template<class MultiVecViewType, class NormsViewType>
struct MultiVecNormInfFunctor {
  typedef typename MultiVecViewType::execution_space execution_space;
  typedef typename MultiVecViewType::size_type size_type;
  typedef typename MultiVecViewType::value_type mv_value_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<mv_value_type> IPT;
  typedef typename IPT::mag_type value_type[];

  typedef MultiVecViewType mv_view_type;
  typedef typename MultiVecViewType::const_type mv_const_view_type;
  typedef NormsViewType norms_view_type;

  mv_const_view_type X_;
  norms_view_type norms_;
  // Kokkos::parallel_reduce wants this, so it needs to be public.
  size_type value_count;

  MultiVecNormInfFunctor (const mv_const_view_type& X,
                          const norms_view_type& norms) :
    X_ (X), norms_ (norms), value_count (X.dimension_1 ())
  {
    if (value_count != norms.dimension_0 ()) {
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
      Kokkos::abort("Kokkos::MultiVecNormInfFunctor: value_count does not match the length of 'norms'");
#else
      std::ostringstream os;
      os << "Kokkos::MultiVecNormInfFunctor: value_count does not match the "
         << "length of 'norms'.  X is " << X.dimension_0 () << " x "
         << X.dimension_1 () << ", norms has length " << norms.dimension_0 ()
         << ", and value_count = " << value_count << ".";
      throw std::invalid_argument (os.str ());
#endif
    }
  }

  KOKKOS_INLINE_FUNCTION void
  operator() (const size_type i, value_type maxes) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      const typename IPT::mag_type curVal = maxes[k];
      const typename IPT::mag_type newVal = IPT::norm (X_(i,k));
      // max(curVal, newVal).  Any comparison predicate involving NaN
      // evaluates to false.  Thus, this will never assign NaN to
      // update[k], unless it contains NaN already.  The initial value
      // is zero, so NaNs won't propagate.  (This definition makes NaN
      // into an "invalid value," which is useful for statistics and
      // other applications that use NaN to indicate a "hole.")
      if (curVal < newVal) {
        maxes[k] = newVal;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION void
  init (value_type update) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      // Zero is a good default value for magnitudes (which are
      // nonnegative by definition).  That way, MPI processes with
      // zero rows won't affect the global maximum.
      update[k] = Kokkos::Details::ArithTraits<typename IPT::mag_type>::zero ();
    }
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type update,
        const volatile value_type source) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      const typename IPT::mag_type curVal = update[k];
      const typename IPT::mag_type newVal = source[k];
      // max(curVal, newVal).  Any comparison predicate involving NaN
      // evaluates to false.  Thus, this will never assign NaN to
      // update[k], unless it contains NaN already.  The initial value
      // is zero, so NaNs won't propagate.  (This definition makes NaN
      // into an "invalid value," which is useful for statistics and
      // other applications that use NaN to indicate a "hole.")
      if (curVal < newVal) {
        update[k] = newVal;
      }
    }
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void
  final (const value_type dst) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      norms_(k) = dst[k];
    }
  }
};


// Implementation detail of Tpetra::MultiVector::dot, for single
// vectors (columns).
template<class VecViewType>
struct VecDotFunctor {
  typedef typename VecViewType::execution_space execution_space;
  typedef typename VecViewType::size_type size_type;
  typedef typename VecViewType::value_type mv_value_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<mv_value_type> IPT;
  typedef typename IPT::dot_type value_type;
  typedef typename VecViewType::const_type vec_const_view_type;
  // This is a nonconst scalar view.  It holds one dot_type instance.
  typedef Kokkos::View<typename IPT::dot_type, execution_space> dot_view_type;
  typedef Kokkos::View<typename IPT::dot_type*, execution_space> dots_view_type;

  vec_const_view_type x_, y_;
  dot_view_type dot_;

  VecDotFunctor (const vec_const_view_type& x,
                 const vec_const_view_type& y,
                 const dot_view_type& dot) :
    x_ (x), y_ (y), dot_ (dot)
  {
    if (x.dimension_0 () != y.dimension_0 ()) {
      std::ostringstream os;
      os << "Kokkos::VecDotFunctor: The dimensions of x and y do not match.  "
        "x.dimension_0() = " << x.dimension_0 ()
         << " != y.dimension_0() = " << y.dimension_0 () << ".";
      throw std::invalid_argument (os.str ());
    }
  }

  KOKKOS_INLINE_FUNCTION void
  operator() (const size_type i, value_type& sum) const {
    sum += IPT::dot (x_(i), y_(i));
  }

  KOKKOS_INLINE_FUNCTION void
  init (value_type& update) const {
    update = Details::ArithTraits<value_type>::zero ();
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& update,
        const volatile value_type& source) const {
    update += source;
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void final (value_type& dst) const {
    // BADNESS HERE
    dot_() = dst;
  }
};


// Compute the square of the two-norm of a single vector.
template<class VecViewType, class NormViewType>
struct VecNorm2SquaredFunctor {
  typedef typename VecViewType::execution_space execution_space;
  typedef typename VecViewType::size_type size_type;
  typedef typename VecViewType::value_type mv_value_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<mv_value_type> IPT;
  typedef typename IPT::mag_type value_type;
  typedef typename VecViewType::const_type vec_const_view_type;
  // This is a nonconst scalar view.  It holds one mag_type instance.
  typedef NormViewType norm_view_type;

  vec_const_view_type x_;
  norm_view_type norm_;

  // Constructor
  //
  // x: the vector for which to compute the square of the two-norm.
  // norm: scalar View into which to put the result.
  VecNorm2SquaredFunctor (const vec_const_view_type& x,
                          const norm_view_type& norm) :
    x_ (x), norm_ (norm)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator() (const size_type i, value_type& sum) const {
    const typename IPT::mag_type tmp = IPT::norm (x_(i));
    sum += tmp * tmp;
  }

  KOKKOS_INLINE_FUNCTION void
  init (value_type& update) const {
    update = Details::ArithTraits<value_type>::zero ();
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& update,
        const volatile value_type& source) const {
    update += source;
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void final (value_type& dst) const {
    norm_ () = dst;
  }
};


// Compute the square of the one-norm of a single vector.
template<class VecViewType, class NormViewType>
struct VecNorm1Functor {
  typedef typename VecViewType::execution_space execution_space;
  typedef typename VecViewType::size_type size_type;
  typedef typename VecViewType::value_type mv_value_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<mv_value_type> IPT;
  typedef typename IPT::mag_type value_type;
  typedef typename VecViewType::const_type vec_const_view_type;
  // This is a nonconst scalar view.  It holds one mag_type instance.
  typedef NormViewType norm_view_type;

  vec_const_view_type x_;
  norm_view_type norm_;

  // Constructor
  //
  // x: the vector for which to compute the one-norm.
  // norm: scalar View into which to put the result.
  VecNorm1Functor (const vec_const_view_type& x,
                   const norm_view_type& norm) :
    x_ (x), norm_ (norm)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator() (const size_type i, value_type& sum) const {
    sum += IPT::norm (x_(i)); // absolute value
  }

  KOKKOS_INLINE_FUNCTION void
  init (value_type& update) const {
    update = Details::ArithTraits<value_type>::zero ();
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& update,
        const volatile value_type& source) const {
    update += source;
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void final (value_type& dst) const {
    norm_ () = dst;
  }
};


// Compute the square of the infinity-norm of a single vector.
template<class VecViewType, class NormViewType>
struct VecNormInfFunctor {
  typedef typename VecViewType::execution_space execution_space;
  typedef typename VecViewType::size_type size_type;
  typedef typename VecViewType::value_type mv_value_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<mv_value_type> IPT;
  typedef typename IPT::mag_type value_type;
  typedef typename VecViewType::const_type vec_const_view_type;
  // This is a nonconst scalar view.  It holds one mag_type instance.
  typedef NormViewType norm_view_type;

  vec_const_view_type x_;
  norm_view_type norm_;

  // Constructor
  //
  // x: the vector for which to compute the infinity-norm.
  // norm: scalar View into which to put the result.
  VecNormInfFunctor (const vec_const_view_type& x,
                     const norm_view_type& norm) :
    x_ (x), norm_ (norm)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator() (const size_type i, value_type& curVal) const {
    const typename IPT::mag_type newVal = IPT::norm (x_(i));
    // max(curVal, newVal).  Any comparison predicate involving NaN
    // evaluates to false.  Thus, this will never assign NaN to
    // update[k], unless it contains NaN already.  The initial value
    // is zero, so NaNs won't propagate.  (This definition makes NaN
    // into an "invalid value," which is useful for statistics and
    // other applications that use NaN to indicate a "hole.")
    if (curVal < newVal) {
      curVal = newVal;
    }
  }

  KOKKOS_INLINE_FUNCTION void
  init (value_type& update) const {
    // Zero is a good default value for magnitudes (which are
    // nonnegative by definition).  That way, MPI processes with
    // zero rows won't affect the global maximum.
    update = Details::ArithTraits<value_type>::zero ();
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& update,
        const volatile value_type& source) const {
    // max(update, source).  Any comparison predicate involving NaN
    // evaluates to false.  Thus, this will never assign NaN to
    // update, unless it contains NaN already.  The initial value is
    // zero, so NaNs won't propagate.  (This definition makes NaN into
    // an "invalid value," which is useful for statistics and other
    // applications that use NaN to indicate a "hole.")
    if (update < source) {
      update = source;
    }
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void final (value_type& dst) const {
    norm_ () = dst;
  }
};



// parallel_for functor for computing the square root, in place, of a
// one-dimensional View.  This is useful for following the MPI
// all-reduce that computes the square of the two-norms of the local
// columns of a Tpetra::MultiVector.
//
// mfh 14 Jul 2014: Carter says that, for now, the standard idiom for
// operating on a single scalar value on the device, is to run in a
// parallel_for with N = 1.
//
// FIXME (mfh 14 Jul 2014): If we could assume C++11, this functor
// would go away.
template<class ViewType>
class SquareRootFunctor {
public:
  typedef typename ViewType::execution_space execution_space;
  typedef typename ViewType::size_type size_type;

  SquareRootFunctor (const ViewType& theView) : theView_ (theView) {}

  KOKKOS_INLINE_FUNCTION void operator() (const size_type i) const {
    typedef typename ViewType::value_type value_type;
    theView_(i) = Kokkos::Details::ArithTraits<value_type>::sqrt (theView_(i));
  }

private:
  ViewType theView_;
};


template<class XVector,class YVector,int UNROLL>
struct MV_DotProduct_Right_FunctorUnroll
{
  typedef typename XVector::execution_space         execution_space;
  typedef typename XVector::size_type             size_type;
  typedef typename XVector::non_const_value_type          xvalue_type;
  typedef Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef typename IPT::dot_type               value_type[];
  size_type value_count;

  typedef typename XVector::const_type        x_const_type;
  typedef typename YVector::const_type        y_const_type;

  x_const_type  m_x ;
  y_const_type  m_y ;

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i, value_type sum ) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
    for(size_type k=0;k<UNROLL;k++)
      sum[k]+= IPT::dot( m_x(i,k), m_y(i,k) );  // m_x(i,k) * m_y(i,k)
  }
  KOKKOS_INLINE_FUNCTION void init( volatile value_type update) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
    for(size_type k=0;k<UNROLL;k++)
      update[k] = 0;
  }
  KOKKOS_INLINE_FUNCTION void join( volatile value_type update ,
                                    const volatile value_type source) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
    for(size_type k=0;k<UNROLL;k++)
      update[k] += source[k] ;
  }
};

/// \brief Compute the dot product(s) of the column(s) of the
///   multivectors (2-D arrays) x and y, and store the result(s) in r.
template<class rVector, class XVector, class YVector>
rVector
MV_Dot (const rVector& r,
        const XVector& x,
        const YVector& y)
{
  const typename XVector::size_type numRows = x.dimension_0 ();
  MV_Dot<rVector, XVector, YVector> (r, x, y, numRows);
}

/// \brief Compute the dot product(s) of the column(s) of the
///   multivectors (2-D arrays) x and y, and store the result(s) in r.
///   Only use the first numRows rows of x and y.
template<class rVector, class XVector, class YVector>
rVector
MV_Dot (const rVector& r,
        const XVector& x,
        const YVector& y,
        typename XVector::size_type numRows)
{
  typedef typename XVector::size_type size_type;
  const size_type numVecs = x.dimension(1);

  if (numVecs > 16) {
    MV_DotProduct_Right_FunctorVector<XVector,YVector> op;
    op.m_x = x;
    op.m_y = y;
    op.value_count = numVecs;

    Kokkos::parallel_reduce (numRows, op, r);
    return r;
  }
  else {
    switch (numVecs) {
    case 16: {
      MV_DotProduct_Right_FunctorUnroll<XVector,YVector,16> op;
      op.m_x = x;
      op.m_y = y;
      op.value_count = numVecs;
      Kokkos::parallel_reduce (numRows, op, r);
      break;
    }
    case 15: {
      MV_DotProduct_Right_FunctorUnroll<XVector,YVector,15> op;
      op.m_x = x;
      op.m_y = y;
      op.value_count = numVecs;
      Kokkos::parallel_reduce (numRows, op, r);
      break;
    }
    case 14: {
      MV_DotProduct_Right_FunctorUnroll<XVector,YVector,14> op;
      op.m_x = x;
      op.m_y = y;
      op.value_count = numVecs;
      Kokkos::parallel_reduce (numRows, op, r);
      break;
    }
    case 13: {
      MV_DotProduct_Right_FunctorUnroll<XVector,YVector,13> op;
      op.m_x = x;
      op.m_y = y;
      op.value_count = numVecs;
      Kokkos::parallel_reduce (numRows, op, r);
      break;
    }
    case 12: {
      MV_DotProduct_Right_FunctorUnroll<XVector,YVector,12> op;
      op.m_x = x;
      op.m_y = y;
      op.value_count = numVecs;
      Kokkos::parallel_reduce (numRows, op, r);
      break;
    }
    case 11: {
      MV_DotProduct_Right_FunctorUnroll<XVector,YVector,11> op;
      op.m_x = x;
      op.m_y = y;
      op.value_count = numVecs;
      Kokkos::parallel_reduce (numRows, op, r);
      break;
    }
    case 10: {
      MV_DotProduct_Right_FunctorUnroll<XVector,YVector,10> op;
      op.m_x = x;
      op.m_y = y;
      op.value_count = numVecs;
      Kokkos::parallel_reduce (numRows, op, r);
      break;
    }
    case 9: {
      MV_DotProduct_Right_FunctorUnroll<XVector,YVector,9> op;
      op.m_x = x;
      op.m_y = y;
      op.value_count = numVecs;
      Kokkos::parallel_reduce (numRows, op, r);
      break;
    }
    case 8: {
      MV_DotProduct_Right_FunctorUnroll<XVector,YVector,8> op;
      op.m_x = x;
      op.m_y = y;
      op.value_count = numVecs;
      Kokkos::parallel_reduce (numRows, op, r);
      break;
    }
    case 7: {
      MV_DotProduct_Right_FunctorUnroll<XVector,YVector,7> op;
      op.m_x = x;
      op.m_y = y;
      op.value_count = numVecs;
      Kokkos::parallel_reduce (numRows, op, r);
      break;
    }
    case 6: {
      MV_DotProduct_Right_FunctorUnroll<XVector,YVector,6> op;
      op.m_x = x;
      op.m_y = y;
      op.value_count = numVecs;
      Kokkos::parallel_reduce (numRows, op, r);
      break;
    }
    case 5: {
      MV_DotProduct_Right_FunctorUnroll<XVector,YVector,5> op;
      op.m_x = x;
      op.m_y = y;
      op.value_count = numVecs;
      Kokkos::parallel_reduce (numRows, op, r);
      break;
    }
    case 4: {
      MV_DotProduct_Right_FunctorUnroll<XVector,YVector,4> op;
      op.m_x = x;
      op.m_y = y;
      op.value_count = numVecs;
      Kokkos::parallel_reduce (numRows, op, r);

      break;
    }
    case 3: {
      MV_DotProduct_Right_FunctorUnroll<XVector,YVector,3> op;
      op.m_x = x;
      op.m_y = y;
      op.value_count = numVecs;
      Kokkos::parallel_reduce (numRows, op, r);
      break;
    }
    case 2: {
      MV_DotProduct_Right_FunctorUnroll<XVector,YVector,2> op;
      op.m_x = x;
      op.m_y = y;
      op.value_count = numVecs;
      Kokkos::parallel_reduce (numRows, op, r);
      break;
    }
    case 1: {
      // For the single-vector case, use the single-vector kernel for
      // better performance.  (It's usually faster to use fewer
      // indices when looking up View entries.)
      typedef View<typename XVector::const_value_type*,
                   typename XVector::execution_space> XVector1D;
      typedef View<typename YVector::const_value_type*,
                   typename YVector::execution_space> YVector1D;

      XVector1D x_1d = Kokkos::subview (x, ALL (), 0);
      YVector1D y_1d = Kokkos::subview (y, ALL (), 0);
      r[0] = V_Dot (x_1d, y_1d, numRows);
      break;
    }
    } // switch
  } // if-else
  return r;
}

/*------------------------------------------------------------------------------------------
 *-------------------------- Compute Sum -------------------------------------------------
 *------------------------------------------------------------------------------------------*/
template<class XVector>
struct MV_Sum_Functor
{
  typedef typename XVector::execution_space        execution_space;
  typedef typename XVector::size_type            size_type;
  typedef typename XVector::non_const_value_type          xvalue_type;
  typedef xvalue_type                         value_type[];

  typename XVector::const_type m_x ;
  size_type value_count;

  MV_Sum_Functor(XVector x):m_x(x),value_count(x.dimension_1()) {}
  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i, value_type sum ) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for(size_type k=0;k<value_count;k++){
      sum[k] += m_x(i,k);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type update) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for(size_type k=0;k<value_count;k++)
      update[k] = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type  update ,
                                    const volatile value_type  source ) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for(size_type k=0;k<value_count;k++){
      update[k] += source[k];
    }
  }
};


template<class normVector, class VectorType>
normVector MV_Sum(const normVector &r, const VectorType & x, int n = -1)
{
  if (n < 0) {
    n = x.dimension_0 ();
  }

  Kokkos::parallel_reduce (n , MV_Sum_Functor<VectorType> (x), r);
  return r;
}

/*------------------------------------------------------------------------------------------
 *-------------------------- Compute Norm1--------------------------------------------------
 *------------------------------------------------------------------------------------------*/
template<class XVector>
struct MV_Norm1_Functor
{
  typedef typename XVector::execution_space        execution_space;
  typedef typename XVector::size_type            size_type;
  typedef typename XVector::non_const_value_type          xvalue_type;
  typedef Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef typename IPT::dot_type               value_type[];

  typename XVector::const_type m_x ;
  size_type value_count;

  MV_Norm1_Functor(XVector x):m_x(x),value_count(x.dimension_1()) {}
  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i, value_type sum ) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for(size_type k=0;k<value_count;k++){
      sum[k] += Kokkos::Details::ArithTraits<typename XVector::non_const_value_type>::abs(m_x(i,k));
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type update) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for(size_type k=0;k<value_count;k++)
      update[k] = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type  update ,
                                    const volatile value_type  source ) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for(size_type k=0;k<value_count;k++){
      update[k] += source[k];
    }
  }
};

template<class normVector, class VectorType>
normVector MV_Norm1(const normVector &r, const VectorType & x, int n = -1)
{
  if (n < 0) {
    n = x.dimension_0 ();
  }

  Kokkos::parallel_reduce (n , MV_Norm1_Functor<VectorType> (x), r);
  return r;
}

/*------------------------------------------------------------------------------------------
 *-------------------------- Compute NormInf--------------------------------------------------
 *------------------------------------------------------------------------------------------*/
template<class XVector>
struct MV_NormInf_Functor
{
  typedef typename XVector::execution_space             execution_space;
  typedef typename XVector::size_type                   size_type;
  typedef typename XVector::non_const_value_type        xvalue_type;
  typedef Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef typename IPT::mag_type                        mag_type;
  typedef mag_type                                      value_type[];

  typename XVector::const_type m_x;
  size_type value_count;

  MV_NormInf_Functor (const XVector& x) :
    m_x(x), value_count (x.dimension_1 ())
  {}

  KOKKOS_INLINE_FUNCTION mag_type
  max (const mag_type& x, const mag_type& y) const {
    return (x > y) ? x : y;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i, value_type sum) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < value_count; ++k) {
      sum[k] = max (sum[k], IPT::norm (m_x(i,k)));
    }
  }

  KOKKOS_INLINE_FUNCTION void init (value_type update) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < value_count; ++k) {
      update[k] = Details::ArithTraits<mag_type>::zero ();
    }
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type update, const volatile value_type source) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < value_count; ++k) {
      update[k] = max (update[k], source[k]);
    }
  }
};


template<class normVector, class VectorType>
normVector MV_NormInf(const normVector &r, const VectorType & x, int n = -1)
{
  if (n < 0) {
    n = x.dimension_0 ();
  }

  Kokkos::parallel_reduce (n , MV_NormInf_Functor<VectorType> (x), r);
  return r;
}

/*------------------------------------------------------------------------------------------
 *-------------------------- Compute Weighted Dot-product (sum(x_i/w_i)^2)----------------------------------
 *------------------------------------------------------------------------------------------*/
template<class WeightVector, class XVector,int WeightsRanks>
struct MV_DotWeighted_Functor{};

template<class WeightVector, class XVector>
struct MV_DotWeighted_Functor<WeightVector,XVector,1>
{
  typedef typename XVector::execution_space        execution_space;
  typedef typename XVector::size_type            size_type;
  typedef typename XVector::non_const_value_type          xvalue_type;
  typedef typename WeightVector::non_const_value_type     wvalue_type;
  typedef Details::InnerProductSpaceTraits<xvalue_type> XIPT;
  typedef Details::InnerProductSpaceTraits<wvalue_type> WIPT;
  typedef typename XIPT::dot_type               value_type[];

  typename WeightVector::const_type m_w ;
  typename XVector::const_type m_x ;
  size_type value_count;

  MV_DotWeighted_Functor(WeightVector w, XVector x):m_w(w),m_x(x),value_count(x.dimension_1()) {}
  //--------------------------------------------------------------------------
  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i, value_type sum ) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for(size_type k=0;k<value_count;k++){
      sum[k] += XIPT::dot( m_x(i,k), m_x(i,k) ) / WIPT::dot( m_w(i), m_w(i) );
    }
  }

  KOKKOS_INLINE_FUNCTION void init( value_type update) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for(size_type k=0;k<value_count;k++)
      update[k] = 0;
  }
  KOKKOS_INLINE_FUNCTION void join( volatile value_type  update ,
                                    const volatile value_type  source ) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for(size_type k=0;k<value_count;k++){
      update[k] += source[k];
    }
  }
};

template<class WeightVector, class XVector>
struct MV_DotWeighted_Functor<WeightVector,XVector,2>
{
  typedef typename XVector::execution_space        execution_space;
  typedef typename XVector::size_type            size_type;
  typedef typename XVector::non_const_value_type          xvalue_type;
  typedef typename WeightVector::non_const_value_type     wvalue_type;
  typedef Details::InnerProductSpaceTraits<xvalue_type> XIPT;
  typedef Details::InnerProductSpaceTraits<wvalue_type> WIPT;
  typedef typename XIPT::dot_type               value_type[];

  typename WeightVector::const_type m_w ;
  typename XVector::const_type m_x ;
  size_type value_count;

  MV_DotWeighted_Functor(WeightVector w, XVector x):m_w(w),m_x(x),value_count(x.dimension_1()) {}
  //--------------------------------------------------------------------------
  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i, value_type sum ) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for(size_type k=0;k<value_count;k++){
      sum[k] += XIPT::dot( m_x(i,k), m_x(i,k) ) / WIPT::dot( m_w(i,k), m_w(i,k) );
    }
  }

  KOKKOS_INLINE_FUNCTION void init( value_type update) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for(size_type k=0;k<value_count;k++)
      update[k] = 0;
  }
  KOKKOS_INLINE_FUNCTION void join( volatile value_type  update ,
                                    const volatile value_type  source ) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for(size_type k=0;k<value_count;k++){
      update[k] += source[k];
    }
  }
};

template<class rVector, class WeightVector, class XVector>
rVector
MV_DotWeighted (const rVector &r,
                const WeightVector & w,
                const XVector & x,
                int n = -1)
{
  if (n < 0) {
    n = x.dimension_0 ();
  }

  typedef MV_DotWeighted_Functor<WeightVector, XVector, WeightVector::Rank> functor_type;
  Kokkos::parallel_reduce (n , functor_type (w, x), r);
  return r;
}

/*------------------------------------------------------------------------------------------
 *-------------------------- Multiply with scalar: y = a * x -------------------------------
 *------------------------------------------------------------------------------------------*/
template<class RVector, class aVector, class XVector>
struct V_MulScalarFunctor
{
  typedef typename XVector::execution_space        execution_space;
  typedef typename XVector::size_type            size_type;

  RVector m_r;
  typename XVector::const_type m_x ;
  typename aVector::const_type m_a ;
  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i) const
  {
    if (m_a[0] == Kokkos::Details::ArithTraits<typename aVector::non_const_value_type>::zero ()) {
      m_r(i) = Kokkos::Details::ArithTraits<typename XVector::non_const_value_type>::zero ();
    } else {
      m_r(i) = m_a[0]*m_x(i);
    }
  }
};

template<class aVector, class XVector>
struct V_MulScalarFunctorSelf
{
  typedef typename XVector::execution_space execution_space;
  typedef typename XVector::size_type       size_type;

  XVector m_x;
  typename aVector::const_type m_a;
  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i) const
  {
    if (m_a(0) == Kokkos::Details::ArithTraits<typename aVector::non_const_value_type>::zero ()) {
      m_x(i) = Kokkos::Details::ArithTraits<typename XVector::non_const_value_type>::zero ();
    } else {
      m_x(i) *= m_a(0);
    }
  }
};

template<class RVector, class DataType,class Layout,class Device, class MemoryManagement,class Specialisation, class XVector>
RVector
V_MulScalar (const RVector& r,
             const typename Kokkos::View<DataType,Layout,Device,MemoryManagement,Specialisation>& a,
             const XVector& x)
{
  typedef typename Kokkos::View<DataType,Layout,Device,MemoryManagement> aVector;
  if (r == x) {
    V_MulScalarFunctorSelf<aVector,XVector> op;
    op.m_x = x;
    op.m_a = a;
    Kokkos::parallel_for (x.dimension (0) , op);
    return r;
  }

  V_MulScalarFunctor<RVector,aVector,XVector> op;
  op.m_r = r;
  op.m_x = x;
  op.m_a = a;
  Kokkos::parallel_for (x.dimension (0), op);
  return r;
}

template<class RVector, class XVector>
struct V_MulScalarFunctor<RVector,typename XVector::non_const_value_type,XVector>
{
  typedef typename XVector::execution_space execution_space;
  typedef typename XVector::size_type       size_type;

  RVector m_r;
  typename XVector::const_type m_x;
  typename XVector::value_type m_a;
  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const
  {
    if (m_a == Kokkos::Details::ArithTraits<typename XVector::non_const_value_type>::zero ()) {
      m_r(i) = Kokkos::Details::ArithTraits<typename XVector::non_const_value_type>::zero ();
    } else {
      m_r(i) = m_a*m_x(i);
    }
  }
};

template<class XVector>
struct V_MulScalarFunctorSelf<typename XVector::non_const_value_type,XVector>
{
  typedef typename XVector::execution_space execution_space;
  typedef typename XVector::size_type       size_type;

  XVector m_x;
  typename XVector::non_const_value_type m_a;

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const
  {
    if (m_a == Kokkos::Details::ArithTraits<typename XVector::non_const_value_type>::zero ()) {
      m_x(i) = Kokkos::Details::ArithTraits<typename XVector::non_const_value_type>::zero ();
    } else {
      m_x(i) *= m_a;
    }
  }
};


template<class RVector, class XVector>
RVector V_MulScalar( const RVector & r, const typename XVector::non_const_value_type &a, const XVector & x)
{
  if (r == x) {
    V_MulScalarFunctorSelf<typename RVector::value_type,RVector> op;
    op.m_x = r;
    op.m_a = a;
    Kokkos::parallel_for (x.dimension (0), op);
    return r;
  }

  V_MulScalarFunctor<RVector,typename XVector::non_const_value_type,XVector> op;
  op.m_r = r;
  op.m_x = x;
  op.m_a = a;
  Kokkos::parallel_for (x.dimension (0), op);
  return r;
}

template<class RVector, class XVector, class YVector, int scalar_x, int scalar_y>
struct V_AddVectorFunctor
{
  typedef typename RVector::execution_space execution_space;
  typedef typename RVector::size_type       size_type;
  typedef typename XVector::non_const_value_type value_type;
  RVector   m_r;
  typename XVector::const_type  m_x;
  typename YVector::const_type   m_y;
  const value_type m_a;
  const value_type m_b;

  //--------------------------------------------------------------------------
  V_AddVectorFunctor(const RVector& r, const value_type& a,const XVector& x,const value_type& b,const YVector& y):
          m_r(r),m_x(x),m_y(y),m_a(a),m_b(b)
  { }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i ) const
  {
        if((scalar_x==1)&&(scalar_y==1))
            m_r(i) = m_x(i) + m_y(i);
        if((scalar_x==1)&&(scalar_y==-1))
            m_r(i) = m_x(i) - m_y(i);
        if((scalar_x==-1)&&(scalar_y==-1))
            m_r(i) = -m_x(i) - m_y(i);
        if((scalar_x==-1)&&(scalar_y==1))
            m_r(i) = -m_x(i) + m_y(i);
        if((scalar_x==2)&&(scalar_y==1))
            m_r(i) = m_a*m_x(i) + m_y(i);
        if((scalar_x==2)&&(scalar_y==-1))
            m_r(i) = m_a*m_x(i) - m_y(i);
        if((scalar_x==1)&&(scalar_y==2))
            m_r(i) = m_x(i) + m_b*m_y(i);
        if((scalar_x==-1)&&(scalar_y==2))
            m_r(i) = -m_x(i) + m_b*m_y(i);
        if((scalar_x==2)&&(scalar_y==2))
            m_r(i) = m_a*m_x(i) + m_b*m_y(i);
  }
};

template<class RVector, class XVector, int scalar_x>
struct V_AddVectorSelfFunctor
{
  typedef typename RVector::execution_space        execution_space;
  typedef typename RVector::size_type            size_type;
  typedef typename XVector::non_const_value_type      value_type;
  RVector   m_r ;
  typename XVector::const_type  m_x ;
  const value_type m_a;

  V_AddVectorSelfFunctor(const RVector& r, const value_type& a,const XVector& x):
    m_r(r),m_x(x),m_a(a)
  { }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i ) const
  {
  if((scalar_x==1))
      m_r(i) += m_x(i);
  if((scalar_x==-1))
      m_r(i) -= m_x(i);
  if((scalar_x==2))
      m_r(i) += m_a*m_x(i);
  }
};
template<class RVector, class XVector, class YVector, int doalpha, int dobeta>
RVector V_AddVector( const RVector & r,const typename XVector::non_const_value_type &av,const XVector & x,
                const typename XVector::non_const_value_type &bv, const YVector & y,int n=-1)
{
  if(n == -1) n = x.dimension_0();
  if(r.ptr_on_device()==x.ptr_on_device() && doalpha == 1) {
    V_AddVectorSelfFunctor<RVector,YVector,dobeta> f(r,bv,y);
    parallel_for(n,f);
  } else if(r.ptr_on_device()==y.ptr_on_device() && dobeta == 1) {
    V_AddVectorSelfFunctor<RVector,XVector,doalpha> f(r,av,x);
    parallel_for(n,f);
  } else {
    V_AddVectorFunctor<RVector,XVector,YVector,doalpha,dobeta> f(r,av,x,bv,y);
    parallel_for(n,f);
  }
  return r;
}

template<class RVector, class XVector, class YVector>
RVector V_AddVector( const RVector & r,const typename XVector::non_const_value_type &av,const XVector & x,
                const typename YVector::non_const_value_type &bv, const YVector & y, int n = -1,
                int a=2,int b=2)
{
        if(a==-1) {
          if(b==-1)
                  V_AddVector<RVector,XVector,YVector,-1,-1>(r,av,x,bv,y,n);
          else if(b==0)
                  V_AddVector<RVector,XVector,YVector,-1,0>(r,av,x,bv,y,n);
          else if(b==1)
              V_AddVector<RVector,XVector,YVector,-1,1>(r,av,x,bv,y,n);
          else
              V_AddVector<RVector,XVector,YVector,-1,2>(r,av,x,bv,y,n);
        } else if (a==0) {
          if(b==-1)
                  V_AddVector<RVector,XVector,YVector,0,-1>(r,av,x,bv,y,n);
          else if(b==0)
                  V_AddVector<RVector,XVector,YVector,0,0>(r,av,x,bv,y,n);
          else if(b==1)
              V_AddVector<RVector,XVector,YVector,0,1>(r,av,x,bv,y,n);
          else
              V_AddVector<RVector,XVector,YVector,0,2>(r,av,x,bv,y,n);
        } else if (a==1) {
          if(b==-1)
                  V_AddVector<RVector,XVector,YVector,1,-1>(r,av,x,bv,y,n);
          else if(b==0)
                  V_AddVector<RVector,XVector,YVector,1,0>(r,av,x,bv,y,n);
          else if(b==1)
              V_AddVector<RVector,XVector,YVector,1,1>(r,av,x,bv,y,n);
          else
              V_AddVector<RVector,XVector,YVector,1,2>(r,av,x,bv,y,n);
        } else if (a==2) {
          if(b==-1)
                  V_AddVector<RVector,XVector,YVector,2,-1>(r,av,x,bv,y,n);
          else if(b==0)
                  V_AddVector<RVector,XVector,YVector,2,0>(r,av,x,bv,y,n);
          else if(b==1)
              V_AddVector<RVector,XVector,YVector,2,1>(r,av,x,bv,y,n);
          else
              V_AddVector<RVector,XVector,YVector,2,2>(r,av,x,bv,y,n);
        }
        return r;
}

template<class RVector,class XVector,class YVector>
RVector V_Add( const RVector & r, const XVector & x, const YVector & y, int n=-1)
{
        return V_AddVector( r,1,x,1,y,n,1,1);
}

template<class RVector,class XVector,class YVector>
RVector V_Add( const RVector & r, const XVector & x, const typename XVector::non_const_value_type  & bv, const YVector & y,int n=-1 )
{
  int b = 2;
  //if(bv == 0) b = 0;
  //if(bv == 1) b = 1;
  //if(bv == -1) b = -1;
  return V_AddVector(r,bv,x,bv,y,n,1,b);
}

template<class RVector,class XVector,class YVector>
RVector V_Add( const RVector & r, const typename XVector::non_const_value_type  & av, const XVector & x, const typename XVector::non_const_value_type  & bv, const YVector & y,int n=-1 )
{
  int a = 2;
  int b = 2;
  //if(av == 0) a = 0;
  //if(av == 1) a = 1;
  //if(av == -1) a = -1;
  //if(bv == 0) b = 0;
  //if(bv == 1) b = 1;
  //if(bv == -1) b = -1;

  return V_AddVector(r,av,x,bv,y,n,a,b);
}

template<class XVector, class YVector>
struct V_DotFunctor
{
  typedef typename XVector::execution_space          execution_space;
  typedef typename XVector::size_type              size_type;
  typedef typename XVector::non_const_value_type xvalue_type;
  typedef Details::InnerProductSpaceTraits<xvalue_type>  IPT;
  typedef typename IPT::dot_type                  value_type;
  XVector  m_x ;
  YVector  m_y ;

  //--------------------------------------------------------------------------
  V_DotFunctor(const XVector& x,const YVector& y):
    m_x(x),m_y(y)
  { }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type &i, value_type &sum ) const
  {
    sum += IPT::dot( m_x(i), m_y(i) );  // m_x(i) * m_y(i)
  }

  KOKKOS_INLINE_FUNCTION
  void init (volatile value_type& update) const
  {
    update = Details::ArithTraits<value_type>::zero ();
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type &update ,
             const volatile value_type &source ) const
  {
    update += source ;
  }
};

//! Return the dot product of the vectors (1-D arrays) x and y.
template<class XVector, class YVector>
typename Details::InnerProductSpaceTraits<typename XVector::non_const_value_type>::dot_type
V_Dot (const XVector& x,
       const YVector& y)
{
  const typename XVector::size_type numRows = x.dimension_0 ();
  return V_Dot<XVector, YVector> (x, y, numRows);
}

/// \brief Return the dot product of the vectors (1-D arrays) x and y.
///   Only use the first numRows entries of x and y.
template<class XVector, class YVector>
typename Details::InnerProductSpaceTraits<typename XVector::non_const_value_type>::dot_type
V_Dot (const XVector& x,
       const YVector& y,
       const typename XVector::size_type numRows)
{
  typedef V_DotFunctor<XVector, YVector> functor_type;
  functor_type f (x, y);
  typename functor_type::value_type ret_val = typename functor_type::value_type();
  parallel_reduce (numRows, f, ret_val);
  return ret_val;
}

template<class WeightVector, class XVector>
struct V_DotWeighted_Functor
{
  typedef typename XVector::execution_space execution_space;
  typedef typename XVector::size_type size_type;
  typedef typename XVector::non_const_value_type xvalue_type;
  typedef typename WeightVector::non_const_value_type wvalue_type;
  typedef Details::InnerProductSpaceTraits<xvalue_type> XIPT;
  typedef Details::InnerProductSpaceTraits<wvalue_type> WIPT;
  typedef typename XIPT::dot_type value_type;

  typename WeightVector::const_type m_w;
  typename XVector::const_type m_x;

  V_DotWeighted_Functor (WeightVector w, XVector x) :
    m_w (w), m_x (x)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i, value_type& sum) const {
    sum += XIPT::dot (m_x(i), m_x(i)) / WIPT::dot (m_w(i), m_w(i));
  }

  KOKKOS_INLINE_FUNCTION void init (value_type& update) const {
    update = Details::ArithTraits<value_type>::zero ();
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& update,
        const volatile value_type& source) const
  {
    update += source;
  }
};

template<class WeightVector, class XVector>
typename Details::InnerProductSpaceTraits<typename XVector::non_const_value_type>::dot_type
V_DotWeighted (const WeightVector& w,
               const XVector& x)
{
  const typename XVector::size_type numRows = x.dimension_0 ();
  return V_DotWeighted<WeightVector, XVector> (w, x, numRows);
}

template<class WeightVector, class XVector>
typename Details::InnerProductSpaceTraits<typename XVector::non_const_value_type>::dot_type
V_DotWeighted (const WeightVector& w,
               const XVector& x,
               const typename XVector::size_type numRows)
{
  typedef Details::InnerProductSpaceTraits<typename XVector::non_const_value_type> IPT;
  typedef typename IPT::dot_type value_type;
  value_type ret_val;

  typedef V_DotWeighted_Functor<WeightVector, XVector> functor_type;
  Kokkos::parallel_reduce (numRows, functor_type (w, x), ret_val);
  return ret_val;
}

/*------------------------------------------------------------------------------------------
 *-------------------------- Compute Sum -------------------------------------------------
 *------------------------------------------------------------------------------------------*/
template<class XVector>
struct V_Sum_Functor
{
  typedef typename XVector::execution_space        execution_space;
  typedef typename XVector::size_type            size_type;
  typedef typename XVector::non_const_value_type          xvalue_type;
  typedef xvalue_type                           value_type;

  typename XVector::const_type m_x ;

  V_Sum_Functor(XVector x):m_x(x) {}
  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i, value_type& sum ) const
  {
      sum += m_x(i);
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type& update) const
  {
    update = Details::ArithTraits<value_type>::zero ();
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type&  update ,
                                    const volatile value_type&  source ) const
  {
      update += source;
  }
};


template<class VectorType>
typename VectorType::non_const_value_type
V_Sum (const VectorType& x, int n = -1)
{
  if (n < 0) {
    n = x.dimension_0 ();
  }

  typedef typename VectorType::non_const_value_type value_type;
  value_type ret_val;
  Kokkos::parallel_reduce (n, V_Sum_Functor<VectorType> (x), ret_val);
  return ret_val;
}

/*------------------------------------------------------------------------------------------
 *-------------------------- Compute Norm1--------------------------------------------------
 *------------------------------------------------------------------------------------------*/
template<class XVector>
struct V_Norm1_Functor
{
  typedef typename XVector::execution_space        execution_space;
  typedef typename XVector::size_type            size_type;
  typedef typename XVector::non_const_value_type          xvalue_type;
  typedef Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef typename IPT::dot_type               value_type;

  typename XVector::const_type m_x ;

  V_Norm1_Functor(XVector x):m_x(x) {}
  //--------------------------------------------------------------------------
  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i, value_type& sum ) const
  {
    sum += Kokkos::Details::ArithTraits<typename XVector::non_const_value_type>::abs(m_x(i));
  }
  KOKKOS_INLINE_FUNCTION void init (value_type& update) const
  {
    update = Details::ArithTraits<value_type>::zero ();
  }
  KOKKOS_INLINE_FUNCTION void join( volatile value_type&  update ,
                                    const volatile value_type&  source ) const
  {
    update += source;
  }
};

template<class VectorType>
typename Details::InnerProductSpaceTraits<typename VectorType::non_const_value_type>::dot_type
V_Norm1( const VectorType & x, int n = -1)
{
  if (n < 0) {
    n = x.dimension_0 ();
  }

  typedef Details::InnerProductSpaceTraits<typename VectorType::non_const_value_type> IPT;
  typedef typename IPT::dot_type value_type;
  value_type ret_val;
  Kokkos::parallel_reduce (n, V_Norm1_Functor<VectorType> (x), ret_val);
  return ret_val;
}
/*------------------------------------------------------------------------------------------
 *-------------------------- Compute NormInf--------------------------------------------------
 *------------------------------------------------------------------------------------------*/
template<class XVector>
struct V_NormInf_Functor
{
  typedef typename XVector::execution_space             execution_space;
  typedef typename XVector::size_type                   size_type;
  typedef typename XVector::non_const_value_type        xvalue_type;
  typedef Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef typename IPT::mag_type                        mag_type;
  typedef mag_type                                      value_type;

  typename XVector::const_type m_x;

  V_NormInf_Functor (const XVector& x) : m_x (x) {}

  KOKKOS_INLINE_FUNCTION mag_type
  max (const mag_type& x, const mag_type& y) const {
    return (x > y) ? x : y;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i, value_type& sum) const
  {
    sum = max (sum, IPT::norm (m_x(i)));
  }

  KOKKOS_INLINE_FUNCTION void init (value_type& update) const
  {
    update = Details::ArithTraits<mag_type>::zero ();
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& update, const volatile value_type& source) const
  {
    update = max (update, source);
  }
};

template<class VectorType>
typename Details::InnerProductSpaceTraits<typename VectorType::non_const_value_type>::mag_type
V_NormInf (const VectorType& x, int n = -1)
{
  if (n < 0) {
    n = x.dimension_0 ();
  }

  typedef Details::InnerProductSpaceTraits<typename VectorType::non_const_value_type> IPT;
  typedef typename IPT::mag_type value_type;
  value_type ret_val;
  Kokkos::parallel_reduce (n, V_NormInf_Functor<VectorType> (x), ret_val);
  return ret_val;
}

/*------------------------------------------------------------------------------------------
 *-------------------------- Reciprocal element wise: y[i] = 1/x[i] ------------------------
 *------------------------------------------------------------------------------------------*/
template<class RVector, class XVector>
struct V_ReciprocalFunctor
{
  typedef typename XVector::execution_space        execution_space;
  typedef typename XVector::size_type            size_type;

  RVector m_r;
  typename XVector::const_type m_x ;

  V_ReciprocalFunctor(RVector r, XVector x):m_r(r),m_x(x) {}
  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i) const
  {
    m_r(i) = Kokkos::Details::ArithTraits<typename XVector::non_const_value_type>::one() / m_x(i);
  }
};

template<class XVector>
struct V_ReciprocalSelfFunctor
{
  typedef typename XVector::execution_space        execution_space;
  typedef typename XVector::size_type            size_type;

  XVector m_x ;

  V_ReciprocalSelfFunctor(XVector x):m_x(x) {}
  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i) const
  {
     m_x(i) = Kokkos::Details::ArithTraits<typename XVector::non_const_value_type>::one() / m_x(i);
  }
};

template<class RVector, class XVector>
RVector V_Reciprocal( const RVector & r, const XVector & x)
{
  // TODO: Add error check (didn't link for some reason?)
  /*if(r.dimension_0() != x.dimension_0())
    Kokkos::Impl::throw_runtime_exception("Kokkos::MV_Reciprocal -- dimension(0) of r and x don't match");
  */


  if(r==x) {
    V_ReciprocalSelfFunctor<XVector> op(x) ;
    Kokkos::parallel_for( x.dimension_0() , op );
    return r;
  }

  V_ReciprocalFunctor<RVector,XVector> op(r,x) ;
  Kokkos::parallel_for( x.dimension_0() , op );
  return r;
}

/*------------------------------------------------------------------------------------------
 *------------------- Reciprocal element wise with threshold: x[i] = 1/x[i] ----------------
 *------------------------------------------------------------------------------------------*/
template<class XVector>
struct V_ReciprocalThresholdSelfFunctor
{
  typedef typename XVector::execution_space           execution_space;
  typedef typename XVector::size_type               size_type;
  typedef typename XVector::non_const_value_type   value_type;
  typedef Kokkos::Details::ArithTraits<value_type>        KAT;
  typedef typename KAT::mag_type                     mag_type;

  const XVector    m_x;
  const value_type m_min_val;
  const value_type m_min_val_mag;

  V_ReciprocalThresholdSelfFunctor(const XVector& x, const value_type& min_val) :
    m_x(x), m_min_val(min_val), m_min_val_mag(KAT::abs(min_val)) {}
  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i) const
  {
    if (KAT::abs(m_x(i)) < m_min_val_mag)
      m_x(i) = m_min_val;
    else
      m_x(i) = KAT::one() / m_x(i);
  }
};

template<class XVector>
XVector V_ReciprocalThreshold( const XVector & x, const typename XVector::non_const_value_type& min_val )
{
  V_ReciprocalThresholdSelfFunctor<XVector> op(x,min_val) ;
  Kokkos::parallel_for( x.dimension_0() , op );
  return x;
}

//! Functor for element-wise absolute value: r(i) = abs(x(i))
template<class RVector, class XVector>
struct V_AbsFunctor
{
  typedef typename XVector::execution_space execution_space;
  typedef typename XVector::size_type       size_type;

  RVector m_r;
  typename XVector::const_type m_x;

  V_AbsFunctor (const RVector& r, const XVector& x) : m_r (r), m_x (x) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const
  {
    m_r(i) = Kokkos::Details::ArithTraits<typename XVector::non_const_value_type>::abs (m_x(i));
  }
};

//! Functor for element-wise absolute value in place: x(i) = abs(x(i))
template<class XVector>
struct V_AbsSelfFunctor
{
  typedef typename XVector::execution_space execution_space;
  typedef typename XVector::size_type       size_type;

  XVector m_x ;
  V_AbsSelfFunctor (const XVector& x) : m_x (x) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const
  {
     m_x(i) = Kokkos::Details::ArithTraits<typename XVector::non_const_value_type>::abs (m_x(i));
  }
};

/// \brief Compute element-wise absolute value: r(i) = abs(x(i)).
///
/// We allow r to alias x.  In that case, we compute in place.
template<class RVector, class XVector>
RVector V_Abs( const RVector & r, const XVector & x)
{
  if (r == x) {
    V_AbsSelfFunctor<XVector> op (x);
    Kokkos::parallel_for (x.dimension_0 () , op);
  } else {
    V_AbsFunctor<RVector, XVector> op (r, x);
    Kokkos::parallel_for (x.dimension_0 (), op);
  }
  return r;
}

/// \brief Functor for element-wise multiply of vectors.
///
/// This functor implements Tpetra::MultiVector::elementWiseMultiply,
/// for the case where all MultiVector instances in question have only
/// a single column.  Thus, the functor computes C(i) = c*C(i) +
/// ab*A(i)*B(i).
template<class CVector, class AVector, class BVector>
struct V_ElementWiseMultiplyFunctor
{
  typedef typename CVector::execution_space        execution_space;
  typedef typename CVector::size_type            size_type;

  typename CVector::const_value_type m_c;
  CVector m_C;
  typename AVector::const_value_type m_ab;
  typename AVector::const_type m_A ;
  typename BVector::const_type m_B ;

  V_ElementWiseMultiplyFunctor (typename CVector::const_value_type c,
                                CVector C,
                                typename AVector::const_value_type ab,
                                typename AVector::const_type A,
                                typename BVector::const_type B) :
    m_c (c), m_C (C), m_ab (ab), m_A (A), m_B (B)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator () (const size_type i) const
  {
    const typename CVector::non_const_value_type zero_C =
      Kokkos::Details::ArithTraits<typename CVector::non_const_value_type>::zero ();
    const typename AVector::non_const_value_type zero_A =
      Kokkos::Details::ArithTraits<typename AVector::non_const_value_type>::zero ();

    if (m_c == zero_C) {
      if (m_ab == zero_A) {
        // Overwrite m_C with zeros, per BLAS update rules.
        m_C(i) = zero_C;
      }
      else { // m_ab != 0, but m_c == 0
        // BLAS update rules say that if m_c == 0, we must overwrite
        // m_C.  This matters only if m_C has entries that are Inf or
        // NaN.
        m_C(i) = m_ab * m_A(i) * m_B(i);
      }
    }
    else { // m_c != 0
      if (m_ab == zero_A) {
        m_C(i) = m_c * m_C(i);
      }
      else { // m_ab != 0, and m_c != 0
        m_C(i) = m_c * m_C(i) + m_ab * m_A(i) * m_B(i);
      }
    }
  }
};


template<class CVector, class AVector, class BVector>
CVector
V_ElementWiseMultiply (const typename CVector::const_value_type& c,
                       const CVector& C,
                       const typename AVector::const_value_type& ab,
                       const AVector& A,
                       const BVector& B)
{
  const typename CVector::non_const_value_type zero_C =
    Kokkos::Details::ArithTraits<typename CVector::non_const_value_type>::zero ();
  const typename AVector::non_const_value_type zero_A =
    Kokkos::Details::ArithTraits<typename AVector::non_const_value_type>::zero ();

  if (ab == zero_A && c == zero_C) {
    // Overwrite m_C with zeros, per BLAS update rules.
    Kokkos::Impl::ViewFill<CVector> (C, zero_C);
  }
  else {
    V_ElementWiseMultiplyFunctor<CVector, AVector, BVector> op (c, C, ab, A, B);
    Kokkos::parallel_for (C.dimension_0 (), op);
  }
  return C;
}

} // namespace Kokkos
#endif /* KOKKOS_MULTIVECTOR_H_ */
