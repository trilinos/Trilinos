#ifndef KOKKOSARRAY_MULTIVECTOR_H_
#define KOKKOSARRAY_MULTIVECTOR_H_

#include <KokkosArray_View.hpp>
#include <KokkosArray_Host.hpp>
#include <KokkosArray_Cuda.hpp>
#include <KokkosArray_Macros.hpp>
#include <KokkosArray_ParallelReduce.hpp>
#include <ctime>

namespace KokkosArray {


template<typename Scalar, class device>
struct MultiVectorDynamic{
#ifdef KOKKOS_USE_CUSPARSE
  typedef typename KokkosArray::LayoutLeft layout;
#else
#ifdef KOKKOS_USE_MKL
  typedef typename KokkosArray::LayoutRight layout;
#else
  typedef typename device::array_layout layout;
#endif
#endif
  typedef typename KokkosArray::View<Scalar**  , layout, device>  type ;
  typedef typename KokkosArray::View<const Scalar**  , layout, device>  const_type ;
  typedef typename KokkosArray::View<const Scalar**  , layout, device, KokkosArray::MemoryRandomRead>  random_read_type ;
  MultiVectorDynamic() {}
  ~MultiVectorDynamic() {}
};

template<typename Scalar, class device, int n>
struct MultiVectorStatic{
  typedef Scalar scalar;
  typedef typename device::array_layout layout;
  typedef typename KokkosArray::View<Scalar*[n]  , layout, device>  type ;
  typedef typename KokkosArray::View<const Scalar*[n]  , layout, device>  const_type ;
  typedef typename KokkosArray::View<const Scalar*[n]  , layout, device, KokkosArray::MemoryRandomRead>  random_read_type ;
  MultiVectorStatic() {}
  ~MultiVectorStatic() {}
};



/*------------------------------------------------------------------------------------------
 *-------------------------- Multiply with scalar: y = a * x -------------------------------
 *------------------------------------------------------------------------------------------*/
template<class RVector, class aVector, class XVector>
struct MV_MulScalarFunctor
{
  typedef typename XVector::device_type        device_type;
  typedef typename XVector::size_type            size_type;

  RVector m_r;
  typename RVector::const_type m_x ;
  typename aVector::const_type m_a ;
  size_type n;
  MV_MulScalarFunctor() {n=1;}
  //--------------------------------------------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const size_type i) const
  {
    #pragma ivdep
	for(size_type k=0;k<n;k++)
	   m_r(i,k) = m_a[k]*m_x(i,k);
  }
};

template<class aVector, class XVector>
struct MV_MulScalarFunctorSelf
{
  typedef typename XVector::device_type        device_type;
  typedef typename XVector::size_type            size_type;

  XVector m_x;
  typename aVector::const_type   m_a ;
  size_type n;
  //--------------------------------------------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const size_type i) const
  {
    #pragma ivdep
	for(size_type k=0;k<n;k++)
	   m_x(i,k) *= m_a[k];
  }
};

template<class RVector, class DataType,class Layout,class Device, class MemoryManagement,class Specialisation, class XVector>
RVector MV_MulScalar( const RVector & r, const typename KokkosArray::View<DataType,Layout,Device,MemoryManagement,Specialisation> & a, const XVector & x)
{
  typedef	typename KokkosArray::View<DataType,Layout,Device,MemoryManagement> aVector;
  if(r==x) {
    MV_MulScalarFunctorSelf<aVector,XVector> op ;
	op.m_x = x ;
	op.m_a = a ;
	op.n = x.dimension(1);
	KokkosArray::parallel_for( x.dimension(0) , op );
	return r;
  }

  MV_MulScalarFunctor<RVector,aVector,XVector> op ;
  op.m_r = r ;
  op.m_x = x ;
  op.m_a = a ;
  op.n = x.dimension(1);
  KokkosArray::parallel_for( x.dimension(0) , op );
  return r;
}

template<class RVector, class XVector>
struct MV_MulScalarFunctor<RVector,typename XVector::scalar_type,XVector>
{
  typedef typename XVector::device_type        device_type;
  typedef typename XVector::size_type            size_type;

  RVector m_r;
  typename RVector::const_type m_x ;
  typename XVector::scalar_type m_a ;
  size_type n;
  MV_MulScalarFunctor() {n=1;}
  //--------------------------------------------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const size_type i) const
  {
    #pragma ivdep
	for(size_type k=0;k<n;k++)
	   m_r(i,k) = m_a*m_x(i,k);
  }
};

template<class XVector>
struct MV_MulScalarFunctorSelf<typename XVector::scalar_type,XVector>
{
  typedef typename XVector::device_type        device_type;
  typedef typename XVector::size_type            size_type;

  XVector m_x;
  typename XVector::scalar_type   m_a ;
  size_type n;
  //--------------------------------------------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const size_type i) const
  {
    #pragma ivdep
	for(size_type k=0;k<n;k++)
	   m_x(i,k) *= m_a;
  }
};

template<class RVector, class XVector>
RVector MV_MulScalar( const RVector & r, const typename XVector::scalar_type &a, const XVector & x)
{
  if(r==x) {
    MV_MulScalarFunctorSelf<typename XVector::scalar_type,XVector> op ;
	op.m_x = x ;
	op.m_a = a ;
	op.n = x.dimension(1);
	KokkosArray::parallel_for( x.dimension(0) , op );
	return r;
  }

  MV_MulScalarFunctor<RVector,typename XVector::scalar_type,XVector> op ;
  op.m_r = r ;
  op.m_x = x ;
  op.m_a = a ;
  op.n = x.dimension(1);
  KokkosArray::parallel_for( x.dimension(0) , op );
  return r;
}
/*------------------------------------------------------------------------------------------
 *-------------------------- Vector Add: r = a*x + b*y -------------------------------------
 *------------------------------------------------------------------------------------------*/

//Unroll for n<=16
template<class RVector,class aVector, class XVector, class bVector, class YVector, int scalar_x, int scalar_y,int UNROLL>
struct MV_AddUnrollFunctor
{
  typedef typename RVector::device_type        device_type;
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

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const size_type i ) const
  {
	if((scalar_x==1)&&(scalar_y==1)){
	#pragma unroll
    for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = m_x(i,k) + m_y(i,k);
	}
	if((scalar_x==1)&&(scalar_y==-1)){
	  #pragma unroll
	  for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = m_x(i,k) - m_y(i,k);
	}
	if((scalar_x==-1)&&(scalar_y==-1)){
#pragma unroll
for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = -m_x(i,k) - m_y(i,k);
	}
	if((scalar_x==-1)&&(scalar_y==1)){
#pragma unroll
for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = -m_x(i,k) + m_y(i,k);
	}
	if((scalar_x==2)&&(scalar_y==1)){
#pragma unroll
for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = m_a(k)*m_x(i,k) + m_y(i,k);
	}
	if((scalar_x==2)&&(scalar_y==-1)){
#pragma unroll
for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = m_a(k)*m_x(i,k) - m_y(i,k);
	}
	if((scalar_x==1)&&(scalar_y==2)){
#pragma unroll
for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = m_x(i,k) + m_b(k)*m_y(i,k);
	}
	if((scalar_x==-1)&&(scalar_y==2)){
#pragma unroll
for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = -m_x(i,k) + m_b(k)*m_y(i,k);
	}
	if((scalar_x==2)&&(scalar_y==2)){
#pragma unroll
for(size_type k=0;k<UNROLL;k++)
      m_r(i,k) = m_a(k)*m_x(i,k) + m_b(k)*m_y(i,k);
	}
  }
};

template<class RVector,class aVector, class XVector, class bVector, class YVector, int scalar_x, int scalar_y>
struct MV_AddVectorFunctor
{
  typedef typename RVector::device_type        device_type;
  typedef typename RVector::size_type            size_type;

  RVector   m_r ;
  XVector  m_x ;
  YVector   m_y ;
  aVector m_a;
  bVector m_b;
  size_type n;

  MV_AddVectorFunctor() {n=1;}
  //--------------------------------------------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const size_type i ) const
  {
	if((scalar_x==1)&&(scalar_y==1))
      #pragma ivdep
      #pragma vector always
      for(size_type k=0;k<n;k++)
	    m_r(i,k) = m_x(i,k) + m_y(i,k);
	if((scalar_x==1)&&(scalar_y==-1))
      #pragma ivdep
	  #pragma vector always
      for(size_type k=0;k<n;k++)
	    m_r(i,k) = m_x(i,k) - m_y(i,k);
	if((scalar_x==-1)&&(scalar_y==-1))
      #pragma ivdep
      #pragma vector always
      for(size_type k=0;k<n;k++)
	    m_r(i,k) = -m_x(i,k) - m_y(i,k);
	if((scalar_x==-1)&&(scalar_y==1))
      #pragma ivdep
      #pragma vector always
      for(size_type k=0;k<n;k++)
	    m_r(i,k) = -m_x(i,k) + m_y(i,k);
	if((scalar_x==2)&&(scalar_y==1))
      #pragma ivdep
      #pragma vector always
      for(size_type k=0;k<n;k++)
	    m_r(i,k) = m_a(k)*m_x(i,k) + m_y(i,k);
	if((scalar_x==2)&&(scalar_y==-1))
      #pragma ivdep
      #pragma vector always
      for(size_type k=0;k<n;k++)
	    m_r(i,k) = m_a(k)*m_x(i,k) - m_y(i,k);
	if((scalar_x==1)&&(scalar_y==2))
      #pragma ivdep
      #pragma vector always
      for(size_type k=0;k<n;k++)
	    m_r(i,k) = m_x(i,k) + m_b(k)*m_y(i,k);
	if((scalar_x==-1)&&(scalar_y==2))
      #pragma ivdep
      #pragma vector always
      for(size_type k=0;k<n;k++)
	    m_r(i,k) = -m_x(i,k) + m_b(k)*m_y(i,k);
	if((scalar_x==2)&&(scalar_y==2))
      #pragma ivdep
      #pragma vector always
      for(size_type k=0;k<n;k++)
	    m_r(i,k) = m_a(k)*m_x(i,k) + m_b(k)*m_y(i,k);

  }
};

template<class RVector,class aVector, class XVector, class bVector, class YVector,int UNROLL>
RVector MV_AddUnroll( const RVector & r,const aVector &av,const XVector & x,
		const bVector &bv, const YVector & y,
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
     KokkosArray::parallel_for( x.dimension(0) , op );
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
     KokkosArray::parallel_for( x.dimension(0) , op );
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
     KokkosArray::parallel_for( x.dimension(0) , op );
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
     KokkosArray::parallel_for( x.dimension(0) , op );
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
     KokkosArray::parallel_for( x.dimension(0) , op );
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
     KokkosArray::parallel_for( x.dimension(0) , op );
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
     KokkosArray::parallel_for( x.dimension(0) , op );
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
     KokkosArray::parallel_for( x.dimension(0) , op );
     return r;
   }
   MV_AddUnrollFunctor<RVector,aVector,XVector,bVector,YVector,2,2,UNROLL> op ;
   op.m_r = r ;
   op.m_x = x ;
   op.m_y = y ;
   op.m_a = av ;
   op.m_b = bv ;
   op.n = x.dimension(1);
   KokkosArray::parallel_for( x.dimension(0) , op );

   return r;
}

template<class RVector,class aVector, class XVector, class bVector, class YVector>
RVector MV_AddUnroll( const RVector & r,const aVector &av,const XVector & x,
		const bVector &bv, const YVector & y,
		int a=2,int b=2)
{
	switch (x.dimension(1)){
      case 1: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 1>( r,av,x,bv,y,a,b);
	          break;
      case 2: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 2>( r,av,x,bv,y,a,b);
	          break;
      case 3: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 3>( r,av,x,bv,y,a,b);
	          break;
      case 4: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 4>( r,av,x,bv,y,a,b);
	          break;
      case 5: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 5>( r,av,x,bv,y,a,b);
	          break;
      case 6: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 6>( r,av,x,bv,y,a,b);
	          break;
      case 7: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 7>( r,av,x,bv,y,a,b);
	          break;
      case 8: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 8>( r,av,x,bv,y,a,b);
	          break;
      case 9: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 9>( r,av,x,bv,y,a,b);
	          break;
      case 10: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 10>( r,av,x,bv,y,a,b);
	          break;
      case 11: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 11>( r,av,x,bv,y,a,b);
	          break;
      case 12: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 12>( r,av,x,bv,y,a,b);
	          break;
      case 13: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 13>( r,av,x,bv,y,a,b);
	          break;
      case 14: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 14>( r,av,x,bv,y,a,b);
	          break;
      case 15: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 15>( r,av,x,bv,y,a,b);
	          break;
      case 16: MV_AddUnroll<RVector, aVector, XVector, bVector, YVector, 16>( r,av,x,bv,y,a,b);
	          break;
	}
	return r;
}


template<class RVector,class aVector, class XVector, class bVector, class YVector>
RVector MV_AddVector( const RVector & r,const aVector &av,const XVector & x,
		const bVector &bv, const YVector & y,
		int a=2,int b=2)
{
   if(a==1&&b==1) {
     MV_AddVectorFunctor<RVector,aVector,XVector,bVector,YVector,1,1> op ;
     op.m_r = r ;
     op.m_x = x ;
     op.m_y = y ;
     op.m_a = av ;
     op.m_b = bv ;
     op.n = x.dimension(1);
     KokkosArray::parallel_for( x.dimension(0) , op );
     return r;
   }
   if(a==1&&b==-1) {
     MV_AddVectorFunctor<RVector,aVector,XVector,bVector,YVector,1,-1> op ;
     op.m_r = r ;
     op.m_x = x ;
     op.m_y = y ;
     op.m_a = av ;
     op.m_b = bv ;
     op.n = x.dimension(1);
     KokkosArray::parallel_for( x.dimension(0) , op );
     return r;
   }
   if(a==-1&&b==1) {
     MV_AddVectorFunctor<RVector,aVector,XVector,bVector,YVector,-1,1> op ;
     op.m_r = r ;
     op.m_x = x ;
     op.m_y = y ;
     op.m_a = av ;
     op.m_b = bv ;
     op.n = x.dimension(1);
     KokkosArray::parallel_for( x.dimension(0) , op );
     return r;
   }
   if(a==-1&&b==-1) {
     MV_AddVectorFunctor<RVector,aVector,XVector,bVector,YVector,-1,-1> op ;
     op.m_r = r ;
     op.m_x = x ;
     op.m_y = y ;
     op.m_a = av ;
     op.m_b = bv ;
     op.n = x.dimension(1);
     KokkosArray::parallel_for( x.dimension(0) , op );
     return r;
   }
   if(a*a!=1&&b==1) {
     MV_AddVectorFunctor<RVector,aVector,XVector,bVector,YVector,2,1> op ;
     op.m_r = r ;
     op.m_x = x ;
     op.m_y = y ;
     op.m_a = av ;
     op.m_b = bv ;
     op.n = x.dimension(1);
     KokkosArray::parallel_for( x.dimension(0) , op );
     return r;
   }
   if(a*a!=1&&b==-1) {
     MV_AddVectorFunctor<RVector,aVector,XVector,bVector,YVector,2,-1> op ;
     op.m_r = r ;
     op.m_x = x ;
     op.m_y = y ;
     op.m_a = av ;
     op.m_b = bv ;
     op.n = x.dimension(1);
     KokkosArray::parallel_for( x.dimension(0) , op );
     return r;
   }
   if(a==1&&b*b!=1) {
     MV_AddVectorFunctor<RVector,aVector,XVector,bVector,YVector,1,2> op ;
     op.m_r = r ;
     op.m_x = x ;
     op.m_y = y ;
     op.m_a = av ;
     op.m_b = bv ;
     op.n = x.dimension(1);
     KokkosArray::parallel_for( x.dimension(0) , op );
     return r;
   }
   if(a==-1&&b*b!=1) {
     MV_AddVectorFunctor<RVector,aVector,XVector,bVector,YVector,-1,2> op ;
     op.m_r = r ;
     op.m_x = x ;
     op.m_y = y ;
     op.m_a = av ;
     op.m_b = bv ;
     op.n = x.dimension(1);
     KokkosArray::parallel_for( x.dimension(0) , op );
     return r;
   }
   MV_AddVectorFunctor<RVector,aVector,XVector,bVector,YVector,2,2> op ;
   op.m_r = r ;
   op.m_x = x ;
   op.m_y = y ;
   op.m_a = av ;
   op.m_b = bv ;
   op.n = x.dimension(1);
   KokkosArray::parallel_for( x.dimension(0) , op );

   return r;
}

template<class RVector,class aVector, class XVector, class bVector, class YVector>
RVector MV_Add( const RVector & r,const aVector &av,const XVector & x,
		const bVector &bv, const YVector & y,
		int a=2,int b=2)
{
	if(x.dimension(1)>16)
		return MV_AddVector( r,av,x,bv,y,a,b);
	return MV_AddUnroll( r,av,x,bv,y,a,b);
}

template<class RVector,class XVector,class YVector>
RVector MV_Add( const RVector & r, const XVector & x, const YVector & y)
{
	KokkosArray::View<typename XVector::value_type*  , KokkosArray::LayoutRight, typename XVector::device_type> dummy;
	if(x.dimension(1)>16)
		return MV_AddVector( r,dummy,x,dummy,y,1,1);
	return MV_AddUnroll( r,dummy,x,dummy,y,1,1);
}

template<class RVector,class XVector,class bVector, class YVector>
RVector MV_Add( const RVector & r, const XVector & x, const bVector & bv, const YVector & y )
{
  if(x.dimension(1)>16)
	return MV_AddVector(r,bv,x,bv,y,1,2);
  return MV_AddUnroll(r,bv,x,bv,y,1,2);
}


template<class XVector,class YVector>
struct MV_DotProduct_Right_FunctorVector
{
  typedef typename XVector::device_type        device_type;
  typedef typename XVector::size_type            size_type;
  typedef typename XVector::value_type        value_type[];
  size_type value_count;


  typedef typename XVector::const_type        x_const_type;
  typedef typename YVector::const_type 	      y_const_type;
  x_const_type  m_x ;
  y_const_type  m_y ;

  //--------------------------------------------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const size_type i, value_type sum ) const
  {
	const int numVecs=value_count;

	  //printf("A\n");
    #pragma ivdep
    #pragma vector always
	for(int k=0;k<numVecs;k++)
      sum[k]+=m_x(i,k)*m_y(i,k);
  }
  static KOKKOSARRAY_INLINE_FUNCTION void init( value_type update, const size_type numVecs)
  {
    #pragma ivdep
    #pragma vector always
	for(size_type k=0;k<numVecs;k++)
	  update[k] = 0;
  }
  static KOKKOSARRAY_INLINE_FUNCTION void join( volatile value_type  update ,
                    const volatile value_type  source,const size_type numVecs )
  {
    #pragma ivdep
    #pragma vector always
	for(size_type k=0;k<numVecs;k++){
	  update[k] += source[k];
	}
  }
};


template<class XVector,class YVector,int UNROLL>
struct MV_DotProduct_Right_FunctorUnroll
{
  typedef typename XVector::device_type        device_type;
  typedef typename XVector::size_type            size_type;
  typedef typename XVector::value_type        value_type[];
  size_type value_count;

  typedef typename XVector::const_type        x_const_type;
  typedef typename YVector::const_type 	      y_const_type;

  x_const_type  m_x ;
  y_const_type  m_y ;

  //--------------------------------------------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const size_type i, value_type sum ) const
  {
    #pragma unroll
    for(size_type k=0;k<UNROLL;k++)
      sum[k]+=m_x(i,k)*m_y(i,k);
  }
  static KOKKOSARRAY_INLINE_FUNCTION void init( volatile value_type update, const size_type numVecs)
  {
    #pragma unroll
	for(size_type k=0;k<UNROLL;k++)
	  update[k] = 0;
  }
  static KOKKOSARRAY_INLINE_FUNCTION void join( volatile value_type update ,
                    const volatile value_type source, const size_type numVecs )
  {
    #pragma unroll
	for(size_type k=0;k<UNROLL;k++)
	 update[k] += source[k] ;
  }
};

template<class rVector, class XVector, class YVector>
rVector MV_Dot(const rVector &r, const XVector & x, const YVector & y)
{
    typedef typename XVector::size_type            size_type;
	const size_type numVecs = x.dimension(1);

    if(numVecs>16){

        MV_DotProduct_Right_FunctorVector<XVector,YVector> op;
        op.m_x = x;
        op.m_y = y;
        op.value_count = numVecs;

        KokkosArray::parallel_reduce( x.dimension(0) , op, r );
        return r;
     }
     else
     switch(numVecs) {
       case 16: {
    	   MV_DotProduct_Right_FunctorUnroll<XVector,YVector,16> op;
           op.m_x = x;
           op.m_y = y;
           op.value_count = numVecs;
           KokkosArray::parallel_reduce( x.dimension(0) , op, r );
      	   break;
       }
       case 15: {
    	   MV_DotProduct_Right_FunctorUnroll<XVector,YVector,15> op;
           op.m_x = x;
           op.m_y = y;
           op.value_count = numVecs;
           KokkosArray::parallel_reduce( x.dimension(0) , op, r );
      	   break;
       }
       case 14: {
    	   MV_DotProduct_Right_FunctorUnroll<XVector,YVector,14> op;
           op.m_x = x;
           op.m_y = y;
           op.value_count = numVecs;
           KokkosArray::parallel_reduce( x.dimension(0) , op, r );
      	   break;
       }
       case 13: {
    	   MV_DotProduct_Right_FunctorUnroll<XVector,YVector,13> op;
           op.m_x = x;
           op.m_y = y;
           op.value_count = numVecs;
           KokkosArray::parallel_reduce( x.dimension(0) , op, r );
      	   break;
       }
       case 12: {
    	   MV_DotProduct_Right_FunctorUnroll<XVector,YVector,12> op;
           op.m_x = x;
           op.m_y = y;
           op.value_count = numVecs;
           KokkosArray::parallel_reduce( x.dimension(0) , op, r );
      	   break;
       }
       case 11: {
    	   MV_DotProduct_Right_FunctorUnroll<XVector,YVector,11> op;
           op.m_x = x;
           op.m_y = y;
           op.value_count = numVecs;
           KokkosArray::parallel_reduce( x.dimension(0) , op, r );
      	   break;
       }
       case 10: {
    	   MV_DotProduct_Right_FunctorUnroll<XVector,YVector,10> op;
           op.m_x = x;
           op.m_y = y;
           op.value_count = numVecs;
           KokkosArray::parallel_reduce( x.dimension(0) , op, r );
      	   break;
       }
       case 9: {
    	   MV_DotProduct_Right_FunctorUnroll<XVector,YVector,9> op;
           op.m_x = x;
           op.m_y = y;
           op.value_count = numVecs;
           KokkosArray::parallel_reduce( x.dimension(0) , op, r );
      	   break;
       }
       case 8: {
    	   MV_DotProduct_Right_FunctorUnroll<XVector,YVector,8> op;
           op.m_x = x;
           op.m_y = y;
           op.value_count = numVecs;
           KokkosArray::parallel_reduce( x.dimension(0) , op, r );
      	   break;
       }
       case 7: {
    	   MV_DotProduct_Right_FunctorUnroll<XVector,YVector,7> op;
           op.m_x = x;
           op.m_y = y;
           op.value_count = numVecs;
           KokkosArray::parallel_reduce( x.dimension(0) , op, r );
      	   break;
       }
       case 6: {
    	   MV_DotProduct_Right_FunctorUnroll<XVector,YVector,6> op;
           op.m_x = x;
           op.m_y = y;
           op.value_count = numVecs;
           KokkosArray::parallel_reduce( x.dimension(0) , op, r );
      	   break;
       }
       case 5: {
    	   MV_DotProduct_Right_FunctorUnroll<XVector,YVector,5> op;
           op.m_x = x;
           op.m_y = y;
           op.value_count = numVecs;
           KokkosArray::parallel_reduce( x.dimension(0) , op, r );
      	   break;
       }
       case 4: {
    	   MV_DotProduct_Right_FunctorUnroll<XVector,YVector,4> op;
           op.m_x = x;
           op.m_y = y;
           op.value_count = numVecs;
           KokkosArray::parallel_reduce( x.dimension(0) , op, r );

      	   break;
       }
       case 3: {
    	   MV_DotProduct_Right_FunctorUnroll<XVector,YVector,3> op;
           op.m_x = x;
           op.m_y = y;
           op.value_count = numVecs;
           KokkosArray::parallel_reduce( x.dimension(0) , op, r );
      	   break;
       }
       case 2: {
    	   MV_DotProduct_Right_FunctorUnroll<XVector,YVector,2> op;
           op.m_x = x;
           op.m_y = y;
           op.value_count = numVecs;
           KokkosArray::parallel_reduce( x.dimension(0) , op, r );
      	   break;
       }
       case 1: {
    	   MV_DotProduct_Right_FunctorUnroll<XVector,YVector,1> op;
           op.m_x = x;
           op.m_y = y;
           op.value_count = numVecs;
           KokkosArray::parallel_reduce(x.dimension(0) , op, r);
      	   break;
       }
     }

    return r;
}

}//end namespace KokkosArray
#endif /* KOKKOSARRAY_MULTIVECTOR_H_ */
