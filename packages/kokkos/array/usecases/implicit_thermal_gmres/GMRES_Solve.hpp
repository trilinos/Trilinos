/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#include <limits>
#include <WAXSBY.hpp>
#include <CRSMatVec.hpp>
#include <impl/Kokkos_Timer.hpp>

template< typename Scalar , class DeviceType , unsigned N = 1 >
struct InvMultiVectorScale ;

template< typename Scalar >
struct InvMultiVectorScale<Scalar, KOKKOS_MACRO_DEVICE ,1> {
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef Scalar                  value_type ;
  typedef device_type::size_type  size_type ;

  Kokkos::MultiVectorView<Scalar,device_type> Y ;
  Kokkos::ValueView<Scalar,device_type>     S ;

  InvMultiVectorScale(
    const Kokkos::MultiVectorView<Scalar,device_type> & argY ,
    const Kokkos::ValueView<Scalar,device_type>       & argS )
    : Y( argY ), S( argS ) {}

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type iwork ) const
  { Y(iwork) /= *S ; }
};


template< typename Scalar , class DeviceType , unsigned N = 1 >
struct DotSingle;

template< typename Scalar >
struct DotSingle<Scalar, KOKKOS_MACRO_DEVICE , 1 > {
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef Scalar                  value_type ;
  typedef device_type::size_type  size_type ;
  typedef Kokkos::MultiVectorView<value_type,device_type> vector_type ;

  vector_type X ;

  DotSingle( const vector_type & argX ) : X( argX ) {}

  DotSingle( const vector_type & argX ,
             const size_type     argI ) : X( argX , argI ) {}

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type iwork , value_type & update ) const
  {
    const Scalar value = X(iwork);
    update += value * value ;
  }

  KOKKOS_MACRO_DEVICE_FUNCTION
  static void join( volatile value_type & update , const volatile value_type & source )
  { update += source ; }

  KOKKOS_MACRO_DEVICE_FUNCTION
  static void init( value_type & update )
  { update = 0 ; }
};

template< typename Scalar , class DeviceType , unsigned N = 1 >
struct MultiVectorYSAX ;

template< typename Scalar >
struct MultiVectorYSAX<Scalar, KOKKOS_MACRO_DEVICE ,1> {
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef Scalar                  value_type ;
  typedef device_type::size_type  size_type ;

  Kokkos::MultiVectorView<Scalar,device_type> Y , X ;
  Kokkos::ValueView<Scalar,device_type>     A ;

  MultiVectorYSAX( const Kokkos::MultiVectorView<Scalar,device_type> & argY ,
                   const Kokkos::ValueView<Scalar,device_type>     & argA ,
                   const Kokkos::MultiVectorView<Scalar,device_type> & argX )
   : Y( argY ), X( argX ), A( argA ) {}

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type iwork ) const
  { Y(iwork) -= X(iwork) * *A ; }
};

template< typename Scalar , class DeviceType , unsigned N = 1 >
struct Dot ;

template< typename Scalar >
struct Dot< Scalar, KOKKOS_MACRO_DEVICE , 1 > {
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef Scalar                  value_type ;
  typedef device_type::size_type  size_type ;
  typedef Kokkos::MultiVectorView<value_type,device_type> vector_type ;

  vector_type X , Y ;

  Dot( const vector_type & argX ,
       const vector_type & argY )
  : X( argX ) , Y( argY ) {}

  Dot( const vector_type argX , const size_type argI ,
       const vector_type argY , const size_type argJ )
  : X( argX , argI ), Y( argY , argJ ) {}

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type iwork , value_type & update ) const
  { update += X(iwork) * Y(iwork); }

  KOKKOS_MACRO_DEVICE_FUNCTION
  static void join( volatile value_type & update , const volatile value_type & source )
  { update += source ; }

  KOKKOS_MACRO_DEVICE_FUNCTION
  static void init( value_type & update )
  { update = 0 ; }
};


template <class Scalar , class DeviceType >
struct GMRES_Solve;


template<class Scalar>
struct GMRES_Solve<Scalar , KOKKOS_MACRO_DEVICE>
{

  typedef KOKKOS_MACRO_DEVICE    device_type;
  typedef device_type::size_type index_type ;
  typedef Kokkos::MultiVectorView<Scalar , device_type>  scalar_vector;
  typedef Kokkos::MultiVectorView<index_type , device_type>    index_vector;
  typedef Kokkos::MultiVectorView<Scalar,device_type> MultiVector ;

  typedef Kokkos::ValueView<Scalar , device_type>     value;
  typedef Kokkos::ValueView<Scalar , Kokkos::DeviceHost> host_val;
  typedef MultiVectorYSAX<  Scalar , device_type , 1 > YSAX ;
  typedef InvMultiVectorScale< Scalar , device_type , 1 > InvScale ;
  typedef Dot< Scalar , device_type , 1 > RealDot ;
  typedef device_type::size_type  size_type ;

  struct Norm2 {
    value norm;

    Norm2(const value & argNorm):
      norm(argNorm)
      {}

    KOKKOS_MACRO_DEVICE_FUNCTION
    void operator()( Scalar & result ) const
    {
      *norm = sqrt( result );
    }
  };

  struct DotM {

    MultiVector R ;
    size_type j , k ;

    DotM( const MultiVector & argR ,
          size_type argJ ,
          size_type argK )
      : R( argR )
      , j( argJ )
      , k( argK )
      {}

    KOKKOS_MACRO_DEVICE_FUNCTION
    void operator()( Scalar & result ) const
    {
       R(j,k) = result ;
    }
  };

  typedef DotSingle< Scalar, device_type > dot_single ;


  
  // Return megaflops / second for iterations

  static double run( scalar_vector & A_value ,
                     index_vector  & A_row ,
                     index_vector & A_offsets ,
                     scalar_vector & b ,
                     scalar_vector & x,
                     const size_t num_iters)
  {

    //Value view used for temp stuff.
    value tmp = Kokkos::create_value<Scalar , device_type>();


    const size_t rows = A_row.length()-1;

    int iteration = 0 ;
    
    value one  = Kokkos::create_value<Scalar , device_type>();
    value zero = Kokkos::create_value<Scalar , device_type>();

    Kokkos::deep_copy( one,  Scalar( 1 ) );
    Kokkos::deep_copy( zero, Scalar( 0 ) );

    // Solvers' working temporaries:
    scalar_vector r = 
      Kokkos::create_labeled_multivector<scalar_vector>("r",rows);
    scalar_vector Ax = 
      Kokkos::create_labeled_multivector<scalar_vector>("Ax",rows);

    //Compute initial residual r = b - A*x

    // create CRSMatVec object for A * x
    CRSMatVec<Scalar,device_type> A_mult_x(A_value, A_row , A_offsets , x , Ax);
    A_mult_x.apply(); // Ap = A * p

    // r = b - Ax
    Kokkos::parallel_for(
      rows , WAXSBY<Scalar , device_type>(one , b , one , Ax , r) );
    
    value beta  = Kokkos::create_value<Scalar , device_type>();
    Kokkos::deep_copy( zero, Scalar( 0 ) );

    // compute norm2 of r
    Kokkos::parallel_reduce(rows, dot_single(r), Norm2(beta) );

    Scalar beta_copy;
    Kokkos::deep_copy(beta_copy,beta);
    if(beta_copy == 0.0){
      //Stuff didn't work
      return 0.0;
    }

    //Allocate space for basis vectors Q. Q has as many rows as b and m+1 
    //columns.
    scalar_vector Q = 
      Kokkos::create_labeled_multivector<scalar_vector>("Q",rows, num_iters+1);

    //Allocate m+1 by m matrix H, and fill it with zeros.
    scalar_vector H = Kokkos::create_labeled_multivector<scalar_vector>(
      "r",
      num_iters, 
      num_iters+1);


    //Q(:, 0) = r ./ beta (elementwise division)
    Kokkos::parallel_for(rows, InvScale(MultiVector(Q,0), beta));

    Kokkos::Impl::Timer wall_clock ;
    for(size_t j =1; j<= num_iters; ++j, ++iteration){
      //Q(:,j) = A * Q(:, j-1)
      CRSMatVec<Scalar,device_type> A_mult_q(
        A_value, A_row , A_offsets , MultiVector(Q, j-1) , MultiVector(Q, j));
      A_mult_q.apply();

      //H(0:j, j-1) = Q(:, 0:j-1)^* Q(:,j)
      for (size_t k = 0; k <= j-1; ++k) { 
        Kokkos::parallel_reduce(rows, 
          RealDot(Q, k, Q, j),
          DotM(H, k, j-1));
      }

      // Q(:, j) = Q(:, j) - Q(:, 0:j-1) * H(0:j, j-1)
      for(size_t i =0; i<j; ++i){
        Kokkos::deep_copy(tmp, H(i,j-1));
        Kokkos::parallel_for(
          rows, YSAX(MultiVector(Q,j), tmp, MultiVector(Q,i)));
      }
     
      //H(j, j-1) = norm(Q(:,j),2) 
      Kokkos::parallel_reduce(rows, dot_single(Q, j), Norm2(tmp) );
      Kokkos::deep_copy(H(j,j-1), tmp);
      //Q(:,j) = Q(:, j) / H(j, j-1)
      Kokkos::parallel_for(rows, InvScale(MultiVector(Q, j), tmp));
       
    }

    const double iter_time = wall_clock.seconds();

    // Compute floating point operations performed during iterations

    size_t iter_dot_flops    = ( 3 * iteration ) * ( 2 * rows );
    size_t iter_matvec_flops = iteration * ( 2 * A_offsets.length() );
    size_t iter_ysax_flops = iteration * (2 * rows);
    size_t iter_scale_flops = iteration * rows;
    size_t iter_total_flops  = iter_dot_flops + iter_matvec_flops + iter_ysax_flops + iter_scale_flops;

    return (double) iter_total_flops / ( iter_time * 1.0e6 );
  }
};


