/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

template< typename Scalar >
struct Dot< Scalar , KOKKOS_MACRO_DEVICE , Impl::unsigned_<2> >
{
  typedef KOKKOS_MACRO_DEVICE              device_type;
  typedef device_type::size_type           size_type;
  typedef MultiVector<Scalar, device_type> scalar_vector;  
  typedef Value < double , device_type >   result_type ;  
  typedef double                           value_type;

private:

  scalar_vector x;
  scalar_vector y;

public:

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( int iwork , value_type & update ) const 
  { update += x(iwork) * y(iwork); }
    
  KOKKOS_MACRO_DEVICE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & source )
  { update += source;    }
    
  KOKKOS_MACRO_DEVICE_FUNCTION
  static void init( value_type & update )
  { update = 0 ; }

  inline static
  void apply( const size_t n ,
              const scalar_vector & x ,
              const scalar_vector & y ,
              const result_type   & result )
  {
    Dot op ; op.x = x ; op.y = y ;
    parallel_reduce( n , op , result );
  }
}; //Dot

template< typename Scalar >
struct Dot< Scalar , KOKKOS_MACRO_DEVICE , Impl::unsigned_<1> >
{
  typedef KOKKOS_MACRO_DEVICE              device_type;
  typedef device_type::size_type           size_type;
  typedef MultiVector<Scalar, device_type> scalar_vector;  
  typedef Value < double , device_type >   result_type ;  
  typedef double                           value_type;

private:

  scalar_vector x;

public:

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( int iwork , value_type & update ) const 
  { const Scalar xi = x(iwork); update += xi * xi ; }
    
  KOKKOS_MACRO_DEVICE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & source )
  { update += source ; }
    
  KOKKOS_MACRO_DEVICE_FUNCTION
  static void init( value_type & update )
  { update = 0 ; }

  inline static
  void apply( const size_t n ,
              const scalar_vector & x ,
              const result_type   & result )
  {
    Dot op ; op.x = x ;
    parallel_reduce( n , op , result );
  }
}; //Dot

//----------------------------------------------------------------------------

template < typename Scalar >
struct FILL<Scalar , KOKKOS_MACRO_DEVICE >
{
  typedef KOKKOS_MACRO_DEVICE               device_type ;
  typedef device_type::size_type            size_type ;
  typedef MultiVector<Scalar, device_type>  scalar_vector ;

private:

  scalar_vector w ;
  Scalar alpha ;

public:

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()(int inode) const
  {
    w(inode) = alpha ;
  }

  inline static
  void apply( const size_t n ,
              const double          alpha ,
              const scalar_vector & w )
  {
    FILL op ;
    op.w = w ;
    op.alpha = alpha ;
    parallel_for( n , op );
  }
};

template < typename Scalar >
struct WAXPBY<Scalar , KOKKOS_MACRO_DEVICE >
{
  typedef KOKKOS_MACRO_DEVICE               device_type ;
  typedef device_type::size_type            size_type ;
  typedef MultiVector<Scalar, device_type>  scalar_vector ;

private:

  scalar_vector w ;
  scalar_vector x ;
  scalar_vector y ;

  Scalar alpha ;
  Scalar beta ;

public:

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()(int inode) const
  {
    w(inode) = alpha * x(inode) + beta * y(inode);
  }

  inline static
  void apply( const size_t n ,
              const double          alpha ,
              const scalar_vector & x ,
              const double          beta ,
              const scalar_vector & y ,
              const scalar_vector & w )
  {
    WAXPBY op ;
    op.w = w ;
    op.x = x ;
    op.y = y ;
    op.alpha = alpha ;
    op.beta  = beta ;
    parallel_for( n , op );
  }
}; // WAXPBY

template < typename Scalar >
struct AXPBY<Scalar , KOKKOS_MACRO_DEVICE >
{
  typedef KOKKOS_MACRO_DEVICE               device_type ;
  typedef device_type::size_type            size_type ;
  typedef MultiVector<Scalar, device_type>  scalar_vector ;

private:

  scalar_vector x ;
  scalar_vector y ;

  Scalar alpha ;
  Scalar beta ;

public:

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()(int inode) const
  {
    Scalar & val = y(inode); val = alpha * x(inode) + beta * val ;
  }

  inline static
  void apply( const size_t n ,
              const double          alpha ,
              const scalar_vector & x ,
              const double          beta ,
              const scalar_vector & y )
  {
    AXPBY op ;
    op.x = x ;
    op.y = y ;
    op.alpha = alpha ;
    op.beta  = beta ;
    parallel_for( n , op );
  }
}; // WAXPBY

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

