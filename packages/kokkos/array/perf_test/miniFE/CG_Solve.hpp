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
#include <WAXPBY.hpp>
#include <Dot.hpp>
#include <Divide.hpp>
#include <CRSMatVec.hpp>
#include <impl/Kokkos_Timer.hpp>


template <class Scalar , class DeviceType >
struct CG_Solve;


template<class Scalar>
struct CG_Solve<Scalar , KOKKOS_MACRO_DEVICE>
{

  typedef KOKKOS_MACRO_DEVICE    device_type;
  typedef device_type::size_type index_type ;
  typedef Kokkos::MultiVectorView<Scalar , device_type>  scalar_vector;
  typedef Kokkos::MultiVectorView<index_type , device_type>    index_vector;

  typedef Kokkos::ValueView<Scalar , device_type>     value;
  typedef Kokkos::ValueView<Scalar , Kokkos::DeviceHost> host_val;

  
  // Return megaflops / second for iterations

  static double run( scalar_vector & A_value ,
                     index_vector  & A_row ,
                     index_vector & A_col ,
                     scalar_vector & b ,
                     scalar_vector & x ,
                     double* times)
  {
    const Scalar tolerance = std::numeric_limits<Scalar>::epsilon();
    const int maximum_iteration = 200 ;
    const size_t rows = A_row.length()-1;

    value one  = Kokkos::create_value<Scalar , device_type>();
    value zero = Kokkos::create_value<Scalar , device_type>();

    Kokkos::deep_copy( one,  Scalar( 1 ) );
    Kokkos::deep_copy( zero, Scalar( 0 ) );

    // Solvers' working temporaries:

    scalar_vector r = Kokkos::create_labeled_multivector<scalar_vector>("r",rows);
    scalar_vector p = Kokkos::create_labeled_multivector<scalar_vector>("p",rows);
    scalar_vector Ap = Kokkos::create_labeled_multivector<scalar_vector>("Ap",rows);

    value rtrans    = Kokkos::create_value<Scalar , device_type>();
    value ptrans    = Kokkos::create_value<Scalar , device_type>();
    value oldrtrans = Kokkos::create_value<Scalar , device_type>();  
    value alpha     = Kokkos::create_value<Scalar , device_type>();          
    value beta      = Kokkos::create_value<Scalar , device_type>();  

    double normr = 1000; 
    
    Kokkos::deep_copy( p , x );

    // create CRSMatVec object for A * p
    CRSMatVec<Scalar,device_type> A_mult_p(A_value, A_row , A_col , p , Ap);

    A_mult_p.MatVec( 1.0 , 0.0 ); // Ap = 1.0 * A * p + 0.0 * Ap
  
    // r = b - Ap
    Kokkos::parallel_for(rows , WAXSBY<Scalar , device_type>(one , b , one , Ap , r) );
    
    // p = r
    Kokkos::deep_copy( p , r );

    int iteration = 0 ;

    Kokkos::Impl::Timer wall_clock ;

    while ( normr > tolerance  && iteration < maximum_iteration )
    {
      int k;

      // Iterate 25 times before checking residual to
      // avoid the device-host synchronization
      // required by a copy-back of residual data.

      for(k = 0; k < 25 ; k++) 
      {
        // Ap = A*p
        A_mult_p.MatVec(1.0 , 0.0); 
        
        // ptrans = p • Ap
        Kokkos::parallel_reduce(rows , Dot<Scalar , device_type>(p , Ap) , ptrans );

        // oldrtrans = r • r
        // alpha = rtrans / ptrans
        Kokkos::parallel_reduce(rows , Dot<Scalar , device_type>(r , r) ,  
                        Divide<Scalar , device_type>(ptrans,alpha,oldrtrans));

        // x = x + alpha * p 
        Kokkos::parallel_for(rows , WAXPBY<Scalar , device_type>(one , x , alpha , p , x) );
    
        // r = rk - alpha * Ap
        Kokkos::parallel_for(rows , WAXSBY<Scalar , device_type>(one , r , alpha , Ap , r) );

        // rtrans = r • r
        // beta = rtrans / oldrtrans
        Kokkos::parallel_reduce(rows , Dot<Scalar , device_type>(r , r) , 
                        Divide<Scalar , device_type>(oldrtrans , beta, rtrans) );
                                  
        // p = r + beta * p
        Kokkos::parallel_for(rows , WAXPBY<Scalar , device_type>(one , r , beta , p , p) );
      }

      Scalar check = 0 ;
      Kokkos::deep_copy(check,oldrtrans);

      iteration += k ;
      normr = sqrt(check);

//      std::cout<<"Iteration: "<<iteration<<" Residual "<<normr<<std::endl;
    }

    const double iter_time = wall_clock.seconds();

    // Compute floating point operations performed during iterations

    size_t iter_waxpby_flops = ( 3 * iteration ) * ( 3 * rows ); // y = a * x + b * y
    size_t iter_dot_flops    = ( 3 * iteration ) * ( 2 * rows );
    size_t iter_matvec_flops = iteration * ( 2 * A_col.length() );
    size_t iter_total_flops  = iter_waxpby_flops + iter_dot_flops + iter_matvec_flops ;

    return (double) iter_total_flops / ( iter_time * 1.0e6 );
  }
};

