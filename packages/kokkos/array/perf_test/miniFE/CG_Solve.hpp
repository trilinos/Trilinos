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

  
  static double run(scalar_vector & A_value , index_vector & A_row , index_vector & A_col , scalar_vector & b , scalar_vector & x , double* times)
  {
    double time = 0.0;
    double total_time = 0.0;
    double waxpby_time = 0.0;
    double matvec_time = 0.0;
    double dot_time = 0.0;
    double iter_time = 0.0;
    host_val one_host = Kokkos::create_value<Scalar , Kokkos::DeviceHost>();
    host_val zero_host = Kokkos::create_value<Scalar , Kokkos::DeviceHost>();
    *one_host = 1.0; *zero_host = 0.0; 
    value one = Kokkos::create_value<Scalar , device_type>();
    value zero = Kokkos::create_value<Scalar , device_type>();
    Kokkos::deep_copy(one,one_host);
    Kokkos::deep_copy(zero,zero_host);
    Scalar tolerance =  std::numeric_limits<Scalar>::epsilon();

    const size_t rows = A_row.length()-1;
    scalar_vector r = Kokkos::create_labeled_multivector<scalar_vector>("r",rows);
    scalar_vector p = Kokkos::create_labeled_multivector<scalar_vector>("p",rows);
    scalar_vector Ap = Kokkos::create_labeled_multivector<scalar_vector>("Ap",rows);


    value rtrans = Kokkos::create_value<Scalar , device_type>();
    value ptrans = Kokkos::create_value<Scalar , device_type>();
    value oldrtrans = Kokkos::create_value<Scalar , device_type>();  
    host_val check = Kokkos::create_value<Scalar , Kokkos::DeviceHost>();
    double normr = 1000; 
    
    timeval start,stop,result;  
    value alpha = Kokkos::create_value<Scalar , device_type>();          
    value beta = Kokkos::create_value<Scalar , device_type>();  

    //create CRSMatVec object for A * x
    CRSMatVec<Scalar,device_type> U(A_value, A_row , A_col , x , Ap);
    
    gettimeofday(&start, NULL);
    U.MatVec(1.0 , 0.0);  //Compute Ap = A * x
    device_type::wait_functor_completion();  //Synchronize threads
    gettimeofday(&stop, NULL);  
    
    timersub(&stop, &start, &result);
    time = (result.tv_sec + result.tv_usec/1000000.0);
    matvec_time += time;
    
    //create CRSMatVec object for A * p
    CRSMatVec<Scalar,device_type> V(A_value, A_row , A_col , p , Ap);
  
    //r = b - Ap
    Kokkos::parallel_for(rows , WAXSBY<Scalar , device_type>(one , b , one , Ap , r),time);
    waxpby_time += time;
    
    //p = r
    Kokkos::parallel_for(rows , WAXPBY<Scalar , device_type>(one , r , zero , r , p), time);
    waxpby_time += time;
            
       double iteration = 0;
    int k;
    while ( normr > tolerance )
    {

      gettimeofday(&start , NULL);
      //Run 25 times before checking residual to allow pipelining
      for(k = 0; k < 25 ; k++) 
      {
        //Ap = A*p
        V.MatVec(1.0 , 0.0); 
        
        //ptrans = p • Ap
        Kokkos::parallel_reduce(rows , Dot<Scalar , device_type>(p , Ap) , ptrans );

        //oldrtrans = r • r
        //alpha = rtrans / ptrans
        Kokkos::parallel_reduce(rows , Dot<Scalar , device_type>(r , r) ,  
                        Divide<Scalar , device_type>(ptrans,alpha,oldrtrans));

        //x = x + alpha*p 
        Kokkos::parallel_for(rows , WAXPBY<Scalar , device_type>(one , x , alpha , p , x) );
    
        //r = rk - alpha*Ap
        Kokkos::parallel_for(rows , WAXSBY<Scalar , device_type>(one , r , alpha , Ap , r) );

        //rtrans = r • r
        //beta = rtrans / oldrtrans
        Kokkos::parallel_reduce(rows , Dot<Scalar , device_type>(r , r) , 
                        Divide<Scalar , device_type>(oldrtrans , beta, rtrans) );
                                  
        // p = r + beta*p
        Kokkos::parallel_for(rows , WAXPBY<Scalar , device_type>(one , r , beta , p , p) );

      }
      
      device_type::wait_functor_completion();
      gettimeofday(&stop, NULL);
      timersub(&stop, &start, &result);
      time = (result.tv_sec + result.tv_usec/1000000.0);
      iter_time += time;

      iteration+=k;
      Kokkos::deep_copy(check,oldrtrans);
      normr = sqrt(*check);

//      std::cout<<"Iteration: "<<iteration<<" Residual "<<normr<<std::endl;
    }
#if 0
    //Compute floating point operations
    double waxpby_flOps =   (3 * (3 * iteration + 2) * rows);
    double dot_flOps =   (2 * (3* iteration) * rows );
    double matvec_flOps = ( (2 * A_col.dimension(0)  + rows) * (iteration + 1) );

    double total_flOps = waxpby_flOps + dot_flOps + matvec_flOps;

    times[2] = (total_flOps/(1e6 * total_time));  //Mega floating point operations/second
#endif

    total_time = waxpby_time + dot_time + matvec_time + iter_time;
    times[4] = iter_time/iteration; //iteration time/# of iterations
    times[5] = total_time; //Total CG Solve Time

    return total_time;
  }
};
