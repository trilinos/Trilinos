/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

//Do not use this it has not been tested at all!
#ifndef KOKKOS_BLAS3_IMPL_HPP_
#define KOKKOS_BLAS3_IMPL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>
#ifdef KOKKOS_USE_MKL
#include <mkl.h>
#endif
#ifdef KOKKOS_ENABLE_CUDA
#include <cublas.h>
#endif
namespace KokkosBlas {
namespace Impl {
//Define block size

size_t block_size=32;

template<class AMat, class BMat, class CMat>
struct blas3_right_2_N_N
{
  typedef typename AMat::execution_space        execution_space;
  typedef AMat::size_type                             size_type;

  AMat  a;
  BMat  b;
  CMat  c;
  AMat::non_const_value_type alpha;
  CMat::non_const_value_type beta;
   (AMat::const_value_type alpha_,const AMat &a_, const BMat &b_,CMat::const_value_type beta, CMat &c_) : alpha(alpha_), a (a_), b (b_), beta(beta_), c(c_) {}

  // Prefer const size_type& to const size_type or size_type,
  // since the compiler has an easier time inlining the former.
  KOKKOS_FORCEINLINE_FUNCTION void
  operator() (const size_type &jjblocknum) const
  {
size_type jj=jjblocknum*block_size;

for(size_type i=0;i<c.extent(0);i++){
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif

	for(size_type j = jj; j<((jj+block_size)>c.extent(1)?c.extent(1):(jj+block_size)); j++){
                                c(i,j) = beta*c(i,j);
                                }
}


        for(size_type kk=0;kk<a.extent(1);kk+= block_size){
          
	      for(size_type i=0;i<c.extent(0);i++){
                               
                     for(size_type k = kk; k<((kk+block_size)>a.extent(1)?a.extent(1):(kk+block_size)); k++){
AMat::const_value_type alpha_a=alpha*a(i,k);
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
                                
                        for(size_type j = jj; j<((jj+block_size)>c.extent(1)?c.extent(1):(jj+block_size)); j++){
			        
				c(i,j) +=alpha_a*b(k,j);
                               
				 }
                               
                        }
                }
        }

  }
};


template<class AMat, class BMat, class CMat>
struct blas3_right_2_N_T
{
  typedef typename AMat::execution_space        execution_space;
  typedef AMat::size_type                             size_type;

  AMat  a;
  BMat  b;
  CMat  c;
  AMat::non_const_value_type alpha;
  CMat::non_const_value_type beta;
   (AMat::const_value_type alpha_,const AMat &a_, const BMat &b_,CMat::const_value_type beta, CMat &c_) : alpha(alpha_), a (a_), b (b_), beta(beta_), c(c_) {}

  // Prefer const size_type& to const size_type or size_type,
  // since the compiler has an easier time inlining the former.
  KOKKOS_FORCEINLINE_FUNCTION void
  operator() (const size_type &jjblocknum) const
  {
size_type jj=jjblocknum*block_size;

for(size_type i=0;i<c.extent(0);i++){
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif

	for(size_type j = jj; j<((jj+block_size)>c.extent(1)?c.extent(1):(jj+block_size)); j++){
                                c(i,j) = beta*c(i,j);
                                }
}


        for(size_type kk=0;kk<a.extent(1);kk+= block_size){
          
	      for(size_type i=0;i<c.extent(0);i++){
                               
                                
                        for(size_type j = jj; j<((jj+block_size)>c.extent(1)?c.extent(1):(jj+block_size)); j++){
CMat::non_const_value_type temp=0;
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif

			        for(size_type k = kk; k<((kk+block_size)>a.extent(1)?a.extent(1):(kk+block_size)); k++){

				temp +=alpha*a(i,k)*b(j,k);
                               
				 }
                               c(i,j) += temp;                              
                        }
                }
        }

  }
};



template<class AMat, class BMat, class CMat>
struct blas3_right_2_T_N
{
  typedef typename AMat::execution_space        execution_space;
  typedef AMat::size_type                             size_type;

  AMat  a;
  BMat  b;
  CMat  c;
  AMat::non_const_value_type alpha;
  CMat::non_const_value_type beta;
   (AMat::const_value_type alpha_,const AMat &a_, const BMat &b_,CMat::const_value_type beta, CMat &c_) : alpha(alpha_), a (a_), b (b_), beta(beta_), c(c_) {}

  // Prefer const size_type& to const size_type or size_type,
  // since the compiler has an easier time inlining the former.
  KOKKOS_FORCEINLINE_FUNCTION void
  operator() (const size_type &jjblocknum) const
  {
size_type jj=jjblocknum*block_size;

for(size_type i=0;i<c.extent(0);i++){
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif

	for(size_type j = jj; j<((jj+block_size)>c.extent(1)?c.extent(1):(jj+block_size)); j++){
                                c(i,j) = beta*c(i,j);
                                }
}


        for(size_type kk=0;kk<a.extent(0);kk+= block_size){
          
                               
                     for(size_type k = kk; k<((kk+block_size)>a.extent(0)?a.extent(0):(kk+block_size)); k++){

                          for(size_type i=0;i<c.extent(0);i++){
AMat::const_value_type alpha_a=alpha*a(k,i);
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
                                
                        for(size_type j = jj; j<((jj+block_size)>c.extent(1)?c.extent(1):(jj+block_size)); j++){
			        
				c(i,j) +=alpha_a*b(k,j);
                               
				 }
                               
                        }
                }
        }

  }
};


template<class AMat, class BMat, class CMat>
struct blas3_right_2_T_T
{
  typedef typename AMat::execution_space        execution_space;
  typedef AMat::size_type                             size_type;

  AMat  a;
  BMat  b;
  CMat  c;
  AMat::non_const_value_type alpha;
  CMat::non_const_value_type beta;
   (AMat::const_value_type alpha_,const AMat &a_, const BMat &b_,CMat::const_value_type beta, CMat &c_) : alpha(alpha_), a (a_), b (b_), beta(beta_), c(c_) {}

  // Prefer const size_type& to const size_type or size_type,
  // since the compiler has an easier time inlining the former.
  KOKKOS_FORCEINLINE_FUNCTION void
  operator() (const size_type &jjblocknum) const
  {
size_type jj=jjblocknum*block_size;

for(size_type i=0;i<c.extent(0);i++){
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif

	for(size_type j = jj; j<((jj+block_size)>c.extent(1)?c.extent(1):(jj+block_size)); j++){
                                c(i,j) = beta*c(i,j);
                                }
}


        for(size_type kk=0;kk<a.extent(0);kk+= block_size){
          
	      for(size_type i=0;i<c.extent(0);i++){
                               
                     for(size_type k = kk; k<((kk+block_size)>a.extent(0)?a.extent(0):(kk+block_size)); k++){
 AMat::const_value_type alpha_a=alpha*a(k,i);
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
                                
                        for(size_type j = jj; j<((jj+block_size)>c.extent(1)?c.extent(1):(jj+block_size)); j++){
			        
				c(i,j) +=alpha_a*b(j,k);
                               
				 }
                               
                        }
                }
        }

  }
};


template<class AMat, class BMat, class CMat>
struct blas3_left_2_N_N
{
  typedef typename AMat::execution_space        execution_space;
  typedef AMat::size_type                             size_type;

  AMat  a;
  BMat  b;
  CMat  c;
  AMat::non_const_value_type alpha;
  CMat::non_const_value_type beta;
   (AMat::const_value_type alpha_,const AMat &a_, const BMat &b_,CMat::const_value_type beta, CMat &c_) : alpha(alpha_), a (a_), b (b_), beta(beta_), c(c_) {}

  // Prefer const size_type& to const size_type or size_type,
  // since the compiler has an easier time inlining the former.
  KOKKOS_FORCEINLINE_FUNCTION void
  operator() (const size_type &jjblocknum) const
  {
size_type jj=jjblocknum*block_size;


	for(size_type j = jj; j<((jj+block_size)>c.extent(1)?c.extent(1):(jj+block_size)); j++){
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif 
                 for(size_type i=0;i<c.extent(0);i++){                
                     c(i,j) = beta*c(i,j);
                          
                 }
          }


        for(size_type kk=0;kk<a.extent(1);kk+= block_size){
          
                               
                        for(size_type j = jj; j<((jj+block_size)>c.extent(1)?c.extent(1):(jj+block_size)); j++){

                                for(size_type k = kk; k<((kk+block_size)>a.extent(1)?a.extent(1):(kk+block_size)); k++){
 BMat::const_value_type alpha_b=alpha*b(k,j);
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
                                      
                                        for(size_type i=0;i<c.extent(0);i++){		        
                                     c(i,j) +=alpha_b*a(i,k);
                                
				 }
                               
                        }
                }
        }

  }
};


template<class AMat, class BMat, class CMat>
struct blas3_left_2_N_T
{
  typedef typename AMat::execution_space        execution_space;
  typedef AMat::size_type                             size_type;

  AMat  a;
  BMat  b;
  CMat  c;
  AMat::non_const_value_type alpha;
  CMat::non_const_value_type beta;
   (AMat::const_value_type alpha_,const AMat &a_, const BMat &b_,CMat::const_value_type beta, CMat &c_) : alpha(alpha_), a (a_), b (b_), beta(beta_), c(c_) {}

  // Prefer const size_type& to const size_type or size_type,
  // since the compiler has an easier time inlining the former.
  KOKKOS_FORCEINLINE_FUNCTION void
  operator() (const size_type &jjblocknum) const
  {
size_type jj=jjblocknum*block_size;


        for(size_type j = jj; j<((jj+block_size)>c.extent(1)?c.extent(1):(jj+block_size)); j++){
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
                 for(size_type i=0;i<c.extent(0);i++){
                     c(i,j) = beta*c(i,j);

                 }
          }


        for(size_type kk=0;kk<a.extent(1);kk+= block_size){


                                for(size_type k = kk; k<((kk+block_size)>a.extent(1)?a.extent(1):(kk+block_size)); k++){

                                        for(size_type j = jj; j<((jj+block_size)>c.extent(1)?c.extent(1):(jj+block_size)); j++){

 BMat::const_value_type alpha_b=alpha*b(j,k);
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif

                                        for(size_type i=0;i<c.extent(0);i++){
                                     c(i,j) +=alpha_b*a(i,k);

                                 }

                        }
                }
        }

  }
};



template<class AMat, class BMat, class CMat>
struct blas3_left_2_T_N
{
  typedef typename AMat::execution_space        execution_space;
  typedef AMat::size_type                             size_type;

  AMat  a;
  BMat  b;
  CMat  c;
  AMat::non_const_value_type alpha;
  CMat::non_const_value_type beta;
   (AMat::const_value_type alpha_,const AMat &a_, const BMat &b_,CMat::const_value_type beta, CMat &c_) : alpha(alpha_), a (a_), b (b_), beta(beta_), c(c_) {}

  // Prefer const size_type& to const size_type or size_type,
  // since the compiler has an easier time inlining the former.
  KOKKOS_FORCEINLINE_FUNCTION void
  operator() (const size_type &jjblocknum) const
  {
size_type jj=jjblocknum*block_size;


	for(size_type j = jj; j<((jj+block_size)>c.extent(1)?c.extent(1):(jj+block_size)); j++){
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
                       for(size_type i=0;i<c.extent(0);i++){
    
                            c(i,j) = beta*c(i,j);
                                }
}


        for(size_type kk=0;kk<a.extent(0);kk+= block_size){

              for(size_type i=0;i<c.extent(0);i++){


                        for(size_type j = jj; j<((jj+block_size)>c.extent(1)?c.extent(1):(jj+block_size)); j++){
CMat::non_const_value_type temp=0;
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif

                                for(size_type k = kk; k<((kk+block_size)>a.extent(0)?a.extent(0):(kk+block_size)); k++){

                                temp +=alpha*a(k,i)*b(k,j);

                                 }
                               c(i,j) += temp;
                        }
                }
        }

  }
};


template<class AMat, class BMat, class CMat>
struct blas3_left_2_T_T
{
  typedef typename AMat::execution_space        execution_space;
  typedef AMat::size_type                             size_type;

  AMat  a;
  BMat  b;
  CMat  c;
  AMat::non_const_value_type alpha;
  CMat::non_const_value_type beta;
   (AMat::const_value_type alpha_,const AMat &a_, const BMat &b_,CMat::const_value_type beta, CMat &c_) : alpha(alpha_), a (a_), b (b_), beta(beta_), c(c_) {}

  // Prefer const size_type& to const size_type or size_type,
  // since the compiler has an easier time inlining the former.
  KOKKOS_FORCEINLINE_FUNCTION void
  operator() (const size_type &jjblocknum) const
  {

size_type jj=jjblocknum*block_size;


        for(size_type j = jj; j<((jj+block_size)>c.extent(1)?c.extent(1):(jj+block_size)); j++){
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
                 for(size_type i=0;i<c.extent(0);i++){
                     c(i,j) = beta*c(i,j);

                 }
          }


        for(size_type kk=0;kk<a.extent(0);kk+= block_size){


                                for(size_type k = kk; k<((kk+block_size)>a.extent(0)?a.extent(0):(kk+block_size)); k++){

                                     for(size_type j = jj; j<((jj+block_size)>c.extent(1)?c.extent(1):(jj+block_size)); j++){
 BMat::const_value_type alpha_b=alpha*b(j,k);
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif

                                        for(size_type i=0;i<c.extent(0);i++){
                                     c(i,j) +=alpha_b*a(k,i);

                                 }

                        }
                }
        }


  }
};

template<class AMat, class BMat, class CMat>
struct blas3_right_3_N_N
{
  typedef typename AMat::execution_space        execution_space;
  typedef AMat::size_type                             size_type;

  AMat  a;
  BMat  b;
  CMat  c;
  AMat::non_const_value_type alpha;
  CMat::non_const_value_type beta;
   (AMat::const_value_type alpha_,const AMat &a_, const BMat &b_,CMat::const_value_type beta, CMat &c_) : alpha(alpha_), a (a_), b (b_), beta(beta_), c(c_) {}

  // Prefer const size_type& to const size_type or size_type,
  // since the compiler has an easier time inlining the former.
  KOKKOS_FORCEINLINE_FUNCTION void
  operator() (const size_type &matnum) const
  {

for(size_type jj=0;jj<c.extent(2);jj+= block_size){
	for(size_type i=0;i<c.extent(1);i++){
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif

	for(size_type j = jj; j<((jj+block_size)>c.extent(2)?c.extent(2):(jj+block_size)); j++){
                                c(matnum,i,j) = beta*c(matnum,i,j);
                                }
}


        for(size_type kk=0;kk<a.extent(2);kk+= block_size){
          
	      for(size_type i=0;i<c.extent(1);i++){
                               
                     for(size_type k = kk; k<((kk+block_size)>a.extent(2)?a.extent(2):(kk+block_size)); k++){
AMat::const_value_type alpha_a=alpha*a(matnum,i,k);
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
                                
                        for(size_type j = jj; j<((jj+block_size)>c.extent(2)?c.extent(2):(jj+block_size)); j++){
			        
				c(matnum,i,j) +=alpha_a*b(matnum,k,j);
                               
				 }
                               
                        }
                }
        }
}
  }
};


template<class AMat, class BMat, class CMat>
struct blas3_right_3_N_T
{
  typedef typename AMat::execution_space        execution_space;
  typedef AMat::size_type                             size_type;

  AMat  a;
  BMat  b;
  CMat  c;
  AMat::non_const_value_type alpha;
  CMat::non_const_value_type beta;
   (AMat::const_value_type alpha_,const AMat &a_, const BMat &b_,CMat::const_value_type beta, CMat &c_) : alpha(alpha_), a (a_), b (b_), beta(beta_), c(c_) {}

  // Prefer const size_type& to const size_type or size_type,
  // since the compiler has an easier time inlining the former.
  KOKKOS_FORCEINLINE_FUNCTION void
  operator() (const size_type &matnum) const
  {
for(size_type jj=0;jj<c.extent(2);jj+= block_size){

for(size_type i=0;i<c.extent(1);i++){
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif

        for(size_type j = jj; j<((jj+block_size)>c.extent(2)?c.extent(2):(jj+block_size)); j++){
                                c(matnum,i,j) = beta*c(matnum,i,j);
                                }
}


        for(size_type kk=0;kk<a.extent(2);kk+= block_size){

              for(size_type i=0;i<c.extent(1);i++){


                        for(size_type j = jj; j<((jj+block_size)>c.extent(2)?c.extent(2):(jj+block_size)); j++){
CMat::non_const_value_type temp=0;
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif

                                for(size_type k = kk; k<((kk+block_size)>a.extent(2)?a.extent(2):(kk+block_size)); k++){

                                temp +=alpha*a(matnum,i,k)*b(matnum,j,k);

                                 }
                               c(matnum,i,j) += temp;
                        }
                }
        }
}
  }
};



template<class AMat, class BMat, class CMat>
struct blas3_right_3_T_N
{
  typedef typename AMat::execution_space        execution_space;
  typedef AMat::size_type                             size_type;

  AMat  a;
  BMat  b;
  CMat  c;
  AMat::non_const_value_type alpha;
  CMat::non_const_value_type beta;
   (AMat::const_value_type alpha_,const AMat &a_, const BMat &b_,CMat::const_value_type beta, CMat &c_) : alpha(alpha_), a (a_), b (b_), beta(beta_), c(c_) {}

  // Prefer const size_type& to const size_type or size_type,
  // since the compiler has an easier time inlining the former.
  KOKKOS_FORCEINLINE_FUNCTION void
  operator() (const size_type &matnum) const
  {
for(size_type jj=0;jj<c.extent(2);jj+= block_size){

for(size_type i=0;i<c.extent(1);i++){
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif

	for(size_type j = jj; j<((jj+block_size)>c.extent(2)?c.extent(2):(jj+block_size)); j++){
                                c(matnum,i,j) = beta*c(matnum,i,j);
                                }
}


        for(size_type kk=0;kk<a.extent(1);kk+= block_size){
          
                               
                     for(size_type k = kk; k<((kk+block_size)>a.extent(1)?a.extent(1):(kk+block_size)); k++){

                          for(size_type i=0;i<c.extent(1);i++){
AMat::const_value_type alpha_a=alpha*a(matnum,k,i);
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
                                
                        for(size_type j = jj; j<((jj+block_size)>c.extent(2)?c.extent(2):(jj+block_size)); j++){
			        
				c(matnum,i,j) +=alpha_a*b(matnum,k,j);
                               
				 }
                               
                        }
                }
        }
}
  }
};


template<class AMat, class BMat, class CMat>
struct blas3_right_3_T_T
{
  typedef typename AMat::execution_space        execution_space;
  typedef AMat::size_type                             size_type;

  AMat  a;
  BMat  b;
  CMat  c;
  AMat::non_const_value_type alpha;
  CMat::non_const_value_type beta;
   (AMat::const_value_type alpha_,const AMat &a_, const BMat &b_,CMat::const_value_type beta, CMat &c_) : alpha(alpha_), a (a_), b (b_), beta(beta_), c(c_) {}

  // Prefer const size_type& to const size_type or size_type,
  // since the compiler has an easier time inlining the former.
  KOKKOS_FORCEINLINE_FUNCTION void
  operator() (const size_type &matnum) const
  {
for(size_type jj=0;jj<c.extent(2);jj+= block_size){

for(size_type i=0;i<c.extent(1);i++){
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif

	for(size_type j = jj; j<((jj+block_size)>c.extent(2)?c.extent(2):(jj+block_size)); j++){
                                c(matnum,i,j) = beta*c(matnum,i,j);
                                }
}


        for(size_type kk=0;kk<a.extent(1);kk+= block_size){
          
	      for(size_type i=0;i<c.extent(1);i++){
                               
                     for(size_type k = kk; k<((kk+block_size)>a.extent(1)?a.extent(1):(kk+block_size)); k++){
 AMat::const_value_type alpha_a=alpha*a(k,i);
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
                                
                        for(size_type j = jj; j<((jj+block_size)>c.extent(2)?c.extent(2):(jj+block_size)); j++){
			        
				c(i,j) +=alpha_a*b(j,k);
                               
				 }
                               
                        }
                }
        }
}
  }
};


template<class AMat, class BMat, class CMat>
struct blas3_left_3_N_N
{
  typedef typename AMat::execution_space        execution_space;
  typedef AMat::size_type                             size_type;

  AMat  a;
  BMat  b;
  CMat  c;
  AMat::non_const_value_type alpha;
  CMat::non_const_value_type beta;
   (AMat::const_value_type alpha_,const AMat &a_, const BMat &b_,CMat::const_value_type beta, CMat &c_) : alpha(alpha_), a (a_), b (b_), beta(beta_), c(c_) {}

  // Prefer const size_type& to const size_type or size_type,
  // since the compiler has an easier time inlining the former.
  KOKKOS_FORCEINLINE_FUNCTION void
  operator() (const size_type &jjblocknum) const
  {
size_type jj=jjblocknum*block_size;


	for(size_type j = jj; j<((jj+block_size)>c.extent(2)?c.extent(2):(jj+block_size)); j++){

                 for(size_type i=0;i<c.extent(1);i++){            
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
                    for(size_type matnum=0;matnum<c.extent(0);matnum++){
    
                     c(matnum,i,j) = beta*c(matnum,i,j);
                        }  
                 }
          }


        for(size_type kk=0;kk<a.extent(2);kk+= block_size){
          
                               
                        for(size_type j = jj; j<((jj+block_size)>c.extent(2)?c.extent(2):(jj+block_size)); j++){

                                for(size_type k = kk; k<((kk+block_size)>a.extent(2)?a.extent(2):(kk+block_size)); k++){
                                      
                                        for(size_type i=0;i<c.extent(1);i++){		        
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
                    for(size_type matnum=0;matnum<c.extent(0);matnum++){

                                     c(matnum,i,j) +=alpha*b(matnum,k,j)*a(matnum,i,k);
                         }       
				 }
                               
                        }
                }
        }

  }
};


template<class AMat, class BMat, class CMat>
struct blas3_left_3_N_T
{
  typedef typename AMat::execution_space        execution_space;
  typedef AMat::size_type                             size_type;

  AMat  a;
  BMat  b;
  CMat  c;
  AMat::non_const_value_type alpha;
  CMat::non_const_value_type beta;
   (AMat::const_value_type alpha_,const AMat &a_, const BMat &b_,CMat::const_value_type beta, CMat &c_) : alpha(alpha_), a (a_), b (b_), beta(beta_), c(c_) {}

  // Prefer const size_type& to const size_type or size_type,
  // since the compiler has an easier time inlining the former.
  KOKKOS_FORCEINLINE_FUNCTION void
  operator() (const size_type &jjblocknum) const
  {
size_type jj=jjblocknum*block_size;


        for(size_type j = jj; j<((jj+block_size)>c.extent(2)?c.extent(2):(jj+block_size)); j++){

                 for(size_type i=0;i<c.extent(1);i++){
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
                    for(size_type matnum=0;matnum<c.extent(0);matnum++){   
                  c(matnum,i,j) = beta*c(matnum,i,j);

                 }
          }
}

        for(size_type kk=0;kk<a.extent(2);kk+= block_size){


                                for(size_type k = kk; k<((kk+block_size)>a.extent(2)?a.extent(2):(kk+block_size)); k++){

                                        for(size_type j = jj; j<((jj+block_size)>c.extent(2)?c.extent(2):(jj+block_size)); j++){


                                        for(size_type i=0;i<c.extent(1);i++){
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif

                                               for(size_type matnum=0;matnum<c.extent(0);matnum++){ 
                                    c(matnum,i,j) +=alpha*b(matnum,j,k)*a(matnum,i,k);

                                 }

                        }
                }
        }

  }
};



template<class AMat, class BMat, class CMat>
struct blas3_left_3_T_N
{
  typedef typename AMat::execution_space        execution_space;
  typedef AMat::size_type                             size_type;

  AMat  a;
  BMat  b;
  CMat  c;
  AMat::non_const_value_type alpha;
  CMat::non_const_value_type beta;
   (AMat::const_value_type alpha_,const AMat &a_, const BMat &b_,CMat::const_value_type beta, CMat &c_) : alpha(alpha_), a (a_), b (b_), beta(beta_), c(c_) {}

  // Prefer const size_type& to const size_type or size_type,
  // since the compiler has an easier time inlining the former.
  KOKKOS_FORCEINLINE_FUNCTION void
  operator() (const size_type &jjblocknum) const
  {
size_type jj=jjblocknum*block_size;


	for(size_type j = jj; j<((jj+block_size)>c.extent(2)?c.extent(2):(jj+block_size)); j++){
                       for(size_type i=0;i<c.extent(1);i++){
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif

                           for(size_type matnum=0;matnum<c.extent(0);matnum++){    
                            c(matnum,i,j) = beta*c(matnum,i,j);
                                }
                }
}

        for(size_type kk=0;kk<a.extent(1);kk+= block_size){

              for(size_type i=0;i<c.extent(1);i++){


                        for(size_type j = jj; j<((jj+block_size)>c.extent(2)?c.extent(2):(jj+block_size)); j++){

                                for(size_type k = kk; k<((kk+block_size)>a.extent(1)?a.extent(1):(kk+block_size)); k++){
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif

                                    for(size_type matnum=0;matnum<c.extent(0);matnum++){
                                c(matnum,i,j) +=alpha*a(matnum,k,i)*b(matnum,k,j);

                                 }
                              
                        }
                }
        }

  }
};


template<class AMat, class BMat, class CMat>
struct blas3_left_3_T_T
{
  typedef typename AMat::execution_space        execution_space;
  typedef AMat::size_type                             size_type;

  AMat  a;
  BMat  b;
  CMat  c;
  AMat::non_const_value_type alpha;
  CMat::non_const_value_type beta;
   (AMat::const_value_type alpha_,const AMat &a_, const BMat &b_,CMat::const_value_type beta, CMat &c_) : alpha(alpha_), a (a_), b (b_), beta(beta_), c(c_) {}

  // Prefer const size_type& to const size_type or size_type,
  // since the compiler has an easier time inlining the former.
  KOKKOS_FORCEINLINE_FUNCTION void
  operator() (const size_type &jjblocknum) const
  {

size_type jj=jjblocknum*block_size;


        for(size_type j = jj; j<((jj+block_size)>c.extent(2)?c.extent(2):(jj+block_size)); j++){
                 for(size_type i=0;i<c.extent(1);i++){
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif

                                    for(size_type matnum=0;matnum<c.extent(0);matnum++){
                     c(matnum,i,j) = beta*c(matnum,i,j);

                 }
          }


        for(size_type kk=0;kk<a.extent(1);kk+= block_size){


                                for(size_type k = kk; k<((kk+block_size)>a.extent(1)?a.extent(1):(kk+block_size)); k++){

                                     for(size_type j = jj; j<((jj+block_size)>c.extent(2)?c.extent(2):(jj+block_size)); j++){


                                        for(size_type i=0;i<c.extent(1);i++){
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif

                                    for(size_type matnum=0;matnum<c.extent(0);matnum++){ 
                                     c(matnum,i,j) +=alpha*b(matnum,j,k)*a(matnum,k,i);

                                 }

                        }
                }
        }


  }
};



  template<typename ScalarA,typename ScalarB, typename ScalarC, typename eSpace,typename LayoutA,typename LayoutB,typename LayoutC,typename sizeType,int rank>
        struct MultiGemm;

template <typename ScalarA,typename ScalarB, typename ScalarC,typename eSpace, typename LayoutA,typename LayoutB,typename LayoutC,typename sizeType>
struct MultiGemm<ScalarA,ScalarB,ScalarC,eSpace,LayoutA,LayoutB,LayoutC,sizetype,2>{
        static void GEMM(const char transA, const char transB, ScalarA alpha,
          Kokkos::View<ScalarA**,LayoutA,eSpace> A, Kokkos::View<ScalarB**,LayoutB,eSpace> B,
          ScalarC beta, Kokkos::View<ScalarB**,LayoutC,eSpace> C){
if((transA=='n'||transA=='N') && (transB=='n'||transB=='N')){
sizetype jblocks=0;
if(C.extent(1)%block_size==0){
jblocks=C.extent(1)/block_size;
}else{
jblocks=C.extent(1)/block_size+1;
}
Kokkos::parallel_for(jblocks,blas3_right_2_N_N<Kokkos::View<ScalarA**,LayoutA,eSpace>,Kokkos::View<ScalarB**,LayoutB,eSpace>,Kokkos::View<ScalarC**,LayoutC,eSpace> >(alpha,A,B,beta,C));
}else if((transA=='n'||transA=='N') && (transB=='t'||transB=='T')){
sizetype jblocks=0;
if(C.extent(1)%block_size==0){
jblocks=C.extent(1)/block_size;
}else{
jblocks=C.extent(1)/block_size+1;
}
Kokkos::parallel_for(jblocks,blas3_right_2_N_T<Kokkos::View<ScalarA**,LayoutA,eSpace>,Kokkos::View<ScalarB**,LayoutB,eSpace>,Kokkos::View<ScalarC**,LayoutC,eSpace> >(alpha,A,B,beta,C));
}else if((transA=='t'||transA=='T') && (transB=='n'||transB=='N')){
sizetype jblocks=0;
if(C.extent(1)%block_size==0){
jblocks=C.extent(1)/block_size;
}else{
jblocks=C.extent(1)/block_size+1;
}
Kokkos::parallel_for(jblocks,blas3_right_2_T_N<Kokkos::View<ScalarA**,LayoutA,eSpace>,Kokkos::View<ScalarB**,LayoutB,eSpace>,Kokkos::View<ScalarC**,LayoutC,eSpace> >(alpha,A,B,beta,C));
}else if((transA=='t'||transA=='T') && (transB=='t'||transB=='T')){
sizetype jblocks=0;
if(C.extent(1)%block_size==0){
jblocks=C.extent(1)/block_size;
}else{
jblocks=C.extent(1)/block_size+1;
}
Kokkos::parallel_for(jblocks,blas3_right_2_T_T<Kokkos::View<ScalarA**,LayoutA,eSpace>,Kokkos::View<ScalarB**,LayoutB,eSpace>,Kokkos::View<ScalarC**,LayoutC,eSpace> >(alpha,A,B,beta,C));
}



}
};

template <typename ScalarA,typename ScalarB, typename ScalarC,typename eSpace, typename LayoutA,typename LayoutB,typename LayoutC,typename sizeType>
struct MultiGemm<ScalarA,ScalarB,ScalarC,eSpace,LayoutA,LayoutB,LayoutC,sizetype,3>{
        static void GEMM(const char transA, const char transB, ScalarA alpha,
          Kokkos::View<ScalarA***,LayoutA,eSpace> A, Kokkos::View<ScalarB***,LayoutB,eSpace> B,
          ScalarC beta, Kokkos::View<ScalarC***,LayoutC,eSpace> C){
if((transA=='n'||transA=='N') && (transB=='n'||transB=='N')){

Kokkos::parallel_for(C.extent(0),blas3_right_3_N_N<Kokkos::View<ScalarA***,LayoutA,eSpace>,Kokkos::View<ScalarB***,LayoutB,eSpace>,Kokkos::View<ScalarC***,LayoutC,eSpace> >(alpha,A,B,beta,C));

}else if((transA=='n'||transA=='N') && (transB=='t'||transB=='T')){

Kokkos::parallel_for(C.extent(0),blas3_right_3_N_T<Kokkos::View<ScalarA***,LayoutA,eSpace>,Kokkos::View<ScalarB***,LayoutB,eSpace>,Kokkos::View<ScalarC***,LayoutC,eSpace> >(alpha,A,B,beta,C));

}else if((transA=='t'||transA=='T') && (transB=='n'||transB=='N')){

Kokkos::parallel_for(C.extent(0),blas3_right_3_T_N<Kokkos::View<ScalarA***,LayoutA,eSpace>,Kokkos::View<ScalarB***,LayoutB,eSpace>,Kokkos::View<ScalarC***,LayoutC,eSpace> >(alpha,A,B,beta,C));

}else if((transA=='t'||transA=='T') && (transB=='t'||transB=='T')){

Kokkos::parallel_for(C.extent(0),blas3_right_3_T_T<Kokkos::View<ScalarA***,LayoutA,eSpace>,Kokkos::View<ScalarB***,LayoutB,eSpace>,Kokkos::View<ScalarC***,LayoutC,eSpace> >(alpha,A,B,beta,C));

}

}
};
template <typename ScalarA,typename ScalarB, typename ScalarC,typename eSpace,typename sizeType>
struct MultiGemm<ScalarA,ScalarB,ScalarC,eSpace,Kokkos::LayoutRight,Kokkos::LayoutRight,Kokkos::LayoutRight,sizetype,2>{
        static void GEMM(const char transA, const char transB, ScalarA alpha,
          Kokkos::View<ScalarA**,Kokkos::LayoutRight,eSpace> A, Kokkos::View<ScalarB**,Kokkos::LayoutRight,eSpace> B,
          ScalarC beta, Kokkos::View<ScalarC**,Kokkos::LayoutRight,eSpace> C){
if((transA=='n'||transA=='N') && (transB=='n'||transB=='N')){
sizetype jblocks=0;
if(C.extent(1)%block_size==0){
jblocks=C.extent(1)/block_size;
}else{
jblocks=C.extent(1)/block_size+1;
}
Kokkos::parallel_for(jblocks,blas3_right_2_N_N<Kokkos::View<ScalarA**,LayoutA,eSpace>,Kokkos::View<ScalarB**,LayoutB,eSpace>,Kokkos::View<ScalarB**,LayoutC,eSpace> >(alpha,A,B,beta,C));
}else if((transA=='n'||transA=='N') && (transB=='t'||transB=='T')){
sizetype jblocks=0;
if(C.extent(1)%block_size==0){
jblocks=C.extent(1)/block_size;
}else{
jblocks=C.extent(1)/block_size+1;
}
Kokkos::parallel_for(jblocks,blas3_right_2_N_T<Kokkos::View<ScalarA**,LayoutA,eSpace>,Kokkos::View<ScalarB**,LayoutB,eSpace>,Kokkos::View<ScalarB**,LayoutC,eSpace> >(alpha,A,B,beta,C));
}else if((transA=='t'||transA=='T') && (transB=='n'||transB=='N')){
sizetype jblocks=0;
if(C.extent(1)%block_size==0){
jblocks=C.extent(1)/block_size;
}else{
jblocks=C.extent(1)/block_size+1;
}
Kokkos::parallel_for(jblocks,blas3_right_2_T_N<Kokkos::View<ScalarA**,LayoutA,eSpace>,Kokkos::View<ScalarB**,LayoutB,eSpace>,Kokkos::View<ScalarB**,LayoutC,eSpace> >(alpha,A,B,beta,C));
}else if((transA=='t'||transA=='T') && (transB=='t'||transB=='T')){
sizetype jblocks=0;
if(C.extent(1)%block_size==0){
jblocks=C.extent(1)/block_size;
}else{
jblocks=C.extent(1)/block_size+1;
}
Kokkos::parallel_for(jblocks,blas3_right_2_T_T<Kokkos::View<ScalarA**,LayoutA,eSpace>,Kokkos::View<ScalarB**,LayoutB,eSpace>,Kokkos::View<ScalarB**,LayoutC,eSpace> >(alpha,A,B,beta,C));
}

}
};

template <typename ScalarA,typename ScalarB, typename ScalarC,typename eSpace,typename sizeType>
struct MultiGemm<ScalarA,ScalarB,ScalarC,eSpace,Kokkos::LayoutRight,Kokkos::LayoutRight,Kokkos::LayoutRight,sizetype,3>{
        static void GEMM(const char transA, const char transB, ScalarA alpha,
          Kokkos::View<ScalarA***,Kokkos::LayoutRight,eSpace> A, Kokkos::View<ScalarB***,Kokkos::LayoutRight,eSpace> B,
          ScalarC beta, Kokkos::View<ScalarC***,Kokkos::LayoutRight,eSpace> C){
if((transA=='n'||transA=='N') && (transB=='n'||transB=='N')){

Kokkos::parallel_for(C.extent(0),blas3_right_3_N_N<Kokkos::View<ScalarA***,Kokkos::LayoutRight,eSpace>,Kokkos::View<ScalarB***,Kokkos::LayoutRight,eSpace>,Kokkos::View<ScalarC***,Kokkos::LayoutRight,eSpace> >(alpha,A,B,beta,C));

}else if((transA=='n'||transA=='N') && (transB=='t'||transB=='T')){

Kokkos::parallel_for(C.extent(0),blas3_right_3_N_T<Kokkos::View<ScalarA***,Kokkos::LayoutRight,eSpace>,Kokkos::View<ScalarB***,Kokkos::LayoutRight,eSpace>,Kokkos::View<ScalarC***,Kokkos::LayoutRight,eSpace> >(alpha,A,B,beta,C));

}else if((transA=='t'||transA=='T') && (transB=='n'||transB=='N')){

Kokkos::parallel_for(C.extent(0),blas3_right_3_T_N<Kokkos::View<ScalarA***,Kokkos::LayoutRight,eSpace>,Kokkos::View<ScalarB***,Kokkos::LayoutRight,eSpace>,Kokkos::View<ScalarC***,Kokkos::LayoutRight,eSpace> >(alpha,A,B,beta,C));

}else if((transA=='t'||transA=='T') && (transB=='t'||transB=='T')){

Kokkos::parallel_for(C.extent(0),blas3_right_3_T_T<Kokkos::View<ScalarA***,Kokkos::LayoutRight,eSpace>,Kokkos::View<ScalarB***,Kokkos::LayoutRight,eSpace>,Kokkos::View<ScalarC***,Kokkos::LayoutRight,eSpace> >(alpha,A,B,beta,C));

}

}
};

template <typename ScalarA,typename ScalarB, typename ScalarC,typename eSpace,typename sizeType>
struct MultiGemm<ScalarA,ScalarB,ScalarC,eSpace,Kokkos::LayoutLeft,Kokkos::LayoutLeft,Kokkos::LayoutLeft,sizetype,2>{
        static void GEMM(const char transA, const char transB, ScalarA alpha,
          Kokkos::View<ScalarA**,Kokkos::LayoutLeft,eSpace> A, Kokkos::View<ScalarB**,Kokkos::LayoutLeft,eSpace> B,
          ScalarC beta, Kokkos::View<ScalarC**,Kokkos::LayoutLeft,eSpace> C){
if((transA=='n'||transA=='N') && (transB=='n'||transB=='N')){
sizetype jblocks=0;
if(C.extent(1)%block_size==0){
jblocks=C.extent(1)/block_size;
}else{
jblocks=C.extent(1)/block_size+1;
}
Kokkos::parallel_for(jblocks,blas3_left_2_N_N<Kokkos::View<ScalarA**,LayoutA,eSpace>,Kokkos::View<ScalarB**,LayoutB,eSpace>,Kokkos::View<ScalarB**,LayoutC,eSpace> >(alpha,A,B,beta,C));
}else if((transA=='n'||transA=='N') && (transB=='t'||transB=='T')){
sizetype jblocks=0;
if(C.extent(1)%block_size==0){
jblocks=C.extent(1)/block_size;
}else{
jblocks=C.extent(1)/block_size+1;
}
Kokkos::parallel_for(jblocks,blas3_left_2_N_T<Kokkos::View<ScalarA**,LayoutA,eSpace>,Kokkos::View<ScalarB**,LayoutB,eSpace>,Kokkos::View<ScalarB**,LayoutC,eSpace> >(alpha,A,B,beta,C));
}else if((transA=='t'||transA=='T') && (transB=='n'||transB=='N')){
sizetype jblocks=0;
if(C.extent(1)%block_size==0){
jblocks=C.extent(1)/block_size;
}else{
jblocks=C.extent(1)/block_size+1;
}
Kokkos::parallel_for(jblocks,blas3_left_2_T_N<Kokkos::View<ScalarA**,LayoutA,eSpace>,Kokkos::View<ScalarB**,LayoutB,eSpace>,Kokkos::View<ScalarB**,LayoutC,eSpace> >(alpha,A,B,beta,C));
}else if((transA=='t'||transA=='T') && (transB=='t'||transB=='T')){
sizetype jblocks=0;
if(C.extent(1)%block_size==0){
jblocks=C.extent(1)/block_size;
}else{
jblocks=C.extent(1)/block_size+1;
}
Kokkos::parallel_for(jblocks,blas3_left_2_T_T<Kokkos::View<ScalarA**,LayoutA,eSpace>,Kokkos::View<ScalarB**,LayoutB,eSpace>,Kokkos::View<ScalarB**,LayoutC,eSpace> >(alpha,A,B,beta,C));
}

}
};

template <typename ScalarA,typename ScalarB, typename ScalarC,typename eSpace,typename sizeType>
struct MultiGemm<ScalarA,ScalarB,ScalarC,eSpace,Kokkos::LayoutLeft,Kokkos::LayoutLeft,Kokkos::LayoutLeft,sizetype,3>{
        static void GEMM(const char transA, const char transB, ScalarA alpha,
          Kokkos::View<ScalarA***,Kokkos::LayoutLeft,eSpace> A, Kokkos::View<ScalarB***,Kokkos::LayoutLeft,eSpace> B,
          ScalarC beta, Kokkos::View<ScalarC***,Kokkos::LayoutLeft,eSpace> C){
if((transA=='n'||transA=='N') && (transB=='n'||transB=='N')){

Kokkos::parallel_for(C.extent(0),blas3_left_3_N_N<Kokkos::View<ScalarA***,Kokkos::LayoutLeft,eSpace>,Kokkos::View<ScalarB***,Kokkos::LayoutLeft,eSpace>,Kokkos::View<ScalarC***,Kokkos::LayoutLeft,eSpace> >(alpha,A,B,beta,C));

}else if((transA=='n'||transA=='N') && (transB=='t'||transB=='T')){

Kokkos::parallel_for(C.extent(0),blas3_left_3_N_T<Kokkos::View<ScalarA***,Kokkos::LayoutLeft,eSpace>,Kokkos::View<ScalarB***,Kokkos::LayoutLeft,eSpace>,Kokkos::View<ScalarC***,Kokkos::LayoutLeft,eSpace> >(alpha,A,B,beta,C));

}else if((transA=='t'||transA=='T') && (transB=='n'||transB=='N')){

Kokkos::parallel_for(C.extent(0),blas3_left_3_T_N<Kokkos::View<ScalarA***,Kokkos::LayoutLeft,eSpace>,Kokkos::View<ScalarB***,Kokkos::LayoutLeft,eSpace>,Kokkos::View<ScalarC***,Kokkos::LayoutLeft,eSpace> >(alpha,A,B,beta,C));

}else if((transA=='t'||transA=='T') && (transB=='t'||transB=='T')){

Kokkos::parallel_for(C.extent(0),blas3_left_3_T_T<Kokkos::View<ScalarA***,Kokkos::LayoutLeft,eSpace>,Kokkos::View<ScalarB***,Kokkos::LayoutLeft,eSpace>,Kokkos::View<ScalarC***,Kokkos::LayoutLeft,eSpace> >(alpha,A,B,beta,C));

}

}
};

#ifdef KOKKOS_USE_MKL


struct MultiGemm<double,double,double,Kokkos::Serial,Kokkos::LayoutRight,Kokkos::LayoutRight,Kokkos::LayoutRight,int,2>{
        static void GEMM(const char transA, const char transB, double alpha,
          Kokkos::View<double**,Kokkos::LayoutRight,Kokkos::Serial> A, Kokkos::View<double**,Kokkos::LayoutRight,Kokkos::Serial> B,
          double beta, Kokkos::View<double**,Kokkos::LayoutRight,Kokkos::Serial> C){
const int m=C.extent(0); 
const int n=C.extent(1);
const int k= (transA == 'N'||transA == 'n') ? A.extent(1) : A.extent(0);	  
const int lda=(transA == 'N'||transA == 'n') ? k : m;
const int ldb=(transA == 'N'||transA == 'n') ? n : k;
const int ldc= n; 
cblas_dgemm (CblasRowMajor, transA, transB, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, C.data(), ldc);
}
};

struct MultiGemm<double,double,double,Kokkos::OpenMP,Kokkos::LayoutRight,Kokkos::LayoutRight,Kokkos::LayoutRight,int,2>{
        static void GEMM(const char transA, const char transB, double alpha,
          Kokkos::View<double**,Kokkos::LayoutRight,Kokkos::OpenMP> A, Kokkos::View<double**,Kokkos::LayoutRight,Kokkos::OpenMP> B,
          double beta, Kokkos::View<double**,Kokkos::LayoutRight,Kokkos::OpenMP> C){
const int m=C.extent(0);
const int n=C.extent(1);
const int k= (transA == 'N'||transA == 'n') ? A.extent(1) : A.extent(0);
const int lda=(transA == 'N'||transA == 'n') ? k : m;
const int ldb=(transA == 'N'||transA == 'n') ? n : k;
const int ldc= n;       
cblas_dgemm (CblasRowMajor, transA, transB, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, C.data(), ldc);

}
};

struct MultiGemm<double,double,double,Kokkos::Serial,Kokkos::LayoutLeft,Kokkos::LayoutLeft,Kokkos::LayoutLeft,int,2>{
        static void GEMM(const char transA, const char transB, double alpha,
          Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::Serial> A, Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::Serial> B,
          double beta, Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::Serial> C){
const int m=C.extent(0);
const int n=C.extent(1);
const int k= (transA == 'N'||transA == 'n') ? A.extent(1) : A.extent(0);
const int lda=(transA == 'N'||transA == 'n') ? m : k;
const int ldb=(transA == 'N'||transA == 'n') ? k : n;
const int ldc= m;
cblas_dgemm (CblasColMajor, transA, transB, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, C.data(), ldc);

}
};

struct MultiGemm<double,double,double,Kokkos::OpenMP,Kokkos::LayoutLeft,Kokkos::LayoutLeft,Kokkos::LayoutLeft,int,2>{
        static void GEMM(const char transA, const char transB, double alpha,
          Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::OpenMP> A, Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::OpenMP> B,
          double beta, Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::OpenMP> C){
const int m=C.extent(0);
const int n=C.extent(1);
const int k= (transA == 'N'||transA == 'n') ? A.extent(1) : A.extent(0);
const int lda=(transA == 'N'||transA == 'n') ? m : k;
const int ldb=(transA == 'N'||transA == 'n') ? k : n;
const int ldc= m;
cblas_dgemm (CblasColMajor, transA, transB, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, C.data(), ldc);


}
};


struct MultiGemm<float,float,float,Kokkos::Serial,Kokkos::LayoutRight,Kokkos::LayoutRight,Kokkos::LayoutRight,int,2>{
        static void GEMM(const char transA, const char transB, float alpha,
          Kokkos::View<float**,Kokkos::LayoutRight,Kokkos::Serial> A, Kokkos::View<float**,Kokkos::LayoutRight,Kokkos::Serial> B,
          float beta, Kokkos::View<float**,Kokkos::LayoutRight,Kokkos::Serial> C){
const int m=C.extent(0);
const int n=C.extent(1);
const int k= (transA == 'N'||transA == 'n') ? A.extent(1) : A.extent(0);
const int lda=(transA == 'N'||transA == 'n') ? k : m;
const int ldb=(transA == 'N'||transA == 'n') ? n : k;
const int ldc= n;       
cblas_sgemm (CblasRowMajor, transA, transB, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, C.data(), ldc);

}
};

struct MultiGemm<float,float,float,Kokkos::OpenMP,Kokkos::LayoutRight,Kokkos::LayoutRight,Kokkos::LayoutRight,int,2>{
        static void GEMM(const char transA, const char transB, float alpha,
          Kokkos::View<float**,Kokkos::LayoutRight,Kokkos::OpenMP> A, Kokkos::View<float**,Kokkos::LayoutRight,Kokkos::OpenMP> B,
          float beta, Kokkos::View<float**,Kokkos::LayoutRight,Kokkos::OpenMP> C){
const int m=C.extent(0);
const int n=C.extent(1);
const int k= (transA == 'N'||transA == 'n') ? A.extent(1) : A.extent(0);
const int lda=(transA == 'N'||transA == 'n') ? k : m;
const int ldb=(transA == 'N'||transA == 'n') ? n : k;
const int ldc= n;
cblas_sgemm (CblasRowMajor, transA, transB, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, C.data(), ldc);

}
};

struct MultiGemm<float,float,float,Kokkos::Serial,Kokkos::LayoutLeft,Kokkos::LayoutLeft,Kokkos::LayoutLeft,int,2>{
        static void GEMM(const char transA, const char transB, float alpha,
          Kokkos::View<float**,Kokkos::LayoutLeft,Kokkos::Serial> A, Kokkos::View<float**,Kokkos::LayoutLeft,Kokkos::Serial> B,
          float beta, Kokkos::View<float**,Kokkos::LayoutLeft,Kokkos::Serial> C){
const int m=C.extent(0);
const int n=C.extent(1);
const int k= (transA == 'N'||transA == 'n') ? A.extent(1) : A.extent(0);
const int lda=(transA == 'N'||transA == 'n') ? m : k;
const int ldb=(transA == 'N'||transA == 'n') ? k : n;
const int ldc= m;
cblas_sgemm (CblasColMajor, transA, transB, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, C.data(), ldc);


}
};

struct MultiGemm<float,float,float,Kokkos::OpenMP,Kokkos::LayoutLeft,Kokkos::LayoutLeft,Kokkos::LayoutLeft,int,2>{
        static void GEMM(const char transA, const char transB, float alpha,
          Kokkos::View<float**,Kokkos::LayoutLeft,Kokkos::OpenMP> A, Kokkos::View<float**,Kokkos::LayoutLeft,Kokkos::OpenMP> B,
          float beta, Kokkos::View<float**,Kokkos::LayoutLeft,Kokkos::OpenMP> C){
const int m=C.extent(0);
const int n=C.extent(1);
const int k= (transA == 'N'||transA == 'n') ? A.extent(1) : A.extent(0);
const int lda=(transA == 'N'||transA == 'n') ? m : k;
const int ldb=(transA == 'N'||transA == 'n') ? k : n;
const int ldc= m;
cblas_sgemm (CblasColMajor, transA, transB, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, C.data(), ldc);


}
};

struct MultiGemm<int,int,int,Kokkos::Serial,Kokkos::LayoutRight,Kokkos::LayoutRight,Kokkos::LayoutRight,int,2>{
        static void GEMM(const char transA, const char transB, int alpha,
          Kokkos::View<int**,Kokkos::LayoutRight,Kokkos::Serial> A, Kokkos::View<int**,Kokkos::LayoutRight,Kokkos::Serial> B,
          int beta, Kokkos::View<int**,Kokkos::LayoutRight,Kokkos::Serial> C){
const int m=C.extent(0);
const int n=C.extent(1);
const int k= (transA == 'N'||transA == 'n') ? A.extent(1) : A.extent(0);
const int lda=(transA == 'N'||transA == 'n') ? k : m;
const int ldb=(transA == 'N'||transA == 'n') ? n : k;
const int ldc= n;
cblas_zgemm (CblasRowMajor, transA, transB, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, C.data(), ldc);

}
};

struct MultiGemm<int,int,int,Kokkos::OpenMP,Kokkos::LayoutRight,Kokkos::LayoutRight,Kokkos::LayoutRight,int,2>{
        static void GEMM(const char transA, const char transB, int alpha,
          Kokkos::View<int**,Kokkos::LayoutRight,Kokkos::OpenMP> A, Kokkos::View<int**,Kokkos::LayoutRight,Kokkos::OpenMP> B,
          int beta, Kokkos::View<int**,Kokkos::LayoutRight,Kokkos::OpenMP> C){
const int m=C.extent(0);
const int n=C.extent(1);
const int k= (transA == 'N'||transA == 'n') ? A.extent(1) : A.extent(0);
const int lda=(transA == 'N'||transA == 'n') ? k : m;
const int ldb=(transA == 'N'||transA == 'n') ? n : k;
const int ldc= n;
cblas_zgemm (CblasRowMajor, transA, transB, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, C.data(), ldc);

}
};

struct MultiGemm<int,int,int,Kokkos::Serial,Kokkos::LayoutLeft,Kokkos::LayoutLeft,Kokkos::LayoutLeft,int,2>{
        static void GEMM(const char transA, const char transB, int alpha,
          Kokkos::View<int**,Kokkos::LayoutLeft,Kokkos::Serial> A, Kokkos::View<int**,Kokkos::LayoutLeft,Kokkos::Serial> B,
          int beta, Kokkos::View<int**,Kokkos::LayoutLeft,Kokkos::Serial> C){
const int m=C.extent(0);
const int n=C.extent(1);
const int k= (transA == 'N'||transA == 'n') ? A.extent(1) : A.extent(0);
const int lda=(transA == 'N'||transA == 'n') ? m : k;
const int ldb=(transA == 'N'||transA == 'n') ? k : n;
const int ldc= m;
cblas_zgemm (CblasColMajor, transA, transB, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, C.data(), ldc);


}
};

struct MultiGemm<int,int,int,Kokkos::OpenMP,Kokkos::LayoutLeft,Kokkos::LayoutLeft,Kokkos::LayoutLeft,int,2>{
        static void GEMM(const char transA, const char transB, int alpha,
          Kokkos::View<int**,Kokkos::LayoutLeft,Kokkos::OpenMP> A, Kokkos::View<int**,Kokkos::LayoutLeft,Kokkos::OpenMP> B,
          int beta, Kokkos::View<int**,Kokkos::LayoutLeft,Kokkos::OpenMP> C){
const int m=C.extent(0);
const int n=C.extent(1);
const int k= (transA == 'N'||transA == 'n') ? A.extent(1) : A.extent(0);
const int lda=(transA == 'N'||transA == 'n') ? m : k;
const int ldb=(transA == 'N'||transA == 'n') ? k : n;
const int ldc= m;
cblas_zgemm (CblasColMajor, transA, transB, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, C.data(), ldc);


}
};


#endif

#ifdef KOKKOS_ENABLE_CUDA

struct MultiGemm<double,double,double,Kokkos::Cuda,Kokkos::LayoutRight,Kokkos::LayoutRight,Kokkos::LayoutRight,int,2>{
        static void GEMM(const char transA, const char transB, double alpha,
          Kokkos::View<double**,Kokkos::LayoutRight,Kokkos::Cuda> A, Kokkos::View<double**,Kokkos::LayoutRight,Kokkos::Cuda> B,
          double beta, Kokkos::View<double**,Kokkos::LayoutRight,Kokkos::Cuda> C){
const int m=C.extent(0); 
const int n=C.extent(1);
const int k= (transA == 'N'||transA == 'n') ? A.extent(1) : A.extent(0);	  
const int lda=(transA == 'N'||transA == 'n') ? k : m;
const int ldb=(transA == 'N'||transA == 'n') ? n : k;
const int ldc= n; 
cublasDgemm(transB,transA, n, m, k, alpha, B.data(), ldb, A.data(), lda, beta, C.data(), ldc);

}
};


struct MultiGemm<double,double,double,Kokkos::Cuda,Kokkos::LayoutLeft,Kokkos::LayoutLeft,Kokkos::LayoutLeft,int,2>{
        static void GEMM(const char transA, const char transB, double alpha,
          Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::Cuda> A, Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::Cuda> B,
          double beta, Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::Cuda> C){
const int m=C.extent(0);
const int n=C.extent(1);
const int k= (transA == 'N'||transA == 'n') ? A.extent(1) : A.extent(0);
const int lda=(transA == 'N'||transA == 'n') ? m : k;
const int ldb=(transA == 'N'||transA == 'n') ? k : n;
const int ldc= m;
cublasDgemm (transA, transB, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, C.data(), ldc);

}
};



struct MultiGemm<float,float,float,Kokkos::Cuda,Kokkos::LayoutRight,Kokkos::LayoutRight,Kokkos::LayoutRight,int,2>{
        static void GEMM(const char transA, const char transB, float alpha,
          Kokkos::View<float**,Kokkos::LayoutRight,Kokkos::Cuda> A, Kokkos::View<float**,Kokkos::LayoutRight,Kokkos::Cuda> B,
          float beta, Kokkos::View<float**,Kokkos::LayoutRight,Kokkos::Cuda> C){
const int m=C.extent(0);
const int n=C.extent(1);
const int k= (transA == 'N'||transA == 'n') ? A.extent(1) : A.extent(0);
const int lda=(transA == 'N'||transA == 'n') ? k : m;
const int ldb=(transA == 'N'||transA == 'n') ? n : k;
const int ldc= n;
cublasSgemm(transB,transA, n, m, k, alpha, B.data(), ldb, A.data(), lda, beta, C.data(), ldc);

}
};


struct MultiGemm<float,float,float,Kokkos::Cuda,Kokkos::LayoutLeft,Kokkos::LayoutLeft,Kokkos::LayoutLeft,int,2>{
        static void GEMM(const char transA, const char transB, float alpha,
          Kokkos::View<float**,Kokkos::LayoutLeft,Kokkos::Cuda> A, Kokkos::View<float**,Kokkos::LayoutLeft,Kokkos::Cuda> B,
          float beta, Kokkos::View<float**,Kokkos::LayoutLeft,Kokkos::Cuda> C){
const int m=C.extent(0);
const int n=C.extent(1);
const int k= (transA == 'N'||transA == 'n') ? A.extent(1) : A.extent(0);
const int lda=(transA == 'N'||transA == 'n') ? m : k;
const int ldb=(transA == 'N'||transA == 'n') ? k : n;
const int ldc= m;
cublasSgemm ( transA, transB, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, C.data(), ldc);


}
};


struct MultiGemm<int,int,int,Kokkos::Cuda,Kokkos::LayoutRight,Kokkos::LayoutRight,Kokkos::LayoutRight,int,2>{
        static void GEMM(const char transA, const char transB, int alpha,
          Kokkos::View<int**,Kokkos::LayoutRight,Kokkos::Cuda> A, Kokkos::View<int**,Kokkos::LayoutRight,Kokkos::Cuda> B,
          int beta, Kokkos::View<int**,Kokkos::LayoutRight,Kokkos::Cuda> C){
const int m=C.extent(0);
const int n=C.extent(1);
const int k= (transA == 'N'||transA == 'n') ? A.extent(1) : A.extent(0);
const int lda=(transA == 'N'||transA == 'n') ? k : m;
const int ldb=(transA == 'N'||transA == 'n') ? n : k;
const int ldc= n;
cublasZgemm(transB,transA, n, m, k, alpha, B.data(), ldb, A.data(), lda, beta, C.data(), ldc);

}
};


struct MultiGemm<int,int,int,Kokkos::Cuda,Kokkos::LayoutLeft,Kokkos::LayoutLeft,Kokkos::LayoutLeft,int,2>{
        static void GEMM(const char transA, const char transB, int alpha,
          Kokkos::View<int**,Kokkos::LayoutLeft,Kokkos::Cuda> A, Kokkos::View<int**,Kokkos::LayoutLeft,Kokkos::Cuda> B,
          int beta, Kokkos::View<int**,Kokkos::LayoutLeft,Kokkos::Cuda> C){
const int m=C.extent(0);
const int n=C.extent(1);
const int k= (transA == 'N'||transA == 'n') ? A.extent(1) : A.extent(0);
const int lda=(transA == 'N'||transA == 'n') ? m : k;
const int ldb=(transA == 'N'||transA == 'n') ? k : n;
const int ldc= m;
cublasZgemm ( transA, transB, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, C.data(), ldc);


}
};




#endif

}
}
