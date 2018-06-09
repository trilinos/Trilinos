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
#ifndef KOKKOS_BLAS3_HPP_
#define KOKKOS_BLAS3_HPP_

#include <Kokkos_Blas3_impl.hpp>
#include <type_traits> // requires C++11

namespace KokkosBlas {

template<class AMat,class BMat,class CMat>
void gemm(const char transA,const char transB,AMat::const_value_type alpha ,const AMat& a, const BMat& b,CMat::const_value_type beta, const CMat &c)
{
  static_assert (Kokkos::Impl::is_view<AMat>::value,
                 "KokkosBlas::gemm: AMat must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<BMat>::value,
                 "KokkosBlas::gemm: BMat must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<CMat>::value,
                 "KokkosBlas::gemm: CMat must be a Kokkos::View.");
  static_assert ( AMat::rank !=  BMat::rank && CMat::rank !=  BMat::rank,
                 "KokkosBlas::gemm: Matrix ranks do not match.");
  static_assert ( AMat::rank != 2 || AMat::rank != 3 ,
                 "KokkosBlas::gemm: Matrix ranks do not match.");
  // Check compatibility of dimensions at run time.
  if ((transA=='n'||transA=='N')&&(transB=='n'||transB=='N')) {
if(rank==2){
if(a.extent(0)!=c.extent(0)||a.extent(1)!=b.extent(0)||b.extent(1)!=c.extent(1)){
    std::ostringstream os;
    os << "KokkosBlas::gemm: Matrix Dimensions Dimensions do not match for A notrans and B notrans: "
       << ", A: " << "("<<a.extent(0)<<","<<a.extent(1)<<")"
       << ", B: " << "("<<b.extent(0)<<","<<b.extent(1)<<")"
       << ", C: " << "("<<c.extent(0)<<","<<c.extent(1)<<")";
    Kokkos::Impl::throw_runtime_exception (os.str ());
}
}
else if(rank==3){
if(a.extent(1)!=c.extent(1)||a.extent(2)!=b.extent(1)||b.extent(2)!=c.extent(2)||a.extent(0)!=b.extent(0)||b.extent(0)!=c.extent(0)){
    std::ostringstream os;
    os << "KokkosBlas::gemm: Matrix Dimensions Dimensions do not match for A notrans and B notrans: "
       << ", A: " << "("<<a.extent(0)<<","<<a.extent(1)<<","<<a.extent(2)<<")"
       << ", B: " << "("<<b.extent(0)<<","<<b.extent(1)<<","<<b.extent(2)<<")"
       << ", C: " << "("<<c.extent(0)<<","<<c.extent(1)<<","<<c.extent(2)<<")";
    Kokkos::Impl::throw_runtime_exception (os.str ());
}
}

}else if((transA=='n'||transA=='N')&&(transB=='t'||transB=='T')){
if(rank==2){
if(a.extent(0)!=c.extent(0)||a.extent(1)!=b.extent(1)||b.extent(0)!=c.extent(1)){
    std::ostringstream os;
    os << "KokkosBlas::gemm: Matrix Dimensions Dimensions do not match for A notrans and B notrans: "
       << ", A: " << "("<<a.extent(0)<<","<<a.extent(1)<<")"
       << ", B: " << "("<<b.extent(0)<<","<<b.extent(1)<<")"
       << ", C: " << "("<<c.extent(0)<<","<<c.extent(1)<<")";
    Kokkos::Impl::throw_runtime_exception (os.str ());
}
}
else if(rank==3){
if(a.extent(1)!=c.extent(1)||a.extent(2)!=b.extent(2)||b.extent(1)!=c.extent(2)||a.extent(0)!=b.extent(0)||b.extent(0)!=c.extent(0)){
    std::ostringstream os;
    os << "KokkosBlas::gemm: Matrix Dimensions Dimensions do not match for A notrans and B notrans: "
       << ", A: " << "("<<a.extent(0)<<","<<a.extent(1)<<","<<a.extent(2)<<")"
       << ", B: " << "("<<b.extent(0)<<","<<b.extent(1)<<","<<b.extent(2)<<")"
       << ", C: " << "("<<c.extent(0)<<","<<c.extent(1)<<","<<c.extent(2)<<")";
    Kokkos::Impl::throw_runtime_exception (os.str ());
}
}


}else if ((transA=='t'||transA=='T')&&(transB=='n'||transB=='N')) {
if(rank==2){
if(a.extent(1)!=c.extent(0)||a.extent(0)!=b.extent(0)||b.extent(1)!=c.extent(1)){
    std::ostringstream os;
    os << "KokkosBlas::gemm: Matrix Dimensions Dimensions do not match for A notrans and B notrans: "
       << ", A: " << "("<<a.extent(0)<<","<<a.extent(1)<<")"
       << ", B: " << "("<<b.extent(0)<<","<<b.extent(1)<<")"
       << ", C: " << "("<<c.extent(0)<<","<<c.extent(1)<<")";
    Kokkos::Impl::throw_runtime_exception (os.str ());
}
}
else if(rank==3){
if(a.extent(2)!=c.extent(1)||a.extent(1)!=b.extent(1)||b.extent(2)!=c.extent(2)||a.extent(0)!=b.extent(0)||b.extent(0)!=c.extent(0)){
    std::ostringstream os;
    os << "KokkosBlas::gemm: Matrix Dimensions Dimensions do not match for A notrans and B notrans: "
       << ", A: " << "("<<a.extent(0)<<","<<a.extent(1)<<","<<a.extent(2)<<")"
       << ", B: " << "("<<b.extent(0)<<","<<b.extent(1)<<","<<b.extent(2)<<")"
       << ", C: " << "("<<c.extent(0)<<","<<c.extent(1)<<","<<c.extent(2)<<")";
    Kokkos::Impl::throw_runtime_exception (os.str ());
}
}

}else if ((transA=='t'||transA=='T')&&(transB=='t'||transB=='T')) {
if(rank==2){
if(a.extent(1)!=c.extent(0)||a.extent(0)!=b.extent(1)||b.extent(0)!=c.extent(1)){
    std::ostringstream os;
    os << "KokkosBlas::gemm: Matrix Dimensions Dimensions do not match for A notrans and B notrans: "
       << ", A: " << "("<<a.extent(0)<<","<<a.extent(1)<<")"
       << ", B: " << "("<<b.extent(0)<<","<<b.extent(1)<<")"
       << ", C: " << "("<<c.extent(0)<<","<<c.extent(1)<<")";
    Kokkos::Impl::throw_runtime_exception (os.str ());
}
}
else if(rank==3){
if(a.extent(2)!=c.extent(1)||a.extent(1)!=b.extent(2)||b.extent(1)!=c.extent(2)||a.extent(0)!=b.extent(0)||b.extent(0)!=c.extent(0)){
    std::ostringstream os;
    os << "KokkosBlas::gemm: Matrix Dimensions Dimensions do not match for A notrans and B notrans: "
       << ", A: " << "("<<a.extent(0)<<","<<a.extent(1)<<","<<a.extent(2)<<")"
       << ", B: " << "("<<b.extent(0)<<","<<b.extent(1)<<","<<b.extent(2)<<")"
       << ", C: " << "("<<c.extent(0)<<","<<c.extent(1)<<","<<c.extent(2)<<")";
    Kokkos::Impl::throw_runtime_exception (os.str ());
}
}

}else{
    std::ostringstream os;
    os << "KokkosBlas::gemm: values for transA or transB should be T t or N n you have input: "
    <<"TransA:"<<transA<<" TransB:"<<transB;
    Kokkos::Impl::throw_runtime_exception (os.str ());

}
  


  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      AMat::rank == 2,
      typename AMat::const_value_type**,
      typename AMat::const_value_type*** >::type,
    typename AMat::array_layout,
    typename AMat::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename AMat::specialize> AMat_Internal;

  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      BMat::rank == 2,
      typename BMat::const_value_type**,
      typename BMat::const_value_type*** >::type,
    typename BMat::array_layout,
    typename BMat::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename BMat::specialize> BMat_Internal;

  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      CMat::rank == 2,
      typename CMat::const_value_type**,
      typename CMat::const_value_type*** >::type,
    typename CMat::array_layout,
    typename CMat::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename CMat::specialize> CMat_Internal;



  AMat_Internal a_i = a;
  BMat_Internal b_i = b;
  CMat_Internal c_i = c;

Impl::MultiGemm<AMat::non_const_value_type,
                BMat::non_const_value_type,
                CMat::non_const_value_type,
                CMat::execution_space,
                AMat::array_layout,
                BMat::array_layout,
                CMat::array_layout,
                CMat::size_type,
                CMat::rank>::GEMM(transA, transB, alpha, A, B, beta, C);

}



} // namespace KokkosBlas

#endif // KOKKOS_BLAS3_HPP_
