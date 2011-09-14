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

template<class Scalar , class DeviceType >
struct WAXPBY;

template<class Scalar >
struct WAXPBY<Scalar , KOKKOS_MACRO_DEVICE >
{
  
  typedef KOKKOS_MACRO_DEVICE               device_type;
  typedef device_type::size_type              size_type;
  typedef Kokkos::MultiVectorView<Scalar, device_type>   scalar_vector;  
  typedef Kokkos::ValueView<Scalar , device_type>      value;


  value alpha;
  scalar_vector x;
  
  value beta;
  scalar_vector y;
  
  scalar_vector w;
  
  WAXPBY(value arg_alpha ,  scalar_vector & arg_x , value arg_beta , scalar_vector & arg_y, 
       scalar_vector & arg_w )
    :  alpha(arg_alpha) , x(arg_x), beta(arg_beta) , y(arg_y) , w(arg_w) { }
    
    
    
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()(int inode) const 
  {
    w(inode) = (*alpha)*x(inode) + (*beta)*y(inode);
  }

}; //WAXPBY

