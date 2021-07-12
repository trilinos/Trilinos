// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid2_ArgExtractor.hpp
    \brief  Header file with various static argument-extractor classes.  These are useful for writing efficient, templated code in terms of a subset of the arguments passed into a specified functor.  See Intrepid2::Data, and specifically its storeInPlaceCombination() implementation, for an example.
    \author Created by Nate Roberts.
*/

#ifndef __Intrepid2_ArgExtractor_HPP__
#define __Intrepid2_ArgExtractor_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_DeviceAssert.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Kokkos_Core.hpp"

namespace Intrepid2 {
  /** \class  Intrepid2::ConstantArgExtractor
      \brief Argument extractor class which ignores the input arguments in favor of passing a single 0 argument to the provided container.
  */
  template<class reference_type>
  struct ConstantArgExtractor
  {
    template<class ViewType, class ...IntArgs>
    static KOKKOS_INLINE_FUNCTION reference_type get(const ViewType &view, const IntArgs&... intArgs)
    {
      return view(0);
    }
  };

  /** \class  Intrepid2::FullArgExtractor
      \brief Argument extractor class which passes all arguments to the provided container.
  */
  template<class reference_type>
  struct FullArgExtractor
  {
    template<class ViewType, class ...IntArgs>
    static KOKKOS_INLINE_FUNCTION reference_type get(const ViewType &view, const IntArgs&... intArgs)
    {
      return view(intArgs...);
    }
  };

  /** \class  Intrepid2::SingleArgExtractor
      \brief Argument extractor class which passes a single argument, indicated by the template parameter whichArg, to the provided container.
  */
  template<class reference_type, int whichArg>
  struct SingleArgExtractor
  {
    template< bool B, class T = reference_type >
    using enable_if_t = typename std::enable_if<B,T>::type;
    
    template<class ViewType, class int_type, int M=whichArg>
    static KOKKOS_INLINE_FUNCTION
    enable_if_t<M == 0>
    get(const ViewType &view, const int_type &i0)
    {
      return view(i0);
    }
    
    template<class ViewType, class int_type, class ...IntArgs, int M=whichArg>
    static KOKKOS_INLINE_FUNCTION
    enable_if_t<M == 0>
    get(const ViewType &view, const int_type &i0, const IntArgs&... intArgs)
    {
      return view(i0);
    }
    
    template<class ViewType, class int_type, int M=whichArg>
    static KOKKOS_INLINE_FUNCTION
    enable_if_t<M == 1>
    get(const ViewType &view, const int_type &i0, const int_type &i1)
    {
      return view(i1);
    }
    
    template<class ViewType, class int_type, class ...IntArgs, int M=whichArg>
    static KOKKOS_INLINE_FUNCTION
    enable_if_t<M == 1>
    get(const ViewType &view, const int_type &i0, const int_type &i1, const IntArgs&... intArgs)
    {
      return view(i1);
    }
    
    template<class ViewType, class int_type, int M=whichArg>
    static KOKKOS_INLINE_FUNCTION
    enable_if_t<M == 2>
    get(const ViewType &view, const int_type &i0, const int_type &i1, const int_type &i2)
    {
      return view(i2);
    }
    
    template<class ViewType, class int_type, class ...IntArgs, int M=whichArg>
    static KOKKOS_INLINE_FUNCTION
    enable_if_t<M == 2>
    get(const ViewType &view, const int_type &i0, const int_type &i1, const int_type &i2, const IntArgs&... intArgs)
    {
      return view(i2);
    }
    
    template<class ViewType, class int_type, int M=whichArg>
    static KOKKOS_INLINE_FUNCTION
    enable_if_t<M == 3>
    get(const ViewType &view, const int_type &i0, const int_type &i1, const int_type &i2, const int_type &i3)
    {
      return view(i3);
    }
    
    template<class ViewType, class int_type, class ...IntArgs, int M=whichArg>
    static KOKKOS_INLINE_FUNCTION
    enable_if_t<M == 3>
    get(const ViewType &view, const int_type &i0, const int_type &i1, const int_type &i2, const int_type &i3, const IntArgs&... intArgs)
    {
      return view(i3);
    }
    
    template<class ViewType, class int_type, int M=whichArg>
    static KOKKOS_INLINE_FUNCTION
    enable_if_t<M == 4>
    get(const ViewType &view, const int_type &i0, const int_type &i1, const int_type &i2, const int_type &i3, const int_type &i4)
    {
      return view(i4);
    }
    
    template<class ViewType, class int_type, class ...IntArgs, int M=whichArg>
    static KOKKOS_INLINE_FUNCTION
    enable_if_t<M == 4>
    get(const ViewType &view, const int_type &i0, const int_type &i1, const int_type &i2, const int_type &i3, const int_type &i4, const IntArgs&... intArgs)
    {
      return view(i4);
    }
    
    template<class ViewType, class int_type, int M=whichArg>
    static KOKKOS_INLINE_FUNCTION
    enable_if_t<M == 5>
    get(const ViewType &view, const int_type &i0, const int_type &i1, const int_type &i2, const int_type &i3, const int_type &i4, const int_type &i5)
    {
      return view(i5);
    }
    
    template<class ViewType, class int_type, class ...IntArgs, int M=whichArg>
    static KOKKOS_INLINE_FUNCTION
    enable_if_t<M == 5>
    get(const ViewType &view, const int_type &i0, const int_type &i1, const int_type &i2, const int_type &i3, const int_type &i4, const int_type &i5, const IntArgs&... intArgs)
    {
      return view(i5);
    }
    
    // the commented-out code below is a cleaner way to implement the above, but we can't support this on CUDA until we can require KOKKOS_ENABLE_CUDA_CONSTEXPR
    /*
    template<class ViewType, class ...IntArgs>
    static KOKKOS_INLINE_FUNCTION
    enable_if_t<whichArg < sizeof...(IntArgs), reference_type>
    get(const ViewType &view, const IntArgs&... intArgs)
    {
      const auto & arg = std::get<whichArg>(std::tuple<IntArgs...>(intArgs...));
      return view(arg);
    }
     */
    
    template<class ViewType, class ...IntArgs>
    static KOKKOS_INLINE_FUNCTION
    enable_if_t<whichArg >= sizeof...(IntArgs), reference_type>
    get(const ViewType &view, const IntArgs&... intArgs)
    {
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true,std::invalid_argument,"calling SingleArgExtractor with out-of-bounds argument");
      Kokkos::abort("Intrepid2::SingleArgExtractor: calling SingleArgExtractor with out-of-bounds argument\n");
      return view(0); // this line added to avoid missing return statement warning under nvcc
    }
  };
}
#endif
