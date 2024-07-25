// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
