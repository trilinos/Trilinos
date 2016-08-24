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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_Utils.hpp
    \brief  Intrepid utilities.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_UTILS_HPP__
#define __INTREPID2_UTILS_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Types.hpp"

namespace Intrepid2 {

  ///
  /// test macros
  ///

#define INTREPID2_TEST_FOR_EXCEPTION(test, x, msg)                      \
  if (test) {                                                           \
    printf("[Intrepid2] Error in file %s, line %d\n",__FILE__,__LINE__); \
    printf("            Test that evaluated to true: %s\n", #test);     \
    printf("            %s \n", msg);                                   \
    throw x(msg);                                                       \
  }

#define INTREPID2_TEST_FOR_ABORT(test, msg)                             \
  if (test) {                                                           \
    printf("[Intrepid2] Error in file %s, line %d\n",__FILE__,__LINE__); \
    printf("            Test that evaluated to true: %s\n", #test);     \
    printf("            %s \n", msg);                                   \
    Kokkos::abort(  "[Intrepid2] Abort\n");                             \
  }

  // check the first error only
#ifdef INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#define INTREPID2_TEST_FOR_DEBUG_ABORT(test, info, msg)                 \
  if (!(info) && (test)) {                                              \
    printf("[Intrepid2] Error in file %s, line %d\n",__FILE__,__LINE__); \
    printf("            Test that evaluated to true: %s\n", #test);     \
    printf("            %s \n", msg);                                   \
    info = true;                                                        \
  }
#else  
#define INTREPID2_TEST_FOR_DEBUG_ABORT(test, info, msg)                 \
  if (!(info) && (test)) {                                              \
    printf("[Intrepid2] Error in file %s, line %d\n",__FILE__,__LINE__); \
    printf("            Test that evaluated to true: %s\n", #test);     \
    printf("            %s \n", msg);                                   \
    info = true ;                                                       \
    Kokkos::abort(  "[Intrepid2] Abort\n");                             \
  }
#endif

  /// 
  /// space overload
  /// 

  template<typename ViewSpaceType, typename UserSpaceType>
  struct ExecSpace {
    typedef UserSpaceType ExecSpaceType;
  };

  template<typename ViewSpaceType>
  struct ExecSpace<ViewSpaceType,void> {
    typedef ViewSpaceType ExecSpaceType;
  };

  ///
  /// utilities device compirable
  ///

  class Util {
  public:

    template<typename IdxType, typename DimType, typename IterType>
    KOKKOS_FORCEINLINE_FUNCTION
    static void 
    unrollIndex(IdxType &i, IdxType &j, 
                const DimType dim0,
                const DimType dim1,
                const IterType iter) {
      // left index
      //j = iter/dim0;
      //i = iter%dim0;

      // right index
      i = iter/dim1;
      j = iter%dim1;
    }
    
    template<typename IdxType, typename DimType, typename IterType>
    KOKKOS_FORCEINLINE_FUNCTION
    static void 
    unrollIndex(IdxType &i, IdxType &j, IdxType &k, 
                const DimType dim0,
                const DimType dim1,
                const DimType dim2,
                const IterType iter) {
      IdxType tmp;

      //unrollIndex(tmp, k, dim0*dim1, dim2, iter);
      //unrollIndex(  i, j, dim0,      dim1,  tmp);

      unrollIndex( i, tmp, dim0, dim1*dim2, iter);
      unrollIndex( j, k,   dim1,      dim2,  tmp);
    }

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static T sign(const T a) {
      return (a < 0 ? -1 : 1);
    }

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static T min(const T a, const T b) {
      return (a < b ? a : b);
    }

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static T max(const T a, const T b) {
      return (a > b ? a : b);
    }

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static T abs(const T a) {
      return (a > 0 ? a : -a);
    }

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static T real(const T a) {
      return a;
    }

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static T imag(const T a) {
      return 0;
    }

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static T conj(const T a) {
      return a;
    }
    
    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static void swap(T &a, T &b) {
      T c(a); a = b; b = c;
    }
  };


} // end namespace Intrepid2

#endif
