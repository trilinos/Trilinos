// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_VECTOR_HPP
#define SHYLUBASKER_VECTOR_HPP

#include "shylubasker_types.hpp"

#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#else
#include <omp.h>
#endif

namespace BaskerNS
{

  template <class Int, class Entry, class Exe_Space>
  struct basker_vector
  {
    template <class TypeOut, class TypeIn>
    BASKER_INLINE
    static
    void init_value(TypeOut a, Int size, TypeIn c)
    {
      //need to add some kind of fast copy
      for(Int i=0; i < size; i++)
      {
        a[i] = c;
      }
    }//end init_value

    template <class TypeOut, class TypeIn >
    BASKER_INLINE
    static
    void init_vector(TypeOut a, Int size, TypeIn c)
    {
      //need to do some fast copy
      for(Int i=0; i< size; i++)
      {
        a[i] = c[i];
      }
    }//end init_vector
  };

}//end BaskerNS

#endif
