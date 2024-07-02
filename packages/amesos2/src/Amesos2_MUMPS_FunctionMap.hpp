// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
   \file   Amesos2_MUMPS_FunctionMap.hpp
   \author Joshua Dennis Booth <jdbooth@sandia.gov>

   \brief  Provides a mechanism to map function calls to the correct Solver
           function based on the scalar type of Matrices and MultiVectors
*/

#ifndef AMESOS2_MUMPS_FUNCTIONMAP_HPP
#define AMESOS2_MUMPS_FUNCTIONMAP_HPP

#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include "Amesos2_FunctionMap.hpp"
#include "Amesos2_MUMPS_TypeMap.hpp"


namespace Amesos2
{
  /*Specializtions*/

  template <>
  struct FunctionMap<MUMPS, float>
  {
    typedef TypeMap<MUMPS, float> typemap;

    static void mumps_c(typemap::MUMPS_STRUC_C *mumps_par)
    {
      MUMPST::smumps_c(mumps_par);
    }   
  };
  template <>
  struct FunctionMap<MUMPS,double>
  {
    typedef TypeMap<MUMPS, double> typemap;
    
    static void mumps_c(typemap::MUMPS_STRUC_C *mumps_par)
    {
      MUMPST::dmumps_c(mumps_par);
    }
  };

#ifdef HAVE_TEUCHOS_COMPLEX
  
  template <>
  struct FunctionMap<MUMPS,MUMPST::CMUMPS_COMPLEX>
  {
    typedef TypeMap<MUMPS,MUMPST::CMUMPS_COMPLEX> typemap;
    
    static void mumps_c(typemap::MUMPS_STRUC_C *mumps_par)
    {
      MUMPST::cmumps_c(mumps_par);
    }
  };
 
  
  template <>
  struct FunctionMap<MUMPS, std::complex<float> >
  {
    typedef TypeMap<MUMPS, std::complex<float> > typemap;
    
    static void mumps_c(typemap::MUMPS_STRUC_C *mumps_par)
    {
      MUMPST::cmumps_c(mumps_par);
    }
  };
   
  template <>
  struct FunctionMap<MUMPS, std::complex<double>  >
  {
    typedef TypeMap<MUMPS, std::complex<double>  > typemap;
    
    static void mumps_c(typemap::MUMPS_STRUC_C *mumps_par)
    {
      MUMPST::zmumps_c(mumps_par);
    }
  };

#endif //complex
} //end namespace Amesos2

#endif  // AMESOS2_MUMPS_FUNCTIONMAP_HPP
