// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
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
