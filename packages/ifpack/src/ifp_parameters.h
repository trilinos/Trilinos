/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef _ifp_parameters_h_
#define _ifp_parameters_h_

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "Ifpack_config.h"

#include <Ifpack_ConfigDefs.h>

#include <Teuchos_map.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_CombineMode.h>

namespace Ifpack {

//define enum values to which parameter names will be mapped.
enum parameter {
  //parameters of type double
  absolute_threshold,
  relative_threshold,
  drop_tolerance,
  fill_tolerance,
  relax_value,

  //parameters of type int
  //(if you add or remove int parameters, be sure to
  //update FIRST_INT_PARAM and LAST_INT_PARAM macros below, as
  //they are used below and in ifp_parameters.cpp)
  level_fill,
  level_overlap,
  num_steps,

  //mixed type parameters
  use_reciprocal,
  overlap_mode
};

#define FIRST_INT_PARAM Ifpack::level_fill
#define LAST_INT_PARAM Ifpack::num_steps

//define struct with union of all Ifpack parameters
struct param_struct {
  int int_params[LAST_INT_PARAM-FIRST_INT_PARAM+1];
  double double_params[FIRST_INT_PARAM];
  bool use_reciprocal;
  Epetra_CombineMode overlap_mode;
};

Teuchos::map<std::string,parameter>& key_map();

// This was marked as IFPACK_DEPRECATED for a very long time,
// but the function never went away and was still used in Ifpack.
// Please consider it still marked as deprecated.
void initialize_string_map();

std::string upper_case(const std::string& s);

// This was marked as IFPACK_DEPRECATED for a very long time,
// but the function never went away and was still used in Ifpack.
// Please consider it still marked as deprecated.
void set_parameters(const Teuchos::ParameterList& parameterlist,
                    param_struct& params,
                    bool cerr_warning_if_unused=false);

}//namespace Ifpack

#endif //_ifp_parameters_h_
