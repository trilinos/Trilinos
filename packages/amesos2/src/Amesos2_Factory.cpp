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

#include "Amesos2_Factory.hpp"

namespace Amesos2 {

  /**********************
   *   QUERY function   *
   **********************/

  bool query(const char* solverName){
    std::string solver = solverName;
    return( query(solver) );
  }

  bool query(const std::string solver_name){
    std::string solverName = tolower(solver_name); // for easier string checking
#ifdef HAVE_AMESOS2_KLU2
    if((solverName == "amesos2_klu2") || (solverName == "klu2") ||
       (solverName == "amesos2_klu")  || (solverName == "klu")){
      return( true );
    }
#endif

#ifdef HAVE_AMESOS2_SUPERLUDIST
    if((solverName == "amesos2_superludist") ||
       (solverName == "superludist") ||
       (solverName == "amesos2_superlu_dist") ||
       (solverName == "superlu_dist")){
      return( true );
    }
#endif

#ifdef HAVE_AMESOS2_SUPERLUMT
    if((solverName == "amesos2_superlumt") ||
       (solverName == "superlumt") ||
       (solverName == "amesos2_superlu_mt") ||
       (solverName == "superlu_mt")){
      return( true );
    }
#endif

#ifdef HAVE_AMESOS2_SUPERLU
    if((solverName == "amesos2_superlu") ||
       (solverName == "superlu")){
      return( true );
    }
#endif

#ifdef HAVE_AMESOS2_PARDISO_MKL
    if((solverName == "amesos2_pardiso_mkl") ||
       (solverName == "pardiso_mkl") ||
       (solverName == "amesos2_pardisomkl")  ||
       (solverName == "pardisomkl")){
      return( true );
    }
#endif

#ifdef HAVE_AMESOS2_LAPACK
    if((solverName == "amesos2_lapack") ||
       (solverName == "lapack")){
      return( true );
    }
#endif

    // Otherwise, the solver is not available
    return( false );
  }


  std::string tolower(const std::string& s)
  {
    std::locale loc;
    std::string rtn = s;
    for (size_t i=0; i<rtn.length(); ++i)
      {
        rtn[i] = tolower(rtn[i],loc);
      }
    return rtn;
  }
} // end namespace Amesos2
