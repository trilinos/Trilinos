/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#ifndef _LibraryFactory_hpp_
#define _LibraryFactory_hpp_

#include <fei_macros.hpp>
#include <fei_SharedPtr.hpp>
#include <fei_LibraryWrapper.hpp>

#include <fei_Factory.hpp>

namespace fei {

  /** Create an instance of LibraryWrapper. Throws std::runtime_error
      if the input libraryName is not recognized. (names are case-sensitive)

      @param libraryName Input name of the solver library that is to be used.
      Valid values of this parameter are:<br>
      <ul>
      <li>Aztec
      <li>FETI
      </ul>

      @return shared-pointer holding newly-created LibraryWrapper instance.
   */
  fei::SharedPtr<LibraryWrapper> create_LibraryWrapper(MPI_Comm comm,
						       const char* libraryName);

  /** Create an instance of the fei::Factory interface. Throws std::runtime_error
      if the input libraryName is not recognized. (names are case-sensitive)

      @param libraryName Input name of solver library, same valid values as for
      'create_LibraryWrapper' above, as well as allowing "Trilinos".

      @return shared-pointer holding newly-created fei::Factory instance.
  */
  fei::SharedPtr<fei::Factory> create_fei_Factory(MPI_Comm comm,
						  const char* libraryName);
}//namespace fei

#endif // _LibraryFactory_hpp_

