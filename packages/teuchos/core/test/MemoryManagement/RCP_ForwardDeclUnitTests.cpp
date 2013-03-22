/*
// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// @HEADER
*/

#include "Teuchos_RCP.hpp"
#include "Teuchos_UnitTestHarness.hpp"


/*
 * This test checks that you can use non-owning Teuchos::RCP with pointers to
 * types that are only forward declared and not defined.
 */

namespace DummyNS {class UndefinedType;}

namespace Teuchos {
TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(DummyNS::UndefinedType);
} // namespace Teuchos


namespace {


using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::rcpFromUndefRef;
using Teuchos::RCP;

using DummyNS::UndefinedType;


TEUCHOS_UNIT_TEST( RCP, ForwardDeclaredUndefined )
{
  // This test ensures that you can declare a null RCP object to an undefined
  // type without trouble.
  RCP<UndefinedType> ut_rcp;
}


TEUCHOS_UNIT_TEST( RCP, ForwardDeclaredUndefined_rcp )
{
  // This test ensures that you can set a pointer to an undefined type without
  // trouble.  Note that this has to be a non-owning RCP otherwise there will
  // be issues with the destructor call.
  UndefinedType *ut_ptr = 0;
  RCP<UndefinedType> ut_rcp =
#if defined(HAS_TEUCHOS_GET_BASE_OBJ_VOID_PTR) 
    rcpFromUndefRef(*ut_ptr)
  // In this case, you have to use rcpFromUndefRef(...) in this case instead
  // of rcpFromRef() because the latter requires the object to be defined in
  // order to call dynamic_cast<const void*>(...) in order to get the base
  // object address needed for RCPNode tracing.
#else
    rcpFromRef(*ut_ptr)
    // In this case, you can use rcpFromRef(...) because the object's baseq
    // address will not be looked up using dynamic_cast and no deallocator
    // needing to know the object's will be compiled.
#endif
    ;
}


} // namespace
