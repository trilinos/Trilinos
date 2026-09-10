// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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


using Teuchos::RCP;
using Teuchos::RCP_UNDEFINED_WEAK_NO_DEALLOC;
using Teuchos::RCP_WEAK_NO_DEALLOC;

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
  // Construct a null non-owning RCP directly from the pointer.  Do not
  // dereference ut_ptr to call rcpFromUndefRef()/rcpFromRef() since binding a
  // reference to a null pointer is undefined behavior.
  RCP<UndefinedType> ut_rcp =
#if defined(HAS_TEUCHOS_GET_BASE_OBJ_VOID_PTR)
    RCP<UndefinedType>(ut_ptr, RCP_UNDEFINED_WEAK_NO_DEALLOC)
#else
    RCP<UndefinedType>(ut_ptr, RCP_WEAK_NO_DEALLOC)
#endif
    ;
}


} // namespace
