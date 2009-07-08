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
  RCP<UndefinedType> ut_rcp = rcpFromRef(*ut_ptr);
  // NOTE: You have to use rcpFroRef(...) and not rcp(..., false) because
  // rcp(..., false) will use the DeallocDelete class that has a delete call
  // that the compiler will complain about.
}


} // namespace
