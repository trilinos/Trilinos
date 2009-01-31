#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_getConst.hpp"
#include "TestClasses.hpp"


namespace {


using Teuchos::null;
using Teuchos::Ptr;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromPtr;
using Teuchos::NullReferenceError;
using Teuchos::DanglingReferenceError;
using Teuchos::RCP_STRONG;
using Teuchos::RCP_WEAK;
using Teuchos::RCP_STRENGTH_INVALID;


TEUCHOS_UNIT_TEST( Ptr, rcpFromPtr )
{
  ECHO(RCP<A> a_rcp = rcp(new A));
  ECHO(Ptr<A> a_ptr = a_rcp.ptr());
  ECHO(RCP<A> a_rcp2 = rcpFromPtr(a_ptr));
  TEST_EQUALITY(a_rcp2.getRawPtr(), a_rcp.getRawPtr());
#ifdef TEUCHOS_DEBUG
  TEST_ASSERT(a_rcp2.shares_resource(a_rcp));
#else
  // In an optimized build, the object a_rcp2 has its own RCPNode object that
  // is unrelated to the orgininal a_rcp object.  This cuts down on overhead.
#endif
  ECHO(a_rcp = null);
#ifdef TEUCHOS_DEBUG
  TEST_THROW(a_ptr.getRawPtr(), DanglingReferenceError);
  TEST_THROW(a_rcp2.getRawPtr(), DanglingReferenceError);
#endif  
}


} // namespace
