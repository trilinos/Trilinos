
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"

//#define THYRA_DEFAULT_SPMD_VECTOR_SPACE_UNIT_TESTS_DUMP

#ifdef THYRA_DEFAULT_SPMD_VECTOR_SPACE_UNIT_TESTS_DUMP
#  include "RTOpPack_SPMD_apply_op_decl.hpp"
#  include "Thyra_SpmdVectorBase.hpp"
#endif // THYRA_DEFAULT_SPMD_VECTOR_SPACE_UNIT_TESTS_DUMP


namespace {


//
// Helper code and declarations
//


using Teuchos::as;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::get_extra_data;
using Thyra::VectorSpaceBase;
using Thyra::VectorBase;
using Thyra::MultiVectorBase;
using Thyra::createMember;
using Thyra::createMembers;
using Thyra::DefaultSpmdVectorSpace;
using Thyra::defaultSpmdVectorSpace;
using Thyra::ConstDetachedVectorView;
using Thyra::DetachedVectorView;
using Thyra::ConstDetachedSpmdVectorView;
using Thyra::DetachedSpmdVectorView;
typedef Thyra::Ordinal Ordinal;


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace_Parallel, empty_proc,
  Scalar )
{
  ECHO(const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm());

  Thyra::Ordinal localSize = comm->getRank()==0 ? 0 : 5;
  RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    Thyra::defaultSpmdVectorSpace<Scalar>(comm,localSize,5);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace_Parallel,
  empty_proc )


} // namespace
