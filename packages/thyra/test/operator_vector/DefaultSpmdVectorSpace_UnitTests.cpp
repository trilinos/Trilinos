
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"


namespace {


//
// Helper code and declarations
//


using Teuchos::as;
using Teuchos::null;
using Teuchos::RCP;
using Thyra::VectorBase;
using Thyra::VectorSpaceBase;
using Thyra::createMember;
using Thyra::DefaultSpmdVectorSpace;
using Thyra::defaultSpmdVectorSpace;
using Thyra::ConstDetachedVectorView;
using Thyra::DetachedVectorView;
using Thyra::ConstDetachedSpmdVectorView;
using Thyra::DetachedSpmdVectorView;
typedef Thyra::Ordinal Ordinal;


const int g_localDim = 4; // ToDo: Make variable!


template<class Scalar>
RCP<VectorSpaceBase<Scalar> > 
createSpmdVectorSpace(const Teuchos_Ordinal localDim)
{
  return defaultSpmdVectorSpace<Scalar>(
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm(),
    localDim, -1 );
}


//
// Unit Tests
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace, defaultConstruct,
  Scalar )
{

  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    Teuchos::rcp(new Thyra::DefaultSpmdVectorSpace<Scalar>));
  TEST_EQUALITY(vs->getComm(), null);
  TEST_EQUALITY(vs->localOffset(), as<Ordinal>(-1));
  TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(-1));
  TEST_EQUALITY(vs->mapCode(), as<Ordinal>(-1));
  TEST_EQUALITY(vs->dim(), as<Ordinal>(-1));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace,
  defaultConstruct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace, serialConstruct,
  Scalar )
{

  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    defaultSpmdVectorSpace<Scalar>(g_localDim));
  TEST_EQUALITY(vs->getComm(), null);
  TEST_EQUALITY(vs->localOffset(), as<Ordinal>(0));
  TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(g_localDim));
  TEST_EQUALITY(vs->mapCode(), as<Ordinal>(g_localDim));
  TEST_EQUALITY(vs->dim(), as<Ordinal>(g_localDim));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace,
  serialConstruct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace, parallelConstruct,
  Scalar )
{
  ECHO(const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm());
  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    defaultSpmdVectorSpace<Scalar>(comm, g_localDim, -1));
  TEST_EQUALITY(vs->getComm(), comm);
  TEST_EQUALITY(vs->localOffset(), as<Ordinal>(comm->getRank()*g_localDim));
  TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(g_localDim));
  TEST_EQUALITY(vs->dim(), as<Ordinal>(comm->getSize()*g_localDim));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace,
  parallelConstruct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace, parallelConstructGlobalDim,
  Scalar )
{
  ECHO(const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm());
  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    defaultSpmdVectorSpace<Scalar>(comm, g_localDim, g_localDim * comm->getSize()));
  TEST_EQUALITY(vs->getComm(), comm);
  TEST_EQUALITY(vs->localOffset(), as<Ordinal>(comm->getRank()*g_localDim));
  TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(g_localDim));
  TEST_EQUALITY(vs->dim(), as<Ordinal>(comm->getSize()*g_localDim));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace,
  parallelConstructGlobalDim )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace, locallyReplicatedParallelConstruct,
  Scalar )
{
  ECHO(const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm());
  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    defaultSpmdVectorSpace<Scalar>(comm, g_localDim, g_localDim));
  TEST_EQUALITY(vs->getComm(), comm);
  TEST_EQUALITY(vs->localOffset(), as<Ordinal>(0));
  TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(g_localDim));
  TEST_EQUALITY(vs->dim(), as<Ordinal>(g_localDim));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace,
  locallyReplicatedParallelConstruct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace, parallelFullExtract,
  Scalar )
{

  const RCP<const VectorSpaceBase<Scalar> > vs =
    createSpmdVectorSpace<Scalar>(g_localDim);

  const RCP<VectorBase<Scalar> > v = createMember(vs);
  {
    out << "\nSetting up v[i] = i, i=0...n-1 ...\n";
    const DetachedSpmdVectorView<Scalar> dv(v);
    const Ordinal localOffset = dv.spmdSpace()->localOffset();
    const Ordinal localSubDim = dv.spmdSpace()->localSubDim();
    for (Ordinal i = 0; i < localSubDim; ++i) {
      dv[i] = as<Scalar>(localOffset + i);
    }
  }

  out << "\nv = " << *v;

  {

    const ConstDetachedVectorView<Scalar> dv(v);

    TEST_EQUALITY(dv.subDim(), vs->dim());
    
    out << "\nTest that dv[i] == i, i=0...n-1 ... ";
    bool local_success = true;
    for (Ordinal i = 0; i < dv.subDim(); ++i) {
      TEST_ARRAY_ELE_EQUALITY( dv, i, as<Scalar>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;

  }
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace,
  parallelFullExtract)


} // namespace
