#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "EpetraThyraAdaptersTestHelpers.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace {


//
// Unit Tests
//


TEUCHOS_UNIT_TEST( get_Epetra_MultiVector, singleBlockProductVector )
{
  using Teuchos::Comm;
  typedef Teuchos_Ordinal Ordinal;
  using Thyra::VectorSpaceBase;
  using Thyra::MultiVectorBase;

  const RCP<const Epetra_Comm> epetra_comm = getEpetraComm();
  const RCP<const Comm<Ordinal> > comm = Thyra::create_Comm(epetra_comm);
  
  const RCP<const Epetra_Map> epetra_map = rcp(new Epetra_Map(g_localDim, 0, *epetra_comm));
  const RCP<const VectorSpaceBase<double> > vs =  Thyra::create_VectorSpace(epetra_map);

  const RCP<const VectorSpaceBase<double> > pvs = Thyra::productVectorSpace(vs, 1);

  const RCP<MultiVectorBase<double> > pmv = Thyra::createMembers(pvs, 1);

  const double alpha = 3.5;
  Thyra::assign<double>( pmv.ptr(), alpha );

  const RCP<Epetra_MultiVector> epetra_mv =
    Thyra::get_Epetra_MultiVector(*epetra_map, pmv);

  const RCP<MultiVectorBase<double> > mv2 =
    Thyra::create_MultiVector(epetra_mv, pvs);

  Thyra::testRelNormDiffErr<double>(
    "*pmv->col(0)", *pmv->col(0),
    "*mv2->col(0)", *mv2->col(0),
    "max-error", 0.0,
    "max-warning", 0.0,
    &out
    );
   
}


} // namespace
