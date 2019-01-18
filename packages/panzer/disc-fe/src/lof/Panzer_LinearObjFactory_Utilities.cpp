#include "Panzer_LinearObjFactory_Utilities.hpp"

#include "Teuchos_RCP.hpp"

#include "Panzer_UniqueGlobalIndexer.hpp"

#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_BlockedTpetraLinearObjFactory.hpp"

namespace panzer {

Teuchos::RCP<const LinearObjFactory<panzer::Traits> > cloneWithNewDomain(const LinearObjFactory<panzer::Traits> & lof,
                                                                         const Teuchos::RCP<const UniqueGlobalIndexerBase> & dUgi)
{
  // This just forwards on to the general case. That makes things much easier
  return cloneWithNewRangeAndDomain(lof,Teuchos::null,dUgi);
}

Teuchos::RCP<const LinearObjFactory<panzer::Traits> > cloneWithNewRange(const LinearObjFactory<panzer::Traits> & lof,
                                                                        const Teuchos::RCP<const UniqueGlobalIndexerBase> & rUgi)
{
  // This just forwards on to the general case. That makes things much easier
  return cloneWithNewRangeAndDomain(lof,rUgi,Teuchos::null);
}

Teuchos::RCP<const LinearObjFactory<panzer::Traits> > cloneWithNewRangeAndDomain(
                                                                        const LinearObjFactory<panzer::Traits> & lof,
                                                                        const Teuchos::RCP<const UniqueGlobalIndexerBase> & rUgi,
                                                                        const Teuchos::RCP<const UniqueGlobalIndexerBase> & dUgi)
{
  using Teuchos::null;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::Ptr;
  using Teuchos::ptr_dynamic_cast;
  using Teuchos::ptrFromRef;

  using TpetraUGI = UniqueGlobalIndexer<int, Ordinal64>;
/*
  typedef UniqueGlobalIndexer<int,int>       EpetraUGI;
  typedef BlockedDOFManager<int,int>         BlockedEpetraUGI;
  typedef BlockedDOFManager<int,Ordinal64>   BlockedTpetraUGI;
*/
  typedef TpetraLinearObjFactory<panzer::Traits,double,int,Ordinal64>        TpetraLOF;
  typedef BlockedEpetraLinearObjFactory<panzer::Traits,int>                  BlockedEpetraLOF;
  typedef BlockedTpetraLinearObjFactory<panzer::Traits,double,int,Ordinal64> BlockedTpetraLOF;

  // This proceeds by casting to a number of known LOF types (all explicitly instantiated)
  // then trying to build a new one. Of course for many of these under implemented operation
  // this fails and an error is thrown.
 
/*
  Ptr<const EpetraLOF> epetra_lof = ptr_dynamic_cast<const EpetraLOF>(ptrFromRef(lof));
  if(epetra_lof!=null) {
    RCP<const EpetraUGI> rangeUGI  = rcp_dynamic_cast<const EpetraUGI>(rUgi==null ? epetra_lof->getRangeGlobalIndexer() : rUgi,true);
    RCP<const EpetraUGI> domainUGI = rcp_dynamic_cast<const EpetraUGI>(dUgi==null ? epetra_lof->getDomainGlobalIndexer() : dUgi,true);
    RCP<Teuchos::MpiComm<int> > mpiComm = rcp(new Teuchos::MpiComm<int>(epetra_lof->getComm()));
    return rcp(new EpetraLOF(mpiComm,rangeUGI,domainUGI));
  }
*/

  Ptr<const TpetraLOF> tpetra_lof = ptr_dynamic_cast<const TpetraLOF>(ptrFromRef(lof));
  if(tpetra_lof!=null) {
    auto rangeUGI  = rcp_dynamic_cast<const TpetraUGI>(rUgi==null ? tpetra_lof->getRangeGlobalIndexer() : rUgi,true);
    auto domainUGI = rcp_dynamic_cast<const TpetraUGI>(dUgi==null ? tpetra_lof->getDomainGlobalIndexer() : dUgi,true);
    auto mpiComm = rcp(new Teuchos::MpiComm<int>(tpetra_lof->getComm()));

    return rcp(new TpetraLOF(mpiComm,rangeUGI,domainUGI));
  }

  Ptr<const BlockedEpetraLOF> blk_epetra_lof = ptr_dynamic_cast<const BlockedEpetraLOF>(ptrFromRef(lof));
  if(blk_epetra_lof!=null) {
    auto rangeUGI  = (rUgi==null ? blk_epetra_lof->getRangeGlobalIndexer() : rUgi);
    auto domainUGI = (dUgi==null ? blk_epetra_lof->getDomainGlobalIndexer() : dUgi);
    RCP<Teuchos::MpiComm<int> > mpiComm = rcp(new Teuchos::MpiComm<int>(blk_epetra_lof->getComm()));
    return rcp(new BlockedEpetraLOF(mpiComm,rangeUGI,domainUGI));
  }

  Ptr<const BlockedTpetraLOF> blk_tpetra_lof = ptr_dynamic_cast<const BlockedTpetraLOF>(ptrFromRef(lof));
  if(blk_tpetra_lof!=null) {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                               "panzer::cloneWithNewRangeAndDomain: Blocked Tpetra LOF does not yet support "
                               "different range and domain indexers!");
  }

  TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                             "panzer::cloneWithNewRangeAndDomain: Could not determine the type of LOF, clone not support!");

  return Teuchos::null;
}

}
