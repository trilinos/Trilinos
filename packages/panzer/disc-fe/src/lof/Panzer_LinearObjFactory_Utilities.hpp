#ifndef __Panzer_LinearObjFactory_Utilities_hpp__
#define __Panzer_LinearObjFactory_Utilities_hpp__

#include "Teuchos_RCP.hpp"

#include "Panzer_Traits.hpp"

namespace panzer {

// forward declaration
template <typename> class LinearObjFactory;
class UniqueGlobalIndexerBase;

/** \brief Clone a linear object factory, but using a different domain. 
  *
  * This is the case where you want to make sure the linear algebra
  * abstractions are compatible, but don't really care which one is used.
  *
  * \param [in] lof Linear object factory to be cloned
  * \param [in] dUgi Domain unique global indexer to be used
  *
  * \note As implemented this only functions for the expicitly instantiated linear
  *       object factory types.
  */ 
Teuchos::RCP<const LinearObjFactory<panzer::Traits> > cloneWithNewDomain(const LinearObjFactory<panzer::Traits> & lof,
                                                                         const Teuchos::RCP<const UniqueGlobalIndexerBase> & dUgi);

/** \brief Clone a linear object factory, but using a different range. 
  *
  * This is the case where you want to make sure the linear algebra
  * abstractions are compatible, but don't really care which one is used.
  *
  * \param [in] lof Linear object factory to be cloned
  * \param [in] rUgi Range unique global indexer to be used
  *
  * \note As implemented this only functions for the expicitly instantiated linear
  *       object factory types.
  */ 
Teuchos::RCP<const LinearObjFactory<panzer::Traits> > cloneWithNewRange(const LinearObjFactory<panzer::Traits> & lof,
                                                                        const Teuchos::RCP<const UniqueGlobalIndexerBase> & rUgi);

/** \brief Clone a linear object factory, but using a different range and domain. 
  *
  * This is the case where you want to make sure the linear algebra
  * abstractions are compatible, but don't really care which one is used.
  *
  * \param [in] lof Linear object factory to be cloned
  * \param [in] rUgi Range unique global indexer to be used
  * \param [in] dUgi Domain unique global indexer to be used
  *
  * \note As implemented this only functions for the expicitly instantiated linear
  *       object factory types.
  */ 
Teuchos::RCP<const LinearObjFactory<panzer::Traits> > cloneWithNewRangeAndDomain(
                                                                        const LinearObjFactory<panzer::Traits> & lof,
                                                                        const Teuchos::RCP<const UniqueGlobalIndexerBase> & rUgi,
                                                                        const Teuchos::RCP<const UniqueGlobalIndexerBase> & dUgi);
}

#endif
