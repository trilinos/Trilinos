#include "Thyra_VectorBase.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"

#include "Panzer_GlobalEvaluationData.hpp"
#include "Panzer_BlockedVector_ReadOnly_GlobalEvaluationData.hpp"

namespace panzer {


BlockedVector_ReadOnly_GlobalEvaluationData::
BlockedVector_ReadOnly_GlobalEvaluationData() 
  : isInitialized_(false) 
{
}

BlockedVector_ReadOnly_GlobalEvaluationData::
BlockedVector_ReadOnly_GlobalEvaluationData(const BlockedVector_ReadOnly_GlobalEvaluationData & src) 
  : isInitialized_(false) 
{ 
  initialize(src.ghostedSpace_,Teuchos::null,src.gedBlocks_); 
}

BlockedVector_ReadOnly_GlobalEvaluationData::
BlockedVector_ReadOnly_GlobalEvaluationData(const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > ghostedSpace,
                                            const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > uniqueSpace,
                                            const std::vector<Teuchos::RCP<ReadOnlyVector_GlobalEvaluationData> > & gedBlocks)
  : isInitialized_(false) 
{ 
  initialize(ghostedSpace,uniqueSpace,gedBlocks); 
}

void 
BlockedVector_ReadOnly_GlobalEvaluationData::
initialize(const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & ghostedSpace,
           const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & uniqueSpace,
           const std::vector<Teuchos::RCP<ReadOnlyVector_GlobalEvaluationData> > & gedBlocks)
{
  using Teuchos::rcp_dynamic_cast;

  // assert all gedBlocks are initialized
  for(std::size_t i=0;i<gedBlocks.size();i++) {
    TEUCHOS_TEST_FOR_EXCEPTION(!gedBlocks[i]->isInitialized(),std::logic_error,
                               "BlockedVector_ReadOnly_GlobalEvaluationData::initialize: GED block is " << i << " is not initialized.");
  }

  gedBlocks_    = gedBlocks;

  ghostedSpace_ = rcp_dynamic_cast<const Thyra::DefaultProductVectorSpace<double> >(ghostedSpace);

  TEUCHOS_TEST_FOR_EXCEPTION(ghostedSpace_==Teuchos::null,std::logic_error,
                             "BlockedVector_ReadOnly_GED::initialize: ghosted space must be a Thyra::DefaultProductVectorSpace");

  isInitialized_ = true; 
}

void 
BlockedVector_ReadOnly_GlobalEvaluationData::
globalToGhost(int mem)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_,std::logic_error,
                             "BlockedVector_ReadOnly_GED has not been initialized, cannot call \"globalToGhost\"!");

  for(std::size_t i=0;i<gedBlocks_.size();i++)
    gedBlocks_[i]->globalToGhost(mem);
}

void 
BlockedVector_ReadOnly_GlobalEvaluationData::
initializeData()
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_,std::logic_error,
                             "BlockedVector_ReadOnly_GED has not been initialized, cannot call \"initializeData\"!");

  for(std::size_t i=0;i<gedBlocks_.size();i++)
    gedBlocks_[i]->initializeData();
}

void 
BlockedVector_ReadOnly_GlobalEvaluationData::
setUniqueVector(const Teuchos::RCP<const Thyra::VectorBase<double> > & uniqueVector)
{
  uniqueVector_ = uniqueVector;

  Teuchos::RCP<const Thyra::ProductVectorBase<double> > blocks = Thyra::castOrCreateProductVectorBase(uniqueVector_);

  TEUCHOS_TEST_FOR_EXCEPTION(blocks->productSpace()->numBlocks()!=Teuchos::as<int>(gedBlocks_.size()),std::logic_error,
                             "BlockedVector_ReadOnly_GED unique vector as the wrong number of blocks!");

  for(std::size_t i=0;i<gedBlocks_.size();i++)
    gedBlocks_[i]->setUniqueVector(blocks->getVectorBlock(i));
}

Teuchos::RCP<const Thyra::VectorBase<double> > 
BlockedVector_ReadOnly_GlobalEvaluationData::
getUniqueVector() const
{
  return uniqueVector_;
}

Teuchos::RCP<Thyra::VectorBase<double> > 
BlockedVector_ReadOnly_GlobalEvaluationData::
getGhostedVector() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_,std::logic_error,
                             "BlockedVector_ReadOnly_GED has not been initialized, cannot call \"getGhostedVector\"!");

  // unfortunately you must build the vector here!
  std::vector<Teuchos::RCP<Thyra::VectorBase<double> > > blocks; 
  for(std::size_t i=0;i<gedBlocks_.size();i++)
    blocks.push_back(gedBlocks_[i]->getGhostedVector());

  // why?????
  const std::vector<Teuchos::RCP<Thyra::VectorBase<double> > > & blocks_c = blocks; 
  return Thyra::defaultProductVector(ghostedSpace_,Teuchos::arrayViewFromVector(blocks_c)); 
}

}
