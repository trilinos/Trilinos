#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_Basis.hpp"

#include "Teuchos_FancyOStream.hpp"

template<typename EvalT,typename Traits,typename LO,typename GO>
panzer::GatherOrientation<EvalT, Traits,LO,GO>::
GatherOrientation(
  const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & indexer,
  const Teuchos::ParameterList& p)
  : globalIndexer_(indexer)
{ 
  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("DOF Names"));

  indexerNames_ = p.get< Teuchos::RCP< std::vector<std::string> > >("Indexer Names");

  Teuchos::RCP<panzer::Basis> basis = 
    p.get< Teuchos::RCP<panzer::Basis> >("Basis");

  gatherFieldOrientations_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    gatherFieldOrientations_[fd] = 
      PHX::MDField<ScalarT,Cell,NODE>(names[fd]+" Orientation",basis->functional);
    this->addEvaluatedField(gatherFieldOrientations_[fd]);
  }

  this->setName("Gather Orientation");
}

// **********************************************************************
template<typename EvalT,typename Traits,typename LO,typename GO>
void panzer::GatherOrientation<EvalT, Traits,LO,GO>::
postRegistrationSetup(typename Traits::SetupData d, 
		      PHX::FieldManager<Traits>& fm)
{
  // globalIndexer_ = d.globalIndexer_;
  TEUCHOS_ASSERT(gatherFieldOrientations_.size() == indexerNames_->size());

  fieldIds_.resize(gatherFieldOrientations_.size());

  for (std::size_t fd = 0; fd < gatherFieldOrientations_.size(); ++fd) {
    // get field ID from DOF manager
    //std::string fieldName = gatherFieldOrientations_[fd].fieldTag().name();
    const std::string& fieldName = (*indexerNames_)[fd];
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // setup the field data object
    this->utils.setFieldData(gatherFieldOrientations_[fd],fm);
  }

  indexerNames_ = Teuchos::null;  // Don't need this anymore
}

// **********************************************************************
template<typename EvalT,typename Traits,typename LO,typename GO>
void panzer::GatherOrientation<EvalT, Traits,LO,GO>::
evaluateFields(typename Traits::EvalData workset)
{ 
   std::vector<typename Traits::GlobalOrdinal> GIDs;
   std::vector<int> LIDs;
   std::vector<double> orientation;
 
   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!
 
   // gather operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];
 
      globalIndexer_->getElementOrientation(cellLocalId,orientation); 
 
      // loop over the fields to be gathered
      for (std::size_t fieldIndex=0; fieldIndex<gatherFieldOrientations_.size();fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
         const std::vector<int> & elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);
 
         // loop over basis functions and fill the fields
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            (gatherFieldOrientations_[fieldIndex])(worksetCellIndex,basis) =orientation[offset];
         }
      }
   }
}
