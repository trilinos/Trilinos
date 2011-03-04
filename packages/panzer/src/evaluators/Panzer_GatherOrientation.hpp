#ifndef PANZER_EVALUATOR_GATHER_ORIENTATION_HPP
#define PANZER_EVALUATOR_GATHER_ORIENTATION_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"

namespace panzer {

template <typename LocalOrdinalT,typename GlobalOrdinalT>
class UniqueGlobalIndexer; //forward declaration

/** \brief Gathers orientations per field from the global indexer and
    stores them in the field manager.
*/
template<typename EvalT, typename Traits,typename LO,typename GO> 
class GatherOrientation
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits>,
    public panzer::CloneableEvaluator  {
   
public:
  
  GatherOrientation(const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & indexer) :
     globalIndexer_(indexer) {}

  GatherOrientation(const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & indexer,
                        const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new GatherOrientation<panzer::Traits::Residual,Traits,LO,GO>(globalIndexer_,pl)); }
  
private:

  typedef typename EvalT::ScalarT ScalarT;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > globalIndexer_;
  std::vector<int> fieldIds_; // field IDs needing mapping

  std::vector< PHX::MDField<ScalarT,Cell,NODE> > gatherFieldOrientations_;

  Teuchos::RCP<std::vector<std::string> > indexerNames_;

  GatherOrientation();
};

}

// **************************************************************

#include "Panzer_GatherOrientationT.hpp"

// **************************************************************
#endif
