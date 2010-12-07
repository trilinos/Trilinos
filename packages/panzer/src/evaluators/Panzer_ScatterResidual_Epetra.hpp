#ifndef PANZER_EVALUATOR_SCATTER_RESIDUAL_EPETRA_HPP
#define PANZER_EVALUATOR_SCATTER_RESIDUAL_EPETRA_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"

namespace panzer {

template <typename LocalOrdinalT,typename GlobalOrdinalT>
class DOFManager;

/** \brief Pushes residual values into the residual vector for a 
           Newton-based solve

    Currently makes an assumption that the stride is constant for dofs
    and that the number of dofs is equal to the size of the solution
    names vector.

*/
template<typename EvalT, typename Traits> class ScatterResidual_Epetra;

// **************************************************************
// **************************************************************
// * Specializations
// **************************************************************
// **************************************************************


// **************************************************************
// Residual 
// **************************************************************
template<typename Traits>
class ScatterResidual_Epetra<panzer::Traits::Residual,Traits>
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, Traits>  {
  
public:
  
  ScatterResidual_Epetra(const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData workset);
  
private:
  typedef typename panzer::Traits::Residual::ScalarT ScalarT;

  // dummy field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  // fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<ScalarT,Cell,NODE> > scatterFields_;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<panzer::DOFManager<panzer::Traits::LocalOrdinal,panzer::Traits::GlobalOrdinal> > dofManager_;
  std::vector<int> fieldIds_; // field IDs needing mapping

  // This maps the scattered field names to the DOF manager field
  // For instance a Navier-Stokes map might look like
  //    fieldMap_["RESIDUAL_Velocity"] --> "Velocity"
  //    fieldMap_["RESIDUAL_Pressure"] --> "Pressure"
  Teuchos::RCP<const std::map<std::string,std::string> > fieldMap_;
};

// **************************************************************
// Jacobian
// **************************************************************
template<typename Traits>
class ScatterResidual_Epetra<panzer::Traits::Jacobian,Traits>
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<panzer::Traits::Jacobian, Traits>  {
  
public:
  
  ScatterResidual_Epetra(const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData workset);
  
private:

  typedef typename panzer::Traits::Jacobian::ScalarT ScalarT;

  // dummy field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  // fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<ScalarT,Cell,NODE> > scatterFields_;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<panzer::DOFManager<panzer::Traits::LocalOrdinal,panzer::Traits::GlobalOrdinal> > dofManager_;
  std::vector<int> fieldIds_; // field IDs needing mapping

  // This maps the scattered field names to the DOF manager field
  // For instance a Navier-Stokes map might look like
  //    fieldMap_["RESIDUAL_Velocity"] --> "Velocity"
  //    fieldMap_["RESIDUAL_Pressure"] --> "Pressure"
  Teuchos::RCP<const std::map<std::string,std::string> > fieldMap_;
};

}

// **************************************************************

#include "Panzer_ScatterResidual_EpetraT.hpp"

// **************************************************************
#endif
