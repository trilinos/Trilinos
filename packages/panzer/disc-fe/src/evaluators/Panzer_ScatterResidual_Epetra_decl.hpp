// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_SCATTER_RESIDUAL_EPETRA_DECL_HPP
#define PANZER_EVALUATOR_SCATTER_RESIDUAL_EPETRA_DECL_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

class Epetra_Vector;

namespace panzer {

class EpetraLinearObjContainer;

class GlobalIndexer;

/** \brief Pushes residual values into the residual vector for a 
           Newton-based solve

    Currently makes an assumption that the stride is constant for dofs
    and that the number of dofs is equal to the size of the solution
    names vector.

    \verbatim
      <ParameterList>
        <Parameter name="Scatter Name" type="string" value="<Name to give to the evaluate (dummy) field>"/>
        <Parameter name="Dependent Names" type="RCP<std::vector<std::string> >" value="<Name of fields to be scattered>"/>
        <Parameter name="Dependent Map" type="RCP<std::map<std::string,std::string> >" value="<Map from scattered field name to DOF name>"/>
        <Parameter name="Basis" type="RCP<const PureBasis>" value="<Basis function for all fields scattered (implicitly associated with DOF)>"/>
        <Parameter name="Global Data Key" type="string" value="<Lookup key for global evaluation data containers, defaults to \"Resiudal Scatter Container\">"/>
      </ParameterList>
    \endverbatim
*/
template<typename EvalT, typename TRAITS,typename LO,typename GO> class ScatterResidual_Epetra;

// **************************************************************
// **************************************************************
// * Specializations
// **************************************************************
// **************************************************************


// **************************************************************
// Residual 
// **************************************************************
template<typename TRAITS,typename LO,typename GO>
class ScatterResidual_Epetra<panzer::Traits::Residual,TRAITS,LO,GO>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, TRAITS>,
    public panzer::CloneableEvaluator {
  
public:
  ScatterResidual_Epetra(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
                         const Teuchos::RCP<const panzer::GlobalIndexer> & /* cIndexer=Teuchos::null */,
                         bool useDiscreteAdjoint=false) 
     : globalIndexer_(indexer),useDiscreteAdjoint_(useDiscreteAdjoint)  {}
  
  ScatterResidual_Epetra(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
                         const Teuchos::RCP<const panzer::GlobalIndexer> & cIndexer,
                         const Teuchos::ParameterList& p,bool=false);
  
  void postRegistrationSetup(typename TRAITS::SetupData d,
			     PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);
  
  void evaluateFields(typename TRAITS::EvalData workset);
  
  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterResidual_Epetra<panzer::Traits::Residual,TRAITS,LO,GO>(globalIndexer_,Teuchos::null,pl)); }

private:
  typedef typename panzer::Traits::Residual::ScalarT ScalarT;

  // dummy field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  // fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<const ScalarT,Cell,NODE> > scatterFields_;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer_;
  std::vector<int> fieldIds_; // field IDs needing mapping

  // This maps the scattered field names to the DOF manager field
  // For instance a Navier-Stokes map might look like
  //    fieldMap_["RESIDUAL_Velocity"] --> "Velocity"
  //    fieldMap_["RESIDUAL_Pressure"] --> "Pressure"
  Teuchos::RCP<const std::map<std::string,std::string> > fieldMap_;

  std::string globalDataKey_; // what global data does this fill?

  Teuchos::RCP<const EpetraLinearObjContainer> epetraContainer_;

  bool useDiscreteAdjoint_;
};

// **************************************************************
// Tangent 
// **************************************************************
template<typename TRAITS,typename LO,typename GO>
class ScatterResidual_Epetra<panzer::Traits::Tangent,TRAITS,LO,GO>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Tangent, TRAITS>,
    public panzer::CloneableEvaluator {
  
public:
  ScatterResidual_Epetra(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
                         const Teuchos::RCP<const panzer::GlobalIndexer> & /* cIndexer=Teuchos::null */,
                         bool useDiscreteAdjoint=false) 
     : globalIndexer_(indexer),useDiscreteAdjoint_(useDiscreteAdjoint)  {}
  
  ScatterResidual_Epetra(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
                         const Teuchos::RCP<const panzer::GlobalIndexer> & cIndexer,
                         const Teuchos::ParameterList& p,bool=false);
  
  void postRegistrationSetup(typename TRAITS::SetupData d,
			     PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);
  
  void evaluateFields(typename TRAITS::EvalData workset);
  
  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterResidual_Epetra<panzer::Traits::Tangent,TRAITS,LO,GO>(globalIndexer_,Teuchos::null,pl)); }

private:
  typedef typename panzer::Traits::Tangent::ScalarT ScalarT;

  // dummy field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  // fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<const ScalarT,Cell,NODE> > scatterFields_;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer_;
  std::vector<int> fieldIds_; // field IDs needing mapping

  // This maps the scattered field names to the DOF manager field
  // For instance a Navier-Stokes map might look like
  //    fieldMap_["RESIDUAL_Velocity"] --> "Velocity"
  //    fieldMap_["RESIDUAL_Pressure"] --> "Pressure"
  Teuchos::RCP<const std::map<std::string,std::string> > fieldMap_;

  std::vector<Teuchos::RCP<Epetra_Vector> > dfdp_vectors_;

  bool useDiscreteAdjoint_;
};

// **************************************************************
// Jacobian
// **************************************************************
template<typename TRAITS,typename LO,typename GO>
class ScatterResidual_Epetra<panzer::Traits::Jacobian,TRAITS,LO,GO>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Jacobian, TRAITS>, 
    public panzer::CloneableEvaluator {
  
public:
  
  ScatterResidual_Epetra(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
                         const Teuchos::RCP<const panzer::GlobalIndexer> & cIndexer=Teuchos::null,
                         bool useDiscreteAdjoint=false) 
     : globalIndexer_(indexer), colGlobalIndexer_(cIndexer), useDiscreteAdjoint_(useDiscreteAdjoint)  {}

  ScatterResidual_Epetra(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
                         const Teuchos::RCP<const panzer::GlobalIndexer> & cIndexer,
                         const Teuchos::ParameterList& pl,bool useDiscreteAdjoint=false);
  
  void postRegistrationSetup(typename TRAITS::SetupData d,
			     PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);
  
  void evaluateFields(typename TRAITS::EvalData workset);
  
  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterResidual_Epetra<panzer::Traits::Jacobian,TRAITS,LO,GO>(globalIndexer_,colGlobalIndexer_,pl,useDiscreteAdjoint_)); }

private:

  typedef typename panzer::Traits::Jacobian::ScalarT ScalarT;

  // dummy field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  // fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<const ScalarT,Cell,NODE> > scatterFields_;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer_, colGlobalIndexer_;
  std::vector<int> fieldIds_; // field IDs needing mapping

  // This maps the scattered field names to the DOF manager field
  // For instance a Navier-Stokes map might look like
  //    fieldMap_["RESIDUAL_Velocity"] --> "Velocity"
  //    fieldMap_["RESIDUAL_Pressure"] --> "Pressure"
  Teuchos::RCP<const std::map<std::string,std::string> > fieldMap_;

  std::string globalDataKey_; // what global data does this fill?

  Teuchos::RCP<const EpetraLinearObjContainer> epetraContainer_;

  ScatterResidual_Epetra();

  bool useDiscreteAdjoint_;

};

}

// optionally include hessian support
#ifdef Panzer_BUILD_HESSIAN_SUPPORT
#include "Panzer_ScatterResidual_Epetra_Hessian.hpp"
#endif

// **************************************************************
#endif
