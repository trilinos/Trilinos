// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_SCATTER_RESIDUAL_BLOCKEDEPETRA_DECL_HPP
#define PANZER_EVALUATOR_SCATTER_RESIDUAL_BLOCKEDEPETRA_DECL_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"
#include "Panzer_BlockedEpetraLinearObjContainer.hpp"

class Epetra_Vector;
class Epetra_CrsMatrix;

namespace panzer {

template <typename LocalOrdinalT,typename GlobalOrdinalT>
class UniqueGlobalIndexer;

template <typename LocalOrdinalT,typename GlobalOrdinalT>
class BlockedDOFManager; //forward declaration

/** \brief Pushes residual values into the residual vector for a 
           Newton-based solve

*/
template<typename EvalT, typename Traits,typename LO,typename GO> class ScatterResidual_BlockedEpetra
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, Traits>,
    public panzer::CloneableEvaluator  {
public:
   typedef typename EvalT::ScalarT ScalarT;
 
   ScatterResidual_BlockedEpetra(const Teuchos::RCP<const BlockedDOFManager<LO,int> > & indexer)
   { }
   ScatterResidual_BlockedEpetra(const Teuchos::RCP<const BlockedDOFManager<LO,int> > & gidProviders,
                                const Teuchos::ParameterList& p);

   virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
   { return Teuchos::rcp(new ScatterResidual_BlockedEpetra<EvalT,Traits,LO,GO>(Teuchos::null,pl)); }

  void postRegistrationSetup(typename Traits::SetupData d, PHX::FieldManager<Traits>& vm)
   { }
  void evaluateFields(typename Traits::EvalData d)
   { std::cout << "unspecialized version of \"ScatterResidual_BlockedEpetra::evaluateFields\" on \""+PHX::TypeString<EvalT>::value+"\" should not be used!" << std::endl;
     TEUCHOS_ASSERT(false); }
};

// **************************************************************
// **************************************************************
// * Specializations
// **************************************************************
// **************************************************************


// **************************************************************
// Residual 
// **************************************************************
template<typename Traits,typename LO,typename GO>
class ScatterResidual_BlockedEpetra<panzer::Traits::Residual,Traits,LO,GO>
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, Traits>,
    public panzer::CloneableEvaluator {
  
public:
  ScatterResidual_BlockedEpetra(const Teuchos::RCP<const BlockedDOFManager<LO,int> > & indexer)
     : globalIndexer_(indexer) {}
  
  ScatterResidual_BlockedEpetra(const Teuchos::RCP<const BlockedDOFManager<LO,int> > & indexer,
                                const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);

  void preEvaluate(typename Traits::PreEvalData d);
  
  void evaluateFields(typename Traits::EvalData workset);
  
  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterResidual_BlockedEpetra<panzer::Traits::Residual,Traits,LO,GO>(globalIndexer_,pl)); }

private:
  typedef typename panzer::Traits::Residual::ScalarT ScalarT;

  // dummy field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  // fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<ScalarT,Cell,NODE> > scatterFields_;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const BlockedDOFManager<LO,int> > globalIndexer_;

  std::vector<int> fieldIds_; // field IDs needing mapping

  // This maps the scattered field names to the DOF manager field
  // For instance a Navier-Stokes map might look like
  //    fieldMap_["RESIDUAL_Velocity"] --> "Velocity"
  //    fieldMap_["RESIDUAL_Pressure"] --> "Pressure"
  Teuchos::RCP<const std::map<std::string,std::string> > fieldMap_;

  std::string globalDataKey_; // what global data does this fill?
  Teuchos::RCP<const BlockedEpetraLinearObjContainer> blockedContainer_;

  ScatterResidual_BlockedEpetra();
};

// **************************************************************
// Jacobian
// **************************************************************

template<typename Traits,typename LO,typename GO>
class ScatterResidual_BlockedEpetra<panzer::Traits::Jacobian,Traits,LO,GO>
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<panzer::Traits::Jacobian, Traits>, 
    public panzer::CloneableEvaluator {
  
public:
  
  /** The parameter list passed takes the following values
      \verbatim
      <ParameterList>
         <Parameter name="Scatter Name" type="string" value=(required)/>
         <Parameter name="Dependent Names" type="RCP<vector<string> >" value="(required)"/>
         <Parameter name="Dependent Map" type="RCP<map<string,string> >" value="(required)"/>
         <Parameter name="Basis" type="RCP<const PureBasis>" value=(required)/>
         <Parameter name="Global Data Key" type="string" value="Residual Scatter Container" (default)/>
      </ParameterList>
      \endverbatim

  * The "Scatter Name" is the name for the dummy field that is computed by this evaluator.
  * This field should be required so that the evaluators is guranteed to run. "Dependent Names"
  * specifies the field to be scatter to the operator.  The "Dependent Map" gives a mapping from the
  * dependent field to the field string used in the global indexer. "Basis" is the basis
  * used to define the size of the "Dependent Names" fields. Finally "Global Data Key" is the key
  * used to index into the GlobalDataContainer object, for finding the operator and residual
  * linear algebra data structures that need to be filled. By default this is the simple residual/jacobian
  * with key "Residual Scatter Container".
  */
  ScatterResidual_BlockedEpetra(const Teuchos::RCP<const BlockedDOFManager<LO,int> > & indexer)
     : globalIndexer_(indexer) {}
  
  ScatterResidual_BlockedEpetra(const Teuchos::RCP<const BlockedDOFManager<LO,int> > & indexer,
                                const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);

  void preEvaluate(typename Traits::PreEvalData d);
  
  void evaluateFields(typename Traits::EvalData workset);
  
  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterResidual_BlockedEpetra<panzer::Traits::Jacobian,Traits,LO,GO>(globalIndexer_,pl)); }

private:

  typedef typename panzer::Traits::Jacobian::ScalarT ScalarT;

  // dummy field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  // fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<ScalarT,Cell,NODE> > scatterFields_;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const BlockedDOFManager<LO,int> > globalIndexer_;

  std::vector<int> fieldIds_; // field IDs needing mapping

  // This maps the scattered field names to the DOF manager field
  // For instance a Navier-Stokes map might look like
  //    fieldMap_["RESIDUAL_Velocity"] --> "Velocity"
  //    fieldMap_["RESIDUAL_Pressure"] --> "Pressure"
  Teuchos::RCP<const std::map<std::string,std::string> > fieldMap_;

  std::string globalDataKey_; // what global data does this fill?
  Teuchos::RCP<const BlockedEpetraLinearObjContainer> blockedContainer_;

  ScatterResidual_BlockedEpetra();
};

}

// **************************************************************
#endif
