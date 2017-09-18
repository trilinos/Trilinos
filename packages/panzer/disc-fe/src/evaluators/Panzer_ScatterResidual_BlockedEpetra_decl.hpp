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

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"
#include "Panzer_BlockedEpetraLinearObjContainer.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace Thyra {
  template <typename> class ProductVectorBase;
  template <typename> class BlockedLinearOpBase;
}

namespace panzer {

template <typename LocalOrdinalT,typename GlobalOrdinalT>
class UniqueGlobalIndexer;

/** \brief Pushes residual values into the residual vector for a 
           Newton-based solve

*/
template<typename EvalT, typename TRAITS,typename LO,typename GO> class ScatterResidual_BlockedEpetra
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, TRAITS>,
    public panzer::CloneableEvaluator  {
public:
   typedef typename EvalT::ScalarT ScalarT;
 
   ScatterResidual_BlockedEpetra(const Teuchos::ParameterList& p) {}

   virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
   { return Teuchos::rcp(new ScatterResidual_BlockedEpetra<EvalT,TRAITS,LO,GO>(pl)); }

  void postRegistrationSetup(typename TRAITS::SetupData d, PHX::FieldManager<TRAITS>& vm)
   { }
  void evaluateFields(typename TRAITS::EvalData d)
   { std::cout << "unspecialized version of \"ScatterResidual_BlockedEpetra::evaluateFields\" on "+PHX::typeAsString<EvalT>()+" \" should not be used!" << std::endl;
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
template<typename TRAITS,typename LO,typename GO>
class ScatterResidual_BlockedEpetra<panzer::Traits::Residual,TRAITS,LO,GO>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, TRAITS>,
    public panzer::CloneableEvaluator {
  
public:

  ScatterResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & rIndexers,
                                const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & /* cIndexers */,
                                bool /* useDiscreteAdjoint=false */)
     : rowIndexers_(rIndexers) {}
  
  ScatterResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & rIndexers,
                                const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & cIndexers,
                                const Teuchos::ParameterList& p,
                                bool useDiscreteAdjoint=false);

  void postRegistrationSetup(typename TRAITS::SetupData d,
			     PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);
  
  void evaluateFields(typename TRAITS::EvalData workset);
  
  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterResidual_BlockedEpetra<panzer::Traits::Residual,TRAITS,LO,GO>(rowIndexers_,colIndexers_,pl)); }

private:
  typedef typename panzer::Traits::Residual::ScalarT ScalarT;

  // dummy field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  // fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<const ScalarT,Cell,NODE> > scatterFields_;

  std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > rowIndexers_;
  std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > colIndexers_;

  std::vector<int> indexerIds_;   // block index
  std::vector<int> subFieldIds_; // sub field numbers

  // This maps the scattered field names to the DOF manager field
  // For instance a Navier-Stokes map might look like
  //    fieldMap_["RESIDUAL_Velocity"] --> "Velocity"
  //    fieldMap_["RESIDUAL_Pressure"] --> "Pressure"
  Teuchos::RCP<const std::map<std::string,std::string> > fieldMap_;

  std::string globalDataKey_; // what global data does this fill?

  Teuchos::RCP<Thyra::ProductVectorBase<double> > r_;

  ScatterResidual_BlockedEpetra();
};

// **************************************************************
// Tangent 
// **************************************************************
template<typename TRAITS,typename LO,typename GO>
class ScatterResidual_BlockedEpetra<panzer::Traits::Tangent,TRAITS,LO,GO>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Tangent, TRAITS>,
    public panzer::CloneableEvaluator {
  
public:

  ScatterResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & rIndexers,
                                const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & /* cIndexers */,
                                bool /* useDiscreteAdjoint=false */)
     : rowIndexers_(rIndexers) {}
  
  ScatterResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & rIndexers,
                                const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & cIndexers,
                                const Teuchos::ParameterList& p,
                                bool useDiscreteAdjoint=false);
  
  void postRegistrationSetup(typename TRAITS::SetupData d,
			     PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);
  
  void evaluateFields(typename TRAITS::EvalData workset);
  
  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterResidual_BlockedEpetra<panzer::Traits::Tangent,TRAITS,LO,GO>(rowIndexers_,colIndexers_,pl)); }

private:
  typedef typename panzer::Traits::Tangent::ScalarT ScalarT;

  // dummy field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  // fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<const ScalarT,Cell,NODE> > scatterFields_;

  std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > rowIndexers_;
  std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > colIndexers_;

  std::vector<int> indexerIds_;   // block index
  std::vector<int> subFieldIds_; // sub field numbers

  // This maps the scattered field names to the DOF manager field
  // For instance a Navier-Stokes map might look like
  //    fieldMap_["RESIDUAL_Velocity"] --> "Velocity"
  //    fieldMap_["RESIDUAL_Pressure"] --> "Pressure"
  Teuchos::RCP<const std::map<std::string,std::string> > fieldMap_;

  std::string globalDataKey_; // what global data does this fill?

  Teuchos::RCP<Thyra::ProductVectorBase<double> > r_;

  ScatterResidual_BlockedEpetra();
};

// **************************************************************
// Jacobian
// **************************************************************

template<typename TRAITS,typename LO,typename GO>
class ScatterResidual_BlockedEpetra<panzer::Traits::Jacobian,TRAITS,LO,GO>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Jacobian, TRAITS>, 
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
  ScatterResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & rIndexers,
                                const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & cIndexers,
                                bool useDiscreteAdjoint=false)
     : rowIndexers_(rIndexers), colIndexers_(cIndexers), useDiscreteAdjoint_(useDiscreteAdjoint) {}
  
  ScatterResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & rIndexers,
                                const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & cIndexers,
                                const Teuchos::ParameterList& p,
                                bool useDiscreteAdjoint=false);
  
  void postRegistrationSetup(typename TRAITS::SetupData d,
			     PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);
  
  void evaluateFields(typename TRAITS::EvalData workset);
  
  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterResidual_BlockedEpetra<panzer::Traits::Jacobian,TRAITS,LO,GO>(rowIndexers_,colIndexers_,pl,useDiscreteAdjoint_)); }

private:

  typedef typename panzer::Traits::Jacobian::ScalarT ScalarT;

  // dummy field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  // fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<const ScalarT,Cell,NODE> > scatterFields_;

  std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > rowIndexers_;
  std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > colIndexers_;

  std::vector<int> indexerIds_;   // block index
  std::vector<int> subFieldIds_; // sub field numbers

  // This maps the scattered field names to the DOF manager field
  // For instance a Navier-Stokes map might look like
  //    fieldMap_["RESIDUAL_Velocity"] --> "Velocity"
  //    fieldMap_["RESIDUAL_Pressure"] --> "Pressure"
  Teuchos::RCP<const std::map<std::string,std::string> > fieldMap_;

  std::string globalDataKey_; // what global data does this fill?
  bool useDiscreteAdjoint_;

  Teuchos::RCP<Thyra::ProductVectorBase<double> > r_;
  Teuchos::RCP<Thyra::BlockedLinearOpBase<double> > Jac_;

  ScatterResidual_BlockedEpetra();
};

}

// optionally include hessian support
#ifdef Panzer_BUILD_HESSIAN_SUPPORT
#include "Panzer_ScatterResidual_BlockedEpetra_Hessian.hpp"
#endif

// **************************************************************
#endif
