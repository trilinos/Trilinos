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

#ifndef PANZER_EVALUATOR_SCATTER_RESIDUAL_BLOCKEDTPETRA_DECL_HPP
#define PANZER_EVALUATOR_SCATTER_RESIDUAL_BLOCKEDTPETRA_DECL_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"
#include "Panzer_BlockedTpetraLinearObjContainer.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {

// Forward declarations
template <typename LocalOrdinalT,typename GlobalOrdinalT> class BlockedDOFManager;
template <typename LocalOrdinalT,typename GlobalOrdinalT> class UniqueGlobalIndexer;

/** \brief Pushes residual values into the residual vector for a 
           Newton-based solve

*/
template <typename EvalT,typename TRAITS,typename LO,typename GO,typename NodeT=panzer::TpetraNodeType>
class ScatterResidual_BlockedTpetra
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, TRAITS>,
    public panzer::CloneableEvaluator  {
public:
   typedef typename EvalT::ScalarT ScalarT;
 
   ScatterResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager<LO,GO> > & /* indexer */)
   { }
   ScatterResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager<LO,GO> > & gidProviders,
                                const Teuchos::ParameterList& p);

   virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
   { return Teuchos::rcp(new ScatterResidual_BlockedTpetra<EvalT,TRAITS,LO,GO,NodeT>(Teuchos::null,pl)); }

  void postRegistrationSetup(typename TRAITS::SetupData /* d */, PHX::FieldManager<TRAITS>& /* vm */)
   { }
  void evaluateFields(typename TRAITS::EvalData /* d */)
   { std::cout << "unspecialized version of \"ScatterResidual_BlockedTpetra::evaluateFields\" on \""+PHX::typeAsString<EvalT>()+"\" should not be used!" << std::endl;
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
template <typename TRAITS,typename LO,typename GO,typename NodeT>
class ScatterResidual_BlockedTpetra<panzer::Traits::Residual,TRAITS,LO,GO,NodeT>
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, TRAITS>,
    public panzer::CloneableEvaluator {
  
public:
  ScatterResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager<LO,GO> > & indexer)
     : globalIndexer_(indexer) {}
  
  ScatterResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager<LO,GO> > & indexer,
                                const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename TRAITS::SetupData d,
			     PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);
  
  void evaluateFields(typename TRAITS::EvalData workset);
  
  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterResidual_BlockedTpetra<panzer::Traits::Residual,TRAITS,LO,GO,NodeT>(globalIndexer_,pl)); }

private:
  typedef typename panzer::Traits::Residual::ScalarT ScalarT;
  typedef typename TRAITS::RealType RealType;

  typedef BlockedTpetraLinearObjContainer<RealType,LO,GO,NodeT> ContainerType;
  typedef Tpetra::Vector<RealType,LO,GO,NodeT> VectorType;
  typedef Tpetra::CrsMatrix<RealType,LO,GO,NodeT> CrsMatrixType;
  typedef Tpetra::CrsGraph<LO,GO,NodeT> CrsGraphType;
  typedef Tpetra::Map<LO,GO,NodeT> MapType;
  typedef Tpetra::Import<LO,GO,NodeT> ImportType;
  typedef Tpetra::Export<LO,GO,NodeT> ExportType;

  //! Dummy evalauted field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  //! Fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<const ScalarT,Cell,NODE> > scatterFields_;

  //! Maps the local (field,element,basis) triplet to a global ID for scattering
  Teuchos::RCP<const BlockedDOFManager<LO,GO> > globalIndexer_;

  //! Vector of global indexers, one for each scattered field respectively
  std::vector<Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO>>> fieldGlobalIndexers_;

  //! Field IDs in the local product vector block (not global field id)
  std::vector<int> fieldIds_;

  //! Returns the index into the Thyra ProductVector sub-block. Size
  //! of number of fields to scatter
  std::vector<int> productVectorBlockIndex_;

  // This maps the scattered field names to the DOF manager field
  // For instance a Navier-Stokes map might look like
  //    fieldMap_["RESIDUAL_Velocity"] --> "Velocity"
  //    fieldMap_["RESIDUAL_Pressure"] --> "Pressure"
  Teuchos::RCP<const std::map<std::string,std::string> > fieldMap_;

  std::string globalDataKey_; // what global data does this fill?
  Teuchos::RCP<const BlockedTpetraLinearObjContainer<RealType,LO,GO,NodeT> > blockedContainer_;

  //! Local indices for unknowns
  Kokkos::View<LO**,PHX::Device> worksetLIDs_;

  //! Offset into the cell lids for each field
  std::vector<Kokkos::View<int*,PHX::Device>> fieldOffsets_;

  ScatterResidual_BlockedTpetra();
};

// **************************************************************
// Jacobian
// **************************************************************

template <typename TRAITS,typename LO,typename GO,typename NodeT>
class ScatterResidual_BlockedTpetra<panzer::Traits::Jacobian,TRAITS,LO,GO,NodeT>
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
  ScatterResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager<LO,GO> > & indexer)
     : globalIndexer_(indexer) {}
  
  ScatterResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager<LO,GO> > & indexer,
                                const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename TRAITS::SetupData d,
			     PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);
  
  void evaluateFields(typename TRAITS::EvalData workset);
  
  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new ScatterResidual_BlockedTpetra<panzer::Traits::Jacobian,TRAITS,LO,GO,NodeT>(globalIndexer_,pl)); }

private:
  typedef typename panzer::Traits::Jacobian::ScalarT ScalarT;
  typedef typename TRAITS::RealType RealType;

  typedef BlockedTpetraLinearObjContainer<RealType,LO,GO,NodeT> ContainerType;
  typedef Tpetra::Operator<RealType,LO,GO,NodeT> OperatorType;
  typedef Tpetra::CrsMatrix<RealType,LO,GO,NodeT> CrsMatrixType;
  typedef Tpetra::Map<LO,GO,NodeT> MapType;

  typedef Thyra::TpetraLinearOp<RealType,LO,GO,NodeT> ThyraLinearOp;

  // dummy field so that the evaluator will have something to do
  Teuchos::RCP<PHX::FieldTag> scatterHolder_;

  // fields that need to be scattered will be put in this vector
  std::vector< PHX::MDField<const ScalarT,Cell,NODE> > scatterFields_;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const BlockedDOFManager<LO,GO> > globalIndexer_;

  std::vector<int> fieldIds_; // field IDs needing mapping

  //! Returns the index into the Thyra ProductVector sub-block. Size
  //! of number of fields to scatter
  std::vector<int> productVectorBlockIndex_;

  // This maps the scattered field names to the DOF manager field
  // For instance a Navier-Stokes map might look like
  //    fieldMap_["RESIDUAL_Velocity"] --> "Velocity"
  //    fieldMap_["RESIDUAL_Pressure"] --> "Pressure"
  Teuchos::RCP<const std::map<std::string,std::string> > fieldMap_;

  std::string globalDataKey_; // what global data does this fill?
  Teuchos::RCP<const BlockedTpetraLinearObjContainer<RealType,LO,GO,NodeT> > blockedContainer_;

  //! Local indices for unknowns
  Kokkos::View<LO**,PHX::Device> worksetLIDs_;

  //! Offset into the cell lids for each field. Size of number of fields to scatter.
  std::vector<Kokkos::View<int*,PHX::Device>> fieldOffsets_;

  //! The offset values of the blocked DOFs per element. Size of number of blocks in the product vector + 1. The plus one is a sentinel.
  Kokkos::View<LO*,PHX::Device> blockOffsets_;

  ScatterResidual_BlockedTpetra();
};

}

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
#include "Panzer_ScatterResidual_BlockedTpetra_Hessian.hpp"
#endif

// **************************************************************
#endif
