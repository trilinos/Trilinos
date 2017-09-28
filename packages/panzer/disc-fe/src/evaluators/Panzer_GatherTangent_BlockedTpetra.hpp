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

#ifndef PANZER_EVALUATOR_GATHER_TANGENT_BLOCKEDTPETRA_DECL_HPP
#define PANZER_EVALUATOR_GATHER_TANGENT_BLOCKEDTPETRA_DECL_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
class BlockedTpetraLinearObjContainer;

template <typename LocalOrdinalT,typename GlobalOrdinalT>
class UniqueGlobalIndexer; //forward declaration

template <typename LocalOrdinalT,typename GlobalOrdinalT>
class BlockedDOFManager; //forward declaration

/** \brief Gathers tangent vectors dx/dp for computing df/dx*dx/dp + df/dp into
    the nodal fields of the field manager.

    This evaluator is very similar to GatherSolution, however it always gathers
    into fields of type double, and it is a no-op if the global evaluation data
    container does not exist (which is an error for GatherSolution).

    Currently makes an assumption that the stride is constant for dofs
    and that the nmber of dofs is equal to the size of the solution
    names vector.
*/
template <typename EvalT,typename TRAITS,typename S,typename LO,typename GO,typename NodeT=panzer::TpetraNodeType>
class GatherTangent_BlockedTpetra
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<EvalT, TRAITS>,
    public panzer::CloneableEvaluator  {
public:

   GatherTangent_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager<LO,GO> > & indexer)
     : gidIndexer_(indexer) {}

   GatherTangent_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager<LO,GO> > & indexer,
                                const Teuchos::ParameterList& p);

  void postRegistrationSetup(typename TRAITS::SetupData d,
                             PHX::FieldManager<TRAITS>& vm);

  void preEvaluate(typename TRAITS::PreEvalData d);

  void evaluateFields(typename TRAITS::EvalData d);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new GatherTangent_BlockedTpetra<EvalT,TRAITS,S,LO,GO>(gidIndexer_,pl)); }

private:

  // We always use RealType for gathering as we never compute derivatives for this evaluator
  //typedef typename panzer::Traits::RealType ScalarT;
  typedef typename EvalT::ScalarT ScalarT;

  typedef BlockedTpetraLinearObjContainer<S,LO,GO,NodeT> ContainerType;
  typedef Tpetra::Vector<S,LO,GO,NodeT> VectorType;
  typedef Tpetra::CrsMatrix<S,LO,GO,NodeT> CrsMatrixType;
  typedef Tpetra::CrsGraph<LO,GO,NodeT> CrsGraphType;
  typedef Tpetra::Map<LO,GO,NodeT> MapType;
  typedef Tpetra::Import<LO,GO,NodeT> ImportType;
  typedef Tpetra::Export<LO,GO,NodeT> ExportType;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<const BlockedDOFManager<LO,GO> > gidIndexer_;

  std::vector<int> fieldIds_; // field IDs needing mapping

  std::vector< PHX::MDField<ScalarT,Cell,NODE> > gatherFields_;

  Teuchos::RCP<std::vector<std::string> > indexerNames_;
  bool useTimeDerivativeSolutionVector_;
  std::string globalDataKey_; // what global data does this fill?

  Teuchos::RCP<const BlockedTpetraLinearObjContainer<S,LO,GO,NodeT> > blockedContainer_;

  GatherTangent_BlockedTpetra();
};

}

#include "Panzer_GatherTangent_BlockedTpetra_impl.hpp"

// **************************************************************
#endif
