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

#ifndef __Panzer_ResponseScatterEvaluator_Probe_hpp__
#define __Panzer_ResponseScatterEvaluator_Probe_hpp__

#include <iostream>
#include <string>

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Response_Probe.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {

class ProbeScatterBase {
public:
  virtual ~ProbeScatterBase() {}

  virtual void scatterDerivative(
    const panzer::Traits::Jacobian::ScalarT& probeValue,
    const size_t cell_index,
    const bool has_probe,
    panzer::Traits::EvalData workset,
    WorksetDetailsAccessor& wda,
    Teuchos::ArrayRCP<double> & dgdx) const = 0;
};

template <typename LO,typename GO>
class ProbeScatter : public ProbeScatterBase {
public:
   ProbeScatter(const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & globalIndexer)
     : globalIndexer_(globalIndexer) { }

   void scatterDerivative(
     const panzer::Traits::Jacobian::ScalarT& probeValue,
     const size_t cell_index,
     const bool has_probe,
     panzer::Traits::EvalData workset,
     WorksetDetailsAccessor& wda,
     Teuchos::ArrayRCP<double> & dgdx) const;

private:

   Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > globalIndexer_;
};

/** This class handles calculation of a DOF at a single point in space
 */
template<typename EvalT, typename Traits, typename LO, typename GO>
class ResponseScatterEvaluator_ProbeBase :
    public panzer::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits>  {
public:

  //! A constructor with concrete arguments instead of a parameter list.
  ResponseScatterEvaluator_ProbeBase(
    const std::string & responseName,
    const std::string & fieldName,
    const int fieldComponent,
    const Teuchos::Array<double>& point,
    const IntegrationRule & ir,
    const Teuchos::RCP<const PureBasis>& basis,
    const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> >& indexer,
    const Teuchos::RCP<ProbeScatterBase> & probeScatter);

  void evaluateFields(typename Traits::EvalData d);

  void preEvaluate(typename Traits::PreEvalData d);

protected:
  typedef typename EvalT::ScalarT ScalarT;

  std::string responseName_;
  std::string fieldName_;
  int fieldComponent_;
  Teuchos::Array<double> point_;
  Teuchos::RCP<const panzer::PureBasis> basis_;
  Teuchos::RCP<Response_Probe<EvalT> > responseObj_;
  Teuchos::RCP<const shards::CellTopology> topology_;
  Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > globalIndexer_;

  Teuchos::RCP<PHX::FieldTag> scatterHolder_; // dummy target
  PHX::MDField<const ScalarT,Cell,BASIS> field_; // holds field values
  Teuchos::RCP<ProbeScatterBase> scatterObj_;

  int cellIndex_;
  size_t num_basis, num_dim;
  Kokkos::DynRankView<double,PHX::Device> basis_values_;

  bool computeBasisValues(typename Traits::EvalData d);
};

/** This class handles calculation of a DOF at a single point in space
 */
template<typename EvalT, typename Traits, typename LO, typename GO>
class ResponseScatterEvaluator_Probe :
    public ResponseScatterEvaluator_ProbeBase<EvalT,Traits,LO,GO>  {
public:

  typedef ResponseScatterEvaluator_ProbeBase<EvalT,Traits,LO,GO> Base;

  //! A constructor with concrete arguments instead of a parameter list.
  ResponseScatterEvaluator_Probe(
    const std::string & responseName,
    const std::string & fieldName,
    const int fieldComponent,
    const Teuchos::Array<double>& point,
    const IntegrationRule & ir,
    const Teuchos::RCP<const PureBasis>& basis,
    const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & indexer,
    const Teuchos::RCP<ProbeScatterBase> & probeScatter) :
    Base(responseName, fieldName, fieldComponent, point,
         ir, basis, indexer, probeScatter) {}
};

/** This class handles calculation of a DOF at a single point in space
 */
template<typename LO, typename GO>
class ResponseScatterEvaluator_Probe<panzer::Traits::Jacobian,panzer::Traits,LO,GO> :
    public ResponseScatterEvaluator_ProbeBase<panzer::Traits::Jacobian,panzer::Traits,LO,GO>  {
public:

  typedef ResponseScatterEvaluator_ProbeBase<panzer::Traits::Jacobian,panzer::Traits,LO,GO> Base;

  //! A constructor with concrete arguments instead of a parameter list.
  ResponseScatterEvaluator_Probe(
    const std::string & responseName,
    const std::string & fieldName,
    const int fieldComponent,
    const Teuchos::Array<double>& point,
    const IntegrationRule & ir,
    const Teuchos::RCP<const PureBasis>& basis,
    const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & indexer,
    const Teuchos::RCP<ProbeScatterBase> & probeScatter) :
    Base(responseName, fieldName, fieldComponent, point,
         ir, basis, indexer, probeScatter) {}

  void evaluateFields(typename panzer::Traits::EvalData d);
};

template <typename LO,typename GO>
void ProbeScatter<LO,GO>::scatterDerivative(
  const panzer::Traits::Jacobian::ScalarT& probeValue,
  const size_t cell_index,
  const bool has_probe,
  panzer::Traits::EvalData workset,
  WorksetDetailsAccessor& wda,
  Teuchos::ArrayRCP<double> & dgdx) const
{

  if (has_probe) {
    Kokkos::View<const LO*, PHX::Device> LIDs = globalIndexer_->getElementLIDs(cell_index);

    // loop over basis functions
    for(std::size_t i=0; i<LIDs.size(); ++i) {
      dgdx[LIDs[i]] += probeValue.dx(i);
    }
  }
}

}

#endif
