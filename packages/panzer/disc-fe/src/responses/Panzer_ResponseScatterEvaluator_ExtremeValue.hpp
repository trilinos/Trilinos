// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_ResponseScatterEvaluator_ExtremeValue_hpp__
#define __Panzer_ResponseScatterEvaluator_ExtremeValue_hpp__

#include <iostream>
#include <string>

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_Response_ExtremeValue.hpp"
#include "Panzer_GlobalIndexer.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {

class ExtremeValueScatterBase {
public:
  virtual ~ExtremeValueScatterBase() {}

  virtual void scatterDerivative(const PHX::MDField<const panzer::Traits::Jacobian::ScalarT,panzer::Cell> & cellExtremeValue,
                                 panzer::Traits::EvalData workset, 
                                 WorksetDetailsAccessor& wda,
                                 Teuchos::ArrayRCP<double> & dgdx) const = 0;
};
 
template <typename LO,typename GO>
class ExtremeValueScatter : public ExtremeValueScatterBase {
public:
   ExtremeValueScatter(const Teuchos::RCP<const panzer::GlobalIndexer> & globalIndexer)
     : globalIndexer_(globalIndexer) { }

   void scatterDerivative(const PHX::MDField<const panzer::Traits::Jacobian::ScalarT,panzer::Cell> & cellExtremeValue,
                         panzer::Traits::EvalData workset, 
                         WorksetDetailsAccessor& wda,
                         Teuchos::ArrayRCP<double> & dgdx) const;
private:
 
   Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer_;
};

/** This class handles responses with values aggregated
  * on each finite element cell.
  */
template<typename EvalT, typename Traits>
class ResponseScatterEvaluator_ExtremeValue : public panzer::EvaluatorWithBaseImpl<Traits>,
                                            public PHX::EvaluatorDerived<EvalT, Traits>  { 
public:

  //! A constructor with concrete arguments instead of a parameter list.
  ResponseScatterEvaluator_ExtremeValue(const std::string & name,const CellData & cd,
                                        bool useMax,
                                        const Teuchos::RCP<ExtremeValueScatterBase> & functionalScatter);
  ResponseScatterEvaluator_ExtremeValue(const std::string & extremeName,const std::string & responseName,const CellData & cd,
                                        bool useMax,
                                        const Teuchos::RCP<ExtremeValueScatterBase> & functionalScatter);

  void evaluateFields(typename Traits::EvalData d);

  void preEvaluate(typename Traits::PreEvalData d);

private:
  typedef typename EvalT::ScalarT ScalarT;

  std::string responseName_;
  Teuchos::RCP<Response_ExtremeValue<EvalT> > responseObj_;

  Teuchos::RCP<PHX::FieldTag> scatterHolder_; // dummy target
  PHX::MDField<const ScalarT,panzer::Cell> cellExtremeValue_; // holds cell integrals
  Teuchos::RCP<ExtremeValueScatterBase> scatterObj_;
  bool useMax_;
};

template <typename LO,typename GO>
void ExtremeValueScatter<LO,GO>::scatterDerivative(
                                        const PHX::MDField<const panzer::Traits::Jacobian::ScalarT,panzer::Cell> & /* cellExtremeValue */,
                                        panzer::Traits::EvalData /* workset */, 
                                        WorksetDetailsAccessor& /* wda */,
                                        Teuchos::ArrayRCP<double> & /* dgdx */) const
{
  TEUCHOS_ASSERT(false);
}

}

#endif
