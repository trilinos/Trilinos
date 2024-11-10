// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_Response_Probe_hpp__
#define __Panzer_Response_Probe_hpp__

#include <string>
#include <limits>

#include <mpi.h> // need for comm

#include "Teuchos_RCP.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"

#include "PanzerDiscFE_config.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MpiComm.h"
#endif

#include "Panzer_ResponseMESupport_Default.hpp"
#include "Panzer_GlobalEvaluationData.hpp"
#include "Panzer_ThyraObjFactory.hpp"
#include "Panzer_ThyraObjContainer.hpp"
#include "Panzer_LinearObjFactory.hpp"


namespace panzer {

/** This class provides a response storage for
  * the value of the solution at a point in space/time.
  */
template <typename EvalT>
class Response_Probe : public ResponseMESupport_Default<EvalT>
                          , public GlobalEvaluationData_BCAdjustment {
public:
  typedef typename EvalT::ScalarT ScalarT;

  ScalarT value;
  bool have_probe;  // Point for probe is on our processor

  Response_Probe(const std::string & responseName, MPI_Comm comm,
                 const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & linObjFact=Teuchos::null);

  //! This simply does global summation, then shoves the result into a vector
  virtual void scatterResponse();

  virtual void initializeResponse();

  // from ResponseMESupport_Default

  //! What is the number of values you need locally
  virtual std::size_t localSizeRequired() const { return 1; }

  //! Is the vector distributed (or replicated)
  virtual bool vectorIsDistributed() const { return false; }

  //! Get ghosted responses (this will be filled by the evaluator)
  Teuchos::RCP<Thyra::VectorBase<double> > getGhostedVector() const
  { return Teuchos::rcp_dynamic_cast<const ThyraObjContainer<double> >(ghostedContainer_)->get_x_th(); }

  void adjustForDirichletConditions(const GlobalEvaluationData & localBCRows,const GlobalEvaluationData & globalBCRows);

private:
  //! Set solution vector space
  void setSolnVectorSpace(const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & soln_vs);

  // hide these methods
  Response_Probe();
  Response_Probe(const Response_Probe &);

  Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > linObjFactory_;
  Teuchos::RCP<const panzer::ThyraObjFactory<double> > thyraObjFactory_;

  Teuchos::RCP<LinearObjContainer> uniqueContainer_;
  Teuchos::RCP<LinearObjContainer> ghostedContainer_;
};

}

#endif
