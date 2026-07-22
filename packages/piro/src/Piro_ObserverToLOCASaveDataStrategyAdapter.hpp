// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_OBSERVERTOLOCASAVEDATASTRATEGYADAPTER_HPP
#define PIRO_OBSERVERTOLOCASAVEDATASTRATEGYADAPTER_HPP

#include "LOCA_Thyra_SaveDataStrategy.H"

#include "Piro_ObserverBase.hpp"

#include "Teuchos_RCP.hpp"

namespace Piro {

class ObserverToLOCASaveDataStrategyAdapter : public LOCA::Thyra::SaveDataStrategy {
public:
  explicit ObserverToLOCASaveDataStrategyAdapter(
      const Teuchos::RCP<ObserverBase<double> > &wrappedObserver);

  // Overriden from LOCA::Thyra::SaveDataStrategy

  virtual void saveSolution(const NOX::Abstract::Vector &x, double p);

private:
  Teuchos::RCP<ObserverBase<double> > wrappedObserver_;
};

} // namespace Piro

#endif /* PIRO_OBSERVERTOLOCASAVEDATASTRATEGYADAPTER_HPP */
