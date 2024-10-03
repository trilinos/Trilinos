// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_ObserverToLOCASaveDataStrategyAdapter.hpp"

#include "LOCA_Thyra.H"

#include "Teuchos_Ptr.hpp"

Piro::ObserverToLOCASaveDataStrategyAdapter::ObserverToLOCASaveDataStrategyAdapter(
    const Teuchos::RCP<ObserverBase<double> > &wrappedObserver) :
  wrappedObserver_(wrappedObserver)
{}

void
Piro::ObserverToLOCASaveDataStrategyAdapter::saveSolution(const NOX::Abstract::Vector &x, double p)
{
  const Teuchos::Ptr<const NOX::Thyra::Vector> x_thyra =
    Teuchos::ptr(dynamic_cast<const NOX::Thyra::Vector *>(&x));

  wrappedObserver_->observeSolution(x_thyra->getThyraVector(), p);
}
