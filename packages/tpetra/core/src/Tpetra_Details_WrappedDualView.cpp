// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_WrappedDualView.hpp"

namespace Tpetra {
namespace Details {

bool wdvTrackingEnabled = true;

void enableWDVTracking()
{
  if(wdvTrackingEnabled)
    throw std::runtime_error("WrappedDualView refcount tracking is already enabled!");
  wdvTrackingEnabled = true;
}

void disableWDVTracking()
{
  if(!wdvTrackingEnabled)
    throw std::runtime_error("WrappedDualView refcount tracking is already disabled!");
  wdvTrackingEnabled = false;
}

} // namespace Details
} // namespace Tpetra
