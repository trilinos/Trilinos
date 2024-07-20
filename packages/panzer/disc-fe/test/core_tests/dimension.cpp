// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Panzer_Dimension.hpp"
#include "Phalanx_MDField.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(dimension, default)
  {
    PHX::MDField<double,Dim,IP,BASIS,NODE,Point,Cell,Edge,PHX::Device> a;
    PHX::MDField<double,Dim,IP,BASIS,NODE,Point,Cell,Edge,PHX::Device> b;
    a = b;
  }

}
