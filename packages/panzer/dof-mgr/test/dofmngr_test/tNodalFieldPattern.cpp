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

#include <string>
#include <iostream>

#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_NodalFieldPattern.hpp"

// include some intrepid basis functions
// 2D basis 
#include "Intrepid2_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_C2_FEM.hpp"

// 3D basis 
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"

#include "Intrepid2_HDIV_TRI_I1_FEM.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {

/////////////////////////////////////////////
// 2D tests
/////////////////////////////////////////////

// triangle tests
TEUCHOS_UNIT_TEST(tNodalFieldPattern, test2d_tri_c1)
{
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis1, basis2;

   basis1 = rcp(new Intrepid2::Basis_HGRAD_TRI_C1_FEM<PHX::exec_space,double,double>);
   basis2 = rcp(new Intrepid2::Basis_HGRAD_TRI_C2_FEM<PHX::exec_space,double,double>);

   Teuchos::RCP<FieldPattern> pattern1 = rcp(new Intrepid2FieldPattern(basis1));
   Teuchos::RCP<FieldPattern> pattern2 = rcp(new Intrepid2FieldPattern(basis2));
   Teuchos::RCP<FieldPattern> nodalPattern = rcp(new NodalFieldPattern(pattern1->getCellTopology()));
  
   TEST_ASSERT(nodalPattern->equals(*pattern1))
   TEST_ASSERT(!nodalPattern->equals(*pattern2))
}

TEUCHOS_UNIT_TEST(tNodalFieldPattern, test3d_HEX_c1)
{
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis1, basis2;

   basis1 = rcp(new Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::exec_space,double,double>);
   basis2 = rcp(new Intrepid2::Basis_HGRAD_HEX_C2_FEM<PHX::exec_space,double,double>);

   Teuchos::RCP<FieldPattern> pattern1 = rcp(new Intrepid2FieldPattern(basis1));
   Teuchos::RCP<FieldPattern> pattern2 = rcp(new Intrepid2FieldPattern(basis2));
   Teuchos::RCP<FieldPattern> nodalPattern = rcp(new NodalFieldPattern(pattern1->getCellTopology()));
  
   TEST_ASSERT(nodalPattern->equals(*pattern1))
   TEST_ASSERT(!nodalPattern->equals(*pattern2))
}

}
