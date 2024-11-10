// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_RCPNode.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace {


using Teuchos::RCPNodeTracer;


TEUCHOS_UNIT_TEST( RCPNodeTracer, defaults )
{
#if defined(TEUCHOS_DEBUG) && defined(HAVE_TEUCHOS_DEBUG_RCP_NODE_TRACING)
  TEST_EQUALITY_CONST(RCPNodeTracer::isTracingActiveRCPNodes(), true);
#else
  TEST_EQUALITY_CONST(RCPNodeTracer::isTracingActiveRCPNodes(), false);
#endif
  TEST_EQUALITY_CONST(RCPNodeTracer::getPrintRCPNodeStatisticsOnExit(), false);
  TEST_EQUALITY_CONST(RCPNodeTracer::getPrintActiveRcpNodesOnExit(), true);
}


TEUCHOS_UNIT_TEST( RCPNodeTracer, changeDefaults )
{
  TEST_EQUALITY_CONST(RCPNodeTracer::getPrintRCPNodeStatisticsOnExit(), false);
  ECHO(RCPNodeTracer::setPrintRCPNodeStatisticsOnExit(true));
  TEST_EQUALITY_CONST(RCPNodeTracer::getPrintRCPNodeStatisticsOnExit(), true);
  TEST_EQUALITY_CONST(RCPNodeTracer::getPrintActiveRcpNodesOnExit(), true);
  (RCPNodeTracer::setPrintActiveRcpNodesOnExit(false));
  TEST_EQUALITY_CONST(RCPNodeTracer::getPrintActiveRcpNodesOnExit(), false);
}


} // namespace
