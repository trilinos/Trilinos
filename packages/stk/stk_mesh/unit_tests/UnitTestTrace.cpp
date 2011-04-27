/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

//#ifndef STK_MESH_TRACE_ENABLED
//#define STK_MESH_TRACE_ENABLED
//#endif

// Replace #ifdef below with above 3 lines once the UseCase dependency
// issue is sorted out
#ifdef STK_MESH_TRACE_ENABLED

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_mesh/base/Trace.hpp>
#include <stk_mesh/base/EntityKey.hpp>

using stk::mesh::EntityKey;

namespace {

///////////////////////////////////////////////////////////////////////////////
STKUNIT_UNIT_TEST(UnitTestTrace, testTrace)
///////////////////////////////////////////////////////////////////////////////
{
  // A simple unit-test for the stk_mesh's tracing infrastructure. Note that
  // none of the Trace macros in the .cpp files in stk_mesh will be active
  // because the will have already been compiled without the macros defined.
  // This limits this unit test to very basic checks.

  EntityKey watch_key(0, 1);
  EntityKey not_watch_key(0, 2);
  
  // Set up a dummy trace configuration. Here, we're telling the tracing
  // system that we want to trace BulkData calls related to entities,
  // specifically Node[1].
  meshlog.setPrintMask(stk::mesh::LOG_ENTITY | stk::mesh::LOG_TRACE);
  stk::mesh::watch(watch_key);
  stk::diag::Trace::addTraceFunction("stk::mesh::BulkData");

  {
    Trace_("stk::mesh::MetaData::should_not_see");
  }

  {
    Trace_("stk::mesh::BulkData::should_see1");
  }

  {
    TraceIf("stk::mesh::MetaData::should_not_see", stk::mesh::LOG_TRACE);
  }

  {
    TraceIf("stk::mesh::BulkData::should_not_see", stk::mesh::LOG_BUCKET);
    DiagIf(stk::mesh::LOG_BUCKET, "should not see");
  }

  {
    TraceIf("stk::mesh::BulkData::should_see2", stk::mesh::LOG_ENTITY);
    DiagIf(stk::mesh::LOG_ENTITY, "should see2");
  }

  {
    TraceIfWatching("stk::mesh::MetaData::should_not_see", stk::mesh::LOG_TRACE, watch_key);
  }

  {
    TraceIfWatching("stk::mesh::BulkData::should_not_see", stk::mesh::LOG_BUCKET, watch_key);
    DiagIfWatching(stk::mesh::LOG_BUCKET, watch_key, "should not see");
  }

  {
    TraceIfWatching("stk::mesh::BulkData::should_not_see", stk::mesh::LOG_TRACE, not_watch_key);
    DiagIfWatching(stk::mesh::LOG_TRACE, not_watch_key, "should not see");
  }

  {
    TraceIfWatching("stk::mesh::BulkData::should_see3", stk::mesh::LOG_TRACE, watch_key);
    DiagIfWatching(stk::mesh::LOG_TRACE, watch_key, "should see3");
  }
}

#undef STK_MESH_TRACE_ENABLED

} // empty namespace

#endif
