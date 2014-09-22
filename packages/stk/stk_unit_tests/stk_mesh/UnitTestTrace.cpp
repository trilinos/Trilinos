/*------------------------------------------------------------------------*/
/*                 Copyright (c) 2013, Sandia Corporation.
/*                 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*                 the U.S. Governement retains certain rights in this software.
/*                 
/*                 Redistribution and use in source and binary forms, with or without
/*                 modification, are permitted provided that the following conditions are
/*                 met:
/*                 
/*                     * Redistributions of source code must retain the above copyright
/*                       notice, this list of conditions and the following disclaimer.
/*                 
/*                     * Redistributions in binary form must reproduce the above
/*                       copyright notice, this list of conditions and the following
/*                       disclaimer in the documentation and/or other materials provided
/*                       with the distribution.
/*                 
/*                     * Neither the name of Sandia Corporation nor the names of its
/*                       contributors may be used to endorse or promote products derived
/*                       from this software without specific prior written permission.
/*                 
/*                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*                 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*                 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*                 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*                 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*                 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*                 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*                 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*                 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*                 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*                 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*                 
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifdef STK_MESH_TRACE_ENABLED

#include <gtest/gtest.h>

#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Trace.hpp>

#include <sstream>

using stk::mesh::EntityKey;

namespace {

// Test globals
const std::string SHOULD_SEE_STR = "should_see";
const std::string SHOULD_NOT_SEE_STR = "should_not_see";
unsigned SHOULD_SEE_COUNTER = 0;

std::string create_should_see_str(unsigned counter)
{
  std::ostringstream oss;
  oss << SHOULD_SEE_STR << counter;
  return oss.str();
}

const char* create_func(const std::string& func, bool should_see)
{
  std::ostringstream& oss = *(new std::ostringstream());
  if (should_see) {
    oss << func << "::" << create_should_see_str(SHOULD_SEE_COUNTER++);
  }
  else {
    oss << func << "::" << SHOULD_NOT_SEE_STR;
  }
  return oss.str().c_str();
}

std::string create_string(bool should_see)
{
  if (should_see) {
    return create_should_see_str(SHOULD_SEE_COUNTER++);
  }
  else {
    return SHOULD_NOT_SEE_STR;
  }
}

///////////////////////////////////////////////////////////////////////////////
TEST(UnitTestTrace, testTrace)
///////////////////////////////////////////////////////////////////////////////
{
  // A simple unit-test for the stk_mesh's tracing infrastructure. Note that
  // none of the Trace macros in the .cpp files in stk_mesh will be active
  // because the will have already been compiled without the macros defined.
  // This limits this unit test to very basic checks.

  // Local constants
  const EntityKey watch_key(0, 1);
  const EntityKey not_watch_key(0, 2);
  const std::string tracing_func = "stk::mesh::BulkData";
  const std::string not_tracing_func = "stk::mesh::MetaData";
  const stk::mesh::LogMask active_mask = stk::mesh::LOG_ENTITY;
  const stk::mesh::LogMask inactive_mask = stk::mesh::LOG_BUCKET;

  // Set up a dummy trace configuration. Here, we're telling the tracing
  // system that we want to trace BulkData calls related to entities,
  // specifically Node[1].
  std::ostringstream trace_output;
  stk::mesh::setStream(trace_output);
  meshlog.setPrintMask(active_mask | stk::mesh::LOG_TRACE);
  stk::mesh::watch(watch_key);
  stk::diag::Trace::addTraceFunction(tracing_func);

  //
  // Make calls to Trace API, some of which should generate trace output. We tag
  // output that should / should-not be in the trace so we can validate later.
  //

  {
    bool should_see = false; // not tracing func
    Trace_(create_func(not_tracing_func, should_see));
  }

  {
    bool should_see = true;
    Trace_(create_func(tracing_func, should_see));
  }

  {
    bool should_see = false; // not tracing func
    TraceIf(create_func(not_tracing_func, should_see), active_mask);
  }

  {
    bool should_see = false; // inactive mask
    TraceIf(create_func(tracing_func, should_see), inactive_mask);
    DiagIf(inactive_mask, create_string(should_see));
  }

  {
    bool should_see = true;
    TraceIf(create_func(tracing_func, should_see), active_mask);
    DiagIf(active_mask, create_string(should_see));
  }

  {
    bool should_see = false; // not tracing func
    TraceIfWatching(create_func(not_tracing_func, should_see), active_mask, watch_key);
  }

  {
    bool should_see = false; // inactive mask
    TraceIfWatching(create_func(tracing_func, should_see), inactive_mask, watch_key);
    DiagIfWatching(inactive_mask, watch_key, create_string(should_see));
  }

  {
    bool should_see = false; // not watching key
    TraceIfWatching(create_func(tracing_func, should_see), active_mask, not_watch_key);
    DiagIfWatching(active_mask, not_watch_key, create_string(should_see));
  }

  {
    bool should_see = true;
    TraceIfWatching(create_func(tracing_func, should_see), active_mask, watch_key);
    DiagIfWatching(active_mask, watch_key, create_string(should_see));
  }

  //
  // Check validity of output
  //

  const std::string trace_output_str = trace_output.str();

  // The not-trace tagged output should not be in the trace output
  ASSERT_EQ(trace_output_str.find(SHOULD_NOT_SEE_STR), std::string::npos);

  // Each occurance of should-see output should be in the trace output
  for (unsigned i = 0; i < SHOULD_SEE_COUNTER; ++i) {
    ASSERT_NE(trace_output_str.find(create_should_see_str(i)), std::string::npos);
  }
}


} // empty namespace

#endif
