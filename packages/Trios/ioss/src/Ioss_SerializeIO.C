// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
//         
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#include <Ioss_SerializeIO.h>

#include <string>
#include <Ioss_ParallelUtils.h>
#include <Ioss_Utils.h>
#include <Ioss_DatabaseIO.h>

namespace Ioss {

int
SerializeIO::s_owner = -1;

int
SerializeIO::s_rank = -1;

int
SerializeIO::s_size = -1;

int
SerializeIO::s_groupSize = -1;

int
SerializeIO::s_groupRank = -1;

int
SerializeIO::s_groupFactor = 0;

SerializeIO::SerializeIO(
  const DatabaseIO *database_io,
  int		    manual_owner_processor)
  : m_databaseIO(database_io),
    m_activeFallThru(s_owner != -1)

{

  const Ioss::ParallelUtils util = m_databaseIO->util();
  if (s_rank == -1) {
    s_rank = util.parallel_rank();
    s_size = util.parallel_size();
    if (s_groupFactor) {
      s_groupRank = s_rank/s_groupFactor;
      s_groupSize = (s_size - 1)/s_groupFactor + 1;
    }
  }

  m_manualOwner = (manual_owner_processor == -1 || s_groupFactor == 0) ? -1 : manual_owner_processor/s_groupFactor;

  if (m_activeFallThru) {
    if (m_manualOwner != -1 && m_manualOwner != s_owner) {
      std::ostringstream errmsg;
      errmsg << "Attempting to replace manual ownership from " << s_owner << " to " << m_manualOwner;
      IOSS_ERROR(errmsg);
    }
  }

  else if (s_groupFactor > 0) {
    if (m_manualOwner == -1) {
      m_databaseIO->openDatabase();
    }
    else {
      if (s_owner != -1 && m_manualOwner != s_owner) {
	std::ostringstream errmsg;
	errmsg << "Attempting to replace manual ownership from " << s_owner << " to " << m_manualOwner;
	IOSS_ERROR(errmsg);
      }
      s_owner = m_manualOwner;
    }
  }
  else
    s_owner = s_groupRank;
}


SerializeIO::~SerializeIO()
{
  try {
  if (m_activeFallThru)
    ;
  else if (s_groupFactor > 0) {
    if (m_manualOwner == -1) {
      m_databaseIO->closeDatabase();
      s_owner = -1;
    }
    else {
      if (s_owner == s_groupRank) {
	m_databaseIO->closeDatabase();
      }
      s_owner = -1;
    }
  }

  else
    s_owner = -1;
  } catch (...) {
  }
}


void
SerializeIO::setGroupFactor(
  int			factor)
{
  if (s_rank != -1)
    IOSS_WARNING << "Mesh I/O serialization group factor cannot be changed once serialized I/O has begun";
  else
    s_groupFactor = factor;
}

} // namespace Ioss
