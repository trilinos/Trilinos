/*
 * Copyright(C) 2010-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include <cstdlib>

#include <string>

#include <EP_Internals.h>
#include <EP_ParallelDisks.h>

#ifdef _WIN32
#include <Shlwapi.h>
#endif

/*****************************************************************************/
Excn::ParallelDisks::ParallelDisks() : number_of_raids(0), raid_offset(0) {}

/*****************************************************************************/
Excn::ParallelDisks::~ParallelDisks() = default;

/*****************************************************************************/
void Excn::ParallelDisks::Number_of_Raids(int i)
{
  number_of_raids = i;
  create_disk_names();
}

/*****************************************************************************/
void Excn::ParallelDisks::Raid_Offset(int i)
{
  raid_offset = i;
  create_disk_names();
}

/*****************************************************************************/
int Excn::ParallelDisks::Number_of_Raids() const { return number_of_raids; }

/*****************************************************************************/
int Excn::ParallelDisks::Raid_Offset() const { return raid_offset; }

/*****************************************************************************/

/*****************************************************************************/
void Excn::ParallelDisks::rename_file_for_mp(const std::string &rootdir, const std::string &subdir,
                                             std::string &name, int node, int numproc) const
{
  // Possible to have node layout without parallel disks

  std::string prepend;
  if (rootdir.length() != 0u) {
    prepend = rootdir + "/";
  }
  else if (Excn::is_path_absolute(name)) {
    prepend = "";
  }
  else {
    prepend = "./";
  }

  int lnn = node;
  if (number_of_raids != 0) {
    int diskn = lnn % number_of_raids;
    Create_IO_Filename(name, lnn, numproc);
    name = disk_names[diskn] + "/" + subdir + "/" + name;
  }
  else {
    Create_IO_Filename(name, lnn, numproc);
    if (subdir.length() != 0u) {
      name = subdir + "/" + name;
    }
  }
  name = prepend + name;
}

/*****************************************************************************/
void Excn::ParallelDisks::create_disk_names()
{

  if (number_of_raids == 0) {
    return;
  }

  disk_names.resize(number_of_raids);
  for (int i = 0; i < number_of_raids; i++) {
    int num = i + raid_offset;
    if (num < 10) {
#ifdef COUGAR
      disk_names[i] = std::to_string(num);
#else
      disk_names[i] = "0" + std::to_string(num);
#endif
    }
    else {
      disk_names[i] = std::to_string(num);
    }
  }
}

/*****************************************************************************/
void Excn::ParallelDisks::Create_IO_Filename(std::string &name, int processor, int num_processors)
{
  // Current format for per-processor file names is:
  // PREFIX/basename.num_proc.cur_proc
  // the 'cur_proc' field is padded to be the same width as
  // the 'num_proc' field
  // Examples: basename.8.1, basename.64.03, basename.128.001

  // Create a std::string containing the total number of processors
  std::string num_proc   = std::to_string(num_processors);
  size_t      proc_width = num_proc.length();

  // Create a std::string containing the current processor number
  std::string cur_proc  = std::to_string(processor);
  size_t      cur_width = cur_proc.length();

  // Build the filename
  name += ".";
  name += num_proc;
  name += ".";

  // Now, pad with zeros so that 'cur_proc' portion is same
  // width as 'num_proc' portion.
  while (cur_width++ < proc_width) {
    name += "0";
  }

  name += cur_proc;
}
