/*
 * Copyright(C) 1999-2021, 2023 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <cstdlib>

#include <string>

#include <EP_Internals.h>
#include <EP_ParallelDisks.h>

#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(_MSC_VER) ||                \
    defined(__MINGW32__) || defined(_WIN64) || defined(__MINGW64__)
#include <Shlwapi.h>
#endif

/*****************************************************************************/
void Excn::ParallelDisks::rename_file_for_mp(const std::string &rootdir, const std::string &subdir,
                                             std::string &name, int node, int numproc) const
{
  // Possible to have node layout without parallel disks

  std::string prepend;
  if (!rootdir.empty()) {
    prepend = rootdir + "/";
  }
  else if (Excn::is_path_absolute(name)) {
    prepend = "";
  }
  else {
    prepend = "./";
  }

  int lnn = node;
  Create_IO_Filename(name, lnn, numproc);
  if (!subdir.empty()) {
    name = subdir + "/" + name;
  }
  name = prepend + name;
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
