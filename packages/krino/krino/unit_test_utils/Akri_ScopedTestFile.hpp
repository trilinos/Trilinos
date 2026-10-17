// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_UNIT_TEST_UTILS_AKRI_SCOPEDTESTFILE_HPP_
#define KRINO_KRINO_UNIT_TEST_UTILS_AKRI_SCOPEDTESTFILE_HPP_

#include <mpi.h>

#include <string>

namespace krino {

// Removes a file written by a unit test when it goes out of scope, so tests leave no files behind,
// even if an assertion fails partway through.  Declare it before the file is written.  Every rank
// waits until all ranks are done with the file, then each removes what it may have written: the
// file itself, and its per-rank "<name>.<numProcs>.<rank>" piece when run in parallel.
class ScopedTestFile
{
public:
  ScopedTestFile(const std::string & fileName, MPI_Comm comm = MPI_COMM_WORLD);
  ~ScopedTestFile();
  ScopedTestFile(const ScopedTestFile &) = delete;
  ScopedTestFile & operator=(const ScopedTestFile &) = delete;

  const std::string & name() const { return myFileName; }

private:
  std::string myFileName;
  MPI_Comm myComm;
};

} // namespace krino

#endif /* KRINO_KRINO_UNIT_TEST_UTILS_AKRI_SCOPEDTESTFILE_HPP_ */
