// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_Util.cpp
 *  \brief Useful namespace methods.
 */

#include <Zoltan2_Util.hpp>
#include <Teuchos_OpaqueWrapper.hpp>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include <unistd.h>

namespace Zoltan2{

/* On a linux node, find the total memory currently allocated
 * to this process.
 * Return the number of kilobytes allocated to this process.
 * Return 0 if it is not possible to determine this.
 */
long getProcessKilobytes()
{
long pageSize;

#ifdef _SC_PAGESIZE
  pageSize = sysconf(_SC_PAGESIZE);
#else
#warning "Page size query is not possible.  No per-process memory stats."
  return 0;
#endif

  pid_t pid = getpid();
  std::ostringstream fname;
  fname << "/proc/" << pid << "/statm";
  std::ifstream memFile;

  try{
    memFile.open(fname.str().c_str());
  }
  catch (...){
    return 0;
  }

  char buf[128];
  memset(buf, 0, 128);
  while (memFile.good()){
    memFile.getline(buf, 128);
    break;
  }

  memFile.close();

  std::istringstream sbuf(buf);
  long totalPages;
  sbuf >> totalPages;

  long pageKBytes = pageSize / 1024;
  totalPages = atol(buf);

  return totalPages * pageKBytes;
}



#ifdef HAVE_ZOLTAN2_MPI

/*! \brief Convert an MPI communicator to a MpiComm object.
 */

RCP<Teuchos::MpiComm<int> >
  MPI2Teuchos(const MPI_Comm &comm)
{
  typedef Teuchos::OpaqueWrapper<MPI_Comm> mpiWrapper_t;
  RCP<mpiWrapper_t>handle = Teuchos::opaqueWrapper<MPI_Comm>(comm);
  RCP<Teuchos::MpiComm<int> > tcommPtr(
    new Teuchos::MpiComm<int>(handle));

  return tcommPtr;
}

RCP<const Teuchos::MpiComm<int> >
  MPI2TeuchosConst(const MPI_Comm &comm)
{
  typedef Teuchos::OpaqueWrapper<MPI_Comm> mpiWrapper_t;
  RCP<mpiWrapper_t>handle = Teuchos::opaqueWrapper<MPI_Comm>(comm);
  RCP<const Teuchos::MpiComm<int> > tcommPtr(
    new Teuchos::MpiComm<int>(handle));

  return tcommPtr;
}

/*! \brief Convert a Teuchos::MpiComm object to the underlying MPI_Comm.
 */

MPI_Comm
  Teuchos2MPI(const RCP<Comm<int> > &comm)
{
  MPI_Comm mpiComm;
  typedef Teuchos::OpaqueWrapper<MPI_Comm> mpiWrapper_t;

  Comm<int> *c = comm.get();
  Teuchos::MpiComm<int> *mc = dynamic_cast<Teuchos::MpiComm<int> *>(c);
  if (mc){
    RCP<const mpiWrapper_t> wrappedComm = mc->getRawMpiComm();
    mpiComm = (*wrappedComm.getRawPtr())();
  }
  else{
    mpiComm = MPI_COMM_SELF;   // or would this be an error?
  }

  return mpiComm;
}

MPI_Comm
  TeuchosConst2MPI(const RCP<const Comm<int> > &comm)
{
  MPI_Comm mpiComm;
  typedef Teuchos::OpaqueWrapper<MPI_Comm> mpiWrapper_t;

  const Comm<int> *cConst = comm.get();
  Comm<int> *c = const_cast<Comm<int> *>(cConst);
  Teuchos::MpiComm<int> *mc = dynamic_cast<Teuchos::MpiComm<int> *>(c);
  if (mc){
    RCP<const mpiWrapper_t> wrappedComm = mc->getRawMpiComm();
    mpiComm = (*wrappedComm.getRawPtr())();
  }
  else{
    mpiComm = MPI_COMM_SELF;   // or would this be an error?
  }

  return mpiComm;
}

#endif

} // namespace Zoltan2
