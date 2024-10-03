// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_Details_gathervPrint.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <algorithm>

namespace Tpetra {
namespace Details {

void
gathervPrint (std::ostream& out,
              const std::string& s,
              const Teuchos::Comm<int>& comm)
{
  using Teuchos::ArrayRCP;
  using Teuchos::CommRequest;
  using Teuchos::ireceive;
  using Teuchos::isend;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::wait;

  const int myRank = comm.getRank ();
  const int rootRank = 0;
  if (myRank == rootRank) {
    out << s; // Proc 0 prints its buffer first
  }

  const int numProcs = comm.getSize ();
  const int sizeTag = 42;
  const int msgTag = 43;

  ArrayRCP<size_t> sizeBuf (1);
  ArrayRCP<char> msgBuf; // to be resized later
  RCP<CommRequest<int> > req;

  for (int p = 1; p < numProcs; ++p) {
    if (myRank == p) {
      sizeBuf[0] = s.size ();
      req = isend<int, size_t> (sizeBuf, rootRank, sizeTag, comm);
      (void) wait<int> (comm, outArg (req));

      const size_t msgSize = s.size ();
      msgBuf.resize (msgSize + 1); // for the '\0'
      std::copy (s.begin (), s.end (), msgBuf.begin ());
      msgBuf[msgSize] = '\0';

      req = isend<int, char> (msgBuf, rootRank, msgTag, comm);
      (void) wait<int> (comm, outArg (req));
    }
    else if (myRank == rootRank) {
      sizeBuf[0] = 0; // just a precaution
      req = ireceive<int, size_t> (sizeBuf, p, sizeTag, comm);
      (void) wait<int> (comm, outArg (req));

      const size_t msgSize = sizeBuf[0];
      msgBuf.resize (msgSize + 1); // for the '\0'
      req = ireceive<int, char> (msgBuf, p, msgTag, comm);
      (void) wait<int> (comm, outArg (req));

      std::string msg (msgBuf.getRawPtr ());
      out << msg;
    }
  }
}

} // namespace Details
} // namespace Tpetra
