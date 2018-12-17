// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_TRANSFER_DEF_HPP
#define TPETRA_DETAILS_TRANSFER_DEF_HPP

#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_Distributor.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include <sstream>

namespace Tpetra {
namespace Details {

template <class LO, class GO, class NT>
void
Transfer<LO, GO, NT>::
describe (Teuchos::FancyOStream& out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  this->describeImpl (out, "Tpetra::Details::Transfer", verbLevel);
}

template<class LO, class GO, class NT>
void
Transfer<LO, GO, NT>::
describeImpl (Teuchos::FancyOStream& out,
              const std::string& className,
              const Teuchos::EVerbosityLevel verbLevel) const
{
  using Teuchos::TypeNameTraits;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using std::endl;
  const Teuchos::EVerbosityLevel vl =
    (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

  if (vl == VERB_NONE) {
    return; // don't print anything
  }
  // If this Transfer's source Map or Comm is null, then the Transfer
  // does not participate in collective operations with the other
  // processes.  In that case, it is not even legal to call this
  // method.  The reasonable thing to do in that case is nothing.
  auto srcMap = this->getSourceMap ();
  if (srcMap.is_null ()) {
    return;
  }
  auto comm = srcMap->getComm ();
  if (comm.is_null ()) {
    return;
  }
  if (this->getTargetMap ().is_null () ||
      this->getTargetMap ()->getComm ().is_null ()) {
    return;
  }

  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  // Only Process 0 should touch the output stream, but this method in
  // general may need to do communication.  Thus, we may need to
  // preserve the current tab level across multiple "if (myRank == 0)
  // { ... }" inner scopes.  This is why we sometimes create OSTab
  // instances by pointer, instead of by value.
  Teuchos::RCP<Teuchos::OSTab> tab0, tab1;

  if (myRank == 0) {
    // At every verbosity level but VERB_NONE, Process 0 prints.
    // By convention, describe() always begins with a tab before
    // printing.
    tab0 = Teuchos::rcp (new Teuchos::OSTab (out));

    out << "\"" << className << "\":" << endl;
    tab1 = Teuchos::rcp (new Teuchos::OSTab (out));

    {
      out << "Template parameters:" << endl;
      Teuchos::OSTab tab2 (out);
      out << "LocalOrdinal: " << TypeNameTraits<LO>::name () << endl
          << "GlobalOrdinal: " << TypeNameTraits<GO>::name () << endl
          << "Node: " << TypeNameTraits<NT>::name () << endl;
    }

    const std::string label = this->getObjectLabel ();
    if (label != "") {
      out << "Label: " << label << endl;
    }
    out << "Number of processes: " << numProcs << endl;
  }

  if (vl > VERB_LOW) {
    // At higher verbosity levels, describe() is allowed to
    // communicate in order to print information from other
    // processes in the object's communicator.
    this->globalDescribe (out, vl);
  }

  // It's illegal to call describe() on a process where either Map is
  // null.  (That implies the process in question is not participating
  // in collective operations with either Map, and describe is
  // collective over the Maps' communicator.)  Thus, we don't have to
  // define behavior when either Map is NULL on any process.  Thus,
  // it's OK that the code below isn't quite right (that is, won't
  // print anything) if either Map is NULL on Process 0.

  if (myRank == 0) {
    out << "Source Map:" << endl;
  }
  // This is collective over the Map's communicator.
  this->getSourceMap ()->describe (out, vl);

  if (myRank == 0) {
    out << "Target Map:" << endl;
  }
  // This is collective over the Map's communicator.
  this->getTargetMap ()->describe (out, vl);

  if (myRank == 0) {
    out << "Distributor:" << endl;
  }
  this->getDistributor ().describe (out, vl);
}

template<class LO, class GO, class NT>
void
Transfer<LO, GO, NT>::
globalDescribe (Teuchos::FancyOStream& out,
                const Teuchos::EVerbosityLevel vl) const
{
  using Teuchos::Comm;
  using Teuchos::OSTab;
  using Teuchos::RCP;
  using Teuchos::toString;
  using std::endl;

  // If this Transfer's source Map or Comm is null, then the Transfer
  // does not participate in collective operations with the other
  // processes.  In that case, it is not even legal to call this
  // method.  The reasonable thing to do in that case is nothing.
  auto srcMap = this->getSourceMap ();
  if (srcMap.is_null ()) {
    return;
  }
  RCP<const Teuchos::Comm<int> > comm = srcMap->getComm ();
  if (comm.is_null ()) {
    return;
  }

  const std::string myStr = localDescribeToString (vl);
  ::Tpetra::Details::gathervPrint (out, myStr, *comm);
}

template<class LO, class GO, class NT>
std::string
Transfer<LO, GO, NT>::
localDescribeToString (const Teuchos::EVerbosityLevel vl) const
{
  using Teuchos::OSTab;
  using Teuchos::RCP;
  using std::endl;

  RCP<std::ostringstream> outString (new std::ostringstream);
  RCP<Teuchos::FancyOStream> outp = Teuchos::getFancyOStream (outString);
  Teuchos::FancyOStream& out = *outp; // only valid during this scope

  RCP<const Teuchos::Comm<int> > comm = this->getSourceMap ()->getComm ();
  if (this->getSourceMap ().is_null () ||
      this->getSourceMap ()->getComm ().is_null ()) {
    // If this Transfer does not participate in the communicator,
    // it's not even legal to call this method.  However, we need to
    // do something in this case.  The reasonable thing to do is not
    // to print anything.
    return std::string ("");
  }
  else {
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    out << "Process " << myRank << " of " << numProcs << ":" << endl;
    OSTab tab1 (out);

    out << "numSameIDs: " << getNumSameIDs () << endl;
    out << "numPermuteIDs: " << getNumPermuteIDs () << endl;
    out << "numRemoteIDs: " << getNumRemoteIDs () << endl;
    out << "numExportIDs: " << getNumExportIDs () << endl;

    // Only print the actual contents of these arrays at the two
    // highest verbosity levels.  Otherwise, just print their counts.
    if (vl <= Teuchos::VERB_MEDIUM) {
      out << "permuteFromLIDs count: " << getPermuteFromLIDs ().size () << endl
          << "permuteToLIDs count: " << getPermuteToLIDs ().size () << endl
          << "remoteLIDs count: " << getRemoteLIDs ().size () << endl
          << "exportLIDs count: " << getExportLIDs ().size () << endl
          << "exportPIDs count: " << getExportPIDs () << endl;
    }
    else { // vl = VERB_HIGH or VERB_EXTREME
      // Build RemoteGIDs
      RCP<const Map<LO,GO,NT> > tmap = getTargetMap();
      RCP<const Map<LO,GO,NT> > smap = getSourceMap();
      Teuchos::Array<GO>  RemoteGIDs(getRemoteLIDs().size());
      Teuchos::Array<int> RemotePIDs(getRemoteLIDs().size());
      for(size_t i=0; i<(size_t)getRemoteLIDs().size(); i++)
        RemoteGIDs[i] = tmap->getGlobalElement(getRemoteLIDs()[i]);

      Teuchos::Array<int> ExportGIDs(getExportLIDs().size());
      for(size_t i=0; i<(size_t)getExportLIDs().size(); i++)
        ExportGIDs[i] = smap->getGlobalElement(getExportLIDs()[i]);

      // Build RemotePIDs (taken from Tpetra_Import_Util.hpp)
      const Tpetra::Distributor & D=getDistributor();
      size_t NumReceives                           = D.getNumReceives();
      Teuchos::ArrayView<const int> ProcsFrom      = D.getProcsFrom();
      Teuchos::ArrayView<const size_t> LengthsFrom = D.getLengthsFrom();
      for (size_t i = 0, j = 0; i < NumReceives; ++i) {
        const int pid = ProcsFrom[i];
        for (size_t k = 0; k < LengthsFrom[i]; ++k) {
          RemotePIDs[j] = pid;
          j++;
        }
      }

      out << "distor.NumRecvs   : "<<NumReceives<<endl
          << "distor.ProcsFrom  : "<<toString(ProcsFrom)<<endl
          << "distor.LengthsFrom: "<<toString(LengthsFrom)<<endl;

      out << "distor.NumSends   : "<<D.getNumSends()<<endl
          << "distor.ProcsTo    : "<<toString(D.getProcsTo())<<endl
          << "distor.LengthsTo  : "<<toString(D.getLengthsTo())<<endl;

      out << "distor.hasSelfMsg : "<<D.hasSelfMessage()<<endl;

      out << "permuteFromLIDs: " << toString (getPermuteFromLIDs ()) << endl
          << "permuteToLIDs: " << toString (getPermuteToLIDs ()) << endl
          << "remoteLIDs: " << toString (getRemoteLIDs ()) << endl
          << "remoteGIDs: " << toString (RemoteGIDs ()) << endl
          << "remotePIDs: " << toString (RemotePIDs ()) << endl
          << "exportLIDs: " << toString (getExportLIDs ()) << endl
          << "exportGIDs: " << toString (ExportGIDs ()) << endl
          << "exportPIDs: " << toString (getExportPIDs ()) << endl;
    }

    out.flush (); // make sure the ostringstream got everything
    return outString->str ();
  }
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_TRANSFER_DEF_HPP
