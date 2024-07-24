// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_TRANSFER_DEF_HPP
#define TPETRA_DETAILS_TRANSFER_DEF_HPP

#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_Distributor.hpp"
#include "Tpetra_ImportExportData.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include <sstream>

namespace { // (anonymous)

  // Assume that dv is sync'd.
  template<class ElementType, class DeviceType>
  Teuchos::ArrayView<const ElementType>
  makeConstArrayViewFromDualView (const Kokkos::DualView<ElementType*, DeviceType>& dv)
  {
    TEUCHOS_ASSERT( ! dv.need_sync_host () );
    auto hostView = dv.view_host ();
    const auto size = hostView.extent (0);
    return Teuchos::ArrayView<const ElementType> (size == 0 ? nullptr : hostView.data (), size);
  }

  template<class DeviceType, class LocalOrdinal>
    struct OrderedViewFunctor {
      OrderedViewFunctor (const Kokkos::View<LocalOrdinal*, DeviceType>& viewToCheck) :
        viewToCheck_ (viewToCheck) {}
      KOKKOS_INLINE_FUNCTION void operator() (const size_t i, unsigned int& isUnordered) const {
        isUnordered |= static_cast<unsigned int>(viewToCheck_(i)+1 != viewToCheck_(i+1));
      }
      Kokkos::View<const LocalOrdinal*, DeviceType> viewToCheck_;
    };

    template<class DeviceType, class LocalOrdinal>
    bool
    isViewOrdered (const Kokkos::View<LocalOrdinal*, DeviceType>& viewToCheck)
    {
      using Kokkos::parallel_reduce;
      typedef DeviceType DT;
      typedef typename DT::execution_space DES;
      typedef Kokkos::RangePolicy<DES, size_t> range_type;

      const size_t size = viewToCheck.extent (0);
      unsigned int isUnordered = 0;
      if (size>1)
        parallel_reduce ("isViewOrdered",
                         range_type (0, size-1),
                         OrderedViewFunctor<DeviceType, LocalOrdinal> (viewToCheck),
                         isUnordered);
      return isUnordered == 0;
    }

} // namespace (anonymous)

namespace Tpetra {
namespace Details {

template <class LO, class GO, class NT>
Transfer<LO, GO, NT>::
Transfer (const Teuchos::RCP<const map_type>& source,
	  const Teuchos::RCP<const map_type>& target,
          const Teuchos::RCP<Teuchos::FancyOStream>& out,
          const Teuchos::RCP<Teuchos::ParameterList>& plist,
	  const std::string& className) :
  TransferData_ (new ImportExportData<LO, GO, NT> (source, target, out, plist))
{
  TEUCHOS_ASSERT( ! TransferData_->out_.is_null () );
  this->setParameterList (plist, className);
}

template <class LO, class GO, class NT>
Transfer<LO, GO, NT>::
Transfer (const Transfer<LO, GO, NT>& rhs, reverse_tag)
{
  TEUCHOS_ASSERT( ! (rhs.TransferData_).is_null () );
  this->TransferData_ = rhs.TransferData_->reverseClone ();
  TEUCHOS_ASSERT( ! this->TransferData_->out_.is_null () );  
}

template <class LO, class GO, class NT>
void
Transfer<LO, GO, NT>::
setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist,
		  const std::string& className)
{
  using ::Tpetra::Details::Behavior;

  const bool verboseEnv = Behavior::verbose (className.c_str ()) ||
    Behavior::verbose ((std::string ("Tpetra::") + className).c_str ());
  
  bool verboseParam = false;
  if (! plist.is_null ()) {
    // FIXME (mfh 03 Feb 2019) Phase out these parameters in favor of
    // TPETRA_VERBOSE.
    if (plist->isType<bool> ("Verbose")) {
      verboseParam = plist->get<bool> ("Verbose");
    }
    else if (plist->isType<bool> ("Debug")) { // backwards compat
      verboseParam = plist->get<bool> ("Debug");
    }
  }
  this->TransferData_->verbose_ = verboseEnv || verboseParam;
}

template <class LO, class GO, class NT>
size_t
Transfer<LO, GO, NT>::
getNumSameIDs () const {
  return TransferData_->numSameIDs_;
}

template <class LO, class GO, class NT>
size_t
Transfer<LO, GO, NT>::
getNumPermuteIDs () const {
  return static_cast<size_t> (TransferData_->permuteFromLIDs_.extent (0));
}

template <class LO, class GO, class NT>  
Kokkos::DualView<const LO*, typename Transfer<LO, GO, NT>::device_type>
Transfer<LO, GO, NT>::
getPermuteFromLIDs_dv () const {
  const auto& dv = TransferData_->permuteFromLIDs_;
  TEUCHOS_TEST_FOR_EXCEPTION
    (dv.need_sync_device (), std::logic_error,
     "Tpetra::Details::Transfer::getPermuteFromLIDs_dv: "
     "DualView needs sync to device" );  
  TEUCHOS_TEST_FOR_EXCEPTION
    (dv.need_sync_host (), std::logic_error,
     "Tpetra::Details::Transfer::getPermuteFromLIDs_dv: "
     "DualView needs sync to host" );
  return dv;
}
  
template <class LO, class GO, class NT>  
Teuchos::ArrayView<const LO>
Transfer<LO, GO, NT>::
getPermuteFromLIDs () const {
  return makeConstArrayViewFromDualView (TransferData_->permuteFromLIDs_);    
}

template <class LO, class GO, class NT>  
Kokkos::DualView<const LO*, typename Transfer<LO, GO, NT>::device_type>
Transfer<LO, GO, NT>::
getPermuteToLIDs_dv () const {
  const auto& dv = TransferData_->permuteToLIDs_;
  TEUCHOS_TEST_FOR_EXCEPTION
    (dv.need_sync_device (), std::logic_error,
     "Tpetra::Details::Transfer::getPermuteToLIDs_dv: "
     "DualView needs sync to device" );  
  TEUCHOS_TEST_FOR_EXCEPTION
    (dv.need_sync_host (), std::logic_error,
     "Tpetra::Details::Transfer::getPermuteToLIDs_dv: "
     "DualView needs sync to host" );
  return dv;
}
  
template <class LO, class GO, class NT>  
Teuchos::ArrayView<const LO>
Transfer<LO, GO, NT>::
getPermuteToLIDs () const {
  return makeConstArrayViewFromDualView (TransferData_->permuteToLIDs_);
}

template <class LO, class GO, class NT>  
size_t
Transfer<LO, GO, NT>::
getNumRemoteIDs () const {
  return static_cast<size_t> (TransferData_->remoteLIDs_.extent (0));
}

template <class LO, class GO, class NT>  
Kokkos::DualView<const LO*, typename Transfer<LO, GO, NT>::device_type>
Transfer<LO, GO, NT>::
getRemoteLIDs_dv () const {
  const auto& dv = TransferData_->remoteLIDs_;
  TEUCHOS_TEST_FOR_EXCEPTION
    (dv.need_sync_device (), std::logic_error,
     "Tpetra::Details::Transfer::getRemoteLIDs_dv: "
     "DualView needs sync to device" );  
  TEUCHOS_TEST_FOR_EXCEPTION
    (dv.need_sync_host (), std::logic_error,
     "Tpetra::Details::Transfer::getRemoteLIDs_dv: "
     "DualView needs sync to host" );
  return dv;
}
  
template <class LO, class GO, class NT>  
Teuchos::ArrayView<const LO>
Transfer<LO, GO, NT>::
getRemoteLIDs () const {
  return makeConstArrayViewFromDualView (TransferData_->remoteLIDs_);
}

template <class LO, class GO, class NT>
size_t
Transfer<LO, GO, NT>::
getNumExportIDs () const {
  return static_cast<size_t> (TransferData_->exportLIDs_.extent (0));
}

template <class LO, class GO, class NT>  
Kokkos::DualView<const LO*, typename Transfer<LO, GO, NT>::device_type>
Transfer<LO, GO, NT>::
getExportLIDs_dv () const {
  const auto& dv = TransferData_->exportLIDs_;
  TEUCHOS_TEST_FOR_EXCEPTION
    (dv.need_sync_device (), std::logic_error,
     "Tpetra::Details::Transfer::getExportLIDs_dv: "
     "DualView needs sync to device" );  
  TEUCHOS_TEST_FOR_EXCEPTION
    (dv.need_sync_host (), std::logic_error,
     "Tpetra::Details::Transfer::getExportLIDs_dv: "
     "DualView needs sync to host" );
  return dv;
}
  
template <class LO, class GO, class NT>
Teuchos::ArrayView<const LO>
Transfer<LO, GO, NT>::
getExportLIDs () const {
  return makeConstArrayViewFromDualView (TransferData_->exportLIDs_);
}

template <class LO, class GO, class NT>
Teuchos::ArrayView<const int>
Transfer<LO, GO, NT>::
getExportPIDs () const {
  return TransferData_->exportPIDs_ ();
}

template <class LO, class GO, class NT>
Teuchos::RCP<const typename Transfer<LO, GO, NT>::map_type>
Transfer<LO, GO, NT>::
getSourceMap () const {
  return TransferData_->source_;
}

template <class LO, class GO, class NT>
Teuchos::RCP<const typename Transfer<LO, GO, NT>::map_type>
Transfer<LO, GO, NT>::
getTargetMap () const {
  return TransferData_->target_;
}

template <class LO, class GO, class NT>
::Tpetra::Distributor&
Transfer<LO, GO, NT>::
getDistributor () const {
  return TransferData_->distributor_;
}

template <class LO, class GO, class NT>
bool
Transfer<LO, GO, NT>::
isLocallyComplete () const {
  return TransferData_->isLocallyComplete_;
}

template <class LO, class GO, class NT>
bool
Transfer<LO, GO, NT>::
isLocallyFitted () const {
  return (getNumSameIDs() == std::min(getSourceMap()->getLocalNumElements(),
                                      getTargetMap()->getLocalNumElements()));
}

template <class LO, class GO, class NT>
void
Transfer<LO, GO, NT>::
detectRemoteExportLIDsContiguous () const {

  // Check that maps are locally fitted
  // TODO: We really want to check here that remote LIDs are sorted last.
  //       The current check is too restrictive in special cases.
  bool ordered = (getNumSameIDs() == std::min(getSourceMap()->getLocalNumElements(),
                                              getTargetMap()->getLocalNumElements()));
  ordered &= (getTargetMap()->getLocalNumElements() == getNumSameIDs() + getNumRemoteIDs());
  if (ordered) {
    const auto& dv = TransferData_->remoteLIDs_;
    TEUCHOS_TEST_FOR_EXCEPTION
      (dv.need_sync_device (), std::logic_error,
       "Tpetra::Details::Transfer::getRemoteLIDs_dv: "
       "DualView needs sync to device" );
    auto v_d = dv.view_device ();
    ordered &= isViewOrdered<device_type, LO>(v_d);
  }
  TransferData_->remoteLIDsContiguous_ = ordered;

  ordered = (getNumSameIDs() == std::min(getSourceMap()->getLocalNumElements(),
                                         getTargetMap()->getLocalNumElements()));
  ordered &= (getSourceMap()->getLocalNumElements() == getNumSameIDs() + getNumExportIDs());
  if (ordered) {
    const auto& dv = TransferData_->exportLIDs_;
    TEUCHOS_TEST_FOR_EXCEPTION
      (dv.need_sync_device (), std::logic_error,
       "Tpetra::Details::Transfer::getRemoteLIDs_dv: "
       "DualView needs sync to device" );
    auto v_d = dv.view_device ();
    ordered &= isViewOrdered<device_type, LO>(v_d);
  }
  TransferData_->exportLIDsContiguous_ = ordered;
}

template <class LO, class GO, class NT>
bool
Transfer<LO, GO, NT>::
areRemoteLIDsContiguous () const {
  return TransferData_->remoteLIDsContiguous_;
}

template <class LO, class GO, class NT>
bool
Transfer<LO, GO, NT>::
areExportLIDsContiguous () const {
  return TransferData_->exportLIDsContiguous_;
}


template <class LO, class GO, class NT>
void
Transfer<LO, GO, NT>::
describe (Teuchos::FancyOStream& out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  this->describeImpl (out, "Tpetra::Details::Transfer", verbLevel);
}

template<class LO, class GO, class NT>
Teuchos::FancyOStream&
Transfer<LO, GO, NT>::
verboseOutputStream () const
{
  Teuchos::FancyOStream* outPtr = TransferData_->out_.getRawPtr ();
  TEUCHOS_ASSERT( outPtr != nullptr );
  return *outPtr;
}

template<class LO, class GO, class NT>
bool
Transfer<LO, GO, NT>::
verbose () const {
  return TransferData_->verbose_;
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


template <class LO, class GO, class NT>
void
expertSetRemoteLIDsContiguous(Transfer<LO, GO, NT> transfer, bool contig) {
  transfer.TransferData_->remoteLIDsContiguous_ = contig;
}


template <class LO, class GO, class NT>
void
expertSetExportLIDsContiguous(Transfer<LO, GO, NT> transfer, bool contig) {
  transfer.TransferData_->exportLIDsContiguous_ = contig;
}


} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_TRANSFER_DEF_HPP
