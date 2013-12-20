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

#ifndef TPETRA_IMPORT_UTIL_HPP
#define TPETRA_IMPORT_UTIL_HPP

/*!
  \file Tpetra_Import_Util.hpp
  \brief Utility functions and macros designed for use with Tpetra::Import and Tpetra::Export objects.
*/

#include "Tpetra_ConfigDefs.hpp" // for map, vector, string, and iostream
#include "Tpetra_Import.hpp"
#include "Tpetra_Distributor.hpp"
#include <Teuchos_Array.hpp>
#include <utility>

namespace Tpetra {

  namespace Import_Util {
    //! Tpetra::Import_Util::getPidGidPairs function
    /*!  For each GID in the TargetMap, find who owns the GID in the SourceMap. 
      This works entirely from the Distributor and has no communication at all.  
      
      The routine returns (by reference) a Teuchos::Array of std::pair<int,GlobalOrdinal> which contains (PID,GID) pairs.
      If the use_minus_one_for_local==true, any GIDs owned by this processor get -1 instead of their PID.
    */
    template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
    void getPidGidPairs(const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> & Importer, Teuchos::Array< std::pair<int,GlobalOrdinal> > & gpids, bool use_minus_one_for_local);
    
    //! Tpetra::Import_Util::getPids function
    /*! Like getPidGidPairs, but just gets the PIDs, ordered by the columnmap 
     */
    template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
    void getPids(const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> & Importer, Teuchos::Array<int> &pids, bool use_minus_one_for_local);

    //! Tpetra::Import_Util::getRemotePIDs
    /*! Gets a list of remote PIDs from an importer in the order corresponding to the RemoteLIDs
     */
    template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
    void getRemotePIDs(const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> & Importer, Teuchos::Array<int> &RemotePIDs);


    //! packAndPrepareWithOwningPIDs.  
    /*! Note: The SourcePids vector should contain a list of owning PIDs for each column in the ColMap, as from Epetra_Util::GetPids,
      without the "-1 for local" option being used.  This routine is basically Tpetra::CrsMatrix::packAndPrepare, but it
      packs the owning PIDs as well as the GIDs.      
      
      \warning This method is intended for expert developer use only, and should never be called by user code.
    */
    template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, typename LocalMatOps>
    void packAndPrepareWithOwningPIDs(const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & SourceMatrix,
				      const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
				      Teuchos::Array<char>& exports,
				      const Teuchos::ArrayView<size_t>& numPacketsPerLID,
				      size_t& constantNumPackets,
				      Distributor &distor,
				      const Teuchos::ArrayView<int>& SourcePids);    
  }// end Import_Util
}//end Tpetra




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//----------------------------------------------------------------------------
template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void Tpetra::Import_Util::getPidGidPairs(const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> & Importer, Teuchos::Array< std::pair<int,GlobalOrdinal> > & gpids, bool use_minus_one_for_local) {
  // Put the (PID,GID) pair in member of Importer.TargetMap() in gpids.  If use_minus_one_for_local==true, put in -1 instead of MyPID.
  const Tpetra::Distributor & D=Importer.getDistributor();

  LocalOrdinal ii;
  size_t  i,j,k;
  int mypid = Importer.getTargetMap()->getComm()->getRank();
  size_t N  = Importer.getTargetMap()->getNodeNumElements();

  // Get the importer's data
  Teuchos::ArrayView<const LocalOrdinal> RemoteLIDs  = Importer.getRemoteLIDs();

  // Get the distributor's data
  size_t NumReceives                           = D.getNumReceives();
  Teuchos::ArrayView<const int> ProcsFrom      = D.getImagesFrom();
  Teuchos::ArrayView<const size_t> LengthsFrom = D.getLengthsFrom();

  // Resize the outgoing data structure
  gpids.resize(N);

  // Start by claiming that I own all the data
  LocalOrdinal lzero = Teuchos::ScalarTraits<LocalOrdinal>::zero();
  if(use_minus_one_for_local)
    for(ii=lzero; Teuchos::as<size_t>(ii)<N; ii++) gpids[ii]=std::make_pair(-1,Importer.getTargetMap()->getGlobalElement(ii));
  else
    for(ii=lzero; Teuchos::as<size_t>(ii)<N; ii++) gpids[ii]=std::make_pair(mypid,Importer.getTargetMap()->getGlobalElement(ii));

  // Now, for each remote ID, record who actually owns it.  This loop follows the operation order in the
  // MpiDistributor so it ought to duplicate that effect.
  for(i=0,j=0; i<NumReceives; i++){
    int pid=ProcsFrom[i];
    for(k=0; k<LengthsFrom[i]; k++){
      if(pid!=mypid) gpids[RemoteLIDs[j]].first=pid;
      j++;
    }    
  }
}




//----------------------------------------------------------------------------
template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void Tpetra::Import_Util::getPids(const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> & Importer, Teuchos::Array<int> &pids, bool use_minus_one_for_local) {
  const Tpetra::Distributor & D=Importer.getDistributor();

  LocalOrdinal ii;
  size_t  i,j,k;
  int mypid = Importer.getTargetMap()->getComm()->getRank();
  size_t N  = Importer.getTargetMap()->getNodeNumElements();

  // Get the importer's data
  Teuchos::ArrayView<const LocalOrdinal> RemoteLIDs  = Importer.getRemoteLIDs();

  // Get the distributor's data
  size_t NumReceives                           = D.getNumReceives();
  Teuchos::ArrayView<const int> ProcsFrom      = D.getImagesFrom();
  Teuchos::ArrayView<const size_t> LengthsFrom = D.getLengthsFrom();

  // Resize the outgoing data structure
  pids.resize(N);

  // Start by claiming that I own all the data
  LocalOrdinal lzero = Teuchos::ScalarTraits<LocalOrdinal>::zero();
  if(use_minus_one_for_local)
    for(ii=lzero; Teuchos::as<size_t>(ii)<N; ii++) pids[ii]=-1;
  else
    for(ii=lzero; Teuchos::as<size_t>(ii)<N; ii++) pids[ii]=mypid;

  // Now, for each remote ID, record who actually owns it.  This loop follows the operation order in the
  // MpiDistributor so it ought to duplicate that effect.
  for(i=0,j=0; i<NumReceives; i++){
    int pid=ProcsFrom[i];
    for(k=0; k<LengthsFrom[i]; k++){
      if(pid!=mypid) pids[RemoteLIDs[j]]=pid;
      j++;
    }    
  }
}


//----------------------------------------------------------------------------
template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void Tpetra::Import_Util::getRemotePIDs(const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> & Importer, Teuchos::Array<int> &RemotePIDs) {
  const Tpetra::Distributor & D=Importer.getDistributor();

  size_t  i,j,k;

  // Get the importer's data
  Teuchos::ArrayView<const LocalOrdinal> RemoteLIDs  = Importer.getRemoteLIDs();

  // Get the distributor's data
  size_t NumReceives                           = D.getNumReceives();
  Teuchos::ArrayView<const int> ProcsFrom      = D.getImagesFrom();
  Teuchos::ArrayView<const size_t> LengthsFrom = D.getLengthsFrom();

  // Resize the outgoing data structure
  RemotePIDs.resize(Importer.getNumRemoteIDs());

  // Now, for each remote ID, record who actually owns it.  This loop follows the operation order in the
  // MpiDistributor so it ought to duplicate that effect.
  for(i=0,j=0; i<NumReceives; i++){
    int pid=ProcsFrom[i];    
    for(k=0; k<LengthsFrom[i]; k++){
      RemotePIDs[j]=pid;
      j++;
    }    
  }
}


//----------------------------------------------------------------------------
template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, typename LocalMatOps>
void Tpetra::Import_Util::packAndPrepareWithOwningPIDs(const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> &SourceMatrix,
						       const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
						       Teuchos::Array<char>& exports,
						       const Teuchos::ArrayView<size_t>& numPacketsPerLID,
						       size_t& constantNumPackets,
						       Distributor &distor,
						       const Teuchos::ArrayView<int>& SourcePids) {

  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::av_reinterpret_cast;
  using Teuchos::RCP;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Map<LocalOrdinal,GlobalOrdinal,Node>  map_type;
  typedef typename ArrayView<const LO>::size_type size_type;
  
  TEUCHOS_TEST_FOR_EXCEPTION(exportLIDs.size() != numPacketsPerLID.size(),
			     std::invalid_argument, "packAndPrepareWithOwningPIDs: exportLIDs.size() = " << exportLIDs.size()
			     << "!= numPacketsPerLID.size() = " << numPacketsPerLID.size() << ".");
  
  // Get a reference to the matrix's row Map.
  const map_type& rowMap = * (SourceMatrix.getRowMap());

  // Sanity
  const bool locallyIndexed = SourceMatrix.isLocallyIndexed();
  TEUCHOS_TEST_FOR_EXCEPTION(!locallyIndexed,
			     std::invalid_argument, "packAndPrepareWithOwningPIDs: SourceMatrix must be locally indexed.");

  constantNumPackets = 0;
  
  // Get the GIDs of the rows we want to pack.
  Array<GO> exportGIDs (exportLIDs.size());
  const size_type numExportGIDs = exportGIDs.size ();
  for (size_type i = 0; i < numExportGIDs; ++i) {
    exportGIDs[i] = rowMap.getGlobalElement(exportLIDs[i]);
  }


  // Record initial sizing
  const size_t sizeOfPacket = sizeof(GO) + sizeof(Scalar) + sizeof(int);
  size_t totalNumEntries = 0;
  size_t maxRowLength = 0;
  for (size_type i = 0; i < exportGIDs.size(); ++i) {
    const size_t curNumEntries = SourceMatrix.getNumEntriesInGlobalRow(exportGIDs[i]);
    numPacketsPerLID[i] = curNumEntries * sizeOfPacket;
    totalNumEntries += curNumEntries;
    maxRowLength = std::max(curNumEntries, maxRowLength);
  }

  // Pack export data by interleaving rows' indices, pids and values in
  // the following way:
  //
  // [inds_row0 pids_row0 vals_row0 inds_row1 pids_row1 vals_row1 ... ]
  if (totalNumEntries > 0) {
    const size_t totalNumBytes = totalNumEntries * sizeOfPacket;
    exports.resize(totalNumBytes);

    // Current position in the 'exports' output array.
    size_t curOffsetInBytes = 0;

    // For each row of the matrix owned by the calling process, pack
    // that row's column indices and values into the exports array.

    // Locally indexed matrices always have a column Map.
    const map_type& colMap = * (SourceMatrix.getColMap());
    ArrayView<const LocalOrdinal> lidsView;
    ArrayView<const Scalar> valsView;
    
    // Temporary buffers for a copy of the column gids/pids
    Array<GO>  gids(as<size_type>(maxRowLength));
    Array<int> pids(as<size_type>(maxRowLength));
    
    const size_type numExportLIDs = exportLIDs.size();
    for (size_type i = 0; i < numExportLIDs; i++) {
      // Get a (locally indexed) view of the current row's data.
      SourceMatrix.getLocalRowView(exportLIDs[i], lidsView, valsView);
      
      // Convert column indices as LIDs to column indices as GIDs.
      const size_type curNumEntries = lidsView.size();
      size_t curNumEntriesST = as<size_t>(curNumEntries);
      ArrayView<GO>  gidsView = gids(0, curNumEntries);
      ArrayView<int> pidsView = pids(0, curNumEntries);
      for (size_type k = 0; k < curNumEntries; ++k) {
	gidsView[k] = colMap.getGlobalElement(lidsView[k]);
	pidsView[k] = SourcePids[lidsView[k]];
      }
      
      // Views of the right places in each array so everthing looks like the right data type	
      ArrayView<char> gidsViewOutChar = exports(curOffsetInBytes, curNumEntriesST*sizeof(GO));
      ArrayView<char> pidsViewOutChar = exports(curOffsetInBytes+curNumEntriesST*sizeof(GO), curNumEntriesST*sizeof(int));
      ArrayView<char> valsViewOutChar = exports(curOffsetInBytes+curNumEntriesST*(sizeof(GO)+sizeof(int)), curNumEntriesST*sizeof(Scalar));
      
      ArrayView<GO> gidsViewOut     = av_reinterpret_cast<GO>(gidsViewOutChar);
      ArrayView<int> pidsViewOut    = av_reinterpret_cast<int>(pidsViewOutChar);
      ArrayView<Scalar> valsViewOut = av_reinterpret_cast<Scalar>(valsViewOutChar);
      
      // Copy the row's data into the views of the exports array.
      std::copy(gidsView.begin(), gidsView.begin() + curNumEntriesST, gidsViewOut.begin());
      std::copy(pidsView.begin(), pidsView.begin() + curNumEntriesST, pidsViewOut.begin());
      std::copy(valsView.begin(), valsView.begin() + curNumEntriesST, valsViewOut.begin());
      
      // Keep track of how many bytes we packed.
      curOffsetInBytes += sizeOfPacket * curNumEntries;
    }

#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(curOffsetInBytes != totalNumBytes,
        std::logic_error, ": At end of method, the final offset bytes count "
        "curOffsetInBytes=" << curOffsetInBytes << " does not equal the total "
        "number of bytes packed totalNumBytes=" << totalNumBytes << ".  Please "
        "report this bug to the Tpetra developers.");
#endif //  HAVE_TPETRA_DEBUG
  }
  
}


#endif // TPETRA_IMPORT_UTIL_HPP
