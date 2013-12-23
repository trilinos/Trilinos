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


    //! unpackAndCombineWithOwningPIDsCount.
    /*! \brief Perform the count for unpacking the imported column indices pids, and values, and combining them into matrix.
      Returns (a ceiling on) the number of local non-zeros in the matrix.  If there are no shared rows in the SourceMatrix this count is exact.

      Note: This routine also counts the copyAndPermute nonzeros in addition to those that come in via import.

      \warning This method is intended for expert developer use only, and should never be called by user code.
    */
    template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, typename LocalMatOps>
    size_t unpackAndCombineWithOwningPIDsCount(const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & SourceMatrix,
					       const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
					       const Teuchos::ArrayView<const char> &imports,
					       const Teuchos::ArrayView<size_t> &numPacketsPerLID,
					       size_t constantNumPackets,
					       Distributor &distor,
					       CombineMode combineMode,
					       size_t numSameIDs,
					       const ArrayView<const LocalOrdinal> &permuteToLIDs,
					       const ArrayView<const LocalOrdinal> &permuteFromLIDs);


    // ===================================================================
    //! UnpackAndCombineIntoCrsArrays
    /*! You should call UnpackWithOwningPIDsCount first and allocate all arrays accordingly.

      Note: The SourcePids vector (on input) should contain of owning PIDs for each column in the (source) ColMap, as from Epetra_Util::GetPids,
      with the "-1 for local" option being used.

      Note: The TargetPids vector (on output) will contain of owning PIDs for each column in the (target) ColMap, as from Epetra_Util::GetPids,
      with the "-1 for local" option being used.
    */
    template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, typename LocalMatOps>
    void unpackAndCombineIntoCrsArrays(const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & SourceMatrix,
				       const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
				       const Teuchos::ArrayView<const char> &imports,
				       const Teuchos::ArrayView<size_t> &numPacketsPerLID,
				       size_t constantNumPackets,
				       Distributor &distor,
				       CombineMode combineMode,
				       size_t numSameIDs,
				       const ArrayView<const LocalOrdinal> &permuteToLIDs,
				       const ArrayView<const LocalOrdinal> &permuteFromLIDs,
				       size_t TargetNumRows,
				       size_t TargetNumNonzeros,
				       const ArrayView<size_t> &rowPointers,
				       const ArrayView<GlobalOrdinal> &columnIndices,
				       const ArrayView<Scalar> &values,
				       const Teuchos::ArrayView<const int> &SourcePids,
				       Teuchos::Array<int> &TargetPids);


    // ===================================================================
    //! sortCrsEntries
    /*! sorts the entries of the matrix by colind w/i each row
    */
    template<typename Scalar, typename Ordinal>
    void sortCrsEntries(const Teuchos::ArrayView<size_t> &CRS_rowptr, const Teuchos::ArrayView<Ordinal> & CRS_colind, const Teuchos::ArrayView<Scalar> &CRS_vals);

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
  using Teuchos::RCP;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Map<LocalOrdinal,GlobalOrdinal,Node>  map_type;
  typedef typename ArrayView<const LO>::size_type size_type;

  TEUCHOS_TEST_FOR_EXCEPTION(!SourceMatrix.isLocallyIndexed(),std::invalid_argument, "packAndPrepareWithOwningPIDs: SourceMatrix must be locally indexed.");
  TEUCHOS_TEST_FOR_EXCEPTION(exportLIDs.size() != numPacketsPerLID.size(),
			     std::invalid_argument, "packAndPrepareWithOwningPIDs: exportLIDs.size() = " << exportLIDs.size()
			     << "!= numPacketsPerLID.size() = " << numPacketsPerLID.size() << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(as<size_t>(SourcePids.size()) != SourceMatrix.getColMap()->getNodeNumElements(),
  			     std::invalid_argument, "packAndPrepareWithOwningPIDs: SourcePids.size() = " << SourcePids.size()
 			     << "!= SourceMatrix.getColMap()->getNodeNumElements() = " << SourceMatrix.getColMap()->getNodeNumElements() << ".");

  // Get a reference to the matrix's row Map.
  const map_type& rowMap = * (SourceMatrix.getRowMap());

  // Sanity


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
      TEUCHOS_TEST_FOR_EXCEPTION(curOffsetInBytes != totalNumBytes,
        std::logic_error, "packAndPrepareWithOwningPids: At end of method, the final offset bytes count "
        "curOffsetInBytes=" << curOffsetInBytes << " does not equal the total "
        "number of bytes packed totalNumBytes=" << totalNumBytes << ".  Please "
        "report this bug to the Tpetra developers.");
#endif //  HAVE_TPETRA_DEBUG
  }

}



//----------------------------------------------------------------------------
template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, typename LocalMatOps>
size_t Tpetra::Import_Util::unpackAndCombineWithOwningPIDsCount(const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & SourceMatrix,
								const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
								const Teuchos::ArrayView<const char> &imports,
								const Teuchos::ArrayView<size_t> &numPacketsPerLID,
								size_t constantNumPackets,
								Distributor &distor,
								CombineMode combineMode,
								size_t numSameIDs,
								const ArrayView<const LocalOrdinal> &permuteToLIDs,
								const ArrayView<const LocalOrdinal> &permuteFromLIDs) {
  typedef LocalOrdinal LO;
  typedef Map<LocalOrdinal,GlobalOrdinal,Node>  map_type;
  typedef typename ArrayView<const LO>::size_type size_type;
  size_t nnz = 0;

  // CopyAndPermuteSection
  TEUCHOS_TEST_FOR_EXCEPTION(permuteToLIDs.size() != permuteFromLIDs.size(),
			     std::invalid_argument, "unpackAndCombineWithOwningPIDsCount: permuteToLIDs.size() = " << permuteToLIDs.size()
			     << "!= permuteFromLIDs.size() = " << permuteFromLIDs.size() << ".");
  const bool locallyIndexed = SourceMatrix.isLocallyIndexed();
  TEUCHOS_TEST_FOR_EXCEPTION(!locallyIndexed,std::invalid_argument, "unpackAndCombineWithOwningPIDsCount: SourceMatrix must be locally indexed.");

  // Copy
  const LO numSameIDs_as_LID = Teuchos::as<LO>(numSameIDs);
  for (LO sourceLID = 0; sourceLID < numSameIDs_as_LID; sourceLID++)
    nnz+=SourceMatrix.getNumEntriesInLocalRow(sourceLID);

  // Permute
  const size_t numPermuteToLIDs = Teuchos::as<size_t>(permuteToLIDs.size());
  for (size_t p = 0; p < numPermuteToLIDs; p++)
    nnz+=SourceMatrix.getNumEntriesInLocalRow(permuteFromLIDs[p]);

  // UnpackAndCombine Section
  TEUCHOS_TEST_FOR_EXCEPTION(importLIDs.size() != numPacketsPerLID.size(),
			     std::invalid_argument, "unpackAndCombineWithOwningPIDsCount: importLIDs.size() = " << importLIDs.size()
			     << "!= numPacketsPerLID.size() = " << numPacketsPerLID.size() << ".");

  const size_t sizeOfPacket    = sizeof(GlobalOrdinal)  + sizeof(int) + sizeof(Scalar);

  size_t curOffsetInBytes = 0;
  for (size_type i = 0; i < importLIDs.size(); ++i) {
    const size_t rowSize = numPacketsPerLID[i] / sizeOfPacket;
    curOffsetInBytes += rowSize * sizeOfPacket;
    nnz +=rowSize;
  }
  return nnz;
}



//----------------------------------------------------------------------------
template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, typename LocalMatOps>
void Tpetra::Import_Util::unpackAndCombineIntoCrsArrays(const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & SourceMatrix,
                                                        const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
                                                        const Teuchos::ArrayView<const char> &imports,
                                                        const Teuchos::ArrayView<size_t> &numPacketsPerLID,
                                                        size_t constantNumPackets,
                                                        Distributor &distor,
                                                        CombineMode combineMode,
                                                        size_t numSameIDs,
                                                        const ArrayView<const LocalOrdinal> &permuteToLIDs,
                                                        const ArrayView<const LocalOrdinal> &permuteFromLIDs,
                                                        size_t TargetNumRows,
                                                        size_t TargetNumNonzeros,
                                                        const ArrayView<size_t> &CSR_rowptr,
                                                        const ArrayView<GlobalOrdinal> &CSR_colind,
                                                        const ArrayView<Scalar> &CSR_vals,
                                                        const Teuchos::ArrayView<const int> &SourcePids,
                                                        Teuchos::Array<int> &TargetPids) {
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::av_reinterpret_cast;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Map<LocalOrdinal,GlobalOrdinal,Node>  map_type;
  typedef typename ArrayView<const LO>::size_type size_type;

  size_t i,j;
  size_t N=TargetNumRows;
  size_t mynnz = TargetNumNonzeros;
  int MyPID = SourceMatrix.getComm()->getRank();

  // Zero the rowptr
  TEUCHOS_TEST_FOR_EXCEPTION(N+1 != as<size_t>(CSR_rowptr.size()),
                             std::invalid_argument, "unpackAndCombineIntoCrsArrays: CsR_rowptr.size() = " << CSR_rowptr.size()
                             << "!= TargetNumRows+1 = " << TargetNumRows+1 << ".");
  for(i=0; i<N+1; i++) CSR_rowptr[i]=0;

  // SameIDs: Always first, always in the same place
  for(i=0; i<numSameIDs; i++)
    CSR_rowptr[i]=SourceMatrix.getNumEntriesInLocalRow(as<LO>(i));

  // PermuteIDs: Still local, but reordered
  TEUCHOS_TEST_FOR_EXCEPTION(permuteToLIDs.size() != permuteFromLIDs.size(),
                             std::invalid_argument, "unpackAndCombineIntoCrsArrays: permuteToLIDs.size() = " << permuteToLIDs.size()
                             << "!= permuteFromLIDs.size() = " << permuteFromLIDs.size() << ".");
  size_t numPermuteIDs = permuteToLIDs.size();
  for(i=0; i<numPermuteIDs; i++)
    CSR_rowptr[permuteToLIDs[i]] = SourceMatrix.getNumEntriesInLocalRow(permuteFromLIDs[i]);



  // Setup CSR_rowptr for remotes
  TEUCHOS_TEST_FOR_EXCEPTION(importLIDs.size() != numPacketsPerLID.size(),
                             std::invalid_argument, "unpackAndCombineIntoCrsArrays: importLIDs.size() = " << importLIDs.size()
                             << "!= numPacketsPerLID.size() = " << numPacketsPerLID.size() << ".");

  const size_t sizeOfPacket     = sizeof(GlobalOrdinal)  + sizeof(int) + sizeof(Scalar);
  const size_t totalNumBytes    = imports.size();
  const size_t RemoteNumEntries = totalNumBytes / sizeOfPacket;
  for (size_type k = 0; k < importLIDs.size(); ++k) {
    const size_t rowSize = numPacketsPerLID[k] / sizeOfPacket;
    CSR_rowptr[importLIDs[k]] += rowSize;
  }

  // If multiple procs contribute to a row;
  Teuchos::Array<size_t> NewStartRow(N+1);

  // Turn row length into a real CSR_rowptr
  size_t last_len = CSR_rowptr[0];
  CSR_rowptr[0] = 0;
  for(i=1; i<N+1; i++){
    size_t new_len    = CSR_rowptr[i];
    CSR_rowptr[i]  = last_len + CSR_rowptr[i-1];
    NewStartRow[i] = CSR_rowptr[i];
    last_len       = new_len;
  }

  // Preseed TargetPids with -1 for local
  if(as<size_t>(TargetPids.size())!=mynnz) TargetPids.resize(mynnz);
  TargetPids.assign(mynnz,-1);

  // Grab pointers for SourceMatrix
  ArrayRCP<const size_t> Source_rowptr;
  ArrayRCP<const LO>     Source_colind;
  ArrayRCP<const Scalar> Source_vals;
  SourceMatrix.getAllValues(Source_rowptr,Source_colind,Source_vals);

  const map_type& sourceColMap = * (SourceMatrix.getColMap());

  // SameIDs: Copy the data over
  for(i=0; i<numSameIDs; i++) {
    size_t FromRow = Source_rowptr[i];
    size_t ToRow   = CSR_rowptr[i];
    NewStartRow[i] += Source_rowptr[i+1]-Source_rowptr[i];

    for(j=Source_rowptr[i]; j<Source_rowptr[i+1]; j++) {
      CSR_vals[ToRow + j - FromRow]   = Source_vals[j];
      CSR_colind[ToRow + j - FromRow] = sourceColMap.getGlobalElement(Source_colind[j]);
      TargetPids[ToRow + j - FromRow] = (SourcePids[Source_colind[j]] != MyPID) ? SourcePids[Source_colind[j]] : -1;
    }
  }

  // PermuteIDs: Copy the data over
  for(i=0; i<numPermuteIDs; i++) {
    LO FromLID     = permuteFromLIDs[i];
    size_t FromRow = Source_rowptr[FromLID];
    size_t ToRow   = CSR_rowptr[permuteToLIDs[i]];

    NewStartRow[permuteToLIDs[i]] += Source_rowptr[FromLID+1]-Source_rowptr[FromLID];

    for(j=Source_rowptr[FromLID]; j<Source_rowptr[FromLID+1]; j++) {
      CSR_vals[ToRow + j - FromRow]   = Source_vals[j];
      CSR_colind[ToRow + j - FromRow] = sourceColMap.getGlobalElement(Source_colind[j]);
      TargetPids[ToRow + j - FromRow] = (SourcePids[Source_colind[j]] != MyPID) ? SourcePids[Source_colind[j]] : -1;
    }
  }

  // RemoteIDs: Loop structure following UnpackAndCombine
  if(RemoteNumEntries > 0) {
    // data packed as follows:
    // [inds_row0 pids_row0 vals_row0 inds_row1 pids_row1 vals_row1 ...]
    ArrayView<const char>   avIndsC, avPidsC, avValsC;
    ArrayView<const GO>     avInds;
    ArrayView<const int>    avPids;
    ArrayView<const Scalar> avVals;

    size_t curOffsetInBytes = 0;
    for (i = 0; i < importLIDs.size(); ++i) {
      const size_t rowSize = numPacketsPerLID[i] / sizeOfPacket;
      LO ToLID     = importLIDs[i];
      int StartRow = NewStartRow[ToLID];
      NewStartRow[ToLID]+=rowSize;
      if (rowSize == 0) continue;

      // Get views of the import (incoming data) buffers.
      avIndsC = imports(curOffsetInBytes, rowSize*sizeof(GO));
      avPidsC = imports(curOffsetInBytes+rowSize*sizeof(GO), rowSize*sizeof(int));
      avValsC = imports(curOffsetInBytes+rowSize*(sizeof(GO)+sizeof(int)), rowSize*sizeof(Scalar));

      avInds = av_reinterpret_cast<const GO>(avIndsC);
      avPids = av_reinterpret_cast<const int>(avPidsC);
      avVals = av_reinterpret_cast<const Scalar>(avValsC);

      for(j=0; j<rowSize; j++){
        CSR_vals[StartRow + j]   = avVals[j];
        CSR_colind[StartRow + j] = avInds[j];
        TargetPids[StartRow + j] = (avPids[j] != MyPID) ? avPids[j] : -1;
      }
      curOffsetInBytes += rowSize * sizeOfPacket;
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(curOffsetInBytes != totalNumBytes,
                               std::logic_error, "unpackAndCombineIntoCrsArrays: After unpacking and counting all the imports, the "
                               "final offset in bytes curOffsetInBytes=" << curOffsetInBytes << " != "
                               "total number of bytes totalNumBytes=" << totalNumBytes << ".  Please "
                               "report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
  }
}


//----------------------------------------------------------------------------
// Because I can't make sense out of the Tpetra sort routines.
template<typename Scalar, typename Ordinal>
void Tpetra::Import_Util::sortCrsEntries(const Teuchos::ArrayView<size_t> &CRS_rowptr, const Teuchos::ArrayView<Ordinal> & CRS_colind, const Teuchos::ArrayView<Scalar> &CRS_vals){
  // For each row, sort column entries from smallest to largest.
  // Use shell sort. Stable sort so it is fast if indices are already sorted.
  // Code copied from  Epetra_CrsMatrix::SortEntries()
  int NumRows = CRS_rowptr.size()-1;
  for(int i = 0; i < NumRows; i++){
    size_t start=CRS_rowptr[i];

    Scalar* locValues   = &CRS_vals[start];
    size_t NumEntries   = CRS_rowptr[i+1] - start;
    Ordinal* locIndices = &CRS_colind[start];

    Ordinal n = NumEntries;
    Ordinal m = n/2;

    while(m > 0) {
      Ordinal max = n - m;
      for(Ordinal j = 0; j < max; j++) {
        for(Ordinal k = j; k >= 0; k-=m) {
          if(locIndices[k+m] >= locIndices[k])
            break;
          Scalar dtemp = locValues[k+m];
          locValues[k+m] = locValues[k];
          locValues[k] = dtemp;
          Ordinal itemp = locIndices[k+m];
          locIndices[k+m] = locIndices[k];
          locIndices[k] = itemp;
        }
      }
      m = m/2;
    }
  }
}

#endif // TPETRA_IMPORT_UTIL_HPP
