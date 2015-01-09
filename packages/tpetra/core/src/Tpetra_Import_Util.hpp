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
#include "Tpetra_HashTable.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Util.hpp"
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




#endif // TPETRA_IMPORT_UTIL_HPP
