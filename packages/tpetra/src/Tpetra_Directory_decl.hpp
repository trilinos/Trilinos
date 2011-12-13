// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

#ifndef TPETRA_DIRECTORY_DECL_HPP
#define TPETRA_DIRECTORY_DECL_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Map_decl.hpp"

namespace Tpetra {

  //! A class used for efficienctly accessing global node/index information from a Map.
  
  /*! For Map objects, a Directory object must be created to allow referencing
      of non-local elements. Tpetra::Directory produces and contains a uniform linear
      Map, allowing the directory to be distributed across all nodes.
      
      This class has a single constructor, accepting the Map object for which the directory 
      is created.

      This class is templated on \c LocalOrdinal and \c GlobalOrdinal. 
      The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
      type, if omitted, defaults to the \c LocalOrdinal type.
  */

  template<class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class Directory : public Teuchos::Describable {
  public:

    //! @name Constructors/Destructor.
    //@{ 

    //! Constructor
    explicit Directory(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map);

    //! Destructor.
    ~Directory();

    //@}

    //! @name Query methods.
    //@{ 

    //! \brief Returns node info for non-local Map entries.
    /*! Given a list of global IDs, this function returns the corresponding list of
      owning node IDs.

      \param globalIDs [in] List of global IDs to look up.

      \param nodeIDs [out] On return, contains node IDs for the global IDs in question. 
      -1 corresponds to global entries not present in the directory.

      \returns IDNotPresent indicates that at least one global ID was not present in the directory. 
               Otherwise, returns AllIDsPresent.

      \note If <tt>nodeIDs.size() != globalIDs.size()</tt>, then a \c std::runtime_error exception is thrown.
    */
    LookupStatus getDirectoryEntries(const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs, 
                                     const Teuchos::ArrayView<int> &nodeIDs) const;
    
    //! \brief Returns node info for non-local Map entries.
    /*! Given a list of global IDs, this function returns the corresponding list of
      owning node IDs and local IDs.

      \param globalIDs [in] List of global IDs to look up.

      \param nodeIDs [out] On return, contains node IDs for the global IDs in question. 
      -1 corresponds to global entries not present in the directory.

      \param localIDs [out] On return contains the local ID of the global on the owning node. 
      Teuchos::OrdinalTraits<LocalOrdinal>::invalid() corresponds to global entries not present in the directory.

      \returns IDNotPresent indicates that at least one global ID was not present in the directory. 
               Otherwise, returns AllIDsPresent.

      \note If <tt>nodeIDs.size() != globalIDs.size()</tt> or 
               <tt>localIDs.size() != globalIDs.size()</tt>, then a \c std::runtime_error exception is thrown.
    */
    LookupStatus getDirectoryEntries(const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs, 
                                     const Teuchos::ArrayView<int> &nodeIDs, 
                                     const Teuchos::ArrayView<LocalOrdinal> &localIDs) const;
    //@}
    
  private:
    Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > map_, directoryMap_;
    Teuchos::RCP<const Teuchos::Comm<int> > comm_;
    std::vector<GlobalOrdinal> allMinGIDs_; // size comm_->getSize()+1; entry i contains minGID for ith node; last entry contains maxGID in the directory
    std::vector<int> nodeIDs_;
    std::vector<LocalOrdinal> LIDs_;

    Directory(const Directory<LocalOrdinal,GlobalOrdinal,Node> &directory);

    //! declared but not defined, to prevent default implementation, do not use
    Directory<LocalOrdinal,GlobalOrdinal,Node> & operator = (const Directory<LocalOrdinal,GlobalOrdinal,Node> &source);

    // common code for both versions of getDirectoryEntries
    LookupStatus getEntries(const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs, 
                            const Teuchos::ArrayView<int> &nodeIDs, 
                            const Teuchos::ArrayView<LocalOrdinal> &localIDs, 
                            bool computeLIDs) const;

    // directory setup for non-contiguous ES
    void generateDirectory();

  }; // class Directory

} // namespace Tpetra

#endif // TPETRA_DIRECTORY_DECL_HPP

