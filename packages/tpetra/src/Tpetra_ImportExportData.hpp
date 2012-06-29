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

#ifndef TPETRA_IMPORTEXPORTDATA_HPP
#define TPETRA_IMPORTEXPORTDATA_HPP

#include <Teuchos_Object.hpp>
#include "Tpetra_Distributor.hpp"

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of Import,Export, needed to prevent circular inclusions
  template<class LocalOrdinal, class GlobalOrdinal, class Node> class Import;
  template<class LocalOrdinal, class GlobalOrdinal, class Node> class Export;
#endif

  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  class ImportExportData : public Teuchos::Object {
    friend class Import<LocalOrdinal,GlobalOrdinal,Node>;
    friend class Export<LocalOrdinal,GlobalOrdinal,Node>;
  public:
    typedef LocalOrdinal local_ordinal_type;
    typedef GlobalOrdinal global_ordinal_type;
    typedef Node node_type;
    typedef Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

    /// \brief Constructor 
    ///
    /// \param source [in] Source Map of the Import or Export
    /// \param target [in] Target Map of the Import or Export
    ImportExportData (const Teuchos::RCP<const map_type>& source, 
                      const Teuchos::RCP<const map_type>& target);

    /// \brief Constructor with ParameterList for Distributor
    ///
    /// \param source [in] Source Map of the Import or Export
    /// \param target [in] Target Map of the Import or Export
    /// \param plist [in/out] List of parameters for the Distributor
    ImportExportData (const Teuchos::RCP<const map_type>& source, 
                      const Teuchos::RCP<const map_type>& target,
                      const Teuchos::RCP<Teuchos::ParameterList>& plist);

    ~ImportExportData();

  protected:
    // OT vectors
    Teuchos::Array<LocalOrdinal> permuteToLIDs_;
    Teuchos::Array<LocalOrdinal> permuteFromLIDs_;
    Teuchos::Array<LocalOrdinal> remoteLIDs_;
    Teuchos::Array<GlobalOrdinal> exportGIDs_;
    // These are ArrayRCP because in the construction of an Import object, they are allocated and returned by a call to 
    Teuchos::ArrayRCP<LocalOrdinal> exportLIDs_;
    Teuchos::ArrayRCP<int> exportImageIDs_;

    /// \brief Number of initial identical IDs.
    ///
    /// The number of IDs that are identical between the source and
    /// target Maps, up to the first different ID.
    size_t numSameIDs_;

    //! Source Map of the Import or Export
    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > source_;

    //! Target Map of the Import or Export
    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > target_;

    //! Communicator over which to distribute data.
    Teuchos::RCP<const Teuchos::Comm<int> > comm_;

    //! Distributor that does the work of distributing data.
    Distributor distributor_;

  private:
    //! Copy constructor (declared but not defined, do not use)
    ImportExportData (const ImportExportData<LocalOrdinal,GlobalOrdinal,Node> &rhs);
    //! Assignment operator (declared but not defined, do not use)
    ImportExportData<LocalOrdinal,GlobalOrdinal,Node>& 
    operator= (const ImportExportData<LocalOrdinal,GlobalOrdinal,Node> & rhs);
  }; // class ImportExportData

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ImportExportData<LocalOrdinal,GlobalOrdinal,Node>::
  ImportExportData (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& source,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& target)
  : numSameIDs_ (0)
  , source_ (source)
  , target_ (target)
  , comm_ (source->getComm())
  , distributor_ (comm_)
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ImportExportData<LocalOrdinal,GlobalOrdinal,Node>::
  ImportExportData (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& source,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& target,
                    const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : numSameIDs_ (0)
  , source_ (source)
  , target_ (target)
  , comm_ (source->getComm())
  , distributor_ (comm_, plist)
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ImportExportData<LocalOrdinal,GlobalOrdinal,Node>::~ImportExportData() 
  {}

} // namespace Tpetra

#endif // TPETRA_IMPORTEXPORTDATA_HPP
