// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TPETRA_MAPDATA_DECL_HPP
#define TPETRA_MAPDATA_DECL_HPP

#include <Teuchos_Object.hpp>
#include "Tpetra_Directory.hpp"

namespace Tpetra {

  /* Don't document this class
   * 
   * Like all data classes, Tpetra::MapData provides no functionality. It's sole purpose is to store 
   * the data associated with a Tpetra::Map object.
   */
  template <class LocalOrdinal, class GlobalOrdinal=LocalOrdinal>
  class MapData {
    friend class Map<LocalOrdinal,GlobalOrdinal>;

    private:
    /*! \brief Default constructor
     */
    MapData(Teuchos_Ordinal indexBase, 
            GlobalOrdinal numGlobalEntries,
            LocalOrdinal numMyEntries,
            GlobalOrdinal minAllGID,
            GlobalOrdinal maxAllGID,
            GlobalOrdinal minMyGID,
            GlobalOrdinal maxMyGID,
            const Teuchos::ArrayRCP<GlobalOrdinal> &lgMap,
            const std::map<GlobalOrdinal,LocalOrdinal> &glMap,
            bool contiguous,
            Teuchos::RCP<const Teuchos::Comm<int> > comm);

    MapData(Teuchos_Ordinal indexBase, 
            GlobalOrdinal numGlobalEntries,
            LocalOrdinal numMyEntries,
            GlobalOrdinal minAllGID,
            GlobalOrdinal maxAllGID,
            GlobalOrdinal minMyGID,
            GlobalOrdinal maxMyGID,
            const Teuchos::ArrayRCP<GlobalOrdinal> &lgMap,
            const std::map<GlobalOrdinal,LocalOrdinal>  &glMap,
            bool contiguous,
            Teuchos::RCP<const Teuchos::Comm<int> > comm,
            bool isLocal);

  public:
		//! Destructor.
		~MapData();

  private:
    // some of the following are globally coherent: that is, they have been guaranteed to 
    // match across all images, and may be assumed to do so
		Teuchos::RCP<const Teuchos::Comm<int> > comm_;
		const GlobalOrdinal numGlobalEntries_;
		const Teuchos_Ordinal indexBase_;
		const LocalOrdinal numMyEntries_;
    const GlobalOrdinal minMyGID_;
    const GlobalOrdinal maxMyGID_;
    const GlobalOrdinal minAllGID_;
    const GlobalOrdinal maxAllGID_;
    const bool contiguous_;
    const bool distributed_;
    Teuchos::ArrayRCP<GlobalOrdinal> lgMap_;
    std::map<GlobalOrdinal, LocalOrdinal> glMap_;
    Teuchos::RCP< Directory<LocalOrdinal,GlobalOrdinal> > directory_;

    bool checkIsDist();

		// declared but not defined, do not use
		MapData(const MapData<LocalOrdinal,GlobalOrdinal> & source);
		MapData<LocalOrdinal,GlobalOrdinal>& operator=(const MapData<LocalOrdinal,GlobalOrdinal> & source);
  };

} // namespace Tpetra

#endif // TPETRA_MAPDATA_DECL_HPP
