// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
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
  template <typename Ordinal>
  class MapData : public Teuchos::Object {
    friend class Map<Ordinal>;

  public:
    /*! \brief Default constructor
     */
    MapData(Ordinal indexBase, 
            Ordinal numGlobalEntries,
            Ordinal numMyEntries,
            Ordinal minAllGID,
            Ordinal maxAllGID,
            Ordinal minMyGID,
            Ordinal maxMyGID,
            const Teuchos::ArrayRCP<Ordinal> &lgMap,
            const std::map<Ordinal,Ordinal> &glMap,
            bool contiguous,
            Teuchos::RCP< Platform<Ordinal> > platform,
            Teuchos::RCP< Teuchos::Comm<Ordinal> > comm);

  MapData(Ordinal indexBase, 
          Ordinal numGlobalEntries,
          Ordinal numMyEntries,
          Ordinal minAllGID,
          Ordinal maxAllGID,
          Ordinal minMyGID,
          Ordinal maxMyGID,
          const Teuchos::ArrayRCP<Ordinal> &lgMap,
          const std::map<Ordinal,Ordinal>  &glMap,
          bool contiguous,
          Teuchos::RCP< Platform<Ordinal> > platform,
          Teuchos::RCP< Teuchos::Comm<Ordinal> > comm,
          bool isLocal);

		//! Destructor.
		~MapData();

  private:
    // some of the following are globally coherent: that is, they have been guaranteed to 
    // match across all images, and may be assumed to do so
    Teuchos::RCP< const Platform<Ordinal> > platform_;
		Teuchos::RCP< Teuchos::Comm<Ordinal> > comm_;
		const Ordinal numGlobalEntries_;
		const Ordinal indexBase_;
		const Ordinal numMyEntries_;
    const Ordinal minMyGID_;
    const Ordinal maxMyGID_;
    const Ordinal minAllGID_;
    const Ordinal maxAllGID_;
    const bool contiguous_;
    const bool distributed_;
    Teuchos::ArrayRCP<Ordinal> lgMap_;
    std::map<Ordinal, Ordinal> glMap_;
    Teuchos::RCP< Directory<Ordinal> > directory_;

		//! Copy constructor (declared but not defined, do not use)
		MapData(MapData<Ordinal> const& source);
		//! Assignment operator (declared but not defined, do not use)
		MapData<Ordinal>& operator = (MapData<Ordinal> const& source);

    bool checkIsDist();
    
  };

} // namespace Tpetra

#endif // TPETRA_MAPDATA_DECL_HPP
