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

namespace Tpetra {

	//! Tpetra::MapData:  The data class for Tpetra::Map
  /*! 
   * Like all data classes, Tpetra::MapData provides no functionality. It's sole purpose is to store 
   * the data associated with a Tpetra::Map object.
   */
  template <typename OrdinalType>
  class MapData : public Teuchos::Object {
    friend class Map<OrdinalType>;

  public:
    /*! \brief Default constructor
     */
    MapData(OrdinalType indexBase, 
            OrdinalType numGlobalEntries,
            OrdinalType numMyEntries,
            OrdinalType minAllGID,
            OrdinalType maxAllGID,
            OrdinalType minMyGID,
            OrdinalType maxMyGID,
            const map<OrdinalType, OrdinalType>& lgMap,
            const map<OrdinalType, OrdinalType>& glMap,
            bool contiguous,
            Teuchos::RCP< Platform<OrdinalType, OrdinalType> > platform,
            Teuchos::RCP< Comm<OrdinalType, OrdinalType> > comm);

		//! Destructor.
		~MapData();

  private:
		const OrdinalType indexBase_;
		const OrdinalType numMyEntries_;
		const OrdinalType numGlobalEntries_;
		Teuchos::RCP< const Comm<OrdinalType> > Comm_;
    const OrdinalType numMyEntries_;
    const OrdinalType minLID_;
    const OrdinalType maxLID_;
    const OrdinalType minMyGID_;
    const OrdinalType maxMyGID_;
    const OrdinalType minAllGID_;
    const OrdinalType maxAllGID_;
    const bool contiguous_;
    const bool global_;
    bool haveDirectory_;
    // TODO: why is lgMap_ const but glMap_ non-const?
    map<OrdinalType, OrdinalType> lgMap_;
    const map<OrdinalType, OrdinalType> glMap_;
    std::vector<OrdinalType> myGlobalEntries_;
    Teuchos::RCP< Directory<OrdinalType> > Directory_;

		//! Copy constructor (declared but not defined, do not use)
		MapData(MapData<OrdinalType> const& source);
		//! Assignment operator (declared but not defined, do not use)
		MapData<OrdinalType>& operator = (MapData<OrdinalType> const& source);

  };

} // namespace Tpetra

#endif // TPETRA_MAPDATA_DECL_HPP
