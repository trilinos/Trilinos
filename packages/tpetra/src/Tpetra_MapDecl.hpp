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

#ifndef TPETRA_MAP_DECL_HPP
#define TPETRA_MAP_DECL_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_Object.hpp>
#include "Tpetra_Platform.hpp"

namespace Tpetra {

  template<typename Ordinal> class MapData;

  //! Tpetra::Map

  template<typename Ordinal>
  class Map : public Teuchos::Object {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    /*! \brief Map constructor with Tpetra-defined contiguous uniform distribution.
     *   The entries are distributed among nodes so that the subsets of global entries
     *   are non-overlapping and contiguous and as evenly distributed across the nodes as 
     *   possible.
     */
    Map(Ordinal numGlobalEntries, Ordinal indexBase, const Platform<Ordinal> &platform, bool local=false);

    /*! \brief Map constructor with a user-defined contiguous distribution.
     *  The entries are distributed among the nodes so that the subsets of global entries
     *  are non-overlapping and contiguous 
     *  
     *  If numGlobalEntries == -1, it will be computed via a global communication.
     *  Otherwise, it must be equal to the sum of the local entries across all 
     *  nodes. This will only be verified if Trilinos was compiled with --enable-teuchos-debug.
     *  If this verification fails, a std::invalid_argument exception will be thrown.
     */
    Map(Ordinal numGlobalEntries, Ordinal numMyEntries, Ordinal indexBase, 
        const Platform<Ordinal> &platform);

    //! Map constructor with user-defined non-contiguous (arbitrary) distribution.
    Map(Ordinal numGlobalEntries, const Teuchos::ArrayView<const Ordinal> &entryList, 
        Ordinal indexBase, const Platform<Ordinal> &platform);

    //! Map copy constructor.
    Map(const Map<Ordinal> &Map);

    //! Map destructor. 
    ~Map();

    //@}


    //! @name Map Attribute Methods
    //@{ 

    //! Returns the number of entries in this Map.
    Ordinal getNumGlobalEntries() const;

    //! Returns the number of entries belonging to the calling image.
    Ordinal getNumMyEntries() const;

    //! Returns the index base for this Map.
    Ordinal getIndexBase() const;

    //! Returns minimum local index
    Ordinal getMinLocalIndex() const;

    //! Returns maximum local index
    Ordinal getMaxLocalIndex() const;

    //! Returns minimum global index owned by this image
    Ordinal getMinGlobalIndex() const;

    //! Returns maximum global index owned by this image
    Ordinal getMaxGlobalIndex() const;

    //! Return the minimum global index over all images
    Ordinal getMinAllGlobalIndex() const;

    //! Return the maximum global index over all images
    Ordinal getMaxAllGlobalIndex() const;

    //! Return the local index for a given global index
    Ordinal getLocalIndex(Ordinal globalIndex) const;

    //! Return the global index for a given local index
    Ordinal getGlobalIndex(Ordinal localIndex) const;

    //! Returns the node IDs and corresponding local indices for a given list of global indices.
    /*! 
      \returns \c true signifies at least one specified global entry was not present in the directory.
     */
    bool getRemoteIndexList(const Teuchos::ArrayView<const Ordinal> & GIDList, 
                            const Teuchos::ArrayView<Ordinal> & imageIDList, 
                            const Teuchos::ArrayView<Ordinal> & LIDList) const;

    //! Returns the node IDs for a given list of global indices.
    /*! 
      \returns \c true signifies at least one specified global entry was not present in the directory.
     */
    bool getRemoteIndexList(const Teuchos::ArrayView<const Ordinal> & GIDList, 
                            const Teuchos::ArrayView<Ordinal> & imageIDList) const;

    //! Return a list of the global entries owned by this image
    Teuchos::ArrayView<const Ordinal> getMyGlobalEntries() const;

    //! Returns true if the local index value passed in is found on the calling image, returns false if it doesn't.
    bool isMyLocalIndex(Ordinal localIndex) const;

    //! Returns true if the global index value passed in is found the calling image, returns false if it doesn't.
    bool isMyGlobalIndex(Ordinal globalIndex) const;

    //! Returns true if this Map is distributed contiguously, returns false otherwise.
    bool isContiguous() const;

    //! Returns true if this Map is distributed across more than one image, returns false otherwise.
    bool isDistributed() const;

    //@}

    //! @name Boolean Tests
    //@{ 

    //! Returns true if \c map is compatible with this Map.
    bool isCompatible (const Map< Ordinal> &map) const;

    //! Returns true if \c map is identical to this Map.
    bool isSameAs (const Map<Ordinal> &map) const;

    //! Returns true if \c map is identical to this Map. Implemented in isSameAs().
    bool operator== (const Map< Ordinal > &map) const;

    //! Returns true if \c map is not identical to this Map. Implemented in isSameAs().
    bool operator!= (const Map< Ordinal > &map) const;

    //@}

    //@{ Misc. 

    //! Assignment operator
    Map<Ordinal>& operator = (const Map<Ordinal> & Source);

    //! Get the Platform object for this Map
    Teuchos::RCP< const Platform<Ordinal> > getPlatform() const;

    //! Get the Comm object for this Map
    Teuchos::RCP<const Teuchos::Comm<Ordinal> > getComm() const;

    //@}

    //@{ Implements Teuchos::Object 

    //! Prints the Map object to the output stream.
    /*! An << operator is inherited from Teuchos::Object, which uses the print method.*/
    void print(std::ostream& os) const;

    //@}


  private:

    Teuchos::RCP< MapData<Ordinal> > MapData_;

    // setup the directory
    void directorySetup();

  }; // Map class

} // Tpetra namespace

#endif // TPETRA_MAP_DECL_HPP

