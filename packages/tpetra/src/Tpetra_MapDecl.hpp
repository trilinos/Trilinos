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

#ifndef TPETRA_MAP_DECL_HPP
#define TPETRA_MAP_DECL_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_Comm.hpp>

namespace Tpetra {

  template <class LocalOrdinal, class GlobalOrdinal> class MapData;

  //! A class for partitioning distributed objects.
  template<typename LocalOrdinal, class GlobalOrdinal=LocalOrdinal>
  class Map : public Teuchos::Describable {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    /*! \brief Map constructor with Tpetra-defined contiguous uniform distribution.
     *   The entries are distributed among nodes so that the subsets of global entries
     *   are non-overlapping and contiguous and as evenly distributed across the nodes as 
     *   possible.
     */
    Map(GlobalOrdinal numGlobalEntries, Teuchos_Ordinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, bool local=false);

    /*! \brief Map constructor with a user-defined contiguous distribution.
     *  The entries are distributed among the nodes so that the subsets of global entries
     *  are non-overlapping and contiguous 
     *  
     *  If numGlobalEntries == Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(), it will be computed via a global communication.
     *  Otherwise, it must be equal to the sum of the local entries across all 
     *  nodes. This will only be verified if Trilinos was compiled with --enable-teuchos-debug.
     *  If this verification fails, a std::invalid_argument exception will be thrown.
     */
    Map(GlobalOrdinal numGlobalEntries, LocalOrdinal numMyEntries, Teuchos_Ordinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm);
        

    /*! \brief Map constructor with user-defined non-contiguous (arbitrary) distribution.
     *  
     *  If numGlobalEntries == Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(), it will be computed via a global communication.
     *  Otherwise, it must be equal to the sum of the local entries across all 
     *  nodes. This will only be verified if Trilinos was compiled with --enable-teuchos-debug.
     *  If this verification fails, a std::invalid_argument exception will be thrown.
     */
    Map(GlobalOrdinal numGlobalEntries, const Teuchos::ArrayView<const GlobalOrdinal> &entryList, 
        Teuchos_Ordinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm);

    //! Map copy constructor.
    Map(const Map<LocalOrdinal,GlobalOrdinal> &Map);

    //! Map destructor. 
    ~Map();

    //@}


    //! @name Map Attribute Methods
    //@{ 

    //! Returns the number of entries in this Map.
    inline GlobalOrdinal getNumGlobalEntries() const;

    //! Returns the number of entries belonging to the calling node.
    inline LocalOrdinal getNumMyEntries() const;

    //! Returns the index base for this Map.
    inline Teuchos_Ordinal getIndexBase() const;

    //! Returns minimum local index
    inline LocalOrdinal getMinLocalIndex() const;

    //! Returns maximum local index
    inline LocalOrdinal getMaxLocalIndex() const;

    //! Returns minimum global index owned by this node
    inline GlobalOrdinal getMinGlobalIndex() const;

    //! Returns maximum global index owned by this node
    inline GlobalOrdinal getMaxGlobalIndex() const;

    //! Return the minimum global index over all nodes
    inline GlobalOrdinal getMinAllGlobalIndex() const;

    //! Return the maximum global index over all nodes
    inline GlobalOrdinal getMaxAllGlobalIndex() const;

    //! \brief Return the local index for a given global index
    /*! If the global index is not owned by this node, returns Teuchos::OrdinalTraits<LocalOrdinal>::invalid(). */
    LocalOrdinal getLocalIndex(GlobalOrdinal globalIndex) const;

    //! Return the global index for a given local index
    /*! If the local index is not valid for this node, returns Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(). */
    GlobalOrdinal getGlobalIndex(LocalOrdinal localIndex) const;

    //! Returns the node IDs and corresponding local indices for a given list of global indices.
    /*! 
      \returns \c true signifies at least one specified global entry was not present in the directory.
     */
    bool getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal> & GIDList, 
                            const Teuchos::ArrayView<int> & nodeIDList, 
                            const Teuchos::ArrayView<LocalOrdinal> & LIDList) const;

    //! Returns the node IDs for a given list of global indices.
    /*! 
      \returns \c true signifies at least one specified global entry was not present in the directory.
     */
    bool getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal> & GIDList, 
                            const Teuchos::ArrayView<int> & nodeIDList) const;

    //! Return a list of the global entries owned by this node
    Teuchos::ArrayView<const GlobalOrdinal> getMyGlobalEntries() const;

    //! Returns true if the local index value passed in is found on the calling node, returns false if it doesn't.
    bool isMyLocalIndex(LocalOrdinal localIndex) const;

    //! Returns true if the global index value passed in is found the calling node, returns false if it doesn't.
    bool isMyGlobalIndex(GlobalOrdinal globalIndex) const;

    //! Returns true if this Map is distributed contiguously, returns false otherwise.
    bool isContiguous() const;

    //! Returns true if this Map is distributed across more than one node, returns false otherwise.
    bool isDistributed() const;

    //@}

    //! @name Boolean Tests
    //@{ 

    //! Returns true if \c map is compatible with this Map.
    bool isCompatible (const Map<LocalOrdinal,GlobalOrdinal> &map) const;

    //! Returns true if \c map is identical to this Map.
    bool isSameAs (const Map<LocalOrdinal,GlobalOrdinal> &map) const;

    //! Returns true if \c map is identical to this Map. Implemented in isSameAs().
    bool operator== (const Map<LocalOrdinal,GlobalOrdinal> &map) const;

    //! Returns true if \c map is not identical to this Map. Implemented in isSameAs().
    bool operator!= (const Map<LocalOrdinal,GlobalOrdinal> &map) const;

    //@}

    //@{ Misc. 

    //! Assignment operator
    Map<LocalOrdinal,GlobalOrdinal>& operator = (const Map<LocalOrdinal,GlobalOrdinal> &source);

    //! Get the Comm object for this Map
    Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

    //@}

    //@{ Implements Teuchos::Object 

    //! Prints the Map object to the output stream.
    /*! An << operator is inherited from Teuchos::Object, which uses the print method.*/
    void print(std::ostream& os) const;

    //@}


  private:

    Teuchos::RCP<MapData<LocalOrdinal,GlobalOrdinal> > MapData_;

    // setup the directory
    void directorySetup();

  }; // Map class

} // Tpetra namespace

#endif // TPETRA_MAP_DECL_HPP

