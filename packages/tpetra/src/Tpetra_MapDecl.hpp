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

// enums and defines
#include "Tpetra_ConfigDefs.hpp"

namespace Tpetra {

  //! A class for partitioning distributed objects.
  /*!
   This class is templated on \c LocalOrdinal and \c GlobalOrdinal. 
   The \c GlobalOrdinal type, if omitted, defaults to the \c LocalOrdinal type.
  */
  template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal>
  class Map : public Teuchos::Describable {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    /*! \brief Map constructor with Tpetra-defined contiguous uniform distribution.
     *   The entries are distributed among nodes so that the subsets of global entries
     *   are non-overlapping and contiguous and as evenly distributed across the nodes as 
     *   possible.
     */
    Map(global_size_t numGlobalEntries, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, LocalGlobal lg=GloballyDistributed);

    /*! \brief Map constructor with a user-defined contiguous distribution.
     *  The entries are distributed among the nodes so that the subsets of global entries
     *  are non-overlapping and contiguous 
     *  
     *  If numGlobalEntries == Teuchos::OrdinalTraits<global_size_t>::invalid(), it will be computed via a global communication.
     *  Otherwise, it must be equal to the sum of the local entries across all 
     *  nodes. This will only be verified if Trilinos was compiled with --enable-teuchos-debug.
     *  If this verification fails, a std::invalid_argument exception will be thrown.
     */
    Map(global_size_t numGlobalEntries, size_t numLocalEntries, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm);
        

    /*! \brief Map constructor with user-defined non-contiguous (arbitrary) distribution.
     *  
     *  If numGlobalEntries == Teuchos::OrdinalTraits<global_size_t>::invalid(), it will be computed via a global communication.
     *  Otherwise, it must be equal to the sum of the local entries across all 
     *  nodes. This will only be verified if Trilinos was compiled with --enable-teuchos-debug.
     *  If this verification fails, a std::invalid_argument exception will be thrown.
     */
    Map(global_size_t numGlobalEntries, const Teuchos::ArrayView<const GlobalOrdinal> &entryList, 
        GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm);

    //! Map destructor. 
    ~Map();

    //@}


    //! @name Map Attribute Methods
    //@{ 

    //! Returns the number of entries in this Map.
    inline global_size_t getNumGlobalEntries() const;

    //! Returns the number of entries belonging to the calling node.
    inline size_t getNumLocalEntries() const;

    //! Returns the index base for this Map.
    inline GlobalOrdinal getIndexBase() const;

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
      \returns IDNotPresent indicates that at least one global ID was not present in the directory. 
               Otherwise, returns AllIDsPresent.
     */
    LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal> & GIDList, 
                                    const Teuchos::ArrayView<                int> & nodeIDList, 
                                    const Teuchos::ArrayView<       LocalOrdinal> & LIDList) const;

    //! Returns the node IDs for a given list of global indices.
    /*! 
      \returns IDNotPresent indicates that at least one global ID was not present in the directory. 
               Otherwise, returns AllIDsPresent.
     */
    LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal> & GIDList, 
                                    const Teuchos::ArrayView<                int> & nodeIDList) const;

    //! Return a list of the global indices owned by this node.
    Teuchos::ArrayView<const GlobalOrdinal> getMyGlobalIndices() const;

    //! Returns true if the local index is valid for this Map on this node; returns false if it isn't.
    bool isMyLocalIndex(LocalOrdinal localIndex) const;

    //! Returns true if the global index is found in this Map on this node; returns false if it isn't.
    bool isMyGlobalIndex(GlobalOrdinal globalIndex) const;

    //! Returns true if this Map is distributed contiguously; returns false otherwise.
    bool isContiguous() const;

    //! Returns true if this Map is distributed across more than one node; returns false otherwise.
    bool isDistributed() const;

    //@}

    //! @name Boolean Tests
    //@{ 

    //! Returns true if \c map is compatible with this Map.
    bool isCompatible (const Map<LocalOrdinal,GlobalOrdinal> &map) const;

    //! Returns true if \c map is identical to this Map.
    bool isSameAs (const Map<LocalOrdinal,GlobalOrdinal> &map) const;

    //@}

    //@{ Misc. 

    //! Get the Comm object for this Map
    Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

    //@}

    //@{ Implements Teuchos::Describable 

    //! \brief Return a simple one-line description of this object.
    std::string description() const;

    //! Print the object with some verbosity level to a \c FancyOStream object.
    void describe( FancyOStream &out, const EVerbosityLevel verbLevel = verbLevel_default) const;

    //@}


  private:

    //! Setup the associated Directory.
    void setupDirectory();

    //! Perform communication to determine whether this is globally distributed or locally replicated.
    bool checkIsDist() const;

		//! Declared but not defined; do not use.
		Map(const Map<LocalOrdinal,GlobalOrdinal> & source);

		//! Declared but not defined; do not use.
		Map<LocalOrdinal,GlobalOrdinal>& operator=(const Map<LocalOrdinal,GlobalOrdinal> & source);

    // some of the following are globally coherent: that is, they have been guaranteed to 
    // match across all images, and may be assumed to do so
		Teuchos::RCP<const Teuchos::Comm<int> > comm_;

    // The based for global IDs in this Map.
		GlobalOrdinal indexBase_;
    //! The number of global IDs located in this Map across all nodes.
		GlobalOrdinal numGlobalEntries_;
    //! The number of global IDs located in this Map on this node.
		size_t numLocalEntries_;
    //! The minimum and maximum global IDs located in this Map on this node.
    GlobalOrdinal minMyGID_, maxMyGID_;
    //! The minimum and maximum global IDs located in this Map across all nodes.
    GlobalOrdinal minAllGID_, maxAllGID_;
    //! Indicates that the range of global indices are contiguous and ordered.
    bool contiguous_;
    //! Indicates that global indices of the map are non-identically distributed among different nodes.
    bool distributed_;
    //! A direct mapping from local IDs to global IDs.
    Teuchos::ArrayRCP<GlobalOrdinal> lgMap_;
    //! A mapping from global IDs to local IDs.
    std::map<GlobalOrdinal, LocalOrdinal> glMap_;
    //! A Directory for looking up nodes for this Map.
    Teuchos::RCP< Directory<LocalOrdinal,GlobalOrdinal> > directory_;

  }; // Map class

  //! Returns true if \c map is identical to this Map. Implemented in isSameAs().
  bool operator== (const Map<LocalOrdinal,GlobalOrdinal> &map1, const Map<LocalOrdinal,GlobalOrdinal> &map2);

  //! Returns true if \c map is not identical to this Map. Implemented in isSameAs().
  bool operator!= (const Map<LocalOrdinal,GlobalOrdinal> &map1, const Map<LocalOrdinal,GlobalOrdinal> &map2);


} // Tpetra namespace

#endif // TPETRA_MAP_DECL_HPP

