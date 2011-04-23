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

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>

// enums and defines
#include "Tpetra_ConfigDefs.hpp"

/** \file Tpetra_Map_decl.hpp 

    The declarations for the class Tpetra::Map and related non-member constructors.
 */

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward dec
  template <class LO, class GO, class N> class Directory;
#endif

  /** \brief A class for partitioning distributed objects.

   This class is templated on \c LocalOrdinal and \c GlobalOrdinal. 
   The \c GlobalOrdinal type, if omitted, defaults to the \c LocalOrdinal type.
  */
  template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class Map : public Teuchos::Describable {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    /** \brief Map constructor with Tpetra-defined contiguous uniform distribution.
     *   The elements are distributed among nodes so that the subsets of global elements
     *   are non-overlapping and contiguous and as evenly distributed across the nodes as 
     *   possible.
     */
    Map(global_size_t numGlobalElements, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, 
        LocalGlobal lg=GloballyDistributed, const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode());

    /** \brief Map constructor with a user-defined contiguous distribution.
     *  The elements are distributed among the nodes so that the subsets of global elements
     *  are non-overlapping and contiguous 
     *  
     *  If numGlobalElements == Teuchos::OrdinalTraits<global_size_t>::invalid(), it will be computed via a global communication.
     *  Otherwise, it must be equal to the sum of the local elements across all 
     *  nodes. This will only be verified if Trilinos was compiled with --enable-teuchos-debug.
     *  If this verification fails, a std::invalid_argument exception will be thrown.
     */
    Map(global_size_t numGlobalElements, size_t numLocalElements, GlobalOrdinal indexBase, 
        const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode());
        

    /** \brief Map constructor with user-defined non-contiguous (arbitrary) distribution.
     *  
     *  If numGlobalElements == Teuchos::OrdinalTraits<global_size_t>::invalid(), it will be computed via a global communication.
     *  Otherwise, it must be equal to the sum of the local elements across all 
     *  nodes. This will only be verified if Trilinos was compiled with --enable-teuchos-debug.
     *  If this verification fails, a std::invalid_argument exception will be thrown.
     */
    Map(global_size_t numGlobalElements, const Teuchos::ArrayView<const GlobalOrdinal> &elementList, GlobalOrdinal indexBase, 
        const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode());

    //! Map destructor. 
    ~Map();

    //@}


    //! @name Map Attribute Methods
    //@{ 

    //! Returns the number of elements in this Map.
    inline global_size_t getGlobalNumElements() const { return numGlobalElements_; }

    //! Returns the number of elements belonging to the calling node.
    inline size_t getNodeNumElements() const { return numLocalElements_; }

    //! Returns the index base for this Map.
    inline GlobalOrdinal getIndexBase() const { return indexBase_; }

    //! Returns minimum local index
    inline LocalOrdinal getMinLocalIndex() const { return Teuchos::OrdinalTraits<LocalOrdinal>::zero(); }

    //! Returns maximum local index
    inline LocalOrdinal getMaxLocalIndex() const { return Teuchos::as<LocalOrdinal>(numLocalElements_-1); }

    //! Returns minimum global index owned by this node
    inline GlobalOrdinal getMinGlobalIndex() const { return minMyGID_; }

    //! Returns maximum global index owned by this node
    inline GlobalOrdinal getMaxGlobalIndex() const { return maxMyGID_; }

    //! Return the minimum global index over all nodes
    inline GlobalOrdinal getMinAllGlobalIndex() const { return minAllGID_; }

    //! Return the maximum global index over all nodes
    inline GlobalOrdinal getMaxAllGlobalIndex() const { return maxAllGID_; }

    //! \brief Return the local index for a given global index
    /** If the global index is not owned by this node, returns Teuchos::OrdinalTraits<LocalOrdinal>::invalid(). */
    LocalOrdinal getLocalElement(GlobalOrdinal globalIndex) const;

    //! Return the global index for a given local index
    /** If the local index is not valid for this node, returns Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(). */
    GlobalOrdinal getGlobalElement(LocalOrdinal localIndex) const;

    //! Returns the node IDs and corresponding local indices for a given list of global indices.
    /** 
      \returns IDNotPresent indicates that at least one global ID was not present in the directory. 
               Otherwise, returns AllIDsPresent.
     */
    LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal> & GIDList, 
                                    const Teuchos::ArrayView<                int> & nodeIDList, 
                                    const Teuchos::ArrayView<       LocalOrdinal> & LIDList) const;

    //! Returns the node IDs for a given list of global indices.
    /** 
      \returns IDNotPresent indicates that at least one global ID was not present in the directory. 
               Otherwise, returns AllIDsPresent.
     */
    LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal> & GIDList, 
                                    const Teuchos::ArrayView<                int> & nodeIDList) const;

    //! Return a list of the global indices owned by this node.
    Teuchos::ArrayView<const GlobalOrdinal> getNodeElementList() const;

    //! Returns true if the local index is valid for this Map on this node; returns false if it isn't.
    bool isNodeLocalElement(LocalOrdinal localIndex) const;

    //! Returns true if the global index is found in this Map on this node; returns false if it isn't.
    bool isNodeGlobalElement(GlobalOrdinal globalIndex) const;

    //! Returns true if this Map is distributed contiguously; returns false otherwise.
    bool isContiguous() const;

    //! Returns true if this Map is distributed across more than one node; returns false otherwise.
    bool isDistributed() const;

    //@}

    //! @name Boolean Tests
    //@{ 

    //! Returns true if \c map is compatible with this Map.
    bool isCompatible (const Map<LocalOrdinal,GlobalOrdinal,Node> &map) const;

    //! Returns true if \c map is identical to this Map.
    bool isSameAs (const Map<LocalOrdinal,GlobalOrdinal,Node> &map) const;

    //@}

    //@{ Misc. 

    //! Get the Comm object for this Map
    const Teuchos::RCP<const Teuchos::Comm<int> > & getComm() const;

    //! Get the Node object for this Map
    const Teuchos::RCP<Node> & getNode() const;

    //@}

    //@{ Implements Teuchos::Describable 

    //! \brief Return a simple one-line description of this object.
    std::string description() const;

    //! Print the object with some verbosity level to a \c FancyOStream object.
    void describe( Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = verbLevel_default) const;

    //@}


  private:

    //! Setup the associated Directory.
    void setupDirectory();

    //! Perform communication to determine whether this is globally distributed or locally replicated.
    bool checkIsDist() const;

		//! Declared but not defined; do not use.
		Map(const Map<LocalOrdinal,GlobalOrdinal,Node> & source);

		//! Declared but not defined; do not use.
		Map<LocalOrdinal,GlobalOrdinal,Node>& operator=(const Map<LocalOrdinal,GlobalOrdinal,Node> & source);

    // some of the following are globally coherent: that is, they have been guaranteed to 
    // match across all images, and may be assumed to do so
		Teuchos::RCP<const Teuchos::Comm<int> > comm_;

    // Map doesn't need node yet, but it likely will later. In the meantime, passing a Node to Map means that we don't have to 
    // pass a Node to downstream classes such as MultiVector, Vector, CrsGraph and CrsMatrix
    Teuchos::RCP<Node> node_;

    // The based for global IDs in this Map.
		GlobalOrdinal indexBase_;
    //! The number of global IDs located in this Map across all nodes.
		global_size_t numGlobalElements_;
    //! The number of global IDs located in this Map on this node.
		size_t numLocalElements_;
    //! The minimum and maximum global IDs located in this Map on this node.
    GlobalOrdinal minMyGID_, maxMyGID_;
    //! The minimum and maximum global IDs located in this Map across all nodes.
    GlobalOrdinal minAllGID_, maxAllGID_;
    //! Indicates that the range of global indices are contiguous and ordered.
    bool contiguous_;
    //! Indicates that global indices of the map are non-identically distributed among different nodes.
    bool distributed_;
    //! A direct mapping from local IDs to global IDs.
    mutable Teuchos::ArrayRCP<GlobalOrdinal> lgMap_;
    //! A mapping from global IDs to local IDs.
    std::map<GlobalOrdinal, LocalOrdinal> glMap_;
    //! A Directory for looking up nodes for this Map. This directory has an rcp(this,false) and is therefore not allowed to persist beyond
    //! the lifetime of this Map. Do not under any circumstance pass this outside of the Map.
    Teuchos::RCP< Directory<LocalOrdinal,GlobalOrdinal,Node> > directory_;

  }; // Map class

  /** \brief Non-member function to create a locally replicated Map with the default node.

      This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

      The Map is configured to use zero-based indexing.

      \relates Map
   */
  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Kokkos::DefaultNode::DefaultNodeType> >
  createLocalMap(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm);

  /** \brief Non-member function to create a locally replicated Map with a specified node.

      The Map is configured to use zero-based indexing.

      \relates Map
   */
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> >
  createLocalMapWithNode(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node);

  /** \brief Non-member function to create a uniform, contiguous Map with the default node.

      This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

      The Map is configured to use zero-based indexing.

      \relates Map
   */
  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Kokkos::DefaultNode::DefaultNodeType> >
  createUniformContigMap(global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm);

  /** \brief Non-member function to create a uniform, contiguous Map with a user-specified node.

      The Map is configured to use zero-based indexing.

      \relates Map
   */
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> >
  createUniformContigMapWithNode(global_size_t numElements,
                                 const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node);

  /** \brief Non-member function to create a (potentially) non-uniform, contiguous Map with the default node.

      This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

      The Map is configured to use zero-based indexing.

      \relates Map
   */
  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Kokkos::DefaultNode::DefaultNodeType> >
  createContigMap(global_size_t numElements, size_t localNumElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm);

  /** \brief Non-member function to create a (potentially) non-uniform, contiguous Map with a user-specified node.

      The Map is configured to use zero-based indexing.

      \relates Map
   */
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> >
  createContigMapWithNode(global_size_t numElements, size_t localNumElements, 
                          const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node);

  /** \brief Non-member function to create a non-contiguous Map with the default node.

      This method returns a Map instantiated on the Kokkos default node type, Kokkos::DefaultNode::DefaultNodeType.

      The Map is configured to use zero-based indexing.

      \relates Map
   */
  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Kokkos::DefaultNode::DefaultNodeType> >
  createNonContigMap(const ArrayView<const GlobalOrdinal> &elementList,
                     const RCP<const Teuchos::Comm<int> > &comm);

  /** \brief Non-member function to create a non-contiguous Map with a user-specified node.

      The Map is configured to use zero-based indexing.

      \relates Map
   */
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> >
  createNonContigMapWithNode(const ArrayView<const GlobalOrdinal> &elementList,
                             const RCP<const Teuchos::Comm<int> > &comm, 
                             const RCP<Node> &node);

  /** \brief Non-member function to create a contiguous Map with user-defined weights and a user-specified node.

      The Map is configured to use zero-based indexing.

      \relates Map
   */
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> >
  createWeightedContigMapWithNode(int thisNodeWeight, global_size_t numElements, 
                                  const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node);

} // Tpetra namespace

/** \brief  Returns true if \c map is identical to this map. Implemented in Tpetra::Map::isSameAs().
    \relates Tpetra::Map */
template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool operator== (const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &map1, const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &map2)
{ return map1.isSameAs(map2); }

/** \brief Returns true if \c map is not identical to this map. Implemented in Tpetra::Map::isSameAs().
    \relates Tpetra::Map */
template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool operator!= (const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &map1, const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &map2)
{ return !map1.isSameAs(map2); }

#endif // TPETRA_MAP_DECL_HPP

